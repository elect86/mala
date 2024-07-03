package casus.mala.dataHandling

import ai.djl.ndarray.NDArray
import ai.djl.ndarray.NDManager
import ai.djl.ndarray.types.Shape
import casus.mala.common.DEFAULT_NP_DATA_DTYPE
import casus.mala.common.Parameters
import casus.mala.descriptors.Descriptor
import casus.mala.targets.Target
import java.nio.file.Path

/**
 * Loads and scales data. Can only process numpy arrays at the moment.
 *
 * Data that is not in a numpy array can be converted using the DataConverter class.
 */
class DataHandler(parameters: Parameters,
                  targetCalculator: Target? = null,
                  descriptorCalculator: Descriptor? = null,
                  inputDataScaler: DataScaler? = null,
                  outputDataScaler: DataScaler? = null,
                  clearData: Boolean = true) : DataHandlerBase(parameters, targetCalculator, descriptorCalculator) {

    // Data will be scaled per user specification.
    val inputDataScaler = inputDataScaler ?: DataScaler(this.parameters.inputRescaling, useHorovod)

    val outputDataScaler = outputDataScaler ?: DataScaler(this.parameters.outputRescaling, useHorovod)

    // Actual data points in the different categories.
    var nrTrainingData = 0L
    var nrTestData = 0L
    var nrValidationData = 0L

    // Number of snapshots in these categories.
    var nrTrainingSnapshots = 0
    var nrTestSnapshots = 0
    var nrValidationSnapshots = 0

    // Arrays and data sets containing the actual data.
    val np = NDManager.newBaseManager() // TODO, remember to close it
    lateinit var trainingDataInputs: NDArray
    lateinit var validationDataInputs: NDArray
    lateinit var testDataInputs: NDArray
    lateinit var trainingDataOutputs: NDArray
    lateinit var validationDataOutputs: NDArray
    lateinit var testDataOutputs: NDArray
    val trainingDataSets = ArrayList<NDArray>()
    val validationDataSets = ArrayList<NDArray>()
    val testDataSets = ArrayList<NDArray>()

    // Needed for the fast tensor data sets.
    val miniBatchSize = parameters.running.miniBatchSize

    init {
        if (clearData)
            clearData()
    }

    //##############################
    // Public methods
    //##############################

    // Adding/Deleting data
    // ########################

    /** Reset the entire data pipeline.
     *
     * Useful when doing multiple investigations in the same python file. */
    override fun clearData() {
        //        self.training_data_sets = []
        //        self.validation_data_sets = []
        //        self.test_data_sets = []
        //        self.nr_training_data = 0
        //        self.nr_test_data = 0
        //        self.nr_validation_data = 0
        //        self.nr_training_snapshots = 0
        //        self.nr_test_snapshots = 0
        //        self.nr_validation_snapshots = 0
        super.clearData()
    }

    // Preparing data
    // ######################

    /** Prepare the data to be used in a training process.
     *
     * This includes:
     *
     * - Checking snapshots for consistency
     * - Parametrizing the DataScalers (if desired)
     * - Building DataSet objects.
     *
     * Parameters
     * ----------
     * @reparametrizeScaler: If true (default), the DataScalers are parametrized based on the training data. */
    fun prepareData(reparametrizeScaler: Boolean = true) {
        var reparametrizeScaler = reparametrizeScaler
        // During data loading, there is no need to save target data to calculators.
        // Technically, this would be no issue, but due to technical reasons (i.e. float64 to float32 conversion) saving
        // the data this way may create copies in memory.
        targetCalculator.saveTargetData = false

        // Do a consistency check of the snapshots so that we don't run into an error later.
        // If there is an error, check_snapshots() will raise an exception .
        println(
            "Checking the snapshots and your inputs for consistency.",
            //            min_verbosity = 1,
        )
        checkSnapshots()
        println("Consistency check successful." /*min_verbosity = 0*/)

        // If the DataHandler is used for inference, i.e. no training or validation snapshots have been provided,
        // than we can definitely not reparametrize the DataScalers.
        if (nrTrainingData == 0L) {
            reparametrizeScaler = false
            if (!inputDataScaler.canTransform || !outputDataScaler.canTransform)
                error("In inference mode, the DataHandler needs parametrized DataScalers, while you provided unparametrized DataScalers.")
        }
        // Parametrize the scalers, if needed.
        if (reparametrizeScaler) {
            println("Initializing the data scalers."/*, min_verbosity = 1*/)
            parametrizeScalers()
            println("Data scalers initialized." /*min_verbosity = 0*/)
        } else if (!parameters.useLazyLoading && nrTrainingData != 0L)
            println("Data scalers already initilized, loading data to RAM." /*min_verbosity = 0*/)

        loadData(Snapshot.Function.training, DataType.inputs)
        loadData(Snapshot.Function.training, DataType.outputs)

        // Build Datasets.
        println("Build datasets." /*min_verbosity = 1*/)
        buildDatasets()
        println("Build dataset: Done." /*min_verbosity = 0*/)

        // After the loading is done, target data can safely be saved again.
        targetCalculator.saveTargetData = true

        // Wait until all ranks are finished with data preparation.
        // It is not uncommon that ranks might be asynchronous in their data preparation by a small amount of minutes.
        // If you notice an elongated wait time at this barrier, check that your file system allows for parallel I / O.
//        barrier()
    }

    /** Load data into the appropriate arrays.
     *
     *  Also transforms them into torch tensors.
     *
     *  Parameters
     *  ----------
     *  @param function: [Snapshot.Function]
     *  @param dataType: Can be "input" or "output". */
    fun loadData(function: Snapshot.Function, dataType: DataType) {
        // Extracting all the information pertaining to the data set.
        val array = when (function) {
            Snapshot.Function.training -> if (dataType == DataType.inputs) trainingDataInputs else trainingDataOutputs
            Snapshot.Function.validation -> if (dataType == DataType.inputs) validationDataInputs else validationDataOutputs
            else -> if (dataType == DataType.inputs) testDataInputs else testDataOutputs
        }
        val calculator = if (dataType == DataType.inputs) descriptorCalculator else targetCalculator

        val featureDimension = Shape(if (dataType == DataType.inputs) inputDimension else outputDimension)

        var snapshotCounter = 0
        var gsOld = 0L
        for (snapshot in parameters.snapshotDirectoriesList) {
            // get the snapshot grid size
            val gsNew = snapshot.gridSize

            // Data scaling is only performed on the training data sets .
            if (snapshot.snapshotFunction != function)
                continue
            val file: Path
            val units: String
            if (dataType == DataType.inputs) {
                file = snapshot.inputNpyFile
                units = snapshot.inputUnits!!
            } else {
                file = snapshot.outputNpyFile
                units = snapshot.outputUnits
            }

            if (snapshot.snapshotType == Snapshot.Type.numpy)
                calculator.readFromNumpyFile(file, units, array["$gsOld : ${gsOld + gsNew}, :"], reshape = true)
            else if (snapshot.snapshotType == Snapshot.Type.openpmd)
                TODO()
            //                array["$gsOld : ${gsOld + gsNew}"] = calculator.read_from_openpmd_file(file, units = units).reshape([gsNew, featureDimension])
            snapshotCounter += 1
            gsOld += gsNew
        }
        // The scalers will later operate on torch Tensors so we have to make sure they are fitted on torch Tensors as well.
        // Preprocessing the numpy data as follows does NOT load it into memory, see test / tensor_memory.py
        // Also, the following bit does not work with getattr, so I had to hard code it.
        // If someone has a smart idea to circumvent this, I am all ears .

        /*
        In python training_data_inputs was loaded from file into a np, and then it convert to torch tensor.
        In Java, there is no np in the picture, you can load file directly into PyTorch NDArray
         */
    }

    /** Build the DataSets that are used during training. */
    fun buildDatasets() {
        if (parameters.useLazyLoading && !parameters.useLazyLoadingPrefetch) {
            TODO()
            // Create the lazy loading data sets.
            //            trainingDataSets += LazyLoadDataset(
            //                    self.input_dimension,
            //                    self.output_dimension,
            //                    self.input_data_scaler,
            //                    self.output_data_scaler,
            //                    self.descriptor_calculator,
            //                    self.target_calculator,
            //                    self.use_horovod,
            //                )
            //            )
            //            self.validation_data_sets.append(
            //                LazyLoadDataset(
            //                    self.input_dimension,
            //                    self.output_dimension,
            //                    self.input_data_scaler,
            //                    self.output_data_scaler,
            //                    self.descriptor_calculator,
            //                    self.target_calculator,
            //                    self.use_horovod,
            //                )
            //            )
            //
            //            if self.nr_test_data != 0:
            //            self.test_data_sets.append(
            //                LazyLoadDataset(
            //                    self.input_dimension,
            //                    self.output_dimension,
            //                    self.input_data_scaler,
            //                    self.output_data_scaler,
            //                    self.descriptor_calculator,
            //                    self.target_calculator,
            //                    self.use_horovod,
            //                    input_requires_grad = True,
            //                )
            //            )
            //
            //            # Add snapshots to the lazy loading data sets .
            //            for snapshot in self.parameters.snapshot_directories_list:
            //            if snapshot.snapshot_function == "tr":
            //            self.training_data_sets[0].add_snapshot_to_dataset(
            //                snapshot
            //            )
            //            if snapshot.snapshot_function == "va":
            //            self.validation_data_sets[0].add_snapshot_to_dataset(
            //                snapshot
            //            )
            //            if snapshot.snapshot_function == "te":
            //            self.test_data_sets[0].add_snapshot_to_dataset(snapshot)
        }
        // I don't think we need to mix them here. We can use the standard ordering for the first epoch and mix it up after.
        // self.training_data_set.mix_datasets()
        // self.validation_data_set.mix_datasets()
        // self.test_data_set.mix_datasets()
        else if (parameters.useLazyLoading && parameters.useLazyLoadingPrefetch) {
            println("Using lazy loading pre-fetching." /*min_verbosity = 2*/)
            // Create LazyLoadDatasetSingle instances per snapshot and add to list.
            for (snapshot in parameters.snapshotDirectoriesList)
                TODO()
            //                when (snapshot.snapshotFunction) {
            //                    Snapshot.Function.training -> trainingDataSets += LazyLoadDatasetSingle(
            //                                self.mini_batch_size,
            //                                snapshot,
            //                                self.input_dimension,
            //                                self.output_dimension,
            //                                self.input_data_scaler,
            //                                self.output_data_scaler,
            //                                self.descriptor_calculator,
            //                                self.target_calculator,
            //                                self.use_horovod,
            //                            )
            //                        )
            //                    if snapshot.snapshot_function == "va":
            //                        self.validation_data_sets.append(
            //                            LazyLoadDatasetSingle(
            //                                self.mini_batch_size,
            //                                snapshot,
            //                                self.input_dimension,
            //                                self.output_dimension,
            //                                self.input_data_scaler,
            //                                self.output_data_scaler,
            //                                self.descriptor_calculator,
            //                                self.target_calculator,
            //                                self.use_horovod,
            //                            )
            //                        )
            //                    if snapshot.snapshot_function == "te":
            //                        self.test_data_sets.append(
            //                            LazyLoadDatasetSingle(
            //                                self.mini_batch_size,
            //                                snapshot,
            //                                self.input_dimension,
            //                                self.output_dimension,
            //                                self.input_data_scaler,
            //                                self.output_data_scaler,
            //                                self.descriptor_calculator,
            //                                self.target_calculator,
            //                                self.use_horovod,
            //                                input_requires_grad = True,
            //                            )
            //                        )
            //                }
        } else {
            TODO()
//            if (nrTrainingData != 0L) {
//                inputDataScaler.transform(self.training_data_inputs)
//                output_data_scaler.transform(self.training_data_outputs)
//                if self.parameters.use_fast_tensor_data_set:
//                printout("Using FastTensorDataset.", min_verbosity = 2)
//                self.training_data_sets.append(
//                    FastTensorDataset(
//                        self.mini_batch_size,
//                        self.training_data_inputs,
//                        self.training_data_outputs,
//                    )
//                )
//                else:
//                self.training_data_sets.append(
//                    TensorDataset(
//                        self.training_data_inputs,
//                        self.training_data_outputs,
//                    )
//                )
//            }
//            if self.nr_validation_data != 0:
//            self.__load_data("validation", "inputs")
//            self.input_data_scaler.transform(self.validation_data_inputs)
//
//            self.__load_data("validation", "outputs")
//            self.output_data_scaler.transform(self.validation_data_outputs)
//            if self.parameters.use_fast_tensor_data_set:
//            printout("Using FastTensorDataset.", min_verbosity = 2)
//            self.validation_data_sets.append(
//                FastTensorDataset(
//                    self.mini_batch_size,
//                    self.validation_data_inputs,
//                    self.validation_data_outputs,
//                )
//            )
//            else:
//            self.validation_data_sets.append(
//                TensorDataset(
//                    self.validation_data_inputs,
//                    self.validation_data_outputs,
//                )
//            )
//
//            if self.nr_test_data != 0:
//            self.__load_data("test", "inputs")
//            self.input_data_scaler.transform(self.test_data_inputs)
//            self.test_data_inputs.requires_grad = True
//
//            self.__load_data("test", "outputs")
//            self.output_data_scaler.transform(self.test_data_outputs)
//            self.test_data_sets.append(
//                TensorDataset(
//                    self.test_data_inputs, self.test_data_outputs
//                )
//            )
        }
    }

    // Scaling
    //#####################

    /** Use the training data to parametrize the DataScalers. */
    private fun parametrizeScalers() {
        //#################
        // Inputs.
        //#################

        // If we do lazy loading, we have to iterate over the files one at a time and add them to the fit,
        // i.e. incrementally updating max/min or mean/std. If we DON'T do lazy loading, we can simply load the
        // training data (we will need it later anyway) and perform the scaling. This should save some performance.

        if (parameters.useLazyLoading) {
            TODO()
            //            inputDataScaler.startIncrementalFitting()
            //            // We need to perform the data scaling over the entirety of the training data .
            //            for (snapshot in parameters.snapshotDirectoriesList) {
            //                // Data scaling is only performed on the training data sets .
            //                if (snapshot.snapshotFunction != Snapshot.Function.training)
            //                    continue
            //                tmp = when (snapshot.snapshotType) {
            //                    Snapshot.Type.numpy -> descriptorCalculator.read_from_numpy_file(
            //                        os.path.join(
            //                            snapshot.input_npy_directory,
            //                            snapshot.input_npy_file,
            //                        ),
            //                        units = snapshot.input_units,
            //                    )
            //                            elif snapshot.snapshot_type == "openpmd":
            //                        tmp = (
            //                                self.descriptor_calculator.read_from_openpmd_file(
            //                                    os.path.join(
            //                                        snapshot.input_npy_directory,
            //                                        snapshot.input_npy_file,
            //                                    )
            //                                )
            //                              )
            //                }
            //
            //                # The scalers will later operate on torch Tensors so we
            //                # have to make sure they are fitted on
            //                # torch Tensors as well . Preprocessing the numpy data as
            //                        # follows does NOT load it into memory, see
            //                # test / tensor_memory.py
            //                tmp = np.array(tmp)
            //                if tmp.dtype != DEFAULT_NP_DATA_DTYPE:
            //                tmp = tmp.astype(DEFAULT_NP_DATA_DTYPE)
            //                tmp = tmp.reshape(
            //                    [snapshot.grid_size, self.input_dimension]
            //                )
            //                tmp = torch.from_numpy(tmp).float()
            //                self.input_data_scaler.incremental_fit(tmp)
            //            }
            //            self.input_data_scaler.finish_incremental_fitting()
        } else {
            loadData(Snapshot.Function.training, DataType.inputs)
            inputDataScaler fit trainingDataInputs
        }

        println("Input scaler parametrized." /*min_verbosity = 1*/)

        //#################
        // Output.
        //#################

        // If we do lazy loading, we have to iterate over the files one at a time and add them to the fit,
        // i.e.incrementally updating max / min or mean / std.
        // If we DON'T do lazy loading, we can simply load the training data (we will need it later anyway)
        // and perform the scaling . This should save some performance.

        if (parameters.useLazyLoading) {
            //            var i = 0
            outputDataScaler.startIncrementalFitting()
            // We need to perform the data scaling over the entirety of the training data .
            for (snapshot in parameters.snapshotDirectoriesList) {
                // Data scaling is only performed on the training data sets.
                if (snapshot.snapshotFunction == Snapshot.Function.training) {
                    val tmp = when (snapshot.snapshotType) {
                        Snapshot.Type.numpy -> targetCalculator.readFromNumpyFile(snapshot.outputNpyFile,
                                                                                  units = snapshot.outputUnits)
                        else -> TODO() // targetCalculator.read_from_openpmd_file (
                        //                        os.path.join(
                        //                            snapshot.output_npy_directory,
                        //                            snapshot.output_npy_file,
                        //                        )
                    }

                    // The scalers will later operate on torch Tensors so we have to make sure they are fitted on
                    // torch Tensors as well.
                    // Preprocessing the numpy data as follows does NOT load it into memory, see test / tensor_memory.py
                    TODO()
                    //                    tmp = np.array(tmp)
                    //                    if tmp.dtype != DEFAULT_NP_DATA_DTYPE:
                    //                    tmp = tmp.astype(DEFAULT_NP_DATA_DTYPE)
                    //                    tmp = tmp.reshape(
                    //                        [snapshot.grid_size, self.output_dimension]
                    //                    )
                    //                    tmp = torch.from_numpy(tmp).float()
                    //                    self.output_data_scaler.incremental_fit(tmp)
                }
                //                i += 1
            }
            outputDataScaler.finishIncrementalFitting()
        } else {
            loadData(Snapshot.Function.training, DataType.outputs)
            outputDataScaler fit trainingDataOutputs
        }

        println("Output scaler parametrized." /*min_verbosity = 1*/)
    }

    // #############################
    //  Private methods
    // #############################

    //  Loading data
    // #####################

    /** Check the snapshots for consistency. */
    override fun checkSnapshots() {
        super.checkSnapshots()

        // Now we need to confirm that the snapshot list has some inner consistency.
        check(parameters.dataSplitting == Parameters.Splitting.bySnapshot) { "Wrong parameter for data splitting provided." }
        // As we are not actually interested in the number of snapshots,
        // but in the number of datasets, we also need to multiply by that.
        for (snapshot in parameters.snapshotDirectoriesList)
            when (snapshot.snapshotFunction) {
                Snapshot.Function.training -> {
                    nrTrainingSnapshots += 1
                    nrTrainingData += snapshot.gridSize
                }
                Snapshot.Function.testing -> {
                    nrTestSnapshots += 1
                    nrTestData += snapshot.gridSize
                }
                else -> {
                    nrValidationSnapshots += 1
                    nrValidationData += snapshot.gridSize
                }
            }

        // Now we need to check whether or not this input is believable.
        val nrOfSnapshots = parameters.snapshotDirectoriesList.size
        check(nrOfSnapshots == nrTrainingSnapshots + nrTestSnapshots + nrValidationSnapshots) {
            "Cannot split snapshots with specified splitting scheme, too few or too many options selected"
        }
        // MALA can either be run in training or test - only mode.
        // But it has to be run in either of those!
        // So either training AND validation snapshots can be provided OR only test snapshots .
        if (nrTestSnapshots != 0) {
            if (nrTrainingSnapshots == 0)
                println(
                    "DataHandler prepared for inference. No training possible with this setup. If this is not what you " +
                    "wanted, please revise the input script. Validation snapshots you may have entered will be ignored.",
                    //            min_verbosity = 0,
                )
        } else {
            if (nrTrainingSnapshots == 0) error("No training snapshots provided.")
            if (nrValidationSnapshots == 0) error("No validation snapshots provided.")
        }

        if (!parameters.useLazyLoading)
            allocateArrays()

        // Reordering the lists.snapshot_order = { "tr": 0, "va": 1, "te": 2 }
        parameters.snapshotDirectoriesList.sortBy(Snapshot::snapshotFunction)
    }

    private fun allocateArrays() {
        if (nrTrainingData > 0) {
            trainingDataInputs = np.zeros(Shape(nrTrainingData, inputDimension), DEFAULT_NP_DATA_DTYPE)
            trainingDataOutputs = np.zeros(Shape(nrTrainingData, outputDimension), DEFAULT_NP_DATA_DTYPE)
        }

        if (nrValidationData > 0) {
            validationDataInputs = np.zeros(Shape(nrValidationData, inputDimension), DEFAULT_NP_DATA_DTYPE)
            validationDataOutputs = np.zeros(Shape(nrValidationData, outputDimension), DEFAULT_NP_DATA_DTYPE)
        }

        if (nrTestData > 0) {
            testDataInputs = np.zeros(Shape(nrTestData, inputDimension), DEFAULT_NP_DATA_DTYPE)
            testDataOutputs = np.zeros(Shape(nrTestData, outputDimension), DEFAULT_NP_DATA_DTYPE)
        }
    }

    enum class DataType { inputs, outputs }
}