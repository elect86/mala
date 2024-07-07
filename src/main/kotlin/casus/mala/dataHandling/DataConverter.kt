package casus.mala.dataHandling

import casus.mala.common.Parameters
import casus.mala.common.PhysicalData
import casus.mala.descriptors.Descriptor
import casus.mala.targets.Target
import java.io.File

/**
 * Converts raw snapshots (direct DFT outputs) into numpy arrays for MALA.
 *
 *     These snapshots can be e.g. Quantum Espresso results.
 *
 *     Parameters
 *     ----------
 *     parameters : mala.common.parameters.Parameters
 *         The parameters object used for creating this instance.
 *
 *     descriptor_calculator : mala.descriptors.descriptor.Descriptor
 *         The descriptor calculator used for parsing/converting fingerprint
 *         data. If None, the descriptor calculator will be created by this
 *         object using the parameters provided. Default: None
 *
 *     target_calculator : mala.targets.target.Target
 *         Target calculator used for parsing/converting target data. If None,
 *         the target calculator will be created by this object using the
 *         parameters provided. Default: None
 *
 *     Attributes
 *     ----------
 *     descriptor_calculator : mala.descriptors.descriptor.Descriptor
 *         Descriptor calculator used for parsing/converting fingerprint data.
 *
 *     target_calculator : mala.targets.target.Target
 *         Target calculator used for parsing/converting target data.
 */
class DataConverter(val parametersFull: Parameters,
                    descriptorCalculator: Descriptor? = null,
                    targetCalculator: Parameters.Target? = null) {

    val parameters = parametersFull.data
    val targetCalculator = targetCalculator ?: Target.from(parametersFull)
    val descriptorCalculator = descriptorCalculator ?: Descriptor.from(parametersFull)

    init {
        if (parametersFull.descriptors.useZSplitting) {
            parametersFull.descriptors.useZSplitting = false
            println("Disabling z-splitting for preprocessing." /*min_verbosity = 0*/)
        }
    }

    val snapshotsToConvert = ArrayList<Snap>()
    val snapshotDescription = ArrayList<Desc>()
    val snapshotUnits = ArrayList<Pair<String?, String?>>()

    // Keep track of what has to be done by this data converter.
    var processDescriptors = false
    var processTargets = false
    var processAdditionalInfo = false

    /**
     * Add a snapshot to be processed.
     *
     * Parameters
     * ----------
     * @param descriptorInputType Type of descriptor data to be processed.
     * See mala.datahandling.data_converter.descriptor_input_types for options.
     *
     * @param descriptorInputPath Path of descriptor data to be processed.
     *
     * @param targetInputType Type of target data to be processed.
     * See mala.datahandling.data_converter.target_input_types for options.
     *
     * @param targetInputPath Path of target data to be processed.
     *
     * @param additionalInfoInputType Type of additional info data to be processed.
     * See mala.datahandling.data_converter.additional_info_input_types for options.
     *
     * @param additionalInfoInputPath Path of additional info data to be processed.
     *
     * @param metadataInputType Type of additional metadata to be processed.
     * See mala.datahandling.data_converter.additional_info_input_types for options.
     * This is essentially the same as additional_info_input_type, but will not affect saving; i.e., the data given
     * here will only be saved in OpenPMD files, not saved separately.
     * If additional_info_input_type is set, this argument will be ignored.
     *
     * @param metadataInputPath Path of additional metadata to be processed.
     * See metadata_input_type for extended info on use.
     *
     * @param descriptorUnits Units for descriptor data processing.
     *
     * @param targetUnits Units for target data processing.
     */
    fun addSnapshot(descriptorInputType: DescriptorInputType? = null,
                    descriptorInputPath: File? = null,
                    targetInputType: TargetInputType? = null,
                    targetInputPath: File? = null,
                    additionalInfoInputType: AdditionalInfoInputType? = null,
                    additionalInfoInputPath: File? = null,
                    descriptorUnits: String? = null,
                    metadataInputType: AdditionalInfoInputType? = null,
                    metadataInputPath: File? = null,
                    targetUnits: String? = null) {
        var metadataInputType = metadataInputType
        var metadataInputPath = metadataInputPath
        // Check the input.
        if (descriptorInputType != null) {
            if (descriptorInputPath == null)
                error("Cannot process descriptor data with no path given.")
            processDescriptors = true
        }
        if (targetInputType != null) {
            if (targetInputPath == null)
                error("Cannot process target data with no path given.")
            processTargets = true
        }
        if (additionalInfoInputType != null) {
            metadataInputType = additionalInfoInputType
            if (additionalInfoInputPath == null)
                error("Cannot process additional info data with no path given.")
            processAdditionalInfo = true
        }
        metadataInputPath = additionalInfoInputPath

        if (metadataInputType != null && metadataInputPath == null)
            error("Cannot process additional info data with no path given.")

        // Assign info.
        snapshotsToConvert += Snap(descriptorInputPath!!, targetInputPath!!,
                                   additionalInfoInputPath!!, metadataInputPath!!)
        snapshotDescription += Desc(descriptorInputType!!, targetInputType!!,
                                    additionalInfoInputType!!, metadataInputType!!)
        snapshotUnits += descriptorUnits to targetUnits
    }

    /**
     * Convert the snapshots in the list to numpy arrays.
     *
     * These can then be used by MALA.
     *
     * @param completeSavePath If not null, the directory in which all snapshots will be saved.
     * Overwrites [descriptorSavePath], [targetSavePath] and [additionalInfoSavePath] if set.
     * @param descriptorSavePath Directory in which to save descriptor data.
     * @param targetSavePath Directory in which to save target data.
     * @param additionalInfoSavePath Directory in which to save additional info data.
     * @param namingScheme String detailing the naming scheme for the snapshots. * symbols
     * will be replaced with the snapshot number.
     * @param startsAt Number of the first snapshot generated using this approach.
     * Default is 0, but may be set to any integer. This is to ensure consistency in naming when converting e.g. only a
     * certain portion of all available snapshots. If set to e.g. 4, the first snapshot generated will be called snapshot4.
     * @param fileBasedCommunication If True, the LDOS will be gathered using a file based mechanism.
     * This is drastically less performant then using MPI, but may be necessary when memory is scarce.
     * Default is False, i.e., the faster MPI version will be used.
     * @param descriptorCalculationKwargs Dictionary with additional keyword arguments for the calculation
     * or parsing of the target quantities.
     * @param targetCalculatorKwargs Dictionary with additional keyword arguments for the calculation
     * or parsing of the descriptor quantities.
     * @param useFp64 If True, data is saved with double precision. If False (default), single precision (FP32) is used.
     * This is advantageous, since internally, torch models are constructed with FP32 anyway.
     */
    fun convertSnapshots(completeSavePath: File? = null,
                         descriptorSavePath: File? = null,
                         targetSavePath: File? = null,
                         additionalInfoSavePath: File? = null,
                         namingScheme: String = "ELEM_snapshot*.npy",
                         startsAt: Int = 0,
                         fileBasedCommunication: Boolean = false,
                         descriptorCalculationKwargs: MutableMap<String, Any> = mutableMapOf(),
                         targetCalculatorKwargs: MutableMap<String, String> = mutableMapOf(),
                         useFp64: Boolean = false) {

        var namingScheme = namingScheme
        val fileEnding = when {
            '.' in namingScheme -> {
                val split = namingScheme.split('.')
                val fileEnding = split.last()
                namingScheme = split.first()
                if (fileEnding != "npy")
                    TODO()
                //                import openpmd_api as io
                //                        if file_ending not in io . file_extensions :
                //                raise Exception ("Invalid file ending selected: " + file_ending)
                fileEnding
            }
            else -> "npy"
        }

        var fileBasedCommunication = fileBasedCommunication
        if (fileEnding == "npy")
        // I will leave the deprecation warning out for now, we re-enable
        // it as soon as we have a precise timeline.
        // parallel_warn("NumPy array based file saving will be deprecated"
        //               "starting in MALA v1.3.0.", min_verbosity=0,
        //               category=FutureWarning)
            Unit
        else
            fileBasedCommunication = false

        var descriptorSavePath = descriptorSavePath
        var targetSavePath = targetSavePath
        var additionalInfoSavePath = additionalInfoSavePath
        if (completeSavePath != null) {
            descriptorSavePath = completeSavePath
            targetSavePath = completeSavePath
            additionalInfoSavePath = completeSavePath
        } else {
            if (processTargets && targetSavePath == null)
                error("No target path specified, cannot process data.")
            if (processDescriptors && descriptorSavePath == null)
                error("No descriptor path specified, cannot process data.")
            if (processAdditionalInfo && additionalInfoSavePath == null)
                error("No additional info path specified, cannot process data.")
        }

        if (fileEnding != "npy") {
            val snapshotName = namingScheme
            TODO()
            //            val seriesName = snapshotName.replace("*", str("%01T"))
            //
            //            if self.process_descriptors:
            //            if self.parameters._configuration["mpi"]:
            //            input_series = io.Series(
            //                os.path.join(
            //                    descriptor_save_path,
            //                    seriesName + ".in." + file_ending,
            //                ),
            //                io.Access.create,
            //                get_comm(),
            //                options = json.dumps(
            //                    self.parameters_full.openpmd_configuration
            //                ),
            //            )
            //            else:
            //            input_series = io.Series(
            //                os.path.join(
            //                    descriptor_save_path,
            //                    seriesName + ".in." + file_ending,
            //                ),
            //                io.Access.create,
            //                options = json.dumps(
            //                    self.parameters_full.openpmd_configuration
            //                ),
            //            )
            //            input_series.set_attribute("is_mala_data", 1)
            //            input_series.set_software(name = "MALA", version = "x.x.x")
            //            input_series.author = "..."
            //
            //            if self.process_targets:
            //            if self.parameters._configuration["mpi"]:
            //            output_series = io.Series(
            //                os.path.join(
            //                    target_save_path,
            //                    seriesName + ".out." + file_ending,
            //                ),
            //                io.Access.create,
            //                get_comm(),
            //                options = json.dumps(
            //                    self.parameters_full.openpmd_configuration
            //                ),
            //            )
            //            else:
            //            output_series = io.Series(
            //                os.path.join(
            //                    target_save_path,
            //                    seriesName + ".out." + file_ending,
            //                ),
            //                io.Access.create,
            //                options = json.dumps(
            //                    self.parameters_full.openpmd_configuration
            //                ),
            //            )
            //
            //            output_series.set_attribute("is_mala_data", 1)
            //            output_series.set_software(name = "MALA", version = mala_version)
            //            output_series.author = "..."
        }

        for (i in snapshotsToConvert.indices) {
            val snapshotNumber = i + startsAt
            var snapshotName = namingScheme
            snapshotName = snapshotName.replace("*", snapshotNumber.toString())

            // Create the paths as needed.
            val infoPath = if (processAdditionalInfo) additionalInfoSavePath!! / "$snapshotName.info.json" else null
            //            input_iteration = None
            //            output_iteration = None

            val descriptorPath: File?
            val targetPath: File?
            var memmap: File?
            if (fileEnding == "npy") {
                // Create the actual paths, if needed.
                descriptorPath = if (processDescriptors) descriptorSavePath!! / "$snapshotName.in.$fileEnding" else null

                memmap = null
                if (processTargets) {
                    targetPath = targetSavePath!! / "$snapshotName.out.$fileEnding"
                    // A memory mapped file is used as buffer for distributed cases.
                    if (parameters.mpi && fileBasedCommunication)
                        memmap = targetSavePath / "$snapshotName.out.npy_temp"
                } else
                    targetPath = null
            } else {
                descriptorPath = null
                targetPath = null
                memmap = null
                TODO()
                //                if (processDescriptors) {
                //                    val input_iteration = input_series.write_iterations()[
                //                        i + starts_at
                //                    ]
                //                    input_iteration.dt = i + starts_at
                //                    input_iteration.time = 0
                //                }
                //                if self.process_targets:
                //                output_iteration = output_series.write_iterations()[
                //                    i + starts_at
                //                ]
                //                output_iteration.dt = i + starts_at
                //                output_iteration.time = 0
            }
            convertSingleSnapshot(i,
                                  descriptorCalculationKwargs,
                                  targetCalculatorKwargs,
                                  inputPath = descriptorPath,
                                  outputPath = targetPath,
                                  useMemmap = memmap,
                //                                  input_iteration = input_iteration,
                //                                  output_iteration = output_iteration,
                                  additionalInfoPath = infoPath,
                                  useFp64 = useFp64)
            TODO()
            //            if get_rank() == 0:
            //            if (
            //                self.parameters._configuration["mpi"]
            //                and file_based_communication
            //            ):
            //            os.remove(memmap)
        }
    }

    /**
     * Convert single snapshot from the conversion lists.
     *
     * Returns the preprocessed data as numpy array, which might be beneficial during testing.
     *
     * Parameters
     * ----------
     * @param snapshotNumber Position of the desired snapshot in the snapshot list.
     * @param descriptorCalculationKwargs Dictionary with additional keyword arguments for the calculation
     * or parsing of the descriptor quantities.
     * @param targetCalculatorKwargs Dictionary with additional keyword arguments for the calculation
     * or parsing of the target quantities.
     * @param inputPath If not None, inputs will be saved in this file.
     * @param outputPath If not None, outputs will be saved in this file.
     * @param useMemmap If not None, a memory mapped file with this name will be used to gather the LDOS.
     * If run in MPI parallel mode, such a file MUST be provided.
     *         output_iteration : OpenPMD iteration
     *              OpenPMD iteration to be used to save the output data of the current
     *              snapshot, as part of an OpenPMD Series.
     *
     *          input_iteration : OpenPMD iteration
     *              OpenPMD iteration to be used to save the input data of the current
     *              snapshot, as part of an OpenPMD Series.
     *
     * @param useFp64 If True, data is saved with double precision. If False (default), single precision (FP32) is used.
     * This is advantageous, since internally, torch models are constructed with FP32 anyway.
     *
     *        return_data : bool
     *             If True, inputs and outputs will be returned directly.
     *
     *         Returns
     *         -------
     *         inputs : numpy.array , optional
     *             Numpy array containing the preprocessed inputs.
     *
     *         outputs : numpy.array , optional
     *             Numpy array containing the preprocessed outputs.
     */
    private fun convertSingleSnapshot(snapshotNumber: Int,
                                      descriptorCalculationKwargs: MutableMap<String, Any>,
                                      targetCalculatorKwargs: MutableMap<String, String>,
                                      inputPath: File? = null,
                                      outputPath: File? = null,
                                      additionalInfoPath: File? = null,
                                      useMemmap: File? = null,
        //                                      output_iteration = None,
        //                                      input_iteration = None,
                                      useFp64: Boolean = false) {

        val snapshot = snapshotsToConvert[snapshotNumber]
        val description = snapshotDescription[snapshotNumber]
        val originalUnits = snapshotUnits[snapshotNumber]

        // Parse and/or calculate the input descriptors.
        when (description.input) {
            DescriptorInputType.`espresso-out` -> {
                originalUnits.first?.let { descriptorCalculationKwargs["units"] = it }
                descriptorCalculationKwargs["useFp64"] = useFp64

                /*tmp_input, local_size =*/ descriptorCalculator.calculateFromQeOut(snapshot.input, kwargs = descriptorCalculationKwargs)
            }
            //In this case, only the output is processed.
            null -> Unit
        }
    }
}

class Snap(val input: File, val output: File,
           val additionalInfo: File, val metadata: File)

class Desc(val input: DescriptorInputType?, val output: TargetInputType,
           val additionalInfo: AdditionalInfoInputType, val metadata: AdditionalInfoInputType)

enum class DescriptorInputType { `espresso-out` }
enum class TargetInputType { cube, xsf }
enum class AdditionalInfoInputType { `espresso-out` }
