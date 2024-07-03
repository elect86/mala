package casus.examples.basic

import casus.mala.common.Parameters
import casus.mala.common.parameters
import casus.mala.dataHandling.Snapshot
import casus.mala.network.Network
import casus.mala.network.Trainer
import kotlin.io.path.div
import kotlin.test.Test

class `ex01 train network` {

    @Test
    fun run() {
        /* ==========
        1. PARAMETERS
            The first step of each MALA workflow is to define a parameters object and
            select the necessary parameters for the application one wants to look into.
        ========== */
        parameters {
            // Specify the data scaling. For regular bispectrum and LDOS data, these have proven successful.
            data {
                inputRescaling = Parameters.Scaling.`feature wise standard`
                outputRescaling = Parameters.Scaling.normal
            }
            // Specify the used activation function.
            network.layerActivations = listOf(Parameters.Activation.ReLU)
            // Specify the training parameters.
            // These may be determined via hyperparameter tuning.
            running {
                maxNumberEpochs = 100
                miniBatchSize = 40
                learningRate = 0.00001f
                training = Parameters.Training.Adam
            }
            /* These parameters characterize how the LDOS and bispectrum descriptors were calculated.
               They are _technically_ not needed to train a simple network. However, it is useful to define them prior
               to training. Then, when using the network later in production, all required parameters are already set. */
            targets {
                target = Parameters.Target.LDOS
                ldosGridsize = 11
                ldosGridspacingEv = 2.5f
                ldosGridoffsetEv = -5f
            }
            descriptors {
                type = Parameters.Descriptors.Type.Bispectrum
                bispectrumTwojmax = 10
                bispectrumCutoff = 4.67637f
            }

            /* ==========
            2. DATA
                Data has to be added to the MALA workflow. The central object for this is the DataHandler class, which
                takes care of all data needs. After data has been added, it is loaded and scaled with the prepare_data function.
            ==========   */

            // Add a snapshot we want to use in to the list.
            dataHandler += Snapshot(be2 / "Be_snapshot0.in.npy", be2 / "Be_snapshot0.out.npy",
                                    Snapshot.Function.training, snapshotType = Snapshot.Type.numpy)
            dataHandler += Snapshot(be2 / "Be_snapshot1.in.npy", be2 / "Be_snapshot1.out.npy",
                                    Snapshot.Function.validation, snapshotType = Snapshot.Type.numpy)
            dataHandler.prepareData()

            /* ####################
            3. NETWORK SETUP
                Now we can set up the NN to be used in the ML-DFT model. The layer_sizes
                list determines the number of neurons in the NN. It can be specified before
                loading data, but it is recommended to do that afterwards, since then
                the input_dimension and output_dimension properties of the data handling
                class can be used to correctly define input and output layer of the NN.
            ################### */
            network.layerSizes = intArrayOf(dataHandler.inputDimension.toInt(),
                                            100,
                                            dataHandler.outputDimension.toInt())
            val testNetwork = Network(this)
            /* ####################
             4. TRAINING THE NETWORK
                 Finally, the network can be trained. Afterwards, it can easily be saved
                 into a .zip archive for inference. It is recommended to load a file
                 containing additional calculation data (e.g., from the QE calculations
                 with which the LDOS data was created) so that things like simulated
                 temperature, information about the pseudopotential, etc. are stored along-
                 side the model. This makes inference easier.
            ################### */

            val testTrainer = Trainer(this, testNetwork, dataHandler)
            testTrainer.trainNetwork()
//            additional_calculation_data = os.path.join(data_path, "Be_snapshot0.out")
//            testTrainer.save_run(
//                "be_model", additional_calculation_data=additional_calculation_data
//            )
        }
    }
}