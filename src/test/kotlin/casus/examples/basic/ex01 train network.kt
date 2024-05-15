package casus.examples.basic

import casus.mala.common.Parameters
import casus.mala.common.parameters
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
        }
        /* ==========
        2. DATA
            Data has to be added to the MALA workflow. The central object for this is the DataHandler class, which
            takes care of all data needs. After data has been added, it is loaded and scaled with the prepare_data function.
        ==========   */
    }
}