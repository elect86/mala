package casus.examples.basic

import casus.mala.common.Parameters
import casus.mala.common.parameters
import casus.mala.dataHandling.AdditionalInfoInputType
import casus.mala.dataHandling.DataConverter
import casus.mala.dataHandling.DescriptorInputType
import casus.mala.dataHandling.TargetInputType
import org.junit.jupiter.api.Test
import java.io.File
import java.lang.foreign.Arena
import java.lang.foreign.FunctionDescriptor
import java.lang.foreign.Linker
import java.lang.foreign.ValueLayout
import kotlin.io.path.div


class `ex03 preprocess data` {

    @Test
    fun preprocess() {

        /*
        Shows how this framework can be used to preprocess
        data. Preprocessing here means converting raw DFT calculation output into
        numpy arrays of the correct size. For the input data, this means descriptor
        calculation.

        REQUIRES LAMMPS.
         */

        //###################
        // 1. PARAMETERS
        // Select the correct parameters with which to preprocess the electronic
        // structure simulation data. For the LDOS, these parameters have to be chosen
        // consistent with the DFT calculation (i.e., if you specified an energy grid
        //   # with 250 values lying between -10 and 15 eV, and therefore 0.1 eV between
        // grid points, you will have to set ldos_gridsize to 250, gridspacing_ev to 0.1
        // and grid_offset_ev to -10). Of course, the values are dependent on the
        // physical system simulated.
        // Likewise, the values for the bispectrum descriptors have to be chosen such
        // that the corresponding descriptors accurately represent atomic enviroment.
        // In the advanced example section, it is shown how a newly developed analysis
        // called ACSD can be used to determine these values automatically.
        // For now, we will use standard values.
        //###################

        parameters {
            // Bispectrum parameters.
            descriptors.type = Parameters.Descriptors.Type.Bispectrum
            descriptors.bispectrumTwojmax = 10
            descriptors.bispectrumCutoff = 4.67637f
            // LDOS parameters.
            targets.target = Parameters.Target.LDOS
            targets.ldosGridsize = 11
            targets.ldosGridspacingEv = 2.5f
            targets.ldosGridoffsetEv = -5f

            //###################
            // 2. ADDING DATA FOR DATA CONVERSION
            // Data conversion itself is simple. We select input and output data
            // to be converted, add this data snapshot-wise and tell MALA to
            // convert snapshots. Inputs and outputs can be processed individually.
            // Further, via the additional_info_input_* keywords, calculation output
            // can be processed from the original simulation *.out output files into
            // more convenient *.json files that can be used in their stead. This saves
            // on disk space.
            // To only process parts of the data, omit/add descriptor_input*, target_input_*
            // and additional_info_input_* at your leisure.
            // Make sure to set the correct units - for QE, this should always be
            // 1/(Ry*Bohr^3).
            //###################

            val dataConverter = DataConverter(this)
            val outFile = be2 / "Be_snapshot0.out"
            val ldosFile = be2 / "cubes/tmp.pp*Be_ldos.cube"

            dataConverter.addSnapshot(DescriptorInputType.`espresso-out`, outFile.toFile(),
                                      TargetInputType.cube, ldosFile.toFile(),
                                      AdditionalInfoInputType.`espresso-out`, outFile.toFile(),
                                      targetUnits = "1/(Ry*Bohr^3)")

            //###################
            // 3. Converting the data
            // To convert the data we now simply have to call the convert_snapshot function.
            // Input (descriptor) and output (target) data can be saved in individual
            // locations, as can the additional output data files. Fine-granular access
            // to the calculators is enabled via the descriptor_calculation_kwargs and
            // target_calculation_kwargs arguments, but usually not needed.
            // If all data files should be stored in the same location, the
            // complete_save_path keyword may be used.
            //###################

            dataConverter.convertSnapshots(descriptorSavePath = File("./"),
                                           targetSavePath = File("./"),
                                           additionalInfoSavePath = File("./"),
                                           namingScheme = "Be_snapshot*.npy",
                                           descriptorCalculationKwargs = mutableMapOf("workingDirectory" to be2,
                                                                                      "profileCalculation" to true))


        }
    }
}