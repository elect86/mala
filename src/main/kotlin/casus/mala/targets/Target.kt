package casus.mala.targets

import casus.mala.common.Parameters
import casus.mala.common.ParametersInterface
import casus.mala.common.PhysicalData

/** Base class for all target quantity parser.
 *
 *  Target parsers read the target quantity (i.e. the quantity the NN will learn to predict) from a specified file
 *  format and performs postprocessing calculations on the quantity. */
abstract class Target(paramsInterface: ParametersInterface): PhysicalData(paramsInterface) {

    var parametersFull: Parameters? = null
    val parameters = when (paramsInterface) {
        is Parameters -> {
            parametersFull = paramsInterface
            paramsInterface.targets
        }
        is Parameters.Targets -> paramsInterface
    }
//    self.fermi_energy_dft = None
//    self.temperature = None
//    self.voxel = None
//    self.number_of_electrons_exact = None
//    self.number_of_electrons_from_eigenvals = None
//    self.band_energy_dft_calculation = None
//    self.total_energy_dft_calculation = None
//    self.entropy_contribution_dft_calculation = None
//    self.atoms = None
//    self.electrons_per_atom = None
//    self.qe_input_data = {
//        "occupations": "smearing",
//        "calculation": "scf",
//        "restart_mode": "from_scratch",
//        "prefix": "MALA",
//        "pseudo_dir": self.parameters.pseudopotential_path,
//        "outdir": "./",
//        "ibrav": None,
//        "smearing": "fermi-dirac",
//        "degauss": None,
//        "ecutrho": None,
//        "ecutwfc": None,
//        "nosym": True,
//        "noinv": True,
//    }
//
//    # It has been shown that the number of k-points
//    # does not really affect the QE post-processing calculation.
//    # This is because we only evaulate density-dependent contributions
//    # with QE. However, there were some (very small) inaccuracies when
//    # operating only at the gamma point. Consequently, MALA defaults
//    # to a small k-grid to ensure best accuracy and performance.
//    # UPDATE 23.04.2021: As per discussion bewteen Normand Modine and Lenz
//    # Fiedler, the number of k-points is moved to 1.
//    # The small inaccuracies are neglected for now.
//    self.kpoints = None  # (2, 2, 2)
//    self.qe_pseudopotentials = {}
//
//    # Local grid and parallelization info for distributed inference.
//    self.local_grid = None
//    self.y_planes = None
//
//    # Control whether target data will be saved.
//    # Can be important for I/O applications.
//    self.save_target_data = True

    /** Get dimension of this target if used as feature in ML. */
    abstract val feaureSize: Int

    companion object {
        /** Create a Target instance.
         *
         *  The correct type of target calculator will automatically be instantiated by this class if possible. You can
         *  also instantiate the desired target directly by calling upon the subclass.
         *  @param parameters: Parameters used to create this target calculator. */
        infix fun from(parameters: ParametersInterface): Target {
            val type = when (parameters) {
                is Parameters -> parameters.targets.target
                is Parameters.Targets -> parameters.target
            }
            val target = when(type) {
                Parameters.Target.LDOS -> LDOS(parameters)
                Parameters.Target.DOS -> DOS(parameters)
                Parameters.Target.Density -> Density(parameters)
            }
            return target
        }
    }
}