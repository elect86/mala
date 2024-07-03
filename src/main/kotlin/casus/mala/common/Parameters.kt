package casus.mala.common

import casus.mala.dataHandling.DataHandler
import casus.mala.dataHandling.Snapshot
import java.io.File
import kotlin.properties.Delegates

fun parameters(init: Parameters.() -> Unit) = Parameters().apply(init)

class Parameters: ParametersInterface {

    var comment = ""

    // Parameters subobjects.
    val network = ParametersNetwork()
    val descriptors = Descriptors()
    val targets = Targets()
    val data = Data()
    val running = ParametersRunning()
//    self.hyperparameters = ParametersHyperparameterOptimization()
//    self.datageneration = ParametersDataGeneration()

    // Attributes.
    var manualSeed: Int? = null

    // Properties
    var useGpu = false
    var useHorovod = false
    var useMpi = false
    var verbosity = 1
    var device = "cpu"
//    self.openpmd_configuration = {}
    // TODO: Maybe as a percentage? Feature dimensions can be quite different.
    var openpmdGranularity = 1
    var useLammps = true


    fun data(init: Data.() -> Unit) = data.apply(init)

    /** Parameters necessary for loading and preprocessing data. */
    class Data: ParametersBase() {
        /** A list of all added snapshots. */
        val snapshotDirectoriesList = ArrayList<Snapshot>()

        /** Specify how the data for validation, test and training is split. Currently, the only supported option
         *  is [Splitting.bySnapshot], which splits the data by snapshot boundaries. It is also the default. */
        val dataSplitting = Splitting.bySnapshot

        /** Specifies how input quantities are normalized. */
        var inputRescaling = Scaling.none

        /** Specifies how output quantities are normalized. */
        var outputRescaling = Scaling.none

        /** If True, data is lazily loaded, i.e. only the snapshots that are currently needed will be kept in memory.
         *  This greatly reduces memory demands, but adds additional computational time. */
        var useLazyLoading = false

        /** If True, will use alternative lazy loading path with prefetching for higher performance */
        var useLazyLoadingPrefetch = false

        /** If True, then the new, fast TensorDataSet implemented by Josh Romero will be used. */
        val useFastTensorDataSet = false

        /** If not None, a seed that will be used to make the shuffling of the data in the DataShuffler class deterministic. */
        var shufflingSeed: Int? = null
    }

    enum class Splitting { bySnapshot }
    enum class Scaling {
        /** No normalization is applied */
        none,

        /** Standardization (Scale to mean 0, standard deviation 1) */
        standard,

        /** Min-Max scaling (Scale to be in range 0...1) */
        normal,

        /** Row Standardization (Scale to mean 0, standard deviation 1) */
        `feature wise standard`,

        /** Row Min-Max scaling (Scale to be in range 0...1) */
        `feature wise normal`
    }

    fun network(init: ParametersNetwork.() -> Unit) = network.apply(init)

    enum class Nn { `feed-forward`, transformer, lstm, gru }
    enum class Activation { Sigmoid, ReLU, LeakyReLU }
    enum class LossFunction {
        /** Mean squared error */
        mse
    }

    fun running(init: ParametersRunning.() -> Unit) = running.apply(init)

    enum class Training {
        /** Stochastic gradient descent */
        SGD,

        /** Adam Optimization Algorithm */
        Adam
    }

    enum class RateScheduler {
        /** No learning rate schedule will be used. */
        None,

        /** The learning rate will be reduced when the validation loss is plateauing. */
        ReduceLROnPlateau
    }

    enum class TrainingMetric { ldos, bandEnergy, totalEnergy }

    fun targets(init: Targets.() -> Unit) = targets.apply(init)

    /** Parameters necessary for calculating/parsing output quantites. */
    class Targets: ParametersInterface {
        /** Number of points in the energy grid that is used to calculate the (L)DOS. */
        var target = Target.LDOS

        /** Gridsize of the LDOS. */
        var ldosGridsize = 0

        /** Gridspacing of the energy grid the (L)DOS is evaluated on [eV]. */
        var ldosGridspacingEv = 0f

        /** Lowest energy value on the (L)DOS energy grid [eV]. */
        var ldosGridoffsetEv = 0f
        var restrictTargets = RestrictTargets.zeroOutNegative

        /** Path at which pseudopotentials are located (for TEM). */
        lateinit var pseudopotentialPath: File // = None

        /** Parameters for calculating the radial distribution function(RDF).
         *  The RDF can directly be calculated via a function call, but if it is calculated e.g. during a MD or MC run,
         *  these parameters will control how.
         *      numberOfBins: [Int]
         *          Number of bins used to create the histogram.
         *      rMax: [Radius]
         *          Radius up to which to calculate the RDF. None by default; this is the suggested behavior, as MALA
         *          will then on its own calculate the maximum radius up until which the calculation of the RDF is
         *          indisputably physically meaningful. Larger radii may be specified, e.g. for a Fourier transformation
         *          to calculate the static structure factor. */
        var rdfParameters = Function(500, Radius.mic)

        /** Parameters for calculating the three particle correlation function (TPCF).
         *  The TPCF can directly be calculated via a function call, but if it is calculated e.g. during a MD or MC run,
         *  these parameters will control how. The following keywords are recognized:
         *      numberOfBins: [Int]
         *          Number of bins used to create the histogram.
         *      rMax: [Radius]
         *          Radius up to which to calculate the TPCF. If None, MALA will determine the maximum radius for which
         *          the TPCF is indisputably defined. Be advised - this may come at increased computational cost. */
        val tpcfParameters = Function(20, Radius.mic)

        /** Parameters for calculating the static structure factor (SSF).
         *  The SSF can directly be calculated via a function call, but if it is calculated e.g. during a MD or MC run,
         *  these parameters will control how. The following keywords are recognized:
         *      numberOfBins: [Int]
         *          Number of bins used to create the histogram.
         *      kMax: [Float]
         *          Maximum wave vector up to which to calculate the SSF. */
        val ssfParameters = SsfFunction(100, 12f)
    }

    enum class Target { LDOS, DOS, Density }
    enum class RestrictTargets { zeroOutNegative, absoluteValues }

    data class Function(val numberOfBins: Int, val rMax: Radius)
    data class SsfFunction(val numberOfBins: Int, val rMax: Float)
    enum class Radius { mic, `2mic` }

    fun descriptors(init: Descriptors.() -> Unit) = descriptors.apply(init)

    val dataHandler by lazy { DataHandler(this) }

    open class Base {
        val configuration = Configuration()
        class Configuration(var gpu: Boolean = false,
                            var horovod: Boolean = false,
                            var mpi: Boolean = false,
                            var device: String = "cpu",
                            var openpmdConfiguration: Map<String, Any> = emptyMap(),
                            var openpmdGranularity: Int = 1,
                            var lammps: Boolean = true)
    }
    /** Parameters necessary for calculating/parsing input descriptors. */
    class Descriptors: Base() {

        /** Type of descriptors that is used to represent the atomic fingerprint. */
        var type = Type.Bispectrum

        // These affect all descriptors, at least as long all descriptors use LAMMPS (which they currently do).

        /** Bispectrum calculation: LAMMPS input file that is used to calculate the Bispectrum descriptors. If this
         *  string is empty, the standard LAMMPS input file found in this repository will be used (recommended). */
        lateinit var lammpsComputeFile: File// = ""

        /** Legacy option. If True, it is assumed that the first three entries of the descriptor vector are the xyz
         *  coordinates, and they are cut from the descriptor vector. If False, no such cutting is performed. */
        var descriptorsContainXyz = true

        // TODO: I would rather handle the parallelization info automatically  and more under the hood. At this stage of
        // the project this would probably be overkill and hard to do, since there are many moving parts, so for now
        // let's keep this here, but in the future, this should be adressed.
        var useZSplitting = true
        var useYSplitting = 0

        // Everything pertaining to the bispectrum descriptors.

        /** Bispectrum calculation: 2*jmax-parameter used for calculation of SNAP descriptors. Default value for jmax is
         *  5, so default value for twojmax is 10. */
        var bispectrumTwojmax = 10

        var bispectrumCutoff: Float = 4.67637f
            set(value) {
                field = value
                atomicDensityCutoff = value
            }
        var bispectrumSwitchflag = 1

        // Everything pertaining to the atomic density.
        // Separate cutoff given here because bispectrum descriptors and atomic density may be used at the same time,
        // if e.g. bispectrum descriptors are used for a full inference, which then uses the atomic density for
        // the calculation of the Ewald sum.
        var useAtomicDensityEnergyFormula = false

        /** Sigma used for the calculation of the Gaussian descriptors. */
        var atomicDensitySigma: Float? = null
        var atomicDensityCutoff: Float? = null

        // Everything concerning the minterpy descriptors.
        var minterpyPointList = listOf<String>()
        var minterpyCutoffCubeSize = 0f
        var minterpyPolynomialDegree = 4
        var minterpyLpNorm = 2

        enum class Type { Bispectrum, AtomicDensity, MinterpyDescriptors }
    }
}