package casus.mala.common

import java.io.File

fun parameters(init: Parameters.() -> Unit) = Parameters().apply(init)

class Parameters: ParametersInterface {

    var comment = ""

    // Parameters subobjects.
    val network = Network()
    val descriptors = Descriptors()
    val targets = Targets()
    val data = Data()
    val running = Running()
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
    class Data {
        /** A list of all added snapshots. */
        val snapshotDirectoriesList = listOf<File>()

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

    fun network(init: Network.() -> Unit) = network.apply(init)

    /** Parameters necessary for constructing a neural network. */
    class Network {
        /** Type of the neural network that will be used. */
        var nn = Nn.`feed-forward`

        /**  A list of integers detailing the sizes of the layer of the neural network. Please note that the input layer
         *  is included therein. Default: [10,10,10] */
        val layerSizes = IntArray(3) { 10 }

        /** A list of strings detailing the activation functions to be used by the neural network. If the dimension of
         * [layerActivations] is smaller than the dimension of layer_sizes-1, than the first entry is used for all layers. */
        var layerActivations = listOf(Activation.Sigmoid)

        /** Loss function for the neural network */
        var lossFunction = LossFunction.mse

        /** for LSTM/Gru + Transformer
         *  Number of hidden layers to be used in lstm or gru or transformer nets */
        var numHiddenLayers = 1

        /** for LSTM/Gru
         *  If True hidden and cell state is assigned to zeros for LSTM Network. false will keep the hidden state active */
        var noHiddenState = false

        /** Sets lstm network size based on bidirectional or just one direction */
        var bidirection = false

        /** for transformer net
         *  Dropout rate for transformer net
         *  0.0 ≤ dropout ≤ 1.0 */
        var dropout = 0.1f

        /** Number of heads to be used in Multi head attention network. This should be a divisor of input dimension */
        var numHeads = 10
    }

    enum class Nn { `feed-forward`, transformer, lstm, gru }
    enum class Activation { Sigmoid, ReLU, LeakyReLU }
    enum class LossFunction {
        /** Mean squared error */
        mse
    }

    fun running(init: Running.() -> Unit) = running.apply(init)

    /** Parameters needed for network runs (train, test or inference).
     *
     *  Some of these parameters only apply to either the train or test or inference case. */
    class Running {
        /** Training type to be used */
        var training = Training.SGD

        /** Learning rate for chosen optimization algorithm. Default: 0.5. */
        var learningRate = 0.5f

        /** Maximum number of epochs to train for. Default: 100. */
        var maxNumberEpochs = 100
        var verbosity = true

        /** Size of the mini batch for the optimization algorihm. Default: 10. */
        var miniBatchSize = 10

        /** Weight decay for regularization. Always refers to L2 regularization. Default: 0f. */
        var weightDecay = 0f

        /** Number of epochs the validation accuracy is allowed to not improve by at least [earlyStoppingThreshold],
         *  before we terminate. If 0, no early stopping is performed. Default: 0. */
        var earlyStoppingEpochs = 0

        /** Minimum fractional reduction in validation loss required to avoid early stopping, e.g. a value of 0.05 means
         *  that validation loss must decrease by 5% within early_stopping_epochs epochs or the training will be stopped
         *  early. More explicitly, validation_loss < validation_loss_old * (1-early_stopping_threshold) or the patience
         *  counter goes up.
         *  Default: 0. Numbers bigger than 0 can make early stopping very aggressive,
         *  while numbers less than 0 make the trainer very forgiving of loss increase. */
        var earlyStoppingThreshold = 0

        /** Learning rate scheduler to be used. If not [RateScheduler.None], an instance of the corresponding pytorch
         *  class will be used to manage the learning rate schedule. */
        var learningRateScheduler = RateScheduler.None

        /** Decay rate to be used in the learning rate (if the chosen scheduler supports that). Default: 0.1 */
        var learningRateDecay = 0.1f

        /** Patience parameter used in the learning rate schedule (how long the validation loss has to plateau before
         *  the schedule takes effect). Default: 0. */
        var learningRatePatience = 0

        /** If True and horovod is used, horovod compression will be used for allreduce communication.
         *  This can improve performance. */
        var useCompression = false

        /** Number of workers to be used for data loading. */
        var numWorkers = 0

        /** If True, the training data will be shuffled in between epochs. If lazy loading is selected, then this
         *  shuffling will be done on a "by snapshot" basis. */
        var useShufflingForSamplers = true

        /** If not 0, checkpoint files will be saved after eac checkpoints_each_epoch epoch. */
        var checkpointsEachEpoch = 0

        /** Name used for the checkpoints. Using this, multiple runs can be performed in the same directory. */
        var checkpointName = "checkpoint_mala"

        /** If True then Tensorboard is activated for visualisation
        case 0: No tensorboard activated
        case 1: tensorboard activated with Loss and learning rate
        case 2; additonally weights and biases and gradient */
        var visualisation = 0

        /** Name of the folder that visualization files will be saved to. */
        lateinit var visualisationDir: File // = os.path.join(".", "mala_logging")

        /** If True, then upon creating visualization files, these will be saved in a subfolder of [visualisationDir]
         *  labelled with the starting date of the visualization, to avoid having to change input scripts often. */
        var visualisationDirAppendDate = true

        var duringTrainingMetric = TrainingMetric.ldos

        var afterBeforeTrainingMetric = TrainingMetric.ldos

        /** List holding the grid to be used for inference in the form of [x,y,z]. */
        var inferenceDataGrid = IntArray(3)

        /** If True, mixed precision computation (via AMP) will be used. */
        var useMixedPrecision = false

        var useGraphs = false

        /** Determines how often detailed performance info is printed during training
         *  (only has an effect if the verbosity is high enough). */
        var trainingReportFrequency = 1000

        /** List with two entries determining with which batch/iteration number the CUDA profiler will start and stop
         *  profiling. Please note that this option only holds significance if the nsys profiler is used. */
        var profilerRange: IntArray? = null  // [1000, 2000]
    }

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

    /** Parameters necessary for calculating/parsing input descriptors. */
    class Descriptors {

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

        var bispectrumCutoff = 4.67637f
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