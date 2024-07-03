package casus.mala.common

import java.io.File

/** Parameters needed for network runs (train, test or inference).
 *
 *  Some of these parameters only apply to either the train or test or inference case. */
class ParametersRunning {
    /** Training type to be used */
    var training = Parameters.Training.SGD

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
    var learningRateScheduler = Parameters.RateScheduler.None

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

    var duringTrainingMetric = Parameters.TrainingMetric.ldos

    var afterBeforeTrainingMetric = Parameters.TrainingMetric.ldos

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