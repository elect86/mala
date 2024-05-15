package casus.mala.dataHandling

import casus.mala.common.Parameters
import casus.mala.descriptors.Descriptor
import casus.mala.targets.Target

/**
 * Loads and scales data. Can only process numpy arrays at the moment.
 *
 * Data that is not in a numpy array can be converted using the DataConverter class.
 */
class DataHandler(parameters: Parameters,
                  targetCalculator: Target? = null,
                  descriptorCalculator: Descriptor? = null,
                  inputDataScaler: DataScaler? = null,
                  outputDataScaler: DataScaler? = null) : DataHandlerBase(parameters, targetCalculator, descriptorCalculator) {

    // Data will be scaled per user specification.
    val inputDataScaler = inputDataScaler ?: DataScaler(this.parameters.inputRescaling, useHorovod)

    val outputDataScaler = outputDataScaler ?: DataScaler(this.parameters.outputRescaling, useHorovod)

    // Actual data points in the different categories.
    var nrTrainingData = 0
    var nrTestData = 0
    var nrValidationData = 0

    // Number of snapshots in these categories.
    var nrTrainingSnapshots = 0
    var nrTestSnapshots = 0
    var nrValidationSnapshots = 0

    // Arrays and data sets containing the actual data.
//    self.training_data_inputs = torch.empty(0)
//    self.validation_data_inputs = torch.empty(0)
//    self.test_data_inputs = torch.empty(0)
//    self.training_data_outputs = torch.empty(0)
//    self.validation_data_outputs = torch.empty(0)
//    self.test_data_outputs = torch.empty(0)
    self.training_data_sets = []
    self.validation_data_sets = []
    self.test_data_sets = []

    # Needed for the fast tensor data sets.
    self.mini_batch_size = parameters.running.mini_batch_size
    if clear_data:
    self.clear_data()
}