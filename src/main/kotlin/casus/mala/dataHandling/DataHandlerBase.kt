package casus.mala.dataHandling

import casus.mala.common.Parameters
import casus.mala.descriptors.Descriptor
import casus.mala.targets.Target

/** Base class for all data handling (loading, shuffling, etc.). */
abstract class DataHandlerBase(parameters: Parameters,
                               targetCalculator: Target? = null,
                               descriptorCalculator: Descriptor? = null) {

    val parameters: Parameters.Data = parameters.data
    val useHorovod = parameters.useHorovod

    // Calculators used to parse data from compatible files.
    val targetCalculator = targetCalculator ?: Target.from(parameters)
    val descriptorCalculator = descriptorCalculator ?: Descriptor.from(parameters)

    // Dimensionalities of data.
    var inputDimension = 0
    var outputDimension = 0
    var nrSnapshots = 0
}