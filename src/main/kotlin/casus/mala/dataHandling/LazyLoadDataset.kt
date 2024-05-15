package casus.mala.dataHandling

import casus.mala.descriptors.Descriptor
import casus.mala.targets.Target

class LazyLoadDataset(inputDimension: Int,
                      outputDimension: Int,
                      inputDataScaler: DataScaler,
                      outputDataScaler: DataScaler,
                      descriptorCalculator: Descriptor,
                      targetCalculator: Target,
                      useHorovod: Boolean,
                      inputRequiresGrad: Boolean = false)/*: Dataset*/ {


}