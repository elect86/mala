package casus.mala.dataHandling

import casus.mala.common.Parameters
import casus.mala.verbosity

/** Scales input and output data.
 *
 *  Sort of emulates the functionality of the scikit-learn library, but by
 *  implementing the class by ourselves we have more freedom. */
class DataScaler(val type: Parameters.Scaling,
                 /** If True, the DataScaler will use horovod to check that data is
                  *  only saved on the root process in parallel execution. */
                 val useHorovod: Boolean = false) {

    var scaleStandard = "standard" in type.name
    var scaleNormal = "normal" in type.name
    var featureWise = "feature wise" in type.name
    var canTransform = !scaleStandard && !scaleNormal
    init {
        if (canTransform)
            verbosity("No data rescaling will be performed.", 1)
        if (scaleStandard && scaleNormal)
            error("Invalid input data rescaling.")
    }
//    self.__parse_typestring()
//
//    self.means = torch.empty(0)
//    self.stds = torch.empty(0)
//    self.maxs = torch.empty(0)
//    self.mins = torch.empty(0)
//    self.total_mean = torch.tensor(0)
//    self.total_std = torch.tensor(0)
//    self.total_max = torch.tensor(float("-inf"))
//    self.total_min = torch.tensor(float("inf"))
//
//    self.total_data_count = 0
}