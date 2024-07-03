package casus.mala.dataHandling

import ai.djl.ndarray.NDArray
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

    lateinit var means: NDArray
    lateinit var stds: NDArray
    //    self.maxs = torch.empty(0)
    //    self.mins = torch.empty(0)
    //    self.total_mean = torch.tensor(0)
    //    self.total_std = torch.tensor(0)
    //    self.total_max = torch.tensor(float("-inf"))
    //    self.total_min = torch.tensor(float("inf"))

    var totalDataCount = 0

    //    def __parse_typestring(self):
    //    """Parse the typestring to class attributes."""
    //    self.scale_standard = False
    //    self.scale_normal = False
    //    self.feature_wise = False
    //
    //    if "standard" in self.typestring:
    //    self.scale_standard = True
    //    if "normal" in self.typestring:
    //    self.scale_normal = True
    //    if "feature-wise" in self.typestring:
    //    self.feature_wise = True
    //    if self.scale_standard is False and self.scale_normal is False:
    //    printout("No data rescaling will be performed.", min_verbosity=1)
    //    self.cantransform = True
    //    return
    //    if self.scale_standard is True and self.scale_normal is True:
    //    raise Exception("Invalid input data rescaling.")

    /** Start the incremental calculation of scaling parameters.
     *
     *  This is necessary for lazy loading. */
    fun startIncrementalFitting() {
        totalDataCount = 0
    }

    /** Indicate that all data has been added to the incremental calculation.
     *
     *  This is necessary for lazy loading. */
    fun finishIncrementalFitting() {
        canTransform = true
    }

    /** Compute the quantities necessary for scaling.
     *
     *  Parameters
     *  ----------
     *  @param unscaled: Data that on which the scaling will be calculated. */
    infix fun fit(unscaled: NDArray) {
        if (!scaleStandard && !scaleNormal)
            return
        // TODO
        //        with torch.no_grad():
        if (featureWise) {

            //#########################
            // Feature - wise - scaling
            //#########################

            if (scaleStandard) {
                means = unscaled.mean(intArrayOf(0), true)
                //                stds = unscaled.std(unscaled, 0, keepdim = True)
            }
            if (scaleNormal) {
                TODO()
                //                maxs = unscaled.max(intArrayOf(0), true).values
                //                self.mins = torch.min(unscaled, 0, keepdim = True).values
            }
        } else {
            TODO()
            //            ##########################
            //            # Total scaling
            //            ##########################
            //
            //            if self.scale_standard:
            //            self.total_mean = torch.mean(unscaled)
            //            self.total_std = torch.std(unscaled)
            //
            //            if self.scale_normal:
            //            self.total_max = torch.max(unscaled)
            //            self.total_min = torch.min(unscaled)
        }
        canTransform = true
    }
}