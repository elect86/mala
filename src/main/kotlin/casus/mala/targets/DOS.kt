package casus.mala.targets

import ai.djl.ndarray.NDArray
import casus.mala.common.ParametersInterface

class DOS(paramsInterface: ParametersInterface): Target(paramsInterface) {

    override val dataName = "DOS"
    override val feaureSize = parameters.ldosGridsize

    override val featureMask: Int
        get() = TODO("Not yet implemented")
    override fun processLoadedArray(array: NDArray, units: String?) {
        TODO("Not yet implemented")
    }
}