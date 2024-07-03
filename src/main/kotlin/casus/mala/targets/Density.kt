package casus.mala.targets

import ai.djl.ndarray.NDArray
import casus.mala.common.ParametersInterface

class Density(paramsInterface: ParametersInterface): Target(paramsInterface) {

    override val dataName = "Density"
    override val feaureSize = 1

    override val featureMask: Int
        get() = TODO("Not yet implemented")

    override fun processLoadedArray(array: NDArray, units: String?) {
        TODO("Not yet implemented")
    }
}