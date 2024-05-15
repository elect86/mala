package casus.mala.targets

import casus.mala.common.ParametersInterface

class DOS(paramsInterface: ParametersInterface): Target(paramsInterface) {

    override val dataName = "DOS"
    override val feaureSize = parameters.ldosGridsize
}