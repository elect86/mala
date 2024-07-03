package casus.mala.descriptors

import casus.mala.common.Parameters
import java.io.File

class AtomicDensity(parameters: Parameters) : Descriptor(parameters) {
    override fun convertUnits(array: Number, inUnits: String?): Number {
        TODO("Not yet implemented")
    }

    override fun calculate(outDir: File, kwargs: Map<String, Any?>): Number {
        TODO("Not yet implemented")
    }

    override val dataName = "AtomicDensity"
}
