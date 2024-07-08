package casus.mala.descriptors

import casus.mala.common.Parameters
import java.io.File
import java.nio.file.Path

class AtomicDensity(parameters: Parameters) : Descriptor(parameters) {
    override fun convertUnits(array: Number, inUnits: String?): Number {
        TODO("Not yet implemented")
    }

    override fun calculate(outDir: Path, kwargs: Map<String, Any>): Pair<Array<Array<Array<FloatArray>>>, Int> {
        TODO("Not yet implemented")
    }

    override val dataName = "AtomicDensity"
}
