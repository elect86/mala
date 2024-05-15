package casus.mala.common

import casus.mala.descriptors.AtomicDensity
import casus.mala.descriptors.Bispectrum

/** Base class for physical data.
 *
 *  Implements general framework to read and write such data to and from files. */
abstract class PhysicalData(val paramsInterface: ParametersInterface) {

    val gridDimensions = IntArray(3)

    /** Get a string that describes the data (for e.g. metadata). */
    abstract val dataName: String
}