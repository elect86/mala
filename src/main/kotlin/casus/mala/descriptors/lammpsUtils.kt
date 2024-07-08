package casus.mala.descriptors

import ase.io.prod
import lammps.Lammps
import lammps.library_h
import java.lang.foreign.MemorySegment
import java.lang.foreign.ValueLayout
import java.util.ArrayList

/**
 *     Add a dicitionary of LAMMPS arguments in a command line argument string.
 *
 *     Parameters
 *     ----------
 *     cmdargs : list
 *         Command line argument string. Will be mutated by this function.
 *
 *     argdict : dict
 *         Dictionary to be added to LAMMPS command line argument string.
 *
 *     Returns
 *     -------
 *     cmdargs : list
 *         New command line argument string.
 */
fun setCmdlineVars(cmdargs: List<String>, argDict: Map<String, Any>): ArrayList<String> =
    (cmdargs + argDict.flatMap { (k, v) -> listOf("-var", k, v.toString()) }) as ArrayList<String>

/**
 *     Convert a lammps compute to a numpy array.
 *
 *     Assumes the compute returns floating point numbers.
 *     Note that the result is a view into the original memory.
 *     If the result type is 0 (scalar) then conversion to numpy is
 *     skipped and a python float is returned.
 *
 *     Parameters
 *     ----------
 *     lmp : lammps.lammps
 *         The LAMMPS object from which data is supposed to be extracted.
 *
 *     name : string
 *         Name of the LAMMPS calculation.
 *
 *     compute_type
 *         Compute type of the LAMMPS calculation.
 *
 *     result_type
 *         Result type of the LAMMPS calculation.
 *
 *     array_shape
 *         Array shape of the LAMMPS calculation.
 *
 *     use_fp64 : bool
 *         If True, return the array with double precision. If False (default),
 *         the array will be processed with single precision. Only has an effect
 *         if result_type equals 2.
 */
fun extractComputeNp(lmp: Lammps, name: String, computeType: library_h.Style, resultType: library_h.Type,
                     arrayShape: IntArray? = null, useFp64: Boolean = false): Any? {

    val ptr = lmp.extractCompute(name, computeType, resultType)
    return when(resultType) {
        library_h.Type.scalar ->  ptr // No casting needed, lammps.py already works
        library_h.Type.array -> {
            val totalSize = arrayShape!!.prod()
            val doubles = ptr.get(ValueLayout.ADDRESS, 0).reinterpret(totalSize.toLong() * Double.SIZE_BYTES)
                    .toArray(ValueLayout.JAVA_DOUBLE)
            if (useFp64) doubles else FloatArray(doubles.size) { doubles[it].toFloat() }
        }
        library_h.Type.sizeRows, library_h.Type.sizeCols -> // ptr is an int
            ptr.get(ValueLayout.JAVA_INT, 0L)
        else -> error("unexpected result type $resultType")
    }
}