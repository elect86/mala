package ase

import ase.io.div
import ase.io.dot

/**
 * Convert k-points between scaled and cartesian coordinates.
 *
 * Given the atomic unit cell, and either the scaled or cartesian k-point
 * coordinates, the other is determined.
 *
 * The k-point arrays can be either a single point, or a list of points,
 * i.e. the dimension k can be empty or multidimensional.
 */
fun kpointConvert(cellCv: Cell, skptsKc: FloatArray? = null, ckptsKv: FloatArray? = null) = when {
    ckptsKv == null -> {
        TODO()
        //            val icellCv = 2 * Math.PI.toFloat() * np.linalg.pinv(cellCv).T
        //            return np.dot(skptsKc, icellCv)
    }
    skptsKc == null -> (ckptsKv dot cellCv.array) / (2 * Math.PI.toFloat())
    else -> error("Either scaled or cartesian coordinates must be given.")
}