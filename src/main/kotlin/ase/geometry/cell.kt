package ase.geometry

import ase.*

/**
 * Calculate complete cell with missing lattice vectors.
 *
 * Returns a new 3x3 ndarray.
 */
fun completeCell(cell: Array<FloatArray>): Array<FloatArray> {
    val missing = (!cell.any(axis = 1)).nonZero()

    check(missing.isEmpty())
    //        if len(missing) == 3:
    //        cell.flat[::4] = 1.0
    //        if len(missing) == 2:
    //        # Must decide two vectors :
    //        V, s, WT = np.linalg.svd(cell.T)
    //        sf = [s[0], 1, 1]
    //        cell = (V @ np.diag(sf) @ WT).T
    //        if np.sign(np.linalg.det(cell)) < 0:
    //        cell[missing[0]] = -cell[missing[0]]
    //        elif len (missing) == 1:
    //        i = missing[0]
    //        cell[i] = np.cross(cell[i - 2], cell[i - 1])
    //        cell[i] /= np.linalg.norm(cell[i])

    return cell
}

/** Check that cell only has stuff in the diagonal. */
val Cell.isOrthorhombic: Boolean
    get() = (array.flatNonZero() % 4).all { it == 0 }