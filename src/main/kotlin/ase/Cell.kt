package ase

import ase.geometry.completeCell
import ase.geometry.isOrthorhombic
import org.ejml.simple.SimpleMatrix

/**
 * Parallel epipedal unit cell of up to three dimensions.
 *
 *     This object resembles a 3x3 array whose [i, j]-th element is the jth
 *     Cartesian coordinate of the ith unit vector.
 *
 *     Cells of less than three dimensions are represented by placeholder
 *     unit vectors that are zero.
 */
class Cell(val array: Array<FloatArray> = Array(3) { FloatArray(3) }) {

    init {
        check(array.size == 3 && array[0].size == 3 && array[1].size == 3 && array[2].size == 3)
    }

    /** Convert missing cell vectors into orthogonal unit vectors. */
    fun complete() = Cell(completeCell(array))

    /** Return a copy of this cell. */
    fun copy() = Cell(array.clone())
    operator fun get(i: Int): FloatArray = array[i]
    operator fun get(r: Int, c: Int): Float = array[r][c]

    /** Return whether this cell is represented by a diagonal matrix. */
    val orthorhombic: Boolean
        get() = isOrthorhombic

    /**
     * Calculate scaled positions from Cartesian positions.
     *
     * The scaled positions are the positions given in the basis
     * of the cell vectors.  For the purpose of defining the basis, cell
     * vectors that are zero will be replaced by unit vectors as per
     * :meth:`~ase.cell.Cell.complete`.
     */
    fun scaledPositions(positions: Array<FloatArray>): Array<FloatArray> {
        val a = complete().array.transpose()
        val b = positions.transpose()
        return SimpleMatrix(a).solve(SimpleMatrix(b)).transpose().array
    }

    /** Calculate Cartesian positions from scaled positions. */
    fun cartesianPositions(scaledPositions: Array<FloatArray>) = SimpleMatrix(scaledPositions).mult(SimpleMatrix(complete().array)).array

    companion object {
        /**
         * Create new cell from any parameters.
         *
         * If cell is three numbers, assume three lengths with right angles.
         *
         * If cell is six numbers, assume three lengths, then three angles.
         *
         * If cell is 3x3, assume three cell vectors.
         */
        infix fun from(cell: Any): Cell {
            TODO()
        }
    }
}