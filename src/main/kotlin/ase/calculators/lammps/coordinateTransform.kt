package ase.calculators.lammps

import ase.*
import ase.io.div
import ase.io.dot
import kotlin.math.pow
import kotlin.math.sqrt
import org.ejml.simple.SimpleMatrix

/**
 * Calculate box parameters
 *
 *     https://docs.lammps.org/Howto_triclinic.html
 */
fun calcBoxParameters(cell: Array<FloatArray>): FloatArray {
    val ax = sqrt(cell[0] dot cell[0])
    val bx = cell[0] dot cell[1] / ax
    val by = sqrt(cell[1] dot cell[1] - bx.pow(2))
    val cx = cell[0] dot cell[2] / ax
    val cy = (cell[1] dot cell[2] - bx * cx) / by
    val cz = sqrt(cell[2] dot cell[2] - cx.pow(2) - cy.pow(2))
    return floatArrayOf(ax, by, cz, bx, cx, cy)
}


/**
 * Calculate rotated cell in LAMMPS coordinates
 *
 *     Parameters
 *     ----------
 *     cell : np.ndarray
 *         Cell to be rotated.
 *
 *     Returns
 *     -------
 *     rotated_cell : np.ndarray
 *         Rotated cell represented by a lower triangular matrix.
 */
fun calcRotatedCell(cell: Array<FloatArray>): Array<FloatArray> {
    val (ax, by, cz, bx, cx, cy) = calcBoxParameters(cell)
    return arrayOf(floatArrayOf(ax, 0f, 0f), floatArrayOf(bx, by, 0f), floatArrayOf(cx, cy, cz))
}

/**
 * The representation of the unit cell in LAMMPS
 *
 *     The main purpose of the prism-object is to create suitable string
 *     representations of prism limits and atom positions within the prism.
 *
 *     Parameters
 *     ----------
 *     cell : np.ndarray
 *         Cell in ASE coordinate system.
 *     pbc : one or three bool
 *         Periodic boundary conditions flags.
 *     reduce_cell : bool
 *         If True, the LAMMPS cell is reduced for short lattice basis vectors.
 *         The atomic positions are always wraped into the reduced cell,
 *         regardress of `wrap` in `vector_to_lammps` and `vector_to_ase`.
 *     tolerance : float
 *         Precision for skewness test.
 *
 *     Methods
 *     -------
 *     vector_to_lammps
 *         Rotate vectors from ASE to LAMMPS coordinates.
 *         Positions can be further wrapped into the LAMMPS cell by `wrap=True`.
 *
 *     vector_to_ase
 *         Rotate vectors from LAMMPS to ASE coordinates.
 *         Positions can be further wrapped into the LAMMPS cell by `wrap=True`.
 *
 *     Notes
 *     -----
 *     LAMMPS prefers triangular matrixes without a strong tilt.
 *     Therefore the 'Prism'-object contains three coordinate systems:
 *
 *     - ase_cell (the simulated system in the ASE coordination system)
 *     - lammps_tilt (ase-cell rotated to be an lower triangular matrix)
 *     - lammps_cell (same volume as tilted cell, but reduce edge length)
 *
 *     The translation between 'ase_cell' and 'lammps_tilt' is done with a
 *     rotation matrix 'rot_mat'. (Mathematically this is a QR decomposition.)
 *
 *     The transformation between 'lammps_tilt' and 'lammps_cell' is done by
 *     changing the off-diagonal elements.
 *
 *     Depending on the option `reduce`, vectors in ASE coordinates are
 *     transformed either `lammps_tilt` or `lammps_cell`.
 *
 *     The vector conversion can fail as depending on the simulation run LAMMPS
 *     might have changed the simulation box significantly. This is for example a
 *     problem with hexagonal cells. LAMMPS might also wrap atoms across periodic
 *     boundaries, which can lead to problems for example NEB calculations.
 */
class Prism(val aseCell: Cell,
            pbc: Boolean = true,
            val isReduced: Boolean = false,
            val tolerance: Float = 1e-8f) {

    //    rot_mat * lammps_tilt^T = ase_cell^T
    // => lammps_tilt * rot_mat^T = ase_cell
    // => lammps_tilt             = ase_cell * rot_mat
    // LAMMPS requires positive diagonal elements of the triangular matrix.
    // The diagonals of `lammps_tilt` are always positive by construction.
    val lammpsTilt = calcRotatedCell(aseCell.array)
    val rotMat = SimpleMatrix(lammpsTilt).solve(SimpleMatrix(cell)).transpose().array
//    self.pbc = np.zeros(3, bool) + pbc
//    self.lammps_cell = calc_reduced_cell(self.lammps_tilt, self.pbc)

    val cell: Array<FloatArray>
        get() = if (isReduced) TODO()/* self.lammps_cell*/ else lammpsTilt

    /**
     * Return box parameters of the rotated cell in LAMMPS coordinates
     *
     *         Returns
     *         -------
     *         np.ndarray
     *             xhi - xlo, yhi - ylo, zhi - zlo, xy, xz, yz
     */
    val lammpsPrism: FloatArray
        get() = floatArrayOf(cell[0][0], cell[1][1], cell[2][2], cell[1][0], cell[2][0], cell[2][1])

    /**
     * Rotate vectors from ASE to LAMMPS coordinates
     *
     *         Parameters
     *         ----------
     *         vec : np.ndarray
     *             Vectors in ASE coordinates to be rotated into LAMMPS coordinates
     *         wrap : bool
     *             If True, the vectors are wrapped into the cell
     *
     *         Returns
     *         -------
     *         np.array
     *             Vectors in LAMMPS coordinates
     */
    fun vectorToLammps(vec: Array<FloatArray>, wrap: Boolean = false): Array<FloatArray> =
        // !TODO: right eps-limit
        // lammps might not like atoms outside the cell
        if (wrap /*or self.is_reduced*/) TODO()
//        return wrap_positions(
//            vec @ self.rot_mat,
//            cell = self.cell,
//            pbc = self.pbc,
//            eps = 1e-18,
//                             )
        else vec * rotMat

    /**
     * Test if the lammps cell is skewed, i.e., monoclinic or triclinic.
     *
     *         Returns
     *         -------
     *         bool
     *             True if the lammps cell is skewed.
     */
    val isSkewed: Boolean
        get() {
            val cellSq = cell.pow(2)
            val onDiag = cellSq.diag().sum()
            val offDiag = cellSq.tril(-1).sum()
            return offDiag / onDiag > tolerance
        }
}