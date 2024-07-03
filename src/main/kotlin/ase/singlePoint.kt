package ase

import ase.io.SinglePointKPoint

/**
 * Special calculator for a single configuration.
 *
 * Used to remember the energy, force and stress for a given
 * configuration.  If the positions, atomic numbers, unit cell, or
 * boundary conditions are changed, then asking for
 * energy/forces/stress will raise an exception.
 */
open class SinglePointCalculator(atoms: Atoms,
                                 val energy: Float? = null, val freeEnergy: Float? = null,
                                 val forces: Array<FloatArray>? = null, val stress: FloatArray? = null,
                                 val magmoms: Any? = null, val dipole: Any? = null) : Calculator() {

    //name = 'unknown'
    val atoms = atoms.copy()
}


class SinglePointDFTCalculator(atoms: Atoms,
                               val efermi: Float? = null, val bzkpts: Any? = null,
                               val ibzkpts: Array<FloatArray>? = null, val bz2ibz: Any? = null,
                               val kpts: MutableList<SinglePointKPoint>? = null,
                               energy: Float? = null, freeEnergy: Float?,
                               forces: Array<FloatArray>? = null, stress: FloatArray? = null,
                               magmoms: Any? = null, dipole: Any? = null) : SinglePointCalculator(atoms,
                                                                                                  energy, freeEnergy,
                                                                                                  forces, stress,
                                                                                                  magmoms, dipole)