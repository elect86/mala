package ase

import kotlin.math.pow
import kotlin.math.sqrt

class PhysicalQuantities(
    /** speed of light, m/s */
    val c: Float = 2.9979245E8f,
    /** permeability of vacuum */
    val mu0: Float,
    /** gravitational constant */
    val Grav: Float,
    /** Planck constant, J s */
    val hplanck: Float,
    /** elementary charge */
    val e: Float,
    /** electron mass */
    val me: Float,
    /** proton mass */
    val mp: Float,
    /** Avogadro number */
    val Nav: Float,
    /** Boltzmann constant, J/K */
    val k: Float,
    /** atomic mass unit, kg */
    val amu: Float)

// this is the hard-coded CODATA values
// all other units are dynamically derived from these values upon import of the module
val CODATA = mapOf(
    // the "original" CODATA version ase used ever since
    // Constants from Konrad Hinsen's PhysicalQuantities module (1986 CODATA)
    // Add the constant pi used to define the mu0 and hbar here for reference as well
    1986 to PhysicalQuantities(mu0 = 4e-7f * Math.PI.toFloat(),
                               Grav = 6.67259e-11f,
                               hplanck = 6.6260757E-34f,
                               e = 1.6021774E-19f,
                               me = 9.10939E-31f,
                               mp = 1.6726232E-27f,
                               Nav = 6.0221367e23f,
                               k = 1.380658e-23f,
                               amu = 1.6605403E-27f),

    // CODATA 1998 taken from
    // https://doi.org/10.1103/RevModPhys.72.351
    1998 to PhysicalQuantities(mu0 = 4.0e-7f * Math.PI.toFloat(),
                               Grav = 6.673e-11f,
                               hplanck = 6.626069E-34f,
                               e = 1.6021765E-19f,
                               me = 9.109382E-31f,
                               mp = 1.6726216E-27f,
                               Nav = 6.022142E23f,
                               k = 1.3806503e-23f,
                               amu = 1.6605387E-27f),

    // CODATA 2002 taken from
    // https://doi.org/10.1103/RevModPhys.77.1
    2002 to PhysicalQuantities(mu0 = 4.0e-7f * Math.PI.toFloat(),
                               Grav = 6.6742e-11f,
                               hplanck = 6.6260693e-34f,
                               e = 1.6021765E-19f,
                               me = 9.1093825E-31f,
                               mp = 1.6726216E-27f,
                               Nav = 6.0221414E23f,
                               k = 1.3806505e-23f,
                               amu = 1.66053886e-27f),

    // CODATA 2006 taken from
    // https://doi.org/10.1103/RevModPhys.80.633
    2006 to PhysicalQuantities(mu0 = 4.0e-7f * Math.PI.toFloat(),
                               Grav = 6.67428e-11f,
                               hplanck = 6.626069E-34f,
                               e = 1.6021765E-19f,
                               me = 9.1093825E-31f,
                               mp = 1.6726216E-27f,
                               Nav = 6.0221417E23f,
                               k = 1.3806505E-23f,
                               amu = 1.6605387E-27f),

    // CODATA 2010 taken from
    // https://doi.org/10.1103/RevModPhys.84.1527
    2010 to PhysicalQuantities(mu0 = 4.0e-7f * Math.PI.toFloat(),
                               Grav = 6.67384e-11f,
                               hplanck = 6.6260697E-34f,
                               e = 1.6021766E-19f,
                               me = 9.1093825E-31f,
                               mp = 1.6726218E-27f,
                               Nav = 6.0221414E23f,
                               k = 1.3806487E-23f,
                               amu = 1.660539E-27f),

    // CODATA 2014 taken from
    // http://arxiv.org/pdf/1507.07956.pdf
    2014 to PhysicalQuantities(mu0 = 4.0e-7f * Math.PI.toFloat(),
                               Grav = 6.67408e-11f,
                               hplanck = 6.62607E-34f,
                               e = 1.6021766E-19f,
                               me = 9.109383E-31f,
                               mp = 1.6726218E-27f,
                               Nav = 6.022141E23f,
                               k = 1.3806486E-23f,
                               amu = 1.6605391E-27f),

    // CODATA 2018 taken from
    // https://physics.nist.gov/cuu/Constants/index.html
    2018 to PhysicalQuantities(mu0 = 4.0e-7f * Math.PI.toFloat(),   // Exact
                               Grav = 6.67430e-11f,                 // +/- 0.000_15e-11
                               hplanck = 6.62607E-34f,              // Exact
                               e = 1.6021766E-19f,                  // Exact
                               me = 9.109383E-31f,                  // +/- 0.000_000_0028e-31
                               mp = 1.672622E-27f,                  // +/- 0.000_000_000_51e-27
                               Nav = 6.0221406E23f,                 // Exact
                               k = 1.380649e-23f,                   // Exact
                               amu = 1.6605391E-27f))               // +/- 0.000_000_000_50e-27

// [JVM] we need to cast some values (~pow) to double to avoid approximation to 0
class Units(version: Int) {
    val u = CODATA[version] ?: throw NotImplementedError("CODATA version \"$version\" not implemented")

    // derived from the CODATA values
    val eps0: Float = 1 / u.mu0 / u.c.pow(2) // permittivity of vacuum
    val hbar: Float = u.hplanck / (2 * Math.PI.toFloat()) // Planck constant / 2pi, J s

    val Ang: Float = 1f
    val Angstrom: Float = 1f
    val nm: Float = 10f
    val Bohr: Float = (4e10 * Math.PI * eps0 * hbar.toDouble().pow(2) / u.me / u.e.pow(2)).toFloat() // Bohr radius

    val eV: Float = 1f
    val Hartree: Float = (u.me * u.e.toDouble().pow(3) / 16 / Math.PI.pow(2) / eps0.toDouble().pow(2) / hbar.toDouble().pow(2)).toFloat()
    val kJ: Float = 1000f / u.e
    val kcal: Float = 4.184f * kJ
    val mol: Float = u.Nav
    val Rydberg: Float = 0.5f * Hartree
    val Ry: Float = Rydberg
    val Ha: Float = Hartree

    val second: Float = 1e10f * sqrt(u.e / u.amu)
    val fs: Float = 1e-15f * second

    val kB: Float = u.k / u.e // Boltzmann constant, eV/K

    val Pascal: Float = (1 / u.e) / 1e30f // J/m^3
    val GPa: Float = 1e9f * Pascal
    val bar: Float = 1e5f * Pascal

    val Debye: Float = 1f / 1e11f / u.e / u.c
    val alpha: Float = u.e.pow(2) / (4 * Math.PI.toFloat() * eps0) / hbar / u.c // fine structure constant
    val invcm: Float = 100 * u.c * u.hplanck / u.e // cm^-1 energy unit

    // Derived atomic units that have no assigned name:
    // atomic unit of time, s:
    val aut: Float = hbar / (alpha.pow(2) * u.me * u.c.pow(2))
    // atomic unit of velocity, m/s:
    val auv: Float = u.e.pow(2) / hbar / (4 * Math.PI.toFloat() * eps0)
    // atomic unit of force, N:
    val auf: Float = (alpha.pow(3) * u.me.toDouble().pow(2) * u.c.pow(3) / hbar).toFloat()
    // atomic unit of pressure, Pa:
    val aup: Float = (alpha.pow(5) * u.me.toDouble().pow(4) * u.c.pow(5) / hbar.toDouble().pow(3)).toFloat()

    val AUT: Float = second * aut

    // SI units
    val m: Float = 1e10f * Ang // metre
    val kg: Float = 1f / u.amu // kilogram
    val s: Float = second // second
    val A: Float = 1f / u.e / s // ampere
    // derived
    val J: Float = kJ / 1000 // Joule = kg * m**2 / s**2
    val C: Float = 1f / u.e // Coulomb = A * s
}
