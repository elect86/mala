package ase

import kotlin.math.pow
import kotlin.math.sqrt

class PhysicalQuantities(
    /** speed of light, m/s */
    val c: Double = 299792458.0,
    /** permeability of vacuum */
    val mu0: Double,
    /** gravitational constant */
    val Grav: Double,
    /** Planck constant, J s */
    val hplanck: Double,
    /** elementary charge */
    val e: Double,
    /** electron mass */
    val me: Double,
    /** proton mass */
    val mp: Double,
    /** Avogadro number */
    val Nav: Double,
    /** Boltzmann constant, J/K */
    val k: Double,
    /** atomic mass unit, kg */
    val amu: Double)

// this is the hard-coded CODATA values
// all other units are dynamically derived from these values upon import of the module
val CODATA = mapOf(
    // the "original" CODATA version ase used ever since
    // Constants from Konrad Hinsen's PhysicalQuantities module (1986 CODATA)
    // Add the constant pi used to define the mu0 and hbar here for reference as well
    1986 to PhysicalQuantities(mu0 = 4e-7 * Math.PI,
                               Grav = 6.67259e-11,
                               hplanck = 6.6260755e-34,
                               e = 1.60217733e-19,
                               me = 9.1093897e-31,
                               mp = 1.6726231e-27,
                               Nav = 6.0221367e23,
                               k = 1.380658e-23,
                               amu = 1.6605402e-27),

    //    /# CODATA 1998 taken from
    //# https://doi.org/10.1103/RevModPhys.72.351
    //'1998': {'_c': 299792458.,
    //    '_mu0': 4.0e-7 * pi,
    //    '_Grav': 6.673e-11,
    //    '_hplanck': 6.62606876e-34,
    //    '_e': 1.602176462e-19,
    //    '_me': 9.10938188e-31,
    //    '_mp': 1.67262158e-27,
    //    '_Nav': 6.02214199e23,
    //    '_k': 1.3806503e-23,
    //    '_amu': 1.66053873e-27},
    //
    //# CODATA 2002 taken from
    //# https://doi.org/10.1103/RevModPhys.77.1
    //'2002': {'_c': 299792458.,
    //    '_mu0': 4.0e-7 * pi,
    //    '_Grav': 6.6742e-11,
    //    '_hplanck': 6.6260693e-34,
    //    '_e': 1.60217653e-19,
    //    '_me': 9.1093826e-31,
    //    '_mp': 1.67262171e-27,
    //    '_Nav': 6.0221415e23,
    //    '_k': 1.3806505e-23,
    //    '_amu': 1.66053886e-27},

    // CODATA 2006 taken from
    // https://doi.org/10.1103/RevModPhys.80.633
    2006 to PhysicalQuantities(mu0 = 4.0e-7 * Math.PI.toFloat(),
                               Grav = 6.67428e-11,
                               hplanck = 6.62606896e-34,
                               e = 1.602176487e-19,
                               me = 9.10938215e-31,
                               mp = 1.672621637e-27,
                               Nav = 6.02214179e23,
                               k = 1.3806504e-23,
                               amu = 1.660538782e-27),

    //    # CODATA 2010 taken from
    //# https://doi.org/10.1103/RevModPhys.84.1527
    //'2010': {'_c': 299792458.,
    //    '_mu0': 4.0e-7 * pi,
    //    '_Grav': 6.67384e-11,
    //    '_hplanck': 6.62606957e-34,
    //    '_e': 1.602176565e-19,
    //    '_me': 9.10938291e-31,
    //    '_mp': 1.672621777e-27,
    //    '_Nav': 6.02214129e23,
    //    '_k': 1.3806488e-23,
    //    '_amu': 1.660538921e-27},
    //
    //# CODATA 2014 taken from
    //# http://arxiv.org/pdf/1507.07956.pdf
    //'2014': {'_c': 299792458.,
    //    '_mu0': 4.0e-7 * pi,
    //    '_Grav': 6.67408e-11,
    //    '_hplanck': 6.626070040e-34,
    //    '_e': 1.6021766208e-19,
    //    '_me': 9.10938356e-31,
    //    '_mp': 1.672621898e-27,
    //    '_Nav': 6.022140857e23,
    //    '_k': 1.38064852e-23,
    //    '_amu': 1.660539040e-27},
    //
    //# CODATA 2018 taken from
    //# https://physics.nist.gov/cuu/Constants/index.html
    //'2018': {'_c': 299792458.,            # Exact
    //    '_mu0': 4.0e-7 * pi,         # Exact
    //    '_Grav': 6.67430e-11,        # +/- 0.000_15e-11
    //    '_hplanck': 6.62607015e-34,  # Exact
    //    '_e': 1.602176634e-19,       # Exact
    //    '_me': 9.1093837015e-31,     # +/- 0.000_000_0028e-31
    //    '_mp': 1.67262192369e-27,    # +/- 0.000_000_000_51e-27
    //    '_Nav': 6.02214076e23,       # Exact
    //    '_k': 1.380649e-23,          # Exact
    //    '_amu': 1.66053906660e-27},  # +/- 0.000_000_000_50e-27
                  )
// [JVM] we need to cast some values (~pow) to double to avoid approximation to 0
class Units(version: Int) {
    val u = CODATA[version] ?: throw NotImplementedError("CODATA version \"$version\" not implemented")

    // derived from the CODATA values
    val eps0: Double = 1 / u.mu0 / u.c.pow(2) // permittivity of vacuum
    val hbar: Double = u.hplanck / (2 * Math.PI) // Planck constant / 2pi, J s

    val Ang: Double = 1.0
    val Angstrom: Double = 1.0
    val nm: Double = 10.0
    val Bohr: Double = 4e10 * Math.PI * eps0 * hbar.pow(2) / u.me / u.e.pow(2) // Bohr radius

    val eV: Double = 1.0
    val Hartree: Double = u.me * u.e.pow(3) / 16 / Math.PI.pow(2) / eps0.pow(2) / hbar.pow(2)
    val kJ: Double = 1000.0 / u.e
    val kcal: Double = 4.184 * kJ
    val mol: Double = u.Nav
    val Rydberg: Double = 0.5 * Hartree
    val Ry: Double = Rydberg
    val Ha: Double = Hartree

    val second: Double = 1e10 * sqrt(u.e / u.amu)
    val fs: Double = 1e-15 * second

    val kB: Double = u.k / u.e // Boltzmann constant, eV/K

    val Pascal: Double = (1 / u.e) / 1e30 // J/m^3
    val GPa: Double = 1e9 * Pascal
    val bar: Double = 1e5 * Pascal

    val Debye: Double = 1.0 / 1e11 / u.e / u.c
    val alpha: Double = u.e.pow(2) / (4 * Math.PI * eps0) / hbar / u.c // fine structure constant
    val invcm: Double = 100 * u.c * u.hplanck / u.e // cm^-1 energy unit

    // Derived atomic units that have no assigned name:
    // atomic unit of time, s:
    val aut: Double = hbar / (alpha.pow(2) * u.me * u.c.pow(2))

    // atomic unit of velocity, m/s:
    val auv: Double = u.e.pow(2) / hbar / (4 * Math.PI * eps0)

    // atomic unit of force, N:
    val auf: Double = alpha.pow(3) * u.me.pow(2) * u.c.pow(3) / hbar

    // atomic unit of pressure, Pa:
    val aup: Double = alpha.pow(5) * u.me.pow(4) * u.c.pow(5) / hbar.pow(3)

    val AUT: Double = second * aut

    // SI units
    val m: Double = 1e10 * Ang // metre
    val kg: Double = 1.0 / u.amu // kilogram
    val s: Double = second // second
    val A: Double = 1.0 / u.e / s // ampere

    // derived
    val J: Double = kJ / 1000.0 // Joule = kg * m**2 / s**2
    val C: Double = 1.0 / u.e // Coulomb = A * s
}
