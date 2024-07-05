package ase.calculators.lammps

import ase.io.units
import kotlin.math.pow


sealed interface UnitSet {

    val mass: Double
    val distance: Double
    val time: Double
    val energy: Double
    val velocity: Double
    val force: Double
    val pressure: Double
    val charge: Double

    data object ase : UnitSet {
        override val mass = 1.0 / units.kg
        override val distance = 1.0 / units.m
        override val time = 1.0 / units.second
        override val energy = 1.0 / units.J
        override val velocity = units.second / units.m
        override val force = units.m / units.J
        override val pressure = 1.0 / units.Pascal
        override val charge = 1.0 / units.C
    }

    data object real : UnitSet {
        override val mass = gramPerMoleSI
        override val distance = angstromSI
        override val time = femtosecondSI
        override val energy = kcalPerMoleSI
        override val velocity = angstromPerFemtosecondSI
        override val force = kcalPerMoleAngstromSI
        val torque = kcalPerMoleSI
        val temperature = kelvinSI
        override val pressure = atmosphereSI
        val dynamic_viscosity = poiseSI
        override val charge = eSI
        val dipole = electronAngstromSI
        val electricField = voltPerAngstromSI
        val density = gramSI / centimeterSI.pow(DIM)
    }

    data object metal : UnitSet {
        override val mass: Double = gramPerMoleSI
        override val distance: Double = angstromSI
        override val time: Double = picosecondSI
        override val energy: Double = evSI
        override val velocity: Double = angstromPerPicosecondSI
        override val force: Double = evPerAngstromSI
        val torque: Double = evSI
        val temperature: Double = kelvinSI
        override val pressure: Double = barSi
        val dynamicViscosity: Double = poiseSI
        override val charge: Double = eSI
        val dipole: Double = electronAngstromSI
        val electricField: Double = voltPerAngstromSI
        val density: Double = gramSI / centimeterSI.pow(DIM)
    }

    //    UNITSETS["si"] = dict(
    //    mass=u.kilogram_si,
    //    distance=u.meter_si,
    //    time=u.second_si,
    //    energy=u.joule_si,
    //    velocity=u.meter_per_second_si,
    //    force=u.newton_si,
    //    torque=u.joule_si,
    //    temperature=u.kelvin_si,
    //    pressure=u.pascal_si,
    //    dynamic_viscosity=u.pascal_si * u.second_si,
    //    charge=u.coulomb_si,
    //    dipole=u.coulomb_meter_si,
    //    electric_field=u.volt_per_meter_si,
    //    density=u.kilogram_si / u.meter_si ** DIM,
    //    )
    //
    //    UNITSETS["cgs"] = dict(
    //    mass=u.gram_si,
    //    distance=u.centimeter_si,
    //    time=u.second_si,
    //    energy=u.erg_si,
    //    velocity=u.centimeter_per_second_si,
    //    force=u.dyne_si,
    //    torque=u.dyne_centimeter_si,
    //    temperature=u.kelvin_si,
    //    pressure=u.dyne_per_centimetersq_si,  # or barye =u. 1.0e-6 bars
    //    dynamic_viscosity=u.poise_si,
    //    charge=u.statcoulomb_si,  # or esu (4.8032044e-10 is a proton)
    //    dipole=u.statcoulomb_centimeter_si,  # =u. 10^18 debye,
    //    electric_field=u.statvolt_per_centimeter_si,  # or dyne / esu
    //    density=u.gram_si / (u.centimeter_si ** DIM),
    //    )
    //
    //    UNITSETS["electron"] = dict(
    //    mass=u.amu_si,
    //    distance=u.bohr_si,
    //    time=u.femtosecond_si,
    //    energy=u.hartree_si,
    //    velocity=u.bohr_per_atu_si,
    //    force=u.hartree_per_bohr_si,
    //    temperature=u.kelvin_si,
    //    pressure=u.pascal_si,
    //    charge=u.e_si,  # multiple of electron charge (1.0 is a proton)
    //    dipole=u.debye_si,
    //    electric_field=u.volt_per_centimeter_si,
    //    )
    //
    //    UNITSETS["micro"] = dict(
    //    mass=u.picogram_si,
    //    distance=u.micrometer_si,
    //    time=u.microsecond_si,
    //    energy=u.picogram_micrometersq_per_microsecondsq_si,
    //    velocity=u.micrometer_per_microsecond_si,
    //    force=u.picogram_micrometer_per_microsecondsq_si,
    //    torque=u.picogram_micrometersq_per_microsecondsq_si,
    //    temperature=u.kelvin_si,
    //    pressure=u.picogram_per_micrometer_microsecondsq_si,
    //    dynamic_viscosity=u.picogram_per_micrometer_microsecond_si,
    //    charge=u.picocoulomb_si,  # (1.6021765e-7 is a proton),
    //    dipole=u.picocoulomb_micrometer_si,
    //    electric_field=u.volt_per_micrometer_si,
    //    density=u.picogram_si / (u.micrometer_si) ** DIM,
    //    )
    //
    //    UNITSETS["nano"] = dict(
    //    mass=u.attogram_si,
    //    distance=u.nanometer_si,
    //    time=u.nanosecond_si,
    //    energy=u.attogram_nanometersq_per_nanosecondsq_si,
    //    velocity=u.nanometer_per_nanosecond_si,
    //    force=u.attogram_nanometer_per_nanosecondsq_si,
    //    torque=u.attogram_nanometersq_per_nanosecondsq_si,
    //    temperature=u.kelvin_si,
    //    pressure=u.attogram_per_nanometer_nanosecondsq_si,
    //    dynamic_viscosity=u.attogram_per_nanometer_nanosecond_si,
    //    charge=u.e_si,  # multiple of electron charge (1.0 is a proton)
    //    dipole=u.electron_nanometer_si,
    //    electric_field=u.volt_si / u.nanometer_si,
    //    density=u.attogram_si / u.nanometer_si ** DIM,
    //    )

    companion object {
        // NOTE: We assume a three-dimensional simulation here!
        val DIM = 3.0
    }
}

/**
 * Convert units between LAMMPS and ASE.
 *
 *     :param value: converted value
 *     :param quantity: mass, distance, time, energy, velocity, force, torque,
 *     temperature, pressure, dynamic_viscosity, charge, dipole,
 *     electric_field or density
 *     :param fromunits: ASE, metal, real or other (see lammps docs).
 *     :param tounits: ASE, metal, real or other
 *     :returns: converted value
 *     :rtype:
 */
infix fun FloatArray.convert(units: ClosedFloatingPointRange<Double>) = FloatArray(size) { (units.start / units.endInclusive * this[it]).toFloat() }
infix fun Array<FloatArray>.convert(units: ClosedFloatingPointRange<Double>) = Array(size) { r ->
    FloatArray(this[0].size) { c ->
        (units.start / units.endInclusive * this[r][c]).toFloat()
    }
}
