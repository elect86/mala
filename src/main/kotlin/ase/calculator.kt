package ase


open class BaseCalculator : GetPropertiesMixin() {

    override fun getProperties(name: String, atoms: Any?, allowCalculation: Boolean) {
        TODO("Not yet implemented")
    }
}

/**
 * Base-class for all ASE calculators.
 *
 * A calculator must raise PropertyNotImplementedError if asked for a
 * property that it can't calculate.  So, if calculation of the
 * stress tensor has not been implemented, get_stress(atoms) should
 * raise PropertyNotImplementedError.  This can be achieved simply by not
 * including the string 'stress' in the list implemented_properties
 * which is a class member.  These are the names of the standard
 * properties: 'energy', 'forces', 'stress', 'dipole', 'charges',
 * 'magmom' and 'magmoms'.7
 */
open class Calculator: BaseCalculator()