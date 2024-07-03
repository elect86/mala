package ase

/**
 * Mixin class which provides get_forces(), get_stress() and so on.
 *
 * Inheriting class must implement get_property().
 */
abstract class GetPropertiesMixin {
    /** Get the named property. */
    abstract fun getProperties(name: String, atoms: Any? = null, allowCalculation: Boolean = true)
}
