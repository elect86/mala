package ase

/**
 * Class for representing a single atom.
 *
 *     Parameters:
 *
 *     symbol: str or int
 *         Can be a chemical symbol (str) or an atomic number (int).
 *     position: sequence of 3 floats
 *         Atomic position.
 *     tag: int
 *         Special purpose tag.
 *     momentum: sequence of 3 floats
 *         Momentum for atom.
 *     mass: float
 *         Atomic mass in atomic units.
 *     magmom: float or 3 floats
 *         Magnetic moment.
 *     charge: float
 *         Atomic charge.
 */
class Atom(symbol: Any = "X",
           val position: FloatArray = FloatArray(3),
           val tag: Any? = null, val momentum: FloatArray? = null, val mass: Any? = null,
           val magmom: FloatArray? = null, val charge: Any? = null,
           val atoms: Atoms? = null, val index: Int? = null) {

    var number: Int? = null
    init {
        if (atoms == null)
            // This atom is not part of any Atoms object:
            number = if (symbol is String) atomicNumbers[symbol]!! else TODO()// symbol
    }
}