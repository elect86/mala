package ase

fun symbols2numbers(symbols: Any): IntArray {
    var symbols = symbols
    if (symbols is String)
        symbols = TODO()// string2symbols (symbols)
    val numbers = ArrayList<Int>()
    symbols as List<String>
    for(s in symbols) {
//        if isinstance(s, str):
        numbers += atomicNumbers[s]!!
//        else:
//        numbers.append(int(s))
    }
    return numbers.toIntArray()
}

/**
 *     A sequence of chemical symbols.
 *
 *     ``atoms.symbols`` is a :class:`ase.symbols.Symbols` object.  This
 *     object works like an editable view of ``atoms.numbers``, except
 *     its elements are manipulated as strings.
 *
 *     Examples:
 *
 *     >>> from ase.build import molecule
 *     >>> atoms = molecule('CH3CH2OH')
 *     >>> atoms.symbols
 *     Symbols('C2OH6')
 *     >>> atoms.symbols[:3]
 *     Symbols('C2O')
 *     >>> atoms.symbols == 'H'  # doctest: +ELLIPSIS
 *     array([False, False, False,  True,  True,  True,  True,  True,  True]...)
 *     >>> atoms.symbols[-3:] = 'Pu'
 *     >>> atoms.symbols
 *     Symbols('C2OH3Pu3')
 *     >>> atoms.symbols[3:6] = 'Mo2U'
 *     >>> atoms.symbols
 *     Symbols('C2OMo2UPu3')
 *     >>> atoms.symbols.formula
 *     Formula('C2OMo2UPu3')
 *
 *     The :class:`ase.formula.Formula` object is useful for extended
 *     formatting options and analysis.
 */
class Symbols(val numbers: IntArray) {

    val size get() = numbers.size
}
