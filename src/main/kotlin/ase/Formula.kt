package ase

/**
 * Chemical formula object.
 *
 *         Parameters
 *         ----------
 *         formula: str
 *             Text string representation of formula.  Examples: ``'6CO2'``,
 *             ``'30Cu+2CO'``, ``'Pt(CO)6'``.
 *         strict: bool
 *             Only allow real chemical symbols.
 *         format: str
 *             Reorder according to *format*.  Must be one of hill, metal,
 *             ab2, a2b, periodic or reduce.
 *
 *         Examples
 *         --------
 *         >>> from ase.formula import Formula
 *         >>> w = Formula('H2O')
 *         >>> w.count()
 *         {'H': 2, 'O': 1}
 *         >>> 'H' in w
 *         True
 *         >>> w == 'HOH'
 *         True
 *         >>> f'{w:latex}'
 *         'H$_{2}$O'
 *         >>> w.format('latex')
 *         'H$_{2}$O'
 *         >>> divmod(6 * w + 'Cu', w)
 *         (6, Formula('Cu'))
 *
 *         Raises
 *         ------
 *         ValueError
 *             on malformed formula
 */
class Formula(formula: Any = "",
              strict: Boolean = false,
              format: Format? = null,
              tree: Any? = null,
              count: Any? = null) {

    init {
        check(formula is String || formula is Formula)
        val formula = if (formula is String) formula else (formula as Formula).toString()
        TODO()
    }

    enum class Format { hill, metal, ab2, a2b, periodic, reduce }
}