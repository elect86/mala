package ase

/**
 * Atoms object.
 *
 *     The Atoms object can represent an isolated molecule, or a
 *     periodically repeated structure.  It has a unit cell and
 *     there may be periodic boundary conditions along any of the three
 *     unit cell axes.
 *     Information about the atoms (atomic numbers and position) is
 *     stored in ndarrays.  Optionally, there can be information about
 *     tags, momenta, masses, magnetic moments and charges.
 *
 *     In order to calculate energies, forces and stresses, a calculator
 *     object has to attached to the atoms object.
 *
 *     Parameters:
 *
 *     symbols: str (formula) or list of str
 *         Can be a string formula, a list of symbols or a list of
 *         Atom objects.  Examples: 'H2O', 'COPt12', ['H', 'H', 'O'],
 *         [Atom('Ne', (x, y, z)), ...].
 *     positions: list of xyz-positions
 *         Atomic positions.  Anything that can be converted to an
 *         ndarray of shape (n, 3) will do: [(x1,y1,z1), (x2,y2,z2),
 *         ...].
 *     scaled_positions: list of scaled-positions
 *         Like positions, but given in units of the unit cell.
 *         Can not be set at the same time as positions.
 *     numbers: list of int
 *         Atomic numbers (use only one of symbols/numbers).
 *     tags: list of int
 *         Special purpose tags.
 *     momenta: list of xyz-momenta
 *         Momenta for all atoms.
 *     masses: list of float
 *         Atomic masses in atomic units.
 *     magmoms: list of float or list of xyz-values
 *         Magnetic moments.  Can be either a single value for each atom
 *         for collinear calculations or three numbers for each atom for
 *         non-collinear calculations.
 *     charges: list of float
 *         Initial atomic charges.
 *     cell: 3x3 matrix or length 3 or 6 vector
 *         Unit cell vectors.  Can also be given as just three
 *         numbers for orthorhombic cells, or 6 numbers, where
 *         first three are lengths of unit cell vectors, and the
 *         other three are angles between them (in degrees), in following order:
 *         [len(a), len(b), len(c), angle(b,c), angle(a,c), angle(a,b)].
 *         First vector will lie in x-direction, second in xy-plane,
 *         and the third one in z-positive subspace.
 *         Default value: [0, 0, 0].
 *     celldisp: Vector
 *         Unit cell displacement vector. To visualize a displaced cell
 *         around the center of mass of a Systems of atoms. Default value
 *         = (0,0,0)
 *     pbc: one or three bool
 *         Periodic boundary conditions flags.  Examples: True,
 *         False, 0, 1, (1, 1, 0), (True, False, False).  Default
 *         value: False.
 *     constraint: constraint object(s)
 *         Used for applying one or more constraints during structure
 *         optimization.
 *     calculator: calculator object
 *         Used to attach a calculator for calculating energies and atomic
 *         forces.
 *     info: dict of key-value pairs
 *         Dictionary of key-value pairs with additional information
 *         about the system.  The following keys may be used by ase:
 *
 *           - spacegroup: Spacegroup instance
 *           - unit_cell: 'conventional' | 'primitive' | int | 3 ints
 *           - adsorbate_info: Information about special adsorption sites
 *
 *         Items in the info attribute survives copy and slicing and can
 *         be stored in and retrieved from trajectory files given that the
 *         key is a string, the value is JSON-compatible and, if the value is a
 *         user-defined object, its base class is importable.  One should
 *         not make any assumptions about the existence of keys.
 *
 *     Examples:
 *
 *     These three are equivalent:
 *
 *     >>> from ase import Atom
 *
 *     >>> d = 1.104  # N2 bondlength
 *     >>> a = Atoms('N2', [(0, 0, 0), (0, 0, d)])
 *     >>> a = Atoms(numbers=[7, 7], positions=[(0, 0, 0), (0, 0, d)])
 *     >>> a = Atoms([Atom('N', (0, 0, 0)), Atom('N', (0, 0, d))])
 *
 *     FCC gold:
 *
 *     >>> a = 4.05  # Gold lattice constant
 *     >>> b = a / 2
 *     >>> fcc = Atoms('Au',
 *     ...             cell=[(0, b, b), (b, 0, b), (b, b, 0)],
 *     ...             pbc=True)
 *
 *     Hydrogen wire:
 *
 *     >>> d = 0.9  # H-H distance
 *     >>> h = Atoms('H', positions=[(0, 0, 0)],
 *     ...           cell=(d, 0, 0),
 *     ...           pbc=(1, 0, 0))
 */
class Atoms(symbols: Any? = null,
            positions: Array<FloatArray>? = null, numbers: IntArray? = null,
            tags: Any? = null, momenta: Any? = null, masses: Any? = null,
            magmoms: Any? = null, charges: Any? = null,
            scaledPositions: Any? = null,
            cell: Array<FloatArray>? = null, pbcArg: Boolean? = null, cellDisp: FloatArray? = null,
            constraint: Any? = null,
            calculator: Calculator? = null,
            var info: Map<String, Any>? = null,
            velocities: Any? = null) : Iterable<Atom> {

    internal var cellObj = Cell()
    internal val pbc = BooleanArray(3)

    val arrays = mutableMapOf<String, Any>()

    var numbers: IntArray?
        get() = arrays["numbers"] as? IntArray
        set(value) = arrays.set("numbers", value as Any)
    var cell: Cell? = null
    private var cellDisp: FloatArray? = null
    var positions: Array<FloatArray>?
        get() = arrays["positions"] as? Array<FloatArray>
        set(value) = arrays.set("positions", value as Any)
    var tags: Any? = null
    var momenta: Any? = null
    var initialMagmoms: Any? = null
    var masses: Any? = null
    var initialCharges: Any? = null
    var calc: Calculator? = null
        set(value) {
            field = value
            //            if hasattr(calc, 'set_atoms'):
            //            calc.set_atoms(self)
        }

    init {
        if (symbols == null && positions == null && numbers == null && tags == null && momenta == null && masses == null
            && magmoms == null && charges == null && scaledPositions == null && cell == null && pbcArg == null &&
            cellDisp == null && constraint == null && calculator == null && info == null && velocities == null)
            Unit // it's a copy
        else {

            var atoms: Atoms? = null
            var symbols = symbols

            //            if hasattr(symbols, 'get_positions'):
            //            atoms = symbols
            //            symbols = None
            /*else*/ if (symbols is List<*> /*(list, tuple)*/ && symbols.isNotEmpty() && symbols[0] is Atom) {
                // Get data from a list or tuple of Atom objects:
                //                data = [[atom.get_raw(name) for atom in symbols]
                //                for name in
                //                        ['position', 'number', 'tag', 'momentum',
                //                            'mass', 'magmom', 'charge']]
                symbols = symbols as List<Atom>
                atoms = Atoms(null, symbols.map(Atom::position).toTypedArray(),
                              symbols.map { it.number!! }.toIntArray())
                symbols = null
            }

            var numbers = numbers
            var positions = positions
            var cell = cell
            var cellDisp = cellDisp
            var pbc: BooleanArray? = null
            var constraint = constraint
            var calculator = calculator

            if (atoms != null) {
                // Get data from another Atoms object :
                if (scaledPositions != null) TODO("not implemented")
                if (symbols == null && numbers == null) numbers = atoms.getAtomicNumbers()
                if (positions == null) positions = atoms.getPositions()
                if (tags == null && atoms.tags != null) TODO() // tags = atoms.get_tags ()
                if (momenta == null && atoms.momenta != null) TODO() // momenta = atoms.get_momenta ()
                if (magmoms == null && atoms.initialMagmoms != null) TODO() // magmoms = atoms.get_initial_magnetic_moments ()
                if (masses == null && atoms.masses != null) TODO() // masses = atoms.get_masses ()
                if (charges == null && atoms.initialCharges != null) TODO() // charges = atoms.get_initial_charges ()
                if (cell == null) cell = atoms.getCell().array
                if (cellDisp == null) cellDisp = atoms.getCellDisp()
                if (pbcArg == null) pbc = atoms.getPbc()
                if (constraint == null && atoms.constraints != null) TODO() // constraint = [c.copy() for c in atoms.constraints]
                if (calculator == null) calculator = atoms.calc
                if (info == null) info = atoms.info?.toMap()
            }

            if (symbols == null) {
                numbers = numbers ?: IntArray(when {
                                                  positions != null -> (positions as Array<*>).size
                                                  scaledPositions != null -> TODO() // len (scaled_positions)
                                                  else -> 0
                                              })
                this.numbers = numbers
            } else
                if (numbers != null)
                    error("Use only one of \"symbols\" and \"numbers\".")
                else
                    this.numbers = symbols2numbers(symbols)

            setCell(cell ?: Array(3) { FloatArray(3) })

            this.cellDisp = cellDisp ?: FloatArray(3)

            this.positions = when {
                positions == null -> when {
                    scaledPositions == null -> Array(this.numbers!!.size) { FloatArray(3) }
                    else -> {
                        TODO()
                        //                assert(cell!!.size == 3)
                        //            positions = np.dot(scaled_positions, self.cell)
                    }
                }
                scaledPositions != null -> error("Use only one of \"symbols\" and \"numbers\".")
                else -> positions
            }

            setConstraint(constraint)
            tags?.let { TODO() } //self.set_tags(default(tags, 0))
            masses?.let { TODO() } //self.set_masses(default(masses, None))
            magmoms?.let { TODO() } //self.set_initial_magnetic_moments(default(magmoms, 0.0))
            charges?.let { TODO() } //self.set_initial_charges(default(charges, 0.0))
            setPbc(pbcArg ?: false)
            momenta?.let { TODO() } //self.set_momenta(default(momenta, (0.0, 0.0, 0.0)), apply_constraint=False)

            if (velocities != null) {
                TODO()
                //            if momenta is None:
                //            self.set_velocities(velocities)
                //            else:
                //            raise TypeError ('Use only one of "momenta" and "velocities".')
            }

            if (info == null) {
                //            self.info = {}
            } else
                TODO() //self.info = dict(info)

            calc = calculator
        }
    }

    val constraints: Any? = null

    /** Get the unit cell displacement vectors. */
    fun getCellDisp() = cellDisp!!.clone()

    /** Get periodic boundary condition flags. */
    fun getPbc() = pbc.clone()

    /** Get integer array of atomic numbers. */
    fun getAtomicNumbers() = numbers!!.clone()

    /**
     * Get array of positions.
     *
     *         Parameters:
     *
     *         wrap: bool
     *             wrap atoms back to the cell before returning positions
     *         wrap_kw: (keyword=value) pairs
     *             optional keywords `pbc`, `center`, `pretty_translation`, `eps`,
     *             see :func:`ase.geometry.wrap_positions`
     */
    fun getPositions(wrap: Boolean = false, wrapKw: Map<String, Any> = emptyMap()): Array<FloatArray> = when {
        wrap -> TODO()
        //        if 'pbc' not in wrap_kw :
        //        wrap_kw['pbc'] = self.pbc
        //        return wrap_positions(self.positions, self.cell, ** wrap_kw)
        else -> positions!!.clone()
    }

    /** Return a copy. */
    fun copy() = Atoms().also {
        pbc.copyInto(it.pbc)
        it.numbers = numbers!!.clone()
        it.cell = Cell(cell!!.array)
        it.cellDisp = cellDisp!!.clone()
        it.positions = Array(positions!!.size) { i -> positions!![i] }
    }

    val len
        get() = positions!!.size

    /** Extend atoms object by appending atoms from *other*. */
    fun extend(other: Any) {
        var other = other
        if (other is Atom)
            other = TODO() // self.__class__ ([other])

        other as Atoms
        val n1 = len
        val n2 = other.len

        for ((name, a1) in arrays)
            when (name) {
                "numbers" -> {
                    a1 as IntArray
                    val a = IntArray(n1 + n2) {
                        when {
                            it < a1.size -> a1[it]
                            else -> other.numbers?.get(it - a1.size) ?: 0
                        }
                    }
                    arrays[name] = a
                }
                "positions" -> {
                    a1 as Array<FloatArray>
                    val a = Array(n1 + n2) { x ->
                        FloatArray(a1[0].size) { y ->
                            when {
                                x < a1.size -> a1[x][y]
                                else -> other.positions?.get(x - a1.size)?.get(y) ?: 0f
                            }
                        }
                    }
                    arrays[name] = a
                }
                else -> TODO()
            }
        for ((name, a2) in other.arrays) {
            if (name in arrays)
                continue
            TODO()
            //            a = np.empty((n1 + n2) + a2.shape[1:], a2.dtype)
            //            a[n1:] = a2
            //            if name == 'masses':
            //            a[:n1] = self.get_masses()[:n1]
            //            else:
            //            a[:n1] = 0
            //
            //            self.set_array(name, a)
        }
    }

    /** [Python] append/__iadd__
     *  Append atom to end. */
    operator fun plusAssign(atom: Atom) {
        extend(Atoms(listOf(atom)))
    }

    /**
     * Return a subset of the atoms.
     *
     * i -- scalar integer, list of integers, or slice object
     * describing which atoms to return.
     *
     * If i is a scalar, return an Atom object. If i is a list or a
     * slice, return an Atoms object with the same cell, pbc, and
     * other associated info as the original Atoms object. The
     * indices of the constraints will be shuffled so that they match
     * the indexing in the subset returned.
     */
    operator fun get(i: Int): Atom {
        if (i !in 0..<len)
            error("Index $i out of range. ($len)")
        return Atom(atoms = this, index = i)
    }

    /**
     * Apply one or more constrains.
     *
     * The *constraint* argument must be one constraint object or a list of constraint objects.
     */
    fun setConstraint(constraint: Any? = null) {
        if (constraint == null) {
            //            self._constraints = []
        } else {
            TODO()
            //            if isinstance(constraint, list):
            //            self._constraints = constraint
            //            elif isinstance (constraint, tuple):
            //            self._constraints = list(constraint)
            //            else:
            //            self._constraints = [constraint]
        }
    }

    /**
     * Set unit cell vectors.
     *
     *         Parameters:
     *
     *         cell: 3x3 matrix or length 3 or 6 vector
     *             Unit cell.  A 3x3 matrix (the three unit cell vectors) or
     *             just three numbers for an orthorhombic cell. Another option is
     *             6 numbers, which describes unit cell with lengths of unit cell
     *             vectors and with angles between them (in degrees), in following
     *             order: [len(a), len(b), len(c), angle(b,c), angle(a,c),
     *             angle(a,b)].  First vector will lie in x-direction, second in
     *             xy-plane, and the third one in z-positive subspace.
     *         scale_atoms: bool
     *             Fix atomic positions or move atoms with the unit cell?
     *             Default behavior is to *not* move the atoms (scale_atoms=False).
     *         apply_constraint: bool
     *             Whether to apply constraints to the given cell.
     *
     *         Examples:
     *
     *         Two equivalent ways to define an orthorhombic cell:
     *
     *         >>> atoms = Atoms('He')
     *         >>> a, b, c = 7, 7.5, 8
     *         >>> atoms.set_cell([a, b, c])
     *         >>> atoms.set_cell([(a, 0, 0), (0, b, 0), (0, 0, c)])
     *
     *         FCC unit cell:
     *
     *         >>> atoms.set_cell([(0, b, b), (b, 0, b), (b, b, 0)])
     *
     *         Hexagonal unit cell:
     *
     *         >>> atoms.set_cell([a, a, c, 90, 90, 120])
     *
     *         Rhombohedral unit cell:
     *
     *         >>> alpha = 77
     *         >>> atoms.set_cell([a, a, a, alpha, alpha, alpha])
     */
    fun setCell(cell: Array<FloatArray>, scaleAtoms: Boolean = false, applyConstraint: Boolean = true) {
        // Override pbcs if and only if given a Cell object:
        val cell = Cell(cell)

        // XXX not working well during initialize due to missing _constraints
        if (applyConstraint && false/* and hasattr(self, '_constraints')*/) {
            TODO()
            //        for constraint in self.constraints:
            //        if hasattr(constraint, 'adjust_cell'):
            //        constraint.adjust_cell(self, cell)
        }
        if (scaleAtoms) {
            TODO()
            //            M = np.linalg.solve(self.cell.complete(), cell.complete())
            //            self.positions[:] = np.dot(self.positions, M)
        }

        this.cell = cell
    }

    /**
     * Get the three unit cell vectors as a `class`:ase.cell.Cell` object.
     *
     * The Cell object resembles a 3x3 ndarray, and cell[i, j]
     * is the jth Cartesian coordinate of the ith cell vector.
     */
    fun getCell(complete: Boolean = false): Cell = if (complete) cell!!.complete() else cell!!.copy()

    /** Set periodic boundary condition flags. */
    fun setPbc(pbc: Boolean) = this.pbc.fill(pbc)

    /**
     * Get positions relative to unit cell.
     *
     * If wrap is True, atoms outside the unit cell will be wrapped into
     * the cell in those directions with periodic boundary conditions
     * so that the scaled coordinates are between zero and one.
     *
     * If any cell vectors are zero, the corresponding coordinates
     * are evaluated as if the cell were completed using
     * ``cell.complete()``.  This means coordinates will be Cartesian
     * as long as the non-zero cell vectors span a Cartesian axis or
     * plane.
     */
    fun getScaledPositions(wrap: Boolean = true): Array<FloatArray> {

        val fractional = cell!!.scaledPositions(positions!!)

        if (wrap)
            for ((i, periodic) in pbc.withIndex())
                if (periodic)
                // Yes, we need to do it twice.
                // See the scaled_positions.py test .
                    for (x in fractional.indices)
                        fractional[x][i] %= 1f
        //        fractional[:, i] %= 1.0

        return fractional
    }

    /** Set positions relative to unit cell. */
    fun setScaledPositions(scaled: Array<FloatArray>) {
        positions = cell!!.cartesianPositions(scaled)
    }

    val indices: IntRange
        get() = 0..<len

    override fun iterator(): Iterator<Atom> = AtomsIterator()

    inner class AtomsIterator : Iterator<Atom> {
        var i = 0
        override fun hasNext(): Boolean = i < len
        override fun next(): Atom = get(i++)
    }
}