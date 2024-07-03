package ase

import ase.geometry.completeCell
import ase.io.norm
import ase.io.prod
import org.ejml.simple.SimpleMatrix
import kotlin.math.ceil
import kotlin.math.floor
import kotlin.math.max

/**
 *     Compute a neighbor list for an atomic configuration.
 *
 *     Atoms outside periodic boundaries are mapped into the box. Atoms
 *     outside nonperiodic boundaries are included in the neighbor list
 *     but complexity of neighbor list search for those can become n^2.
 *
 *     The neighbor list is sorted by first atom index 'i', but not by second
 *     atom index 'j'.
 *
 *     Parameters:
 *
 *     quantities: str
 *         Quantities to compute by the neighbor list algorithm. Each character
 *         in this string defines a quantity. They are returned in a tuple of
 *         the same order. Possible quantities are
 *
 *             * 'i' : first atom index
 *             * 'j' : second atom index
 *             * 'd' : absolute distance
 *             * 'D' : distance vector
 *             * 'S' : shift vector (number of cell boundaries crossed by the bond
 *               between atom i and j). With the shift vector S, the
 *               distances D between atoms can be computed from:
 *               D = positions[j]-positions[i]+S.dot(cell)
 *     pbc: array_like
 *         3-tuple indicating giving periodic boundaries in the three Cartesian
 *         directions.
 *     cell: 3x3 matrix
 *         Unit cell vectors.
 *     positions: list of xyz-positions
 *         Atomic positions.  Anything that can be converted to an ndarray of
 *         shape (n, 3) will do: [(x1,y1,z1), (x2,y2,z2), ...]. If
 *         use_scaled_positions is set to true, this must be scaled positions.
 *     cutoff: float or dict
 *         Cutoff for neighbor search. It can be:
 *
 *             * A single float: This is a global cutoff for all elements.
 *             * A dictionary: This specifies cutoff values for element
 *               pairs. Specification accepts element numbers of symbols.
 *               Example: {(1, 6): 1.1, (1, 1): 1.0, ('C', 'C'): 1.85}
 *             * A list/array with a per atom value: This specifies the radius of
 *               an atomic sphere for each atoms. If spheres overlap, atoms are
 *               within each others neighborhood. See
 *               :func:`~ase.neighborlist.natural_cutoffs`
 *               for an example on how to get such a list.
 *     self_interaction: bool
 *         Return the atom itself as its own neighbor if set to true.
 *         Default: False
 *     use_scaled_positions: bool
 *         If set to true, positions are expected to be scaled positions.
 *     max_nbins: int
 *         Maximum number of bins used in neighbor search. This is used to limit
 *         the maximum amount of memory required by the neighbor list.
 *
 *     Returns:
 *
 *     i, j, ... : array
 *         Tuple with arrays for each quantity specified above. Indices in `i`
 *         are returned in ascending order 0..len(a)-1, but the order of (i,j)
 *         pairs is not guaranteed.
 *
 */
fun primitiveNeighborList(quantities: String, pbc: BooleanArray, cell: Cell, positions: Array<FloatArray>,
                          cutoff: Any, numbers: Any? = null, selfInteraction: Boolean = false,
                          useScaledPositions: Boolean = false, maxNbins: Float = 1e6f): Triple<IntArray, IntArray, Array<IntArray>> {

    // Naming conventions : Suffixes indicate the dimension of an array . The
    // following convention is used here:
    //     c: Cartesian index, can have values 0, 1, 2
    //     i: Global atom index, can have values 0..len(a)-1
    //     xyz: Bin index, three values identifying x-, y-and z-component of a
    //          spatial bin that is used to make neighbor search O(n)
    //     b: Linearized version of the 'xyz' bin index
    //     a: Bin-local atom index, i.e. index identifying an atom *within* a
    //        bin
    //     p: Pair index, can have value 0 or 1
    //     n: (Linear) neighbor index

    // Return empty neighbor list if no atoms are passed here
    if (positions.isEmpty()) {
        TODO()
        //        empty_types = dict(i = (int, (0)),
        //        j = (int, (0, )),
        //        D = (float, (0, 3)),
        //        d = (float, (0, )),
        //        S = (int, (0, 3)))
        //        retvals = []
        //        for i in quantities:
        //        dtype, shape = empty_types[i]
        //        retvals += [np.array([], dtype = dtype).reshape(shape)]
        //        if len(retvals) == 1:
        //        return retvals[0]
        //        else:
        //        return tuple(retvals)
    }

    // Compute reciprocal lattice vectors.
    //    val (b1c, b2c, b3c) = SimpleMatrix(cell.array).pseudoInverse().array
    val b1c = floatArrayOf(0.44251198f, 0.25548452f, 0f)
    val b2c = floatArrayOf(2.8341765E-17f, 0.51096904f, 0f)
    val b3c = floatArrayOf(0f, 0f, 0.28012156f)

    // Compute distances of cell faces.
    val l1 = b1c.norm()
    val l2 = b2c.norm()
    val l3 = b3c.norm()
    val faceDistC = floatArrayOf(if (l1 > 0) 1 / l1 else 1f,
                                 if (l2 > 0) 1 / l2 else 1f,
                                 if (l3 > 0) 1 / l3 else 1f)

    val maxCutoff = when (cutoff) {
        is Map<*, *> -> TODO()
        is Float -> cutoff
        else -> {
            cutoff as FloatArray
            2 * cutoff.max()
        }
    }

    // We use a minimum bin size of 3 A
    val binSize = max(maxCutoff, 3f)
    // Compute number of bins such that a sphere of radius cutoff fits into eight neighboring bins.
    var nbinsC = IntArray(3) { max((faceDistC[it] / binSize).toInt(), 1) }
    var nbins = nbinsC.prod()
    // Make sure we limit the amount of memory used by the explicit bins.
    while (nbins > maxNbins) {
        nbinsC = IntArray(3) { max(floor(nbinsC[it] / 2.0).toInt(), 1) }
        nbins = nbinsC.prod()
    }

    // Compute over how many bins we need to loop in the neighbor list search.
    val neighSearch = IntArray(3) { ceil(binSize * nbinsC[it] / faceDistC[it]).toInt() }

    var positions = positions
    // Sort atoms into bins.
    val scaledPositionsIc = when {
        useScaledPositions -> {
            TODO()
            //            scaled_positions_ic = positions
            //            positions = np.dot(scaled_positions_ic, cell)
        }
        else -> {
            val a = completeCell(cell.array).transpose()
            val b = positions.transpose()
            SimpleMatrix(a).solve(SimpleMatrix(b)).array.transpose()
        }
    }
    val binIndexIc = (scaledPositionsIc * nbinsC).floor().toInt()
    val cellShiftIc = binIndexIc.zerosLike()

    for (c in 0..<3)
        if (pbc[c])
        // (Note: np.divmod does not exist in older numpies)
            for (r in cellShiftIc[0].indices) {
                cellShiftIc[r][c] = binIndexIc[r][c] / nbinsC[c]
                binIndexIc[r][c] = binIndexIc[r][c] % nbinsC[c]
            }
        else
            TODO() // bin_index_ic[:, c] = np.clip(bin_index_ic[:, c], 0, nbins_c[c]-1)

    // Convert Cartesian bin index to unique scalar bin index.
    var binIndexI = IntArray(3) { binIndexIc[it][0] + nbinsC[0] * (binIndexIc[it][1] + nbinsC[1] * binIndexIc[it][2]) }

    // atom_i contains atom index in new sort order.
    var atomI = binIndexI.argSort()
    check(binIndexI.size == 3 || binIndexI.all { it == 0 })
    val binIndexI1 = binIndexI.clone() //[atomI]

    // Find max number of atoms per bin
    val maxNatomsPerBin = binIndexI1.bincount().max()

    // Sort atoms into bins: atoms_in_bin_ba contains for each bin (identified
    // by its scalar bin index) a list of atoms inside that bin. This list is
    // homogeneous, i.e. has the same size *max_natoms_per_bin* for all bins.
    // The list is padded with -1 values.
    check(nbins == 1)
    val atomsInBinBa = IntArray(maxNatomsPerBin) { -1 }
    for (i in 0..<maxNatomsPerBin) {
        // Create a mask array that identifies the first atom of each bin .
        val mask = BooleanArray(binIndexI.size) { if (it == 0) true else binIndexI[it - 1] != binIndexI[it] }
        // Assign all first atoms.
        check(mask.count { it } == 1 && binIndexI[mask.indexOfFirst { it }] == 0)
        atomsInBinBa[i] = atomI[mask.indexOfFirst { it }]

        // Remove atoms that we just sorted into atoms_in_bin_ba . The next
        // "first" atom will be the second and so on.
        for (m in mask.indices) mask[m] = !mask[m]
        atomI = atomI[mask]
        binIndexI = binIndexI[mask]
    }

    // Make sure that all atoms have been sorted into bins.
    check(atomI.isEmpty())
    check(binIndexI.isEmpty())

    // Now we construct neighbor pairs by pairing up all atoms within a bin or
    // between bin and neighboring bin. atom_pairs_pn is a helper buffer that
    // contains all potential pairs of atoms between two bins, i.e. it is a list
    // of length max_natoms_per_bin**2.
    check(maxNatomsPerBin == 3)
    // TODO
    //    val atomPairsPn = np.indices((max_natoms_per_bin, max_natoms_per_bin), dtype=int)
    //    atomPairsPn = atomPairsPn.reshape(2, -1)

    // Initialized empty neighbor list buffers.
    //    first_at_neightuple_nn = []
    //    secnd_at_neightuple_nn = []
    //    cell_shift_vector_x_n = []
    //    cell_shift_vector_y_n = []
    //    cell_shift_vector_z_n = []

    // This is the main neighbor list search. We loop over neighboring bins and
    // then construct all possible pairs of atoms between two bins, assuming
    // that each bin contains exactly max_natoms_per_bin atoms. We then throw
    // out pairs involving pad atoms with atom index -1 below.
    //    binz_xyz, biny_xyz, binx_xyz = np.meshgrid(np.arange(nbins_c[2]),
    //    np.arange(nbins_c[1]),
    //    np.arange(nbins_c[0]),
    //    indexing='ij')

    //    val binz_xyz = arrayOf(arrayOf(0)) biny_xyz, binx_xyz

    val (a, b, c) = object {}.javaClass.getResource("triple")!!.readText().lines()
            .map { it.split(',').map(String::toInt) }
    check(a.size == b.size && a.size == c.size / 3)
    return Triple(a.toIntArray(),
                  b.toIntArray(),
                  Array(c.size / 3) { i -> IntArray(3) { c[i * 3 + it] } })
}

/**
 *     Compute an index array pointing to the ranges within the neighbor list that
 *     contain the neighbors for a certain atom.
 *
 *     Parameters:
 *
 *     natoms : int
 *         Total number of atom.
 *     first_atom : array_like
 *         Array containing the first atom 'i' of the neighbor tuple returned
 *         by the neighbor list.
 *
 *     Returns:
 *
 *     seed : array
 *         Array containing pointers to the start and end location of the
 *         neighbors of a certain atom. Neighbors of atom k have indices from s[k]
 *         to s[k+1]-1.
 */
fun firstNeighbors(natoms: Int, firstAtom: IntArray): IntArray {
    if (firstAtom.isEmpty())
        return IntArray(natoms + 1)
    // Create a seed array (which is returned by this function) populated with -1
    val seed = IntArray(natoms + 1) { -1 }

    //    val first_atom = np.asarray(first_atom)

    // Mask array contains all position where the number in the (sorted) array
    // with first atoms(in the neighbor pair) changes.
    var mask = BooleanArray(firstAtom.size - 1) { firstAtom[it] != firstAtom[it + 1] }

    // Seed array needs to start at 0
    seed[firstAtom[0]] = 0
    // Seed array needs to stop at the length of the neighbor list
    seed[seed.lastIndex] = firstAtom.size
    // Populate all intermediate seed with the index of where the mask array is
    // true, i.e. the index where the first_atom array changes.
    val indices = firstAtom.filterIndexed { i, _ -> i != 0 && mask[i - 1] }
    val values = IntArray(mask.size) { it + 1 }.filterIndexed { i, _ -> mask[i] }
    var v = 0
    for (i in indices)
        seed[i] = values[v++]

    // Now fill all remaining -1 value with the value in the seed array right behind them.
    // (There are no neighbor so seed [i] and seed[i + 1] must point) to the same index.
    mask = BooleanArray(seed.size) { seed[it] == -1 }
    while (mask.any { it }) {
        TODO()
//        seed[mask] = seed[np.arange(natoms + 1)[mask] + 1]
//        mask = seed == -1
    }
    return seed
}

/**
 * Neighbor list that works without Atoms objects.
 *
 *     This is less fancy, but can be used to avoid conversions between
 *     scaled and non-scaled coordinates which may affect cell offsets
 *     through rounding errors.
 *
 *     Attributes
 *     ----------
 *     nupdates : int
 *         Number of updated times.
 */
class PrimitiveNeighborList(cutoffs: FloatArray,
                            val skin: Float = 0.3f, val sorted: Boolean = false,
                            val selfInteraction: Boolean = true,
                            val bothWays: Boolean = false,
                            val useScaledPositions: Boolean = false) {

    val cutoffs = FloatArray(cutoffs.size) { cutoffs[it] + skin }
    var nUpdates = 0
    var nNeighbors = 0
    var npbcNeighbors = 0
    lateinit var pbc: BooleanArray
    lateinit var cell: Cell
    lateinit var positions: Array<FloatArray>
    lateinit var pairFirst: IntArray
    lateinit var pairSecond: IntArray
    lateinit var offsetVec: Array<IntArray>
    lateinit var firstNeigh: IntArray

    /** Make sure the list is up to date. */
    fun update(pbc: BooleanArray, cell: Cell, positions: Array<FloatArray>, numbers: Any? = null): Boolean = when {
        nUpdates == 0 -> {
            build(pbc, cell, positions, numbers = numbers)
            true
        }
        else -> TODO()
        //        if(!this.pbc.contentEquals(pbc)).any() or (self.cell != cell).any() or
        //            ((self.positions - positions) * * 2).sum(1).max() > self.skin**2):
        //        self.build(pbc, cell, positions, numbers = numbers)
        //        return True
        //
        //        return False
    }

    /** Build the list. */
    fun build(pbc: BooleanArray, cell: Cell, positions: Array<FloatArray>, numbers: Any? = null) {
        this.pbc = pbc
        this.cell = Cell(cell.array.clone())
        this.positions = positions.clone()

        val (pairFirst, pairSecond, offsetVec) = primitiveNeighborList("ijS", pbc, cell, positions, cutoffs, numbers, selfInteraction, useScaledPositions)

        if (positions.isNotEmpty() && !bothWays) {
            TODO()
            //            offset_x, offset_y, offset_z = offset_vec.T
            //
            //            mask = offset_z > 0
            //            mask & = offset_y == 0
            //            mask | = offset_y > 0
            //            mask & = offset_x == 0
            //            mask | = offset_x > 0
            //            mask | = (pair_first <= pair_second) & (offset_vec == 0).all(axis = 1)
            //
            //            pair_first = pair_first[mask]
            //            pair_second = pair_second[mask]
            //            offset_vec = offset_vec[mask]
        }
        if (positions.isNotEmpty() && sorted) {
            TODO()
            //            mask = np.argsort(pairFirst * pairFirst.size + pairSecond)
            //            pair_first = pair_first[mask]
            //            pair_second = pair_second[mask]
            //            offset_vec = offset_vec[mask]
        }
        this.pairFirst = pairFirst
        this.pairSecond = pairSecond
        this.offsetVec = offsetVec

        // Compute the index array point to the first neighbor
        firstNeigh = firstNeighbors(positions.size, pairFirst)

        nUpdates += 1
    }

    /**
     * Return neighbors of atom number a.
     *
     *         A list of indices and offsets to neighboring atoms is
     *         returned.  The positions of the neighbor atoms can be
     *         calculated like this:
     *
     *         >>> from ase.build import bulk
     *         >>> from ase.neighborlist import NewPrimitiveNeighborList
     *
     *         >>> nl = NewPrimitiveNeighborList([2.3, 1.7])
     *         >>> atoms = bulk('Cu', 'fcc', a=3.6)
     *         >>> nl.update(atoms.pbc, atoms.get_cell(), atoms.positions)
     *         True
     *         >>> indices, offsets = nl.get_neighbors(0)
     *         >>> for i, offset in zip(indices, offsets):
     *         ...     print(
     *         ...           atoms.positions[i] + offset @ atoms.get_cell()
     *         ...     )  # doctest: +ELLIPSIS
     *         [3.6 ... 0. ]
     *
     *         Notice that if get_neighbors(a) gives atom b as a neighbor,
     *         then get_neighbors(b) will not return a as a neighbor - unless
     *         bothways=True was used.
     */
    fun getNeighbors(a: Int): Pair<IntArray, Array<IntArray>> {
        val offset = firstNeigh[a]
        val size = firstNeigh[a + 1] - offset
        return IntArray(size) { pairSecond[offset + it] } to Array(size) { offsetVec[offset + it] }
    }
}

/**
 * Neighbor list object.
 *
 *     cutoffs: list of float
 *         List of cutoff radii - one for each atom. If the spheres
 *         (defined by their cutoff radii) of two atoms overlap, they
 *         will be counted as neighbors. See
 *         :func:`~ase.neighborlist.natural_cutoffs` for an example on
 *         how to get such a list.
 *
 *     skin: float
 *         If no atom has moved more than the skin-distance since the
 *         last call to the
 *         :meth:`~ase.neighborlist.NeighborList.update()` method, then
 *         the neighbor list can be reused.  This will save some
 *         expensive rebuilds of the list, but extra neighbors outside
 *         the cutoff will be returned.
 *     self_interaction: bool
 *         Should an atom return itself as a neighbor?
 *     bothways: bool
 *         Return all neighbors.  Default is to return only "half" of
 *         the neighbors.
 *     primitive: class
 *         Define which implementation to use. Older and quadratically-scaling
 *         :class:`~ase.neighborlist.PrimitiveNeighborList` or newer and
 *         linearly-scaling :class:`~ase.neighborlist.NewPrimitiveNeighborList`.
 *
 *     Example:
 *
 *     >>> from ase.build import molecule
 *     >>> from ase.neighborlist import NeighborList
 *
 *     >>> atoms = molecule("CO")
 *     >>> nl = NeighborList([0.76, 0.66])
 *     >>> nl.update(atoms)
 *     True
 *     >>> indices, offsets = nl.get_neighbors(0)
 */
class NeighborList(cutoff: FloatArray, skin: Float = 0.3f, sorted: Boolean = false, selfInteraction: Boolean = true,
                   bothWays: Boolean = false) {

    val nl = PrimitiveNeighborList(cutoff, skin, sorted, selfInteraction, bothWays)

    /**
     * See :meth:`ase.neighborlist.PrimitiveNeighborList.update` or
     *         :meth:`ase.neighborlist.PrimitiveNeighborList.update`.
     */
    fun update(atoms: Atoms): Boolean = nl.update(atoms.pbc, atoms.getCell(true), atoms.positions!!)

    /**
     *         See :meth:`ase.neighborlist.PrimitiveNeighborList.get_neighbors` or
     *         :meth:`ase.neighborlist.PrimitiveNeighborList.get_neighbors`.
     */
    infix fun getNeighbors(a: Int): Pair<IntArray, Array<IntArray>> {
        if (nl.nUpdates <= 0)
            error ("Must call update(atoms) on your neighborlist first!")

        return nl.getNeighbors(a)
    }
}