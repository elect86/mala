package casus.mala.descriptors

import ai.djl.ndarray.NDArray
import ase.cdist
import ase.factorial
import ase.now
import casus.mala.common.Parameters
import java.io.File
import kotlin.math.max
import kotlin.math.min
import kotlin.math.sqrt

/** Class for calculation and parsing of bispectrum descriptors. */
class Bispectrum(parameters: Parameters) : Descriptor(parameters) {

    override val dataName = "Bispectrum"

    //def __init__(self, parameters):
    //super(Bispectrum, self).__init__(parameters)
    //
    //# Index arrays needed only when computing the bispectrum descriptors
    //# via python.
    //# They are later filled in the __init_index_arrays() function.
    lateinit var indexUBlock: IntArray
    var indexUMax = 0

    private lateinit var cglist: FloatArray

    //self.__index_u_one_initialized = None
    private val indexUFull = ArrayList<Int>()
    private val indexUSymmetryPos = ArrayList<Int>()
    private val indexUSymmetryNeg = ArrayList<Int>()
    private val indexU1Full = ArrayList<Int>()
    private val indexU1SymmetryPos = ArrayList<Int>()
    private val indexU1SymmetryNeg = ArrayList<Int>()
    private val rootpqFull1 = ArrayList<Float>()
    private val rootpqFull2 = ArrayList<Float>()

    private val indexZU1r = ArrayList<Int>()
    private val indexZU1i = ArrayList<Int>()
    private val indexZU2r = ArrayList<Int>()
    private val indexZU2i = ArrayList<Int>()
    private val indexZIcga = ArrayList<Int>()
    private val indexZIcgb = ArrayList<Int>()
    private val indexZJjz = ArrayList<Int>()
    private lateinit var indexZBlock: Array<Array<IntArray>>

    var indexBMax = 0
    //self.__index_b = None
    var rMin0 = 0f
    var rFac0 = 0f
    var bZeroFlag = false
    var wselfallFlag = false
    var bnormFlag = false
    var quadraticFlag = false
    var numberElements = 0
    var wSelf = 0f

    fun convertUnits(array: NDArray): Nothing = TODO() // convertUnits (array, null)

    override fun calculate(outDir: File, kwargs: Map<String, Any?>): Number = when {
        else -> calculate(kwargs)
        //        parameters.configuration.lammps
        //        if find_spec("lammps") is None:
        //        printout(
        //            "No LAMMPS found for descriptor calculation, "
        //            "falling back to python."
        //        )
        //        return self.__calculate_python(** kwargs)
        //        else:
        //        return self.__calculate_lammps(outdir, ** kwargs)
        //        else:
        //        return self.__calculate_python(** kwargs)

    }

    /**
     * Perform bispectrum calculation using python.
     *
     * The code used to this end was adapted from the LAMMPS implementation.
     * It serves as a fallback option whereever LAMMPS is not available.
     * This may be useful, e.g., to students or people getting started with
     * MALA who just want to look around. It is not intended for production
     * calculations.
     * Compared to the LAMMPS implementation, this implementation has quite a
     * few limitations. Namely
     *
     * - It is roughly an order of magnitude slower for small systems
     *      and doesn't scale too great (more information on the optimization below)
     * - It only works for ONE chemical element
     * - It has no MPI or GPU support
     *
     * Some options are hardcoded in the same manner the LAMMPS implementation
     * hard codes them. Compared to the LAMMPS implementation, some
     * essentially never used options are not maintained/optimized.
     */
    private fun calculate(kwargs: Map<String, Any?>): Number {

        println("Using kotlin for descriptor calculation. The resulting calculation will be slow for large systems.")

        // The entire bispectrum calculation may be extensively profiled.
        val profileCalculation = kwargs["profile_calculation"] as Boolean? ?: false

        var timingDistances = 0
        var timingUi = 0
        var timingZi = 0
        var timingBi = 0
        var timingGridpoints = 0

        // Set up the array holding the bispectrum descriptors.
        val ncoeff = parameters.bispectrumTwojmax.let { (it + 2) * (it + 3) * (it + 4) }
        //        ncoeff = ncoeff // 24  # integer division
        fingerprintLength = ncoeff + 3
        val bispectrumNp =
            Array(gridDimensions[0]) {
                Array(gridDimensions[1]) {
                    Array(gridDimensions[2]) { FloatArray(fingerprintLength) }
                }
            }

        // Create a list of all potentially relevant atoms.
        val allAtoms = setupAtomList()

        // These are technically hyperparameters. We currently simply set them to set values for everything.
        rMin0 = 0f
        rFac0 = 0.99363f
        bZeroFlag = false
        wselfallFlag = false
        // Currently not supported
        bnormFlag = false
        // Currently not supported
        quadraticFlag = false
        numberElements = 1
        wSelf = 1f

        // What follows is the python implementation of the
        // bispectrum descriptor calculation.
        //
        // It was developed by first copying the code directly and
        // then optimizing it just enough to be usable. LAMMPS is
        // written in C++, and as such, many for-loops which are
        // optimized by the compiler can be employed. This is
        // drastically inefficient in python, so functions were
        // rewritten to use optimized vector-operations
        // (e.g. via numpy) where possible. This requires the
        // precomputation of quite a few index arrays. Thus,
        // this implementation is memory-intensive, which should
        // not be a problem given the intended use.
        //
        // There is still quite some optimization potential here.
        // I have decided to not optimized this code further just
        // now, since we do not know yet whether the bispectrum
        // descriptors will be used indefinitely, or if, e.g.
        // other types of neural networks will be used.
        // The implementation here is fast enough to be used for
        // tests of small test systems during development,
        // which is the sole purpose. If we eventually decide to
        // stick with bispectrum descriptors and feed-forward
        // neural networks, this code can be further optimized and
        // refined. I will leave some guidance below on what to
        // try/what has already been done, should someone else
        // want to give it a try.
        //
        // Final note: if we want to ship MALA with its own
        // bispectrum descriptor calculation to be used at scale,
        // the best way would potentially be via self-maintained
        // C++-functions.
        //
        //#######
        // Initialize index arrays.
        //
        // This function initializes a couple of lists of indices for
        // matrix multiplication/summation. By doing so, nested for-loops
        // can be avoided.
        //#######

        val tBegin = now
        initIndexArrays()
        val timingIndexInit = now - tBegin

        for (x in 0..<gridDimensions[0])
            for (y in 0..<gridDimensions[1])
                for (z in 0..<gridDimensions[2]) {
                    // Compute the grid point.
                    val tGrid = now
                    val coord = gridToCoord(intArrayOf(x, y, z))
                    for (c in 0..2)
                        bispectrumNp[x][y][z][c] = coord[c]

                    //#######
                    // Distance matrix calculation.
                    //
                    // Here, the distances to all atoms within our
                    // targeted cutoff are calculated.
                    //#######

                    val t0 = now
                    val distances = np.squeeze(
                        distance.cdist(
                            [bispectrumNp[x][y][z].cdist(allAtoms)
                    )
                    )
                }

        TODO()
    }

    /** Convert the units of a bispectrum descriptor.
     *
     *  Since these do not really have units this function does nothing yet.
     *
     *  Parameters
     *  ----------
     *  @param array: Data for which the units should be converted.
     *
     *  @param inUnits: Units of array.
     *
     *  @returns: converted_array, Data in MALA units. */
    override fun convertUnits(array: Number, inUnits: String?): Number {
        if (inUnits == null || inUnits == "None")
            return array
        else
            error("Unsupported unit for bispectrum descriptors.")
    }

    //#######
    // Functions and helper classes for calculating the bispectrum descriptors.
    //
    // The ZIndices and BIndices classes are useful stand-ins for structs used
    // in the original C++ code.
    //#######

    private class ZIndices {
        var j1 = 0
        var j2 = 0
        var j = 0
        var ma1min = 0
        var ma2max = 0
        var mb1min = 0
        var mb2max = 0
        var na = 0
        var nb = 0
        var jju = 0f
    }

    private class BIndices {
        var j1 = 0
        var j2 = 0
        var j = 0
    }

    /**
     *         Initialize index arrays.
     *
     *         This function initializes a couple of lists of indices for
     *         matrix multiplication/summation. By doing so, nested for-loops
     *         can be avoided.
     *
     *         FURTHER OPTIMIZATION: This function relies on nested for-loops.
     *         They may be optimized. I have not done so, because it is non-trivial
     *         in some cases and not really needed. These arrays are the same
     *         for each grid point, so the overall overhead is rather small.
     */
    private fun initIndexArrays() {
        // Needed for the Clebsch -Gordan product matrices(below)

        fun deltaCg(j1: Int, j2: Int, j: Int): Double {
            val sfaccg = ((j1 + j2 + j).floorDiv(2) + 1).factorial()
            return sqrt((j1 + j2 - j).floorDiv(2).factorial() *
                        (j1 - j2 + j).floorDiv(2).factorial() *
                        (-j1 + j2 + j).floorDiv(2).factorial() / sfaccg.toDouble())
        }

        //#######
        // Indices for compute_ui.
        //#######

        // First, the ones also used in LAMMPS.
        var idxuCount = 0
        indexUBlock = IntArray(parameters.bispectrumTwojmax + 1) {
            val res = idxuCount
            for (mb in 0..it)
                for (ma in 0..it)
                    idxuCount++
            res
        }
        indexUMax = idxuCount

        val rootPqArray = Array(parameters.bispectrumTwojmax + 2) { r ->
            FloatArray(parameters.bispectrumTwojmax + 2) { c ->
                if (r == 0 || c == 0) 0f
                else sqrt(r / c.toFloat())
            }
        }

        // These are only for optimization purposes.
        var indexUOneInitialized: IntArray? = null
        for (j in 0..parameters.bispectrumTwojmax) {
            val stop = when {
                j < parameters.bispectrumTwojmax -> indexUBlock[j + 1]
                else -> indexUMax
            }
            val range = indexUBlock[j]..stop step j + 2
            val iterator = range.iterator()
            val array = IntArray(range.count()) { iterator.nextInt() }
            indexUOneInitialized = if (indexUOneInitialized == null) array else indexUOneInitialized + array
        }

        for (j in 1..parameters.bispectrumTwojmax) {
            var jju = indexUBlock[j]
            var jjup = indexUBlock[j - 1]

            for (mb in 0..j.floorDiv(2)) {
                for (ma in 0..<j) {
                    rootpqFull1 += rootPqArray[j - ma][j - mb]
                    rootpqFull2 += rootPqArray[ma + 1][j - mb]
                    indexUFull += jju
                    indexU1Full += jjup
                    jju += 1
                    jjup += 1
                }
                jju += 1
            }

            var mbpar = 1
            jju = indexUBlock[j]
            jjup = jju + (j + 1) * (j + 1) - 1

            for (mb in 0..j.floorDiv(2)) {
                var mapar = mbpar
                for (ma in 0..j) {
                    if (mapar == 1) {
                        indexUSymmetryPos += jju
                        indexU1SymmetryPos += jjup
                    } else {
                        indexUSymmetryNeg += jju
                        indexU1SymmetryNeg += jjup
                    }
                    mapar = -mapar
                    jju += 1
                    jjup -= 1
                }
                mbpar = -mbpar
            }
        }

        //#######
        // Indices for compute_zi.
        //#######

        // First, the ones also used in LAMMPS.
        var idxzCount = 0
        for (j1 in 0..parameters.bispectrumTwojmax)
            for (j2 in 0..j1)
                for (j in (j1 - j2)..min(parameters.bispectrumTwojmax, j1 + j2) step 2)
                    for (mb in 0..j.floorDiv(2))
                        for (ma in 0..j)
                            idxzCount++

        val idxzMax = idxzCount
        val idxz = Array(idxzMax) { ZIndices() }
        val size = parameters.bispectrumTwojmax + 1
        indexZBlock = Array(size) { Array(size) { IntArray(size) } }

        idxzCount = 0
        for (j1 in 0..parameters.bispectrumTwojmax)
            for (j2 in 0..j1)
                for (j in (j1 - j2)..min(parameters.bispectrumTwojmax, j1 + j2) step 2) {
                    indexZBlock[j1][j2][j] = idxzCount

                    for (mb in 0..j.floorDiv(2))
                        for (ma in 0..j) {
                            val id = idxz[idxzCount]
                            id.j1 = j1
                            id.j2 = j2
                            id.j = j
                            id.ma1min = max(0, (2 * ma - j - j2 + j1).floorDiv(2))
                            id.ma2max = (2 * ma - j - (2 * id.ma1min - j1) + j2).floorDiv(2)
                            id.na = min(j1, (2 * ma - j + j2 + j1).floorDiv(2)) - id.ma1min + 1
                            id.mb1min = max(0, (2 * mb - j - j2 + j1).floorDiv(2))
                            id.mb2max = (2 * mb - j - (2 * id.mb1min - j1) + j2).floorDiv(2)
                            id.nb = min(j1, (2 * mb - j + j2 + j1).floorDiv(2)) - id.mb1min + 1

                            id.jju = (indexUBlock[j] + (j + 1) * mb + ma).toFloat()

                            idxzCount++
                        }
                }

        val idxcgBlock = Array(size) { Array(size) { IntArray(size) } }
        var idxcgCount = 0
        for (j1 in 0..parameters.bispectrumTwojmax)
            for (j2 in 0..j1)
                for (j in (j1 - j2)..min(parameters.bispectrumTwojmax, j1 + j2) step 2) {
                    idxcgBlock[j1][j2][j] = idxcgCount
                    for (m1 in 0..j1)
                        for (m2 in 0..j2)
                            idxcgCount++
                }

        cglist = FloatArray(idxcgCount)

        idxcgCount = 0
        for (j1 in 0..parameters.bispectrumTwojmax)
            for (j2 in 0..j1)
                for (j in (j1 - j2)..min(parameters.bispectrumTwojmax, j1 + j2) step 2)
                    for (m1 in 0..j1) {
                        val aa2 = 2 * m1 - j1
                        for (m2 in 0..j2) {
                            val bb2 = 2 * m2 - j2
                            val m = (aa2 + bb2 + j).floorDiv(2)
                            if (m < 0 || m > j) {
                                cglist[idxcgCount] = 0f
                                idxcgCount++
                                continue
                            }
                            var cgsum = 0f
                            val x = max(-(j - j2 + aa2).floorDiv(2), -(j - j1 - bb2).floorDiv(2))
                            val y = min((j1 - aa2).floorDiv(2), (j2 + bb2).floorDiv(2))
                            for (z in max(0, x)..min((j1 + j2 - j).floorDiv(2), y)) {
                                val ifac = if (z % 2 != 0) -1 else 1
                                cgsum += ifac / (z.factorial() *
                                                 ((j1 + j2 - j).floorDiv(2) - z).factorial() *
                                                 ((j1 - aa2).floorDiv(2) - z).factorial() *
                                                 ((j2 + bb2).floorDiv(2) - z).factorial() *
                                                 ((j - j2 + aa2).floorDiv(2) + z).factorial() *
                                                 ((j - j1 - bb2).floorDiv(2) + z).factorial())
                            }
                            val cc2 = 2 * m - j
                            val dcg = deltaCg(j1, j2, j)
                            val sfaccg = sqrt(((j1 + aa2).floorDiv(2).factorial() *
                                               (j1 - aa2).floorDiv(2).factorial() *
                                               (j2 + bb2).floorDiv(2).factorial() *
                                               (j2 - bb2).floorDiv(2).factorial() *
                                               (j + cc2).floorDiv(2).factorial() *
                                               (j - cc2).floorDiv(2).factorial() * (j + 1)).toDouble())
                            cglist[idxcgCount] = (cgsum * dcg * sfaccg).toFloat()
                            idxcgCount++
                        }
                    }

        // These are only for optimization purposes.
        //        self.__index_z_icga = []
        //        self.__index_z_icgb = []
        for (jjz in 0..<idxzMax) {
            val id = idxz[jjz]
            val j1 = id.j1
            val j2 = id.j2
            val j = id.j
            val ma1min = id.ma1min
            val ma2max = id.ma2max
            val na = id.na
            val mb1min = id.mb1min
            val mb2max = id.mb2max
            val nb = id.nb
            var jju1 = indexUBlock[j1] + (j1 + 1) * mb1min
            var jju2 = indexUBlock[j2] + (j2 + 1) * mb2max

            var icgb = mb1min * (j2 + 1) + mb2max
            for (ib in 0..<nb) {
                var ma1 = ma1min
                var ma2 = ma2max
                var icga = ma1min * (j2 + 1) + ma2max
                for (ia in 0..<na) {
                    indexZJjz += jjz
                    indexZIcgb += idxcgBlock[j1][j2][j] + icgb
                    indexZIcga += idxcgBlock[j1][j2][j] + icga
                    indexZU1r += jju1 + ma1
                    indexZU1i += jju1 + ma1
                    indexZU2r += jju2 + ma2
                    indexZU2i += jju2 + ma2
                    ma1 += 1
                    ma2 -= 1
                    icga += j2
                }
                jju1 += j1 + 1
                jju2 -= j2 + 1
                icgb += j2
            }
        }

        //#######
        // Indices for compute_bi.
        //#######

        // These are identical to LAMMPS, because we do not optimize compute_bi.
        var idxbCount = 0
        for (j1 in 0..parameters.bispectrumTwojmax)
            for (j2 in 0..j1)
                for (j in (j1 - j2)..min(parameters.bispectrumTwojmax, j1 + j2) step 2)
                    if (j >= j1)
                        idxbCount++

        indexBMax = idxbCount
        val indexB = Array(indexBMax) { BIndices() }

        idxbCount = 0
        for (j1 in 0..parameters.bispectrumTwojmax)
            for (j2 in 0..j1)
                for (j in (j1 - j2)..min(parameters.bispectrumTwojmax, j1 + j2) step 2)
                    if (j >= j1) {
                        indexB[idxbCount].j1 = j1
                        indexB[idxbCount].j2 = j2
                        indexB[idxbCount].j = j
                        idxbCount++
                    }
    }
}