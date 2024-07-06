package casus.mala.descriptors

import ase.*
import ase.io.prod
import casus.mala.common.Parameters
import casus.mala.dataHandling.div
import java.io.File
import kotlin.math.*

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

    private var indexUOneInitialized: IntArray? = null
    private val indexUFull = ArrayList<Int>()
    private val indexUSymmetryPos = ArrayList<Int>()
    private val indexUSymmetryNeg = ArrayList<Int>()
    private val indexU1Full = ArrayList<Int>()
    private val indexU1SymmetryPos = ArrayList<Int>()
    private val indexU1SymmetryNeg = ArrayList<Int>()
    private val rootpqFull1 = ArrayList<Float>()
    private val rootpqFull2 = ArrayList<Float>()

    private lateinit var indexZU1r: IntArray
    private lateinit var indexZU1i: IntArray
    private lateinit var indexZU2r: IntArray
    private lateinit var indexZU2i: IntArray
    private lateinit var indexZIcga: IntArray
    private lateinit var indexZIcgb: IntArray
    private lateinit var indexZJjz: IntArray
    private lateinit var indexZBlock: Array<Array<IntArray>>

    private var indexBMax = 0

    private lateinit var indexB: Array<BIndices>
    var rMin0 = 0f
    var rFac0 = 0f
    var bZeroFlag = false
    var wselfallFlag = false
    var bnormFlag = false
    var quadraticFlag = false
    var numberElements = 0
    var wSelf = 0f

    //    fun convertUnits(array: NDArray): Nothing = TODO() // convertUnits (array, null)

    override fun calculate(outDir: File, kwargs: Map<String, Any>): Number = when {
        //        else -> calculateKotlin(kwargs)
        else -> calculateLammps(outDir, kwargs) as Number
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
     *         Perform bispectrum calculation using LAMMPS.
     *
     *         Creates a LAMMPS instance with appropriate call parameters and uses
     *         it for the calculation.
     */
    private fun calculateLammps(outdir: File, kwargs: Map<String, Any>): Number {
        // For version compatibility; older lammps versions(the serial version
        // we still use on some machines) have these constants as part of the
        // general LAMMPS import.
        val useFp64 = kwargs["use_fp64"] as Boolean? ?: false
        val keepLogs = kwargs["keep_logs"] as Boolean? ?: false

        val lammpsFormat = "lammps-data"
        lammpsTemporaryInput = File.createTempFile("lammps_input_", ".tmp")

        ase.io.write(lammpsTemporaryInput!!, atoms, lammpsFormat)

        val (nx, ny, nz) = gridDimensions

        // Create LAMMPS instance.
        val lammpsDict: MutableMap<String, Any> = mutableMapOf("twojmax" to parameters.bispectrumTwojmax,
                                                               "rcutfac" to parameters.bispectrumCutoff)

        lammpsTemporaryLog = File.createTempFile("lammps_bgrid_log_", ".tmp")

        val lmp = setupLammps(nx, ny, nz, lammpsDict)

        // An empty string means that the user wants to use the standard input.
        // What that is differs depending on serial/parallel execution.
        if (parameters.lammpsComputeFile == null)
            parameters.lammpsComputeFile = when {
                parameters.configuration.mpi -> TODO()
//            if self.parameters.use_z_splitting:
//            self.parameters.lammps_compute_file = os.path.join(
//                filepath, "in.bgridlocal.python"
//                                                              )
//            else:
//            self.parameters.lammps_compute_file = os.path.join(
//                filepath, "in.bgridlocal_defaultproc.python"
//                                                              )
            else -> File(this::class.java.getResource("in.bgrid.python")!!.toURI())
        }

        // Do the LAMMPS calculation and clean up.
        lmp.file(self.parameters.lammps_compute_file)
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
    private fun calculateKotlin(kwargs: Map<String, Any?>): Number {

        println("Using kotlin for descriptor calculation. The resulting calculation will be slow for large systems.")

        // The entire bispectrum calculation may be extensively profiled.
        val profileCalculation = true //kwargs["profile_calculation"] as Boolean? ?: false

        var timingDistances = 0f
        var timingUi = 0f
        var timingZi = 0f
        var timingBi = 0f
        var timingGridpoints = 0f

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
                    val tGrid = if (profileCalculation) now else 0f
                    val coord = gridToCoord(intArrayOf(x, y, z))
                    for (c in 0..2)
                        bispectrumNp[x][y][z][c] = coord[c]

                    //#######
                    // Distance matrix calculation.
                    //
                    // Here, the distances to all atoms within our
                    // targeted cutoff are calculated.
                    //#######

                    var t0 = if (profileCalculation) now else 0f
                    val distances = bispectrumNp[x][y][z] cdist allAtoms
                    val indices = buildList { for (i in distances.indices) if (distances[i] < parameters.bispectrumCutoff) add(i) }
                    val distancesCutoff = FloatArray(indices.size) { abs(distances[indices[it]]) }
                    val atomsCutoff = Array(indices.size) { allAtoms[indices[it]] }
                    val nrAtoms = atomsCutoff.size
                    if (profileCalculation) timingDistances += now - t0

                    //#######
                    // Compute ui.
                    //
                    // This calculates the expansion coefficients of the
                    // hyperspherical harmonics (usually referred to as ui).
                    //#######

                    if (profileCalculation) t0 = now
                    val (ulisttotR, ulisttotI) = computeUi(nrAtoms,
                                                           atomsCutoff,
                                                           distancesCutoff,
                                                           bispectrumNp[x][y][z])
                    if (profileCalculation) timingUi += now - t0

                    //#######
                    // Compute zi.
                    //
                    // This calculates the bispectrum components through
                    // triple scalar products/Clebsch-Gordan products.
                    //#######

                    if (profileCalculation) t0 = now
                    val (zlistR, zlistI) = computeZi(ulisttotR, ulisttotI)
                    if (profileCalculation) timingZi += now - t0

                    //#######
                    // Compute the bispectrum descriptors itself.
                    //
                    // This essentially just extracts the descriptors from
                    // the expansion coeffcients.
                    //#######
                    if (profileCalculation) t0 = now
                    computeBi(ulisttotR, ulisttotI, zlistR, zlistI).copyInto(bispectrumNp[x][y][z], destinationOffset = 3)
                    if (profileCalculation) {
                        timingGridpoints += now - tGrid
                        timingBi += now - t0
                    }
                }

        if (profileCalculation) {
            val timingTotal = now - tBegin
            val prod = gridDimensions.prod()
            println("""
                Python-based bispectrum descriptor calculation timing: 
                Index matrix initialization $timingIndexInit]
                Overall calculation time $timingTotal
                Calculation time per gridpoint [s/gridpoint] ${timingGridpoints / prod}
                Timing contributions per gridpoint:
                Distance matrix [s/gridpoint] ${timingDistances / prod}
                Compute ui [s/gridpoint] ${timingUi / prod}
                Compute zi [s/gridpoint] ${timingZi / prod}
                Compute bi [s/gridpoint] ${timingBi / prod}""".trimIndent())
        }
        return 3
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
        indexUOneInitialized = null
        for (j in 0..parameters.bispectrumTwojmax) {
            val stop = when {
                j < parameters.bispectrumTwojmax -> indexUBlock[j + 1]
                else -> indexUMax
            }
            val range = indexUBlock[j]..stop step j + 2
            val iterator = range.iterator()
            val array = IntArray(range.count()) { iterator.nextInt() }
            indexUOneInitialized = indexUOneInitialized?.let { it + array } ?: array
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
                    jju++
                    jjup++
                }
                jju++
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
                    jju++
                    jjup--
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
                                val ifac = if (z % 2 != 0) -1f else 1f
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
        val indexZU1r = ArrayList<Int>()
        val indexZU1i = ArrayList<Int>()
        val indexZU2r = ArrayList<Int>()
        val indexZU2i = ArrayList<Int>()
        val indexZIcga = ArrayList<Int>()
        val indexZIcgb = ArrayList<Int>()
        val indexZJjz = ArrayList<Int>()
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
                    ma1++
                    ma2--
                    icga += j2
                }
                jju1 += j1 + 1
                jju2 -= j2 + 1
                icgb += j2
            }
        }
        this.indexZU1r = indexZU1r.toIntArray()
        this.indexZU1i = indexZU1i.toIntArray()
        this.indexZU2r = indexZU2r.toIntArray()
        this.indexZU2i = indexZU2i.toIntArray()
        this.indexZIcga = indexZIcga.toIntArray()
        this.indexZIcgb = indexZIcgb.toIntArray()
        this.indexZJjz = indexZJjz.toIntArray()

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
        indexB = Array(indexBMax) { BIndices() }

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

    /**
     *         Compute ui.
     *
     *         This calculates the expansion coefficients of the
     *         hyperspherical harmonics (usually referred to as ui).
     *
     *         FURTHER OPTIMIZATION: This originally was a huge nested for-loop.
     *         By vectorizing over the atoms and pre-initializing a bunch of arrays,
     *         a  massive amount of time could be saved. There is one principal
     *         for-loop remaining - I have not found an easy way to optimize it out.
     *         Also, I have not tried numba or some other just-in-time compilation,
     *         may help.
     */
    private fun computeUi(nrAtoms: Int, atomsCutoff: Array<FloatArray>,
                          distancesCutoff: FloatArray, grid: FloatArray): Pair<FloatArray, FloatArray> {

        // Precompute and prepare ui stuff
        val theta0 = FloatArray(nrAtoms) {
            (distancesCutoff[it] - rMin0) * rFac0 * Math.PI.toFloat() / (parameters.bispectrumCutoff - rMin0)
        }
        val z0 = FloatArray(nrAtoms) { distancesCutoff[it] / tan(theta0[it]) }

        val ulistRIj = Array(nrAtoms) { FloatArray(indexUMax) { if (it == 0) 1f else 0f } }
        val ulistIIj = Array(nrAtoms) { FloatArray(indexUMax) }
        val ulisttotR = FloatArray(indexUMax)
        val ulisttotI = FloatArray(indexUMax)
        val r0inv = FloatArray(nrAtoms) { 1f / sqrt(distancesCutoff[it] * distancesCutoff[it] + z0[it] * z0[it]) }
        for (i in indexUOneInitialized!!) ulisttotR[i] = 1f
        val distanceVector = Array(nrAtoms) { r -> FloatArray(3) { -1f * (atomsCutoff[r][it] - grid[it]) } }

        // Cayley-Klein parameters for unit quaternion.
        if (nrAtoms > 0) {
            val aR = FloatArray(nrAtoms) { r0inv[it] * z0[it] }
            val aI = FloatArray(nrAtoms) { -r0inv[it] * distanceVector[it][2] }
            val bR = FloatArray(nrAtoms) { r0inv[it] * distanceVector[it][1] }
            val bI = FloatArray(nrAtoms) { -r0inv[it] * distanceVector[it][0] }

            // This encapsulates the compute_uarray function
            var jju1 = 0
            var jju2 = 0
            var jju3 = 0
            for (jjuOuter in 0..<indexUMax) {
                if (jjuOuter in indexUFull) {
                    for (i in 0..<nrAtoms) {
                        var rootpq = rootpqFull1[jju1]
                        val listR = ulistRIj[i]
                        val listI = ulistIIj[i]
                        val u = indexUFull[jju1]
                        val u1 = indexU1Full[jju1]
                        listR[u] += rootpq * (aR[i] * listR[u1] + aI[i] * listI[u1])
                        listI[u] += rootpq * (aR[i] * listI[u1] - aI[i] * listR[u1])

                        rootpq = rootpqFull2[jju1]
                        listR[u + 1] = -1f * rootpq * (bR[i] * listR[u1] + bI[i] * listI[u1])
                        listI[u + 1] = -1f * rootpq * (bR[i] * listI[u1] - bI[i] * listR[u1])
                    }
                    jju1++
                }
                if (jjuOuter in indexU1SymmetryPos) {
                    for (i in 0..<nrAtoms) {
                        val listR = ulistRIj[i]
                        val listI = ulistIIj[i]
                        val u = indexUSymmetryPos[jju2]
                        val u1 = indexU1SymmetryPos[jju2]
                        listR[u1] = listR[u]
                        listI[u1] = -listI[u]
                    }
                    jju2++
                }
                if (jjuOuter in indexU1SymmetryNeg) {
                    for (i in 0..<nrAtoms) {
                        val listR = ulistRIj[i]
                        val listI = ulistIIj[i]
                        val u = indexUSymmetryNeg[jju3]
                        val u1 = indexU1SymmetryNeg[jju3]
                        listR[u1] = -listR[u]
                        listI[u1] = listI[u]
                    }
                    jju3++
                }
            }

            // This emulates add_uarraytot.
            // First, we compute sfac.
            val sfac = when {
                parameters.bispectrumSwitchflag == 0 -> FloatArray(nrAtoms) { 1f }
                else -> when {
                    nrAtoms > 1 -> {
                        val rcutfac = Math.PI.toFloat() / (parameters.bispectrumCutoff - rMin0)
                        FloatArray(nrAtoms) {
                            when {
                                distancesCutoff[it] <= rMin0 -> 1f
                                distancesCutoff[it] > parameters.bispectrumCutoff -> 0f
                                else -> 0.5f * (cos((distancesCutoff[it] - rMin0) * rcutfac) + 1f)
                            }
                        }
                    }
                    else -> FloatArray(nrAtoms)
                }
            }

            // sfac technically has to be weighted according to the chemical
            // species. But this is a minimal implementation only for a single
            // chemical species, so I am ommitting this for now. It would
            // look something like
            // sfac *= weights[a]
            // Further, some things have to be calculated if
            // switch_inner_flag is true. If I understand correctly, it
            // essentially never is in our case. So I am ommitting this
            // (along with some other similar lines) here for now.
            // If this becomes relevant later, we of course have to
            // add it.

            // Now use sfac for computations.
            for (jju in 0..<indexUMax) {
                ulisttotR[jju] += sfac.sumIndexed { i, f -> f * ulistRIj[i][jju] }
                ulisttotI[jju] += sfac.sumIndexed { i, f -> f * ulistIIj[i][jju] }
            }
            println()
        }
        return ulisttotR to ulisttotI
    }

    /**
     *         Compute zi.
     *
     *         This calculates the bispectrum components through
     *         triple scalar products/Clebsch-Gordan products.
     *
     *         FURTHER OPTIMIZATION: In the original code, this is a huge nested
     *         for-loop. Even after optimization, this is the principal
     *         computational cost (for realistic systems). I have found this
     *         implementation to be the most efficient without any major refactoring.
     *         However, due to the usage of np.unique, numba cannot trivially be used.
     *         A different route that then may employ just-in-time compilation
     *         could be fruitful.
     */
    fun computeZi(ulisttotR: FloatArray, ulisttotI: FloatArray): Pair<FloatArray, FloatArray> {
        val tmpReal = FloatArray(indexZIcgb.size) {
            cglist[indexZIcgb[it]] * cglist[indexZIcga[it]] *
                    (ulisttotR[indexZU1r[it]] * ulisttotR[indexZU2r[it]] -
                            ulisttotI[indexZU1i[it]] * ulisttotI[indexZU2i[it]])
        }
        val tmpImag = FloatArray(indexZIcgb.size) {
            cglist[indexZIcgb[it]] * cglist[indexZIcga[it]] *
                    (ulisttotR[indexZU1r[it]] * ulisttotI[indexZU2i[it]] +
                            ulisttotI[indexZU1i[it]] * ulisttotR[indexZU2r[it]])
        }

        // Summation over an array based on indices stored in a different
        // array.
        // Taken from: https://stackoverflow.com/questions/67108215/how-to-get-sum-of-values-in-a-numpy-array-based-on-another-array-with-repetitive
        // Under "much better version".
        val zlistR = tmpReal.sumBy(indexZJjz)
        val zlistI = tmpImag.sumBy(indexZJjz)
        return zlistR to zlistI
    }

    /**
     *         Compute the bispectrum descriptors itself.
     *
     *         This essentially just extracts the descriptors from
     *         the expansion coeffcients.
     *
     *         FURTHER OPTIMIZATION: I have not optimized this function AT ALL.
     *         This is due to the fact that its computational footprint is miniscule
     *         compared to the other parts of the bispectrum descriptor calculation.
     *         It contains multiple for-loops, that may be optimized out.
     */
    fun computeBi(ulisttotR: FloatArray, ulisttotI: FloatArray, zlistR: FloatArray, zlistI: FloatArray): FloatArray {
        // For now set the number of elements to 1.
        // This also has some implications for the rest of the function.
        // This currently really only works for one element .
        val nElements = 1
        val nElementPairs = nElements * nElements
        val nElementTriples = nElementPairs * nElements
        var iElem = 0
        val bList = FloatArray(indexBMax * nElementTriples)
        var iTriple = 0
        var iDouble = 0

        if (bZeroFlag) {
            TODO()
            //            var wself = 1.0
            //            val bZero = np.zeros(self.parameters.bispectrum_twojmax + 1)
            //            www = wself * wself * wself
            //            for j in range(self.parameters.bispectrum_twojmax + 1):
            //            if self.bnorm_flag:
            //            bZero[j] = www
            //            else:
            //            bZero[j] = www * (j + 1)
        }

        for (elem1 in 0..<nElements)
            for (elem2 in 0..<nElements) {
                for (elem3 in 0..<nElements) {
                    for (jjb in 0..<indexBMax) {
                        val b = indexB[jjb]
                        val j1 = b.j1
                        val j2 = b.j2
                        val j = b.j
                        var jjz = indexZBlock[j1][j2][j]
                        var jju = indexUBlock[j]
                        var sumZu = 0f
                        for (mb in 0..<ceil((j / 2).toFloat()).toInt())
                            for (ma in 0..j) {
                                sumZu += ulisttotR[elem3 * indexUMax + jju] * zlistR[jjz] +
                                        ulisttotI[elem3 * indexUMax + jju] * zlistI[jjz]
                                jjz++
                                jju++
                            }
                        if (j % 2 == 0) {
                            val mb = j.floorDiv(2)
                            for (ma in 0..<mb) {
                                sumZu += ulisttotR[elem3 * indexUMax + jju] * zlistR[jjz] +
                                        ulisttotI[elem3 * indexUMax + jju] * zlistI[jjz]
                                jjz++
                                jju++
                            }
                            sumZu += 0.5f * (ulisttotR[elem3 * indexUMax + jju] * zlistR[jjz] +
                                    ulisttotI[elem3 * indexUMax + jju] * zlistI[jjz])
                        }
                        bList[iTriple * indexBMax + jjb] = 2f * sumZu
                    }
                    iTriple++
                }
                iDouble++
            }

        if (bZeroFlag) {
            TODO()
            //            if not self.wselfall_flag:
            //            itriple = (
            //                    ielem * number_elements + ielem
            //                    ) * number_elements + ielem
            //            for jjb in range(self.__index_b_max):
            //            j = self.__index_b[jjb].j
            //            blist[itriple * self.__index_b_max + jjb] -= bzero[j]
            //            else:
            //            itriple = 0
            //            for elem1 in range(number_elements):
            //            for elem2 in range(number_elements):
            //            for elem3 in range(number_elements):
            //            for jjb in range(self.__index_b_max):
            //            j = self.__index_b[jjb].j
            //            blist[
            //                itriple * self.__index_b_max + jjb
            //            ] -= bzero[j]
            //            itriple += 1
        }

        // Untested  & Unoptimized
        if (quadraticFlag) {
            TODO()
            //            xyz_length = 3 if self.parameters.descriptors_contain_xyz else 0
            //            ncount = self.fingerprint_length - xyz_length
            //            for icoeff in range(ncount):
            //            bveci = blist[icoeff]
            //            blist[3 + ncount] = 0.5 * bveci * bveci
            //            ncount += 1
            //            for jcoeff in range(icoeff + 1, ncount):
            //            blist[xyz_length + ncount] = bveci * blist[jcoeff]
            //            ncount += 1
        }

        return bList
    }
}