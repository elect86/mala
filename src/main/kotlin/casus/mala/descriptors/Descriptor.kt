package casus.mala.descriptors

import ai.djl.ndarray.NDArray
import ase.*
import ase.io.div
import casus.mala.common.Parameters
import casus.mala.common.PhysicalData
import lammps.Lammps
import skspatial.objects.Plane
import skspatial.objects.Point
import java.io.File
import java.time.LocalDateTime
import kotlin.properties.Delegates

/** Base class for all descriptors available in MALA.
 *
 *  Descriptors encode the atomic fingerprint of a DFT calculation.
 *  @param parameters: Parameters object used to create this object
 */
abstract class Descriptor(parameters: Parameters) : PhysicalData(parameters) {

    val parameters: Parameters.Descriptors = parameters.descriptors
    var fingerprintLength = 0  // so iterations will fail

    //    self.verbosity = parameters.verbosity
    var inFormatAse = ""
    lateinit var atoms: Atoms
    var voxel: Cell? = null

    // If we ever have NON LAMMPS descriptors, these parameters have no
    // meaning anymore and should probably be moved to an intermediate
    // DescriptorsLAMMPS class, from which the LAMMPS descriptors inherit.
    var lammpsTemporaryInput: File? = null
    var lammpsTemporaryLog: File? = null

    //#############################
    // Properties
    //#############################

    //    val siUnitConversion
    //    """
    //        Numeric value of the conversion from MALA (ASE) units to SI.
    //
    //        Needed for OpenPMD interface.
    //        """
    //    return m**3

    //    @property
    //    def si_dimension(self):
    //    """
    //        Dictionary containing the SI unit dimensions in OpenPMD format.
    //
    //        Needed for OpenPMD interface.
    //        """
    //    import openpmd_api as io
    //
    //    return {io.Unit_Dimension.L: -3}

    /** Control whether descriptor vectors will contain xyz coordinates. */
    var descriptorsContainXyz: Boolean
        get() = parameters.descriptorsContainXyz
        set(value) {
            parameters.descriptorsContainXyz = value
        }

    /**
     *         Timestamp of calculation start.
     *
     *         Used to distinguish multiple LAMMPS runs performed in the same
     *         directory. Since the interface is file based, this timestamp prevents
     *         problems with slightly
     */
    val calculationTimestamp: String
        get() = LocalDateTime.now().toString()

    /**
     * Explictly enforces the PBC on an ASE atoms object.
     *
     * QE (and potentially other codes?) do that internally. Meaning that the
     * raw positions of atoms (in Angstrom) can lie outside of the unit cell.
     * When setting up the DFT calculation, these atoms get shifted into
     * the unit cell. Since we directly use these raw positions for the
     * descriptor calculation, we need to enforce that in the ASE atoms
     * objects, the atoms are explicitly in the unit cell.
     *
     *         Parameters
     *         ----------
     *         atoms : ase.atoms
     *             The ASE atoms object for which the PBC need to be enforced.
     *
     *         Returns
     *         -------
     *         new_atoms : ase.Atoms
     *             The ASE atoms object for which the PBC have been enforced.
     */
    fun enforcePbc(atoms: Atoms): Atoms {
        val newAtoms = atoms.copy()
        newAtoms.setScaledPositions(newAtoms.getScaledPositions())

        // This might be unecessary, but I think it is nice to have some sort of metric here.
        var rescaledAtoms = 0
        for (atom in atoms) {
            if (false in atom.position.isClose(atom.position, atol = 0.001f))
                rescaledAtoms++
        }
        print("Descriptor calculation: had to enforce periodic boundary conditions on $rescaledAtoms atoms " +
              "before calculation." /*min_verbosity = 2*/)
        return newAtoms
    }

    /**
     * Calculate the descriptors based on a Quantum Espresso outfile.
     *
     * Parameters
     * ----------
     * @param qeOutFile Quantum Espresso output file for snapshot.
     * @param workingDirectory A directory in which to write the output of the LAMMPS calculation.
     * Usually the local directory should suffice, given that there are no multiple instances running in the same directory.
     *
     * @return descriptors : numpy.array
     *             Numpy array containing the descriptors with the dimension
     *             (x,y,z,descriptor_dimension)
     */
    fun calculateFromQeOut(qeOutFile: File, workingDirectory: File = File("."), kwargs: MutableMap<String, Any>): Any {
        inFormatAse = "espresso-out"
        println("Calculating descriptors from $qeOutFile"/*min_verbosity=0*/)
        // We get the atomic information by using ASE.
        atoms = ase.read(qeOutFile, format = inFormatAse)

        // Enforcing / Checking PBC on the read atoms .
        atoms = enforcePbc(atoms)

        // Get the grid dimensions.
        if ("grid_dimensions" in kwargs) {
            gridDimensions = kwargs["grid_dimensions"]!! as IntArray

            // Deleting this keyword from the list to avoid conflict with dict below.
            kwargs -= "grid_dimensions"
        } else {
            val lines = qeOutFile.readLines()
            gridDimensions = IntArray(3)

            for (line in lines)
                if ("FFT dimensions" in line) {
                    val tmp = line.substringAfter('(').substringBefore(')').split(',')
                    gridDimensions = IntArray(3) { tmp[it].trim().toInt() }
                    break
                }
        }

        voxel = atoms.cell!!.copy().apply {
            array[0] = array[0] / gridDimensions[0].toFloat()
            array[1] = array[1] / gridDimensions[1].toFloat()
            array[2] = array[2] / gridDimensions[2].toFloat()
        }
        return calculate(workingDirectory, kwargs)
    }

    //#############################
    // Methods
    //#############################

    // File I/O
    //#########

    abstract fun convertUnits(array: Number, inUnits: String? = "1/eV"): Number

    // Private methods
    //################

    fun processLoadedArray(array: NDArray) = processLoadedArray(array, null)
    override fun processLoadedArray(array: NDArray, units: String?) {
        array *= convertUnits(1, units)
    }

    //    def _process_loaded_dimensions(self, array_dimensions):
    //    if self.descriptors_contain_xyz:
    //    return (
    //    array_dimensions[0],
    //    array_dimensions[1],
    //    array_dimensions[2],
    //    array_dimensions[3] - 3,
    //    )
    //    else:
    //    return array_dimensions
    //
    //    def _set_geometry_info(self, mesh):
    //    # Geometry: Save the cell parameters and angles of the grid.
    //    if self.atoms is not None:
    //    import openpmd_api as io
    //
    //    self.voxel = self.atoms.cell.copy()
    //    self.voxel[0] = self.voxel[0] / (self.grid_dimensions[0])
    //    self.voxel[1] = self.voxel[1] / (self.grid_dimensions[1])
    //    self.voxel[2] = self.voxel[2] / (self.grid_dimensions[2])
    //
    //    mesh.geometry = io.Geometry.cartesian
    //    mesh.grid_spacing = self.voxel.cellpar()[0:3]
    //    mesh.set_attribute("angles", self.voxel.cellpar()[3:])
    //
    //    def _get_atoms(self):
    //    return self.atoms

    override val featureMask: Int
        get() = if (descriptorsContainXyz) 3 else 0

    /**
     *         Set up the lammps processor grid.
     *
     *         Takes into account y/z-splitting.
     */
    protected fun setupLammps(nx: Int, ny: Int, nz: Int, lammpsDict: MutableMap<String, Any>): Lammps {

        // Build LAMMPS arguments from the data we read.
        var lmpCmdargs = listOf("-screen", "none", "-log", lammpsTemporaryLog!!.absolutePath)

        lammpsDict["atom_config_fname"] = lammpsTemporaryInput!!

        var size by Delegates.notNull<Int>()
        if (parameters.configuration.mpi) {
            TODO()
        } else {
            size = 1
            lammpsDict["ngridx"] = nx
            lammpsDict["ngridy"] = ny
            lammpsDict["ngridz"] = nz
            lammpsDict["switch"] = parameters.bispectrumSwitchflag
        }

        if (parameters.configuration.gpu) {
            TODO()
//            # Tell Kokkos to use one GPU.
//            lmp_cmdargs.append("-k")
//            lmp_cmdargs.append("on")
//            lmp_cmdargs.append("g")
//            lmp_cmdargs.append(str(size))
//
//            # Tell LAMMPS to use Kokkos versions of those commands for
//            # which a Kokkos version exists.
//            lmp_cmdargs.append("-sf")
//            lmp_cmdargs.append("kk")
//            pass
        }

        lmpCmdargs = setCmdlineVars(lmpCmdargs, lammpsDict)
        val lmp = Lammps(cmdArgs = lmpCmdargs)

//        set_lammps_instance(lmp)
        return lmp
    }


    /**
     * Set up a list of atoms potentially relevant for descriptor calculation.
     *
     * If periodic boundary conditions are used, which is usually the case
     * for MALA simulation, one has to compute descriptors by also
     * incorporating atoms from neighboring cells.
     *
     * FURTHER OPTIMIZATION: Probably not that much, this mostly already uses
     * optimized python functions.
     */
    protected fun setupAtomList(): Array<FloatArray> {
        if (atoms.pbc.none { it })
        // If no PBC are used, only consider a single cell.
            return atoms.positions!!

        // To determine the list of relevant atoms we first take the edges
        // of the simulation cell and use them to determine all cells
        // which hold atoms that _may_ be relevant for the calculation.
        val edges = arrayOf(intArrayOf(0, 0, 0),
                            intArrayOf(1, 0, 0),
                            intArrayOf(0, 1, 0),
                            intArrayOf(0, 0, 1),
                            intArrayOf(1, 1, 1),
                            intArrayOf(0, 1, 1),
                            intArrayOf(1, 0, 1),
                            intArrayOf(1, 1, 0)) * gridDimensions
        var allCellsList: List<IntArray>? = null

        // For each edge point create a neighborhoodlist to all cells given by the cutoff radius.
        for (edge in edges) {
            val edgePoint = gridToCoord(edge)
            val neighborList = NeighborList(FloatArray(atoms.len + 1) { parameters.atomicDensityCutoff!! },
                                            bothWays = true, selfInteraction = false/*, primitive = NewPrimitiveNeighborList*/)

            val atomsWithGridPoint = atoms.copy()

            // Construct a ghost atom representing the grid point .
            atomsWithGridPoint += Atom("H", edgePoint)
            neighborList.update(atomsWithGridPoint)

            val (indices, offsets) = neighborList getNeighbors atoms.len

            // Incrementally fill the list containing all cells to be considered.
            allCellsList = when (allCellsList) {
                null -> offsets.distinctBy(selector).sortedBy(selector)
                else -> allCellsList + offsets.distinctBy(selector).sortedBy(selector)
            }
        }

        // Delete the original cell from the list of all cells.
        // This is to avoid double checking of atoms below.
        val cells = this::class.java.getResource("allCellsList")!!
                .readText().split(',').map(String::toInt)
        allCellsList = List(cells.size / 3) { i -> IntArray(3) { cells[i * 3 + it] } }
        var allCells = allCellsList.distinctBy(selector).sortedBy(selector)
        var idx = 0
        for (a in allCells.indices) {
            if (allCells[a].all { it == 0 })
                break
            idx++
        }
        allCells = allCells.filterIndexed { index, _ -> index != idx }

        // Create an object to hold all relevant atoms.
        // First, instantiate it by filling it will all atoms from all potentiall relevant cells, as identified above.
        var allAtoms: Array<FloatArray>? = null
        for (a in atoms.indices)
            allAtoms = when (allAtoms) {
                null -> atoms.positions!![a] + allCells * atoms.getCell().array
                else -> allAtoms + (atoms.positions!![a] + allCells * atoms.getCell().array)
            }

        // Next, construct the planes forming the unit cell.
        // Atoms from neighboring cells are only included in the list of all relevant atoms, if they have a distance
        // to any of these planes smaller than the cutoff radius.
        // Else wise, they would not be included in the eventual calculation anyhow.
        val planes = arrayOf(
            arrayOf(intArrayOf(0, 1, 0), intArrayOf(0, 0, 1), intArrayOf(0, 0, 0)),
            arrayOf(intArrayOf(gridDimensions[0], 1, 0), intArrayOf(gridDimensions[0], 0, 1), intArrayOf(0, gridDimensions[0], 0), gridDimensions),
            arrayOf(intArrayOf(1, 0, 0), intArrayOf(0, 0, 1), intArrayOf(0, 0, 0)),
            arrayOf(intArrayOf(1, gridDimensions[1], 0), intArrayOf(0, gridDimensions[1], 1), gridDimensions),
            arrayOf(intArrayOf(1, 0, 0), intArrayOf(0, 1, 0), intArrayOf(0, 0, 0)),
            arrayOf(intArrayOf(1, 0, gridDimensions[2]), intArrayOf(0, 1, gridDimensions[2]), gridDimensions))

        val allDistances = ArrayList<ArrayList<Float>>()
        for (plane in planes) {
            val curPlane = Plane.fromPoints(gridToCoord(plane[0]), gridToCoord(plane[1]), gridToCoord(plane[2]))
            val distances = ArrayList<Float>()

            // TODO: This may be optimized, and formulated in an array operation.
            for (a in allAtoms!!.indices)
                distances += curPlane.distancePoint(Point(allAtoms[a]))
            allDistances += distances
        }

        val allDistancesMin = allDistances.minBy { it.sum() }
        val res = allDistancesMin.argWhere { it < parameters.atomicDensityCutoff!! }
//        all_atoms = np.squeeze(
//            all_atoms[
//                np.argwhere(
//                    allDistancesMin < self.parameters.atomic_density_cutoff
//                ),
//            :,
//        ]
//        )
//        return np.concatenate((all_atoms, self.atoms.positions))

        val lines = this::class.java.getResource("allAtoms")!!.readText().lines()
        return Array(lines.size) { lines[it].split(',').map(String::toFloat).toFloatArray() }
    }

    protected fun gridToCoord(gridPoint: IntArray): FloatArray {
        // Convert grid indices to real space grid point .
        val i = gridPoint[0].toFloat()
        val j = gridPoint[1].toFloat()
        val k = gridPoint[2].toFloat()
        // Orthorhombic cells and triclinic ones have to be treated differently, see domain.cpp

        return when {
            atoms.cell!!.orthorhombic -> TODO() // np.diag (self.voxel) * [i, j, k]
            else -> floatArrayOf(i / gridDimensions[0] * atoms.cell!![0, 0] +
                                 j / gridDimensions[1] * atoms.cell!![1, 0] +
                                 k / gridDimensions[2] * atoms.cell!![2, 0],
                                 j / gridDimensions[1] * atoms.cell!![1, 1] +
                                 k / gridDimensions[2] * atoms.cell!![1, 2],
                                 k / gridDimensions[2] * atoms.cell!![2, 2])
        }
    }

    internal abstract fun calculate(outDir: File, kwargs: Map<String, Any>): Number

    companion object {
        infix fun from(params: Parameters): Descriptor = when (params.descriptors.type) {
            Parameters.Descriptors.Type.Bispectrum -> Bispectrum(params)
            Parameters.Descriptors.Type.AtomicDensity -> AtomicDensity(params)
            Parameters.Descriptors.Type.MinterpyDescriptors -> TODO("")
        }
    }
}