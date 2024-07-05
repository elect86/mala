package ase.io

import ase.*
import java.io.File
import kotlin.math.pow
import kotlin.math.sqrt

// Quantum ESPRESSO uses CODATA 2006 internally
val units = Units(2006)

// Section identifiers
private const val PW_START = "Program PWSCF"
private const val PW_END = "End of self-consistent calculation"
private const val PW_CELL = "CELL_PARAMETERS"
private const val PW_POS = "ATOMIC_POSITIONS"
private const val PW_MAGMOM = "Magnetic moment per site"
private const val PW_FORCE = "Forces acting on atoms"
private const val PW_TOTEN = "!    total energy"
private const val PW_STRESS = "total   stress"
private const val PW_FERMI = "the Fermi energy is"
private const val PW_HIGHEST_OCCUPIED = "highest occupied level"
private const val PW_HIGHEST_OCCUPIED_LOWEST_FREE = "highest occupied, lowest unoccupied level"
private const val PW_KPTS = "number of k points="
private val PW_BANDS = PW_END
private const val PW_BANDSTRUCTURE = "End of band structure calculation"
private const val PW_DIPOLE = "Debye"
private const val PW_DIPOLE_DIRECTION = "Computed dipole along edir"

// ibrav error message
const val ibravErrorMessage = "ASE does not support ibrav != 0. Note that with ibrav == 0, Quantum ESPRESSO will " +
                              "still detect the symmetries of your system because the CELL_PARAMETERS are defined to " +
                              "a high level of precision."

/**
 * Reads Quantum ESPRESSO output files.
 *
 * The atomistic configurations as well as results (energy, force, stress, magnetic moments) of the calculation are read
 * for all configurations within the output file.
 *
 * Will probably raise errors for broken or incomplete files.
 *
 * @param fileobj A file like object or filename
 * @param index The index of configurations to extract.
 * @param resultsRequired If True, atomistic configurations that do not have any associated results will not be included.
 * This prevents double printed configurations and incomplete calculations from being returned as the final
 * configuration with no results data.
 */
fun readEspressoOut(fileobj: Any, /*index=slice(None),*/ resultsRequired: Boolean = true): Atoms {

    // work with a copy in memory for faster random access
    fileobj as File
    val pwoLines = fileobj.readLines()
    // TODO: index -1 special case?
    // Index all the interesting points
    val indexes: Map<String, ArrayList<Int>> = listOf(PW_START, PW_END, PW_CELL, PW_POS, PW_MAGMOM, PW_FORCE,
                                                      PW_TOTEN, PW_STRESS, PW_FERMI, PW_HIGHEST_OCCUPIED,
                                                      PW_HIGHEST_OCCUPIED_LOWEST_FREE, PW_KPTS, PW_BANDS,
                                                      PW_BANDSTRUCTURE, PW_DIPOLE, PW_DIPOLE_DIRECTION)
            .associateWith { ArrayList() }

    for ((idx, line) in pwoLines.withIndex())
        for (identifier in indexes.keys)
            if (identifier in line)
                indexes[identifier]!! += idx

    // Configurations are either at the start, or defined in ATOMIC_POSITIONS
    // in a subsequent step. Can deal with concatenated output files.
    val allConfigIndexes = (indexes[PW_START]!! + indexes[PW_POS]!!).sorted()

    // Slice only requested indexes
    // setting results_required argument stops configuration-only structures from being returned.
    // This ensures the [-1] structure is one that has results. Two cases:
    // - SCF of last configuration is not converged, job terminated abnormally.
    // - 'relax' and 'vc-relax' re-prints the final configuration but only 'vc-relax' recalculates.
    val imageIndexes = when {
        resultsRequired -> {
            val resultsIndexes = (indexes[PW_TOTEN]!! + indexes[PW_FORCE]!! +
                                  indexes[PW_STRESS]!! + indexes[PW_MAGMOM]!! +
                                  indexes[PW_BANDS]!! +
                                  indexes[PW_BANDSTRUCTURE]!!).sorted()

            // Prune to only configurations with results data before the next configuration
            val resultsConfigIndexes = ArrayList<Int>()
            for ((configIndex, configIndexNext) in allConfigIndexes.zip(allConfigIndexes.drop(1) + pwoLines.size))
                if (resultsIndexes.any { it in configIndex..<configIndexNext })
                    resultsConfigIndexes += configIndex

            // slice from the subset
            resultsConfigIndexes //[index]
        }
        else -> allConfigIndexes //[index]
    }

    // Extract initialisation information each time PWSCF starts to add to subsequent configurations.
    // Use None so slices know when to fill in the blanks.
    val pwscfStartInfo = mutableMapOf<Int, Info>()

    for (imageIndex in imageIndexes) {
        // Find the nearest calculation start to parse info. Needed in, for example, relaxation where cell is only
        // printed at the start.
        val prevStartIndex = if (imageIndex in indexes[PW_START]!!) imageIndex
        else
        // The greatest start index before this structure
            indexes[PW_START]!!.last { it < imageIndex }

        // add structure to reference if not there
        if (pwscfStartInfo[prevStartIndex] == null)
            pwscfStartInfo[prevStartIndex] = parsePwoStart(pwoLines, prevStartIndex)

        // Get the bounds for information for this structure. Any associated values will be between the image_index and
        // the following one, EXCEPT for cell, which will be 4 lines before if it exists.
        val nextIndex = allConfigIndexes.find { it > imageIndex } ?: pwoLines.size // right to the end of the file

        // Get the structure
        // Use this for any missing data
        val prevStructure = pwscfStartInfo[prevStartIndex]!!.atoms
        val structure = when {
            imageIndex in indexes[PW_START]!! -> prevStructure.copy() // parsed from start info
            else -> {
                TODO()
                //                if (PW_CELL in pwoLines[imageIndex - 5])
                //                // CELL_PARAMETERS would be just before positions if present
                //                val (cell, cell_alat) = getCellParameters(pwoLines.drop(imageIndex - 5:image_index])
                //                else:
                //                cell = prevStructure.cell
                //                cell_alat = pwscf_start_info[prev_start_index]['alat']
                //
                //                # give at least enough lines to parse the positions
                //                # should be same format as input card
                //                n_atoms = len(prevStructure)
                //                positions_card = get_atomic_positions(
                //                    pwo_lines[image_index:image_index +n_atoms + 1],
                //                n_atoms = n_atoms, cell = cell, alat = cell_alat)
                //
                //                # convert to Atoms object
                //                symbols = [label_to_symbol(position[0]) for position in
                //                positions_card]
                //                positions = [position[1] for position in positions_card]
                //                structure = Atoms(symbols = symbols, positions = positions, cell = cell, pbc = True)
            }
        }

        // Extract calculation results
        // Energy
        var energy: Float? = null
        for (energyIndex in indexes[PW_TOTEN]!!)
            if (energyIndex in imageIndex..<nextIndex) {
                val split = pwoLines[energyIndex].splitSpaces()
                energy = split[split.size - 2].toFloat() * units.Ry.toFloat()
            }

        // Forces
        var forces: Array<FloatArray>? = null
        for (forceIndex in indexes[PW_FORCE]!!) {
            if (forceIndex in imageIndex..<nextIndex) {
                // Before QE 5.3 'negative rho' added 2 lines before forces
                // Use exact lines to stop before 'non-local' forces in high verbosity
                val offset = if (pwoLines[forceIndex + 2].trim().isEmpty()) 4 else 2
                // assume contiguous
                val rb = (units.Ry / units.Bohr).toFloat()
                forces = pwoLines.drop(forceIndex + offset).take(structure.len).map { forceLine ->
                    val f = forceLine.splitSpaces()
                    floatArrayOf(f[f.size - 3].toFloat() * rb, f[f.size - 2].toFloat() * rb, f[f.size - 1].toFloat() * rb)
                }.toTypedArray()
            }
        }

        // Stress
        var stress: FloatArray? = null
        for (stressIndex in indexes[PW_STRESS]!!) {
            if (stressIndex in imageIndex..<nextIndex) {
                val (sxx, sxy, sxz) = pwoLines[stressIndex + 1].splitSpaces().take(3).map(String::toFloat)
                val (_, syy, syz) = pwoLines[stressIndex + 2].splitSpaces().take(3).map(String::toFloat)
                val (_, _, szz) = pwoLines[stressIndex + 3].splitSpaces().take(3).map(String::toFloat)
                stress = floatArrayOf(sxx, syy, szz, syz, sxz, sxy)
                // sign convention is opposite of ase
                stress *= (-1 * units.Ry / units.Bohr.pow(3)).toFloat()
            }
        }

        // Magmoms
        var magmoms: Any? = null
        for (magmomsIndex in indexes[PW_MAGMOM]!!) {
            TODO()
            //            if image_index < magmomsIndex < next_index:
            //            magmoms = [
            //                float(mag_line.split()[-1]) for mag_line
            //            in pwo_lines[magmoms_index+1:
            //            magmomsIndex + 1 + len(structure)]]
        }

        // Dipole moment
        var dipole: Any? = null
        if (indexes[PW_DIPOLE]!!.isNotEmpty()) {
            TODO()
            //            for dipole_index in indexes[_PW_DIPOLE]:
            //            if image_index < dipole_index < next_index:
            //            _dipole = float(pwo_lines[dipole_index].split()[-2])
            //
            //            for dipole_index in indexes[_PW_DIPOLE_DIRECTION]:
            //            if image_index < dipole_index < next_index:
            //            _direction = pwo_lines[dipole_index].strip()
            //            prefix = 'Computed dipole along edir('
            //            _direction = _direction[len(prefix):]
            //            _direction = int(_direction[0])
            //
            //            dipole = np.eye(3)[_direction - 1] * _dipole * units['Debye']
        }

        // Fermi level / highest occupied level
        var efermi: Float? = null
        for (fermiIndex in indexes[PW_FERMI]!!) {
            if (fermiIndex in imageIndex..<nextIndex)
                efermi = pwoLines[fermiIndex].splitSpaces().run { get(size - 2).toFloat() }
        }

        if (efermi == null) {
            TODO()
            //            for ho_index in indexes[_PW_HIGHEST_OCCUPIED]:
            //            if image_index < ho_index < next_index:
            //            efermi = float(pwo_lines[ho_index].split()[-1])
        }

        if (efermi == null) {
            TODO()
            //            for holf_index in indexes[_PW_HIGHEST_OCCUPIED_LOWEST_FREE]:
            //            if image_index < holf_index < next_index:
            //            efermi = float(pwo_lines[holf_index].split()[-2])
        }

        // K-points
        lateinit var ibzkpts: Array<FloatArray>
        lateinit var weights: FloatArray
        var kpointsWarning = "Number of k-points >= 100: set verbosity='high' to print them."

        for (i in indexes[PW_KPTS]!!) {
            val nkpts = Regex("\\b\\d+\\b").findAll(pwoLines[i]).first().groupValues[0].toInt()
            val kptsIndex = i + 2

            if (pwoLines[kptsIndex].trim() == kpointsWarning)
                continue

            // QE prints the k -points in units of 2 * pi / alat
            // with alat defined as the length of the first cell vector
            val cell = structure.getCell()
            val alat = cell[0].norm()
            ibzkpts = Array(nkpts) { FloatArray(3) }
            weights = FloatArray(nkpts)
            for (pt in 0..<nkpts) {
                val L = pwoLines[kptsIndex + pt].splitSpaces()
                weights[pt] = L.last().toFloat()
                val coord = floatArrayOf(L[L.size - 6].toFloat(), L[L.size - 5].toFloat(), L[L.size - 4].dropLast(2).toFloat())
                coord *= 2 * Math.PI.toFloat() / alat
                ibzkpts[pt] = kpointConvert(cell, ckptsKv = coord)
            }
        }

        // Bands
        lateinit var kpts: MutableList<SinglePointKPoint>
        kpointsWarning = "Number of k-points >= 100: set verbosity='high' to print the bands."

        for (idx in indexes[PW_BANDS]!! + indexes[PW_BANDSTRUCTURE]!!)
            if (idx in imageIndex..<nextIndex) {
                var bandsIndex = idx + 1
                // skip over the lines with DFT +U occupation matrices
                if ("enter write_ns" in pwoLines[bandsIndex])
                    while ("exit write_ns" !in pwoLines[bandsIndex])
                        bandsIndex++
                bandsIndex++

                if (pwoLines[bandsIndex].trim() == kpointsWarning)
                    continue

                //                assert ibzkpts is not None
                var spin = 0
                val bands = mutableListOf<Float>()
                val eigenvalues = Array(2) { mutableListOf<MutableList<Float>>() }

                while (true) {
                    val L = pwoLines[bandsIndex].replace("-", " -").splitSpaces()
                    when {
                        L.isEmpty() ->
                            if (bands.isNotEmpty()) {
                                eigenvalues[spin].add(bands.toMutableList())
                                bands.clear()
                            }
                        L == listOf("occupation", "numbers") ->
                            // Skip the lines with the occupation numbers
                            bandsIndex += eigenvalues[spin][0].size // 8 + 1
                        L[0] == "k" && L[1].startsWith('=') -> Unit
                        "SPIN" in L -> if ("DOWN" in L) spin += 1
                        else -> try {
                            bands += L.map(String::toFloat)
                        } catch (ex: Exception) {
                            break
                        }
                    }
                    bandsIndex++
                }
                if (spin == 1)
                    check(eigenvalues[0].size == eigenvalues[1].size)
                check(eigenvalues[0].size == ibzkpts.size) { "np.shape(eigenvalues), len(ibzkpts)" }

                kpts = mutableListOf()
                for (s in 0..<spin + 1)
                    for (i in weights.indices)
                        kpts += SinglePointKPoint(weights[i], s, ibzkpts[i], epsN = eigenvalues[s])
            }

        // Put everything together
        //
        // In PW the forces are consistent with the "total energy"; that's why its value must be assigned to free_energy.
        // PW doesn't compute the extrapolation of the energy to 0K smearing the closer thing to this is again the total
        // energy that contains the correct (i.e. variational) form of the band energy is
        //   Eband = \int e N(e) de   for e<Ef , where N(e) is the DOS
        // This differs by the term (-TS)  from the sum of KS eigenvalues:
        //    Eks = \sum wg(n,k) et(n,k)
        // which is non variational. When a Fermi-Dirac function is used for a given T, the variational energy is REALLY
        // the free energy F, and F = E - TS , with E = non variational energy.

        val calc = SinglePointDFTCalculator(structure, energy = energy,
                                            freeEnergy = energy, kpts = kpts,
                                            forces = forces, stress = stress,
                                            magmoms = magmoms, efermi = efermi,
                                            ibzkpts = ibzkpts, dipole = dipole)

        return structure
    }
    error("")
}

/**
 * Parse Quantum ESPRESSO calculation info from lines, starting from index. Return a dictionary containing extracted
 * information.
 *
 * - `celldm(1)`: lattice parameters (alat)
 * - `cell`: unit cell in Angstrom
 * - `symbols`: element symbols for the structure
 * - `positions`: cartesian coordinates of atoms in Angstrom
 * - `atoms`: an `ase.Atoms` object constructed from the extracted data
 *
 * Parameters
 * ----------
 * @param lines Contents of PWSCF output file.
 * @param index Line number to begin parsing. Only first calculation will be read.
 *
 * @return
 * -------
 * info : dict
 *         Dictionary of calculation parameters, including `celldm(1)`, `cell`,
 *         `symbols`, `positions`, `atoms`.
 *
 *     Raises
 *     ------
 *     KeyError
 *         If interdependent values cannot be found (especially celldm(1))
 *         an error will be raised as other quantities cannot then be
 *         calculated (e.g. cell and positions).
 *
 */
fun parsePwoStart(lines: List<String>, index: Int = 0): Info {
    // TODO: extend with extra DFT info?

    val info = Info()

    for ((idx, line) in lines.withIndex()) {
        if (idx < index) continue
        fun last() = line.substringAfterLast(' ').trimEnd('\n')
        fun String.splitFloats(start: Int, end: Int) = splitSpaces().subList(start, end).map(String::toFloat).toFloatArray()
        if ("celldm(1)" in line) {
            // celldm(1) has more digits than alat !!
            info.`celldm(1)` = line.splitSpaces()[1].toFloat() * units.Bohr.toFloat()
            info.alat = info.`celldm(1)`
        } else if ("number of atoms/cell" in line)
            info.nat = last().toInt()
        else if ("number of atomic types" in line)
            info.ntyp = last().toInt()
        else if ("crystal axes:" in line)
            info.cell = info.`celldm(1)`!! * arrayOf(lines[idx + 1].splitFloats(3, 6),
                                                     lines[idx + 2].splitFloats(3, 6),
                                                     lines[idx + 3].splitFloats(3, 6))
        else if ("positions (alat units)" in line) {
            for (atLine in lines.drop(idx + 1).take(info.nat!!)) {
                val (sym, x, y, z) = atLine.parsePositionLine()
                info.symbols += sym.labelToSymbol()
                info.positions += floatArrayOf(x.toFloat() * info.`celldm(1)`!!,
                                               y.toFloat() * info.`celldm(1)`!!,
                                               z.toFloat() * info.`celldm(1)`!!)
            }
            // This should be the end of interesting info.
            // Break here to avoid dealing with large lists of kpoints.
            // Will need to be extended for DFTCalculator info.
            break
        }
    }

    // Make atoms for convenience
    info.atoms = Atoms(symbols = info.symbols,
                       positions = info.positions.toTypedArray(),
                       cell = info.cell, pbcArg = true)

    return info
}

/**
 * Parse a single line from a pw.x output file.
 *
 * The line must contain information about the atomic symbol and the position, e.g.
 *
 * 995           Sb  tau( 995) = (   1.4212023   0.7037863   0.1242640  )
 *
 * @receiver Line to be parsed.
 *
 * @return
 * sym, // Atomic symbol.
 * x,   // float, x-position.
 * y,   // float, y-position.
 * z    // float, z-position.
 */
fun String.parsePositionLine(): List<String> {
    val regex = Regex("""\s*\d+\s*(\S+)\s*tau\(\s*\d+\)\s*=\s*\(\s*(\S+)\s+(\S+)\s+(\S+)\s*\)""")
    return regex.find(this)!!.groupValues.drop(1)
}

/**
 * Convert a valid espresso ATOMIC_SPECIES label to a
 *     chemical symbol.
 *
 *     Parameters
 *     ----------
 *     label : str
 *         chemical symbol X (1 or 2 characters, case-insensitive)
 *         or chemical symbol plus a number or a letter, as in
 *         "Xn" (e.g. Fe1) or "X_*" or "X-*" (e.g. C1, C_h;
 *         max total length cannot exceed 3 characters).
 *
 *     Returns
 *     -------
 *     symbol : str
 *         The best matching species from ase.utils.chemical_symbols
 *
 *     Raises
 *     ------
 *     KeyError
 *         Couldn't find an appropriate species.
 *
 *     Notes
 *     -----
 *         It's impossible to tell whether e.g. He is helium
 *         or hydrogen labelled 'e'.
 */
fun String.labelToSymbol(): String {
    // possibly a two character species
    // ase Atoms need proper case of chemical symbols .
    if (length >= 2) {
        val testSymbol = this[0].uppercase() + this[1].lowercase()
        if (testSymbol in chemicalSymbols)
            return testSymbol
    }
    // finally try with one character
    val testSymbol = this[0].uppercase()
    if (testSymbol in chemicalSymbols)
        return testSymbol
    else error("Could not parse species from label $this.")
}

class Info {
    var `celldm(1)`: Float? = null
    var alat: Float? = null
    var nat: Int? = null
    var ntyp: Int? = null
    lateinit var cell: Array<FloatArray>
    val symbols = ArrayList<String>()
    val positions = ArrayList<FloatArray>()
    lateinit var atoms: Atoms
}

/**
 * Parse unit cell from CELL_PARAMETERS card.
 *
 *     Parameters
 *     ----------
 *     lines : list[str]
 *         A list with lines containing the CELL_PARAMETERS card.
 *     alat : float | None
 *         Unit of lattice vectors in Angstrom. Only used if the card is
 *         given in units of alat. alat must be None if CELL_PARAMETERS card
 *         is in Bohr or Angstrom. For output files, alat will be parsed from
 *         the card header and used in preference to this value.
 *
 *     Returns
 *     -------
 *     cell : np.array | None
 *         Cell parameters as a 3x3 array in Angstrom. If no cell is found
 *         None will be returned instead.
 *     cell_alat : float | None
 *         If a value for alat is given in the card header, this is also
 *         returned, otherwise this will be None.
 *
 *     Raises
 *     ------
 *     ValueError
 *         If CELL_PARAMETERS are given in units of bohr or angstrom
 *         and alat is not
 */
fun getCellParameters(lines: List<String>, alat: Any? = null): Pair<Any, Any> {
    TODO()
    //    cell = None
    //    cell_alat = None
    //    # no blanks or comment lines, can take three lines for cell
    //    trimmed_lines = (line for line in lines if line.strip() and line[0] != '#')
    //
    //    for line in trimmed_lines:
    //    if line.strip().startswith('CELL_PARAMETERS'):
    //    if cell is not None :
    //    # multiple definitions
    //            raise ValueError ('CELL_PARAMETERS specified multiple times')
    //    # Priority and behaviour tested with QE 5.3
    //    if 'bohr' in line.lower():
    //    if alat is not None :
    //    raise ValueError ('Lattice parameters given in '
    //    '&SYSTEM celldm/A and CELL_PARAMETERS '
    //    'bohr')
    //    cell_units = units['Bohr']
    //    elif 'angstrom' in line.lower():
    //    if alat is not None :
    //    raise ValueError ('Lattice parameters given in '
    //    '&SYSTEM celldm/A and CELL_PARAMETERS '
    //    'angstrom')
    //    cell_units = 1.0
    //    elif 'alat' in line.lower():
    //    # Output file has(alat = value)(in Bohrs)
    //    if '=' in line:
    //    alat = float(line.strip(') \n').split()[-1]) * units['Bohr']
    //    cell_alat = alat
    //    elif alat is None :
    //    raise ValueError ('Lattice parameters must be set in '
    //    '&SYSTEM for alat units')
    //    cell_units = alat
    //    elif alat is None :
    //    # may be DEPRECATED in future
    //    cell_units = units['Bohr']
    //    else:
    //    # may be DEPRECATED in future
    //    cell_units = alat
    //    # Grab the parameters; blank lines have been removed
    //    cell = [[ffloat(x) for x in next(trimmed_lines).split()[:3]],
    //    [ffloat(x) for x in next(trimmed_lines).split()[:3]],
    //    [ffloat(x) for x in next(trimmed_lines).split()[:3]]]
    //    cell = np.array(cell) * cell_units
    //
    //    return cell, cell_alat
}

operator fun Float.times(floats: Array<FloatArray>): Array<FloatArray> = Array(floats.size) { i -> FloatArray(floats[i].size) { this * floats[i][it] } }
fun String.splitSpaces() = trim().split(Regex("\\s+")).filter { it.isNotBlank() }
operator fun FloatArray.timesAssign(float: Float) {
    for (i in indices)
        this[i] *= float
}

infix fun FloatArray.dot(floats: FloatArray): Float {
    check(size == floats.size)
    var dot = 0f
    for (i in indices)
        dot += this[i] * floats[i]
    return dot
}

infix fun FloatArray.dot(arrays: Array<FloatArray>): FloatArray {
    check(size == arrays.size)
    return FloatArray(size) { i ->
        var dot = 0f
        for (j in indices)
            dot += this[j] * arrays[j][i]
        dot
    }
}

class SinglePointKPoint(val weight: Float,
                        /** spin index */
                        val s: Int,
                        /** k-point index */
                        val k: FloatArray,
                        val epsN: MutableList<MutableList<Float>> = mutableListOf(),
                        val fN: Any = Unit)

fun FloatArray.norm(): Float = sqrt(this dot this)
operator fun FloatArray.div(float: Float): FloatArray = FloatArray(size) { this[it] / float }
fun IntArray.prod(): Int {
    var res = 1
    for (i in this)
        res *= i
    return res
}
