@file:Suppress("UNCHECKED_CAST")

package ase.io

import ase.Atoms
import ase.calculators.lammps.Prism
import ase.calculators.lammps.UnitSet
import ase.calculators.lammps.convert
import ase.component6
import java.io.File

/**
 * Write atomic structure data to a LAMMPS data file.
 *
 *     Parameters
 *     ----------
 *     fd : file|str
 *         File to which the output will be written.
 *     atoms : Atoms
 *         Atoms to be written.
 *     specorder : list[str], optional
 *         Chemical symbols in the order of LAMMPS atom types, by default None
 *     force_skew : bool, optional
 *         Force to write the cell as a
 *         `triclinic <https://docs.lammps.org/Howto_triclinic.html>`__ box,
 *         by default False
 *     reduce_cell : bool, optional
 *         Whether the cell shape is reduced or not, by default False
 *     prismobj : Prism|None, optional
 *         Prism, by default None
 *     write_image_flags : bool, default False
 *         If True, the image flags, i.e., in which images of the periodic
 *         simulation box the atoms are, are written.
 *     masses : bool, optional
 *         Whether the atomic masses are written or not, by default False
 *     velocities : bool, optional
 *         Whether the atomic velocities are written or not, by default False
 *     units : str, optional
 *         `LAMMPS units <https://docs.lammps.org/units.html>`__,
 *         by default 'metal'
 *     bonds : bool, optional
 *         Whether the bonds are written or not. Bonds can only be written
 *         for atom_style='full', by default True
 *     atom_style : {'atomic', 'charge', 'full'}, optional
 *         `LAMMPS atom style <https://docs.lammps.org/atom_style.html>`__,
 *         by default 'atomic'
 */
fun writeLammpsData(file: File, atoms: Any, args: Map<String, Any>) {

    // FIXME: We should add a check here that the encoding of the file object
    //        is actually ascii once the 'encoding' attribute of IOFormat objects
    //        starts functioning in implementation (currently it doesn't do
    //         anything).

    val specOrder: List<String>? = args["specirder"] as List<String>?
    val forceSkew = args["force_skew"] as Boolean? ?: false
    //    prismobj: Prism? = null,
    val writeImageFlags = args["write_image_flags"] as Boolean? ?: false
    val masses = args["masses"] as Boolean? ?: false
    val velocities = args["velocities"] as Boolean? ?: false
    val bonds = args["bonds"] as Boolean? ?: true
    val atomStyle = args["atomStyle"] as AtomStyle? ?: AtomStyle.atomic

    check(atoms is Atoms) { "Can only write one configuration to a lammps data file!" }

    val fd = StringBuilder("(written by ASE)\n\n")
    operator fun String.unaryPlus() = fd.append(this)

    val symbols = atoms.chemicalSymbols
    val nAtoms = symbols.size
    +"$nAtoms atoms\n"

    val species =
        // To index elements in the LAMMPS data file (indices must correspond to order in the potential file)
        specOrder ?:
        // This way it is assured that LAMMPS atom types are always assigned predictably according to the alphabetic order
        symbols.toSet().sorted()

    val nAtomTypes = species.size
    +"$nAtomTypes atom types\n\n"

    val bondsIn: List<Int> = emptyList()
    if (bonds && atomStyle == AtomStyle.full && atoms.arrays["bonds"] != null) {
        TODO()
        //        n_bonds = 0
        //        n_bond_types = 1
        //        for i, bondsi in enumerate(atoms.arrays['bonds']):
        //        if bondsi != '_':
        //        for bond in bondsi.split(','):
        //        dummy1, dummy2 = bond.split('(')
        //        bond_type = int(dummy2.split(')')[0])
        //        at1 = int(i) + 1
        //        at2 = int(dummy1) + 1
        //        bonds_in.append((bond_type, at1, at2))
        //        n_bonds = n_bonds + 1
        //        if bond_type > n_bond_types:
        //        n_bond_types = bond_type
        //        fd.write(f'{n_bonds} bonds\n')
        //        fd.write(f'{n_bond_types} bond types\n\n')
    }

    val prismObj = args["prismobj"] as Prism? ?: Prism(atoms.getCell(), isReduced = args["reduce_cell"] as Boolean? ?: false)

    // Get cell parameters and convert from ASE units to LAMMPS units
    val units = args["units"] as UnitSet? ?: UnitSet.metal
    val (xhi, yhi, zhi, xy, xz, yz) = prismObj.lammpsPrism convert UnitSet.ase.distance..units.distance

    +"0.0 $xhi  xlo xhi\n"
    +"0.0 $yhi  ylo yhi\n"
    +"0.0 $zhi  zlo zhi\n"

    if (forceSkew || prismObj.isSkewed)
        +"$xy $xz $yz  xy xz yz\n"
    +"\n"

    if (masses) TODO()
    //    _write_masses(fd, atoms, species, units)

    // Write (unwrapped) atomic positions.  If wrapping of atoms back into the
    // cell along periodic directions is desired, this should be done manually
    // on the Atoms object itself beforehand.
    +"Atoms # $atomStyle\n\n"

    if (writeImageFlags) TODO()
    //    scaled_positions = atoms.get_scaled_positions(wrap=False)
    //    image_flags = np.floor(scaled_positions).astype(int)

    // when `write_image_flags` is True, the positions are wrapped while the
    // unwrapped positions can be recovered from the image flags
    var pos = prismObj.vectorToLammps(atoms.getPositions(), writeImageFlags)

    when(atomStyle) {
        AtomStyle.atomic -> {
            // Convert position from ASE units to LAMMPS units
            pos = pos convert UnitSet.ase.distance..units.distance
            for ((i, r) in pos.withIndex()) {
                val s = species.indexOf(symbols[i]) + 1
                +"${i+1} $s ${r[0]} ${r[1]} ${r[2]}"
                if (writeImageFlags) {
                    TODO()
//                    img = image_flags[i]
//                    line += f' {img[0]:6d} {img[1]:6d} {img[2]:6d}'
                }
                +"\n"
            }
        }
        AtomStyle.charge -> TODO()
//        charges = atoms.get_initial_charges()
//            # Convert position and charge from ASE units to LAMMPS units
//                pos = convert(pos, 'distance', 'ASE', units)
//        charges = convert(charges, 'charge', 'ASE', units)
//            for i, (q, r) in enumerate(zip(charges, pos)):
//        s = species.index(symbols[i]) + 1
//        line = (
//                f'{i+1:>6} {s:>3} {q:>5}'
//            f' {r[0]:23.17g} {r[1]:23.17g} {r[2]:23.17g}'
//        )
//        if write_image_flags:
//        img = image_flags[i]
//        line += f' {img[0]:6d} {img[1]:6d} {img[2]:6d}'
//            line += '\n'
//        fd.write(line)
        AtomStyle.full -> TODO()
//        charges = atoms.get_initial_charges()
//            # The label 'mol-id' has apparenlty been introduced in read earlier,
//        # but so far not implemented here. Wouldn't a 'underscored' label
//            # be better, i.e. 'mol_id' or 'molecule_id'?
//        if atoms.has('mol-id'):
//        molecules = atoms.get_array('mol-id')
//        if not np.issubdtype(molecules.dtype, np.integer):
//        raise TypeError(
//                f'If "atoms" object has "mol-id" array, then '
//            f'mol-id dtype must be subtype of np.integer, and '
//            f'not {molecules.dtype!s:s}.')
//        if (len(molecules) != len(atoms)) or (molecules.ndim != 1):
//        raise TypeError(
//                'If "atoms" object has "mol-id" array, then '
//        'each atom must have exactly one mol-id.')
        else -> TODO()
            // Assigning each atom to a distinct molecule id would seem
            // preferableabove assigning all atoms to a single molecule
            // id per default, as done within ase <= v 3.19.1. I.e.,
            // molecules = np.arange(start=1, stop=len(atoms)+1,
            // step=1, dtype=int) However, according to LAMMPS default
            // behavior,
//            molecules = np.zeros(len(atoms), dtype=int)
            // which is what happens if one creates new atoms within LAMMPS
            // without explicitly taking care of the molecule id.
            // Quote from docs at https://lammps.sandia.gov/doc/read_data.html:
            //    The molecule ID is a 2nd identifier attached to an atom.
            //    Normally, it is a number from 1 to N, identifying which
            //    molecule the atom belongs to. It can be 0 if it is a
            //    non-bonded atom or if you don't care to keep track of molecule
            //    assignments.
    }

    if (velocities /*&& atoms.get_velocities() is not None*/) {
        TODO()
//        fd.write('\n\nVelocities\n\n')
//        vel = prismobj.vector_to_lammps(atoms.get_velocities())
//        # Convert velocity from ASE units to LAMMPS units
//                vel = convert(vel, 'velocity', 'ASE', units)
//        for i, v in enumerate(vel):
//        fd.write(f'{i+1:>6} {v[0]:23.17g} {v[1]:23.17g} {v[2]:23.17g}\n')
    }

    file.writeText(fd.toString())
}

enum class AtomStyle { atomic, charge, full }