package ase

import ase.io.readEspressoOut
import java.io.File

class UnknownFileTypeError(name: String): Exception(name)

val ioFormats: Map<String, IOFormat> by lazy {
    buildMap {
        operator fun String.invoke(desc: String, code: String, module: String? = null,
                                   ext: Any? = null, glob: Any? = null, magic: Any? = null,
                                   encoding: String? = null,
                                   magicRegex: String? = null, external: Boolean = false) {
            var module = module ?: replace('-', '_')

            if (!external)
                module = "ase.io.$module"

            fun Any?.normalizePatterns(): List<String> = when {
                this == null -> emptyList()
                this is String -> listOf(this)
                else -> this as List<String>
            }

            val fmt = IOFormat(this, desc, code, module, encoding)
            fmt.extensions += ext.normalizePatterns()
            fmt.globs += glob.normalizePatterns()
            fmt.magic += magic.normalizePatterns()

            if (magicRegex != null)
                fmt.magicRegex = magicRegex

            for (ext in fmt.extensions) {
                if (ext in extension2format)
                    error("extension \"$ext\" already registered")
                extension2format[ext] = fmt
            }

            put(this, fmt)
        }
        "abinit-gsr"("ABINIT GSR file", "1S", "abinit", glob = "*o_GSR.nc")
        "abinit-in"("ABINIT input file", "1F", "abinit", magic = "*znucl *")
        "abinit-out"("ABINIT output file", "1F", "abinit", magic = "*.Version * of ABINIT")
        "aims"("FHI-aims geometry file", "1S", ext = "in")
        "aims-output"("FHI-aims output", "+S", "aims", magic = "*Invoking FHI-aims ...")
        "bundletrajectory"("ASE bundle trajectory", "+S")
        "castep-castep"("CASTEP output file", "+F", "castep", "castep")
        "castep-cell"("CASTEP geom file", "1F", "castep", "cell")
        "castep-geom"("CASTEP trajectory file", "+F", "castep", "geom")
        "castep-md"("CASTEP molecular dynamics file", "+F", "castep", "md")
        "castep-phonon"("CASTEP phonon file", "1F", "castep", "phonon")
        "cfg"("AtomEye configuration", "1F")
        "cif"("CIF-file", "+B", ext = "cif")
        "cmdft"("CMDFT-file", "1F", glob = "*I_info")
        "cjson"("Chemical json file", "1F", ext = "cjson")
        "cp2k-dcd"("CP2K DCD file", "+B", "cp2k", "dcd")
        "cp2k-restart"("CP2K restart file", "1F", "cp2k", "restart")
        "crystal"("Crystal fort.34 format", " '1F'", ext = listOf("f34", "34"), glob = listOf("f34", "34"))
        "cube"("CUBE file", "1F", ext = "cube")
        "dacapo-text"("Dacapo text output", "1F", "dacapo", magic = "*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n")
        "db"("ASE SQLite database file", "+S")
        "dftb"("DftbPlus input file", "1S", magic = "Geometry")
        "dlp4"("DL_POLY_4 CONFIG file", "1F", "dlp4", "config", "*CONFIG*")
        "dlp-history"("DL_POLY HISTORY file", "+F", "dlp4", glob = "HISTORY")
        "dmol-arc"("DMol3 arc file", "+S", "dmol", "arc")
        "dmol-car"("DMol3 structure file", "1S", "dmol", "car")
        "dmol-incoor"("DMol3 structure file", "1S", "dmol")
        "elk"("ELK atoms definition from GEOMETRY.OUT", "1F", glob = "GEOMETRY.OUT")
        "elk-in"("ELK input file", "1F", "elk")
        "eon"("EON CON file", "+F", ext = "con")
        "eps"("Encapsulated Postscript", "1S")
        "espresso-in"("Quantum espresso in file", "1F", "espresso", "pwi", magic = listOf("*\n&system", "*\n&SYSTEM"))
        "espresso-out"("Quantum espresso out file", "+F", "espresso", listOf("pwo", "out"), magic="*Program PWSCF")
        //F('exciting', 'exciting input', '1F', module='exciting', glob='input.xml')
        //F('exciting', 'exciting output', '1F', module='exciting', glob='INFO.out')
        //F('extxyz', 'Extended XYZ file', '+F', ext='xyz')
        //F('findsym', 'FINDSYM-format', '+F')
        //F('gamess-us-out', 'GAMESS-US output file', '1F',
        //module='gamess_us', magic=b'*GAMESS')
        //F('gamess-us-in', 'GAMESS-US input file', '1F',
        //module='gamess_us')
        //F('gamess-us-punch', 'GAMESS-US punchcard file', '1F',
        //module='gamess_us', magic=b' $DATA', ext='dat')
        //F('gaussian-in', 'Gaussian com (input) file', '1F',
        //module='gaussian', ext=['com', 'gjf'])
        //F('gaussian-out', 'Gaussian output file', '+F',
        //module='gaussian', ext='log', magic=b'*Entering Gaussian System')
        //F('acemolecule-out', 'ACE output file', '1S',
        //module='acemolecule')
        //F('acemolecule-input', 'ACE input file', '1S',
        //module='acemolecule')
        //F('gen', 'DFTBPlus GEN format', '1F')
        //F('gif', 'Graphics interchange format', '+S',
        //module='animation')
        //F('gpaw-out', 'GPAW text output', '+F',
        //magic=b'*  ___ ___ ___ _ _ _')
        //F('gpumd', 'GPUMD input file', '1F', glob='xyz.in')
        //F('gpw', 'GPAW restart-file', '1S',
        //magic=[b'- of UlmGPAW', b'AFFormatGPAW'])
        //F('gromacs', 'Gromacs coordinates', '1F',
        //ext='gro')
        //F('gromos', 'Gromos96 geometry file', '1F', ext='g96')
        //F('html', 'X3DOM HTML', '1F', module='x3d')
        //F('json', 'ASE JSON database file', '+F', ext='json', module='db')
        //F('jsv', 'JSV file format', '1F')
        //F('lammps-dump-text', 'LAMMPS text dump file', '+F',
        //module='lammpsrun', magic_regex=b'.*?^ITEM: TIMESTEP$')
        //F('lammps-dump-binary', 'LAMMPS binary dump file', '+B',
        //module='lammpsrun')
        //F('lammps-data', 'LAMMPS data file', '1F', module='lammpsdata',
        //encoding='ascii')
        //F('magres', 'MAGRES ab initio NMR data file', '1F')
        //F('mol', 'MDL Molfile', '1F')
        //F('mp4', 'MP4 animation', '+S',
        //module='animation')
        //F('mustem', 'muSTEM xtl file', '1F',
        //ext='xtl')
        //F('mysql', 'ASE MySQL database file', '+S',
        //module='db')
        //F('netcdftrajectory', 'AMBER NetCDF trajectory file', '+S',
        //magic=b'CDF')
        //F('nomad-json', 'JSON from Nomad archive', '+F',
        //ext='nomad-json')
        //F('nwchem-in', 'NWChem input file', '1F',
        //module='nwchem', ext='nwi')
        //F('nwchem-out', 'NWChem output file', '+F',
        //module='nwchem', ext='nwo',
        //magic=b'*Northwest Computational Chemistry Package')
        //F('octopus-in', 'Octopus input file', '1F',
        //module='octopus', glob='inp')
        //F('onetep-out', 'ONETEP output file', '+F',
        //module='onetep',
        //magic=b'*Linear-Scaling Ab Initio Total Energy Program*')
        //F('onetep-in', 'ONETEP input file', '1F',
        //module='onetep',
        //magic=[b'*lock species ',
        //b'*LOCK SPECIES ',
        //b'*--- INPUT FILE ---*'])
        //F('proteindatabank', 'Protein Data Bank', '+F',
        //ext='pdb')
        //F('png', 'Portable Network Graphics', '1B')
        //F('postgresql', 'ASE PostgreSQL database file', '+S', module='db')
        //F('pov', 'Persistance of Vision', '1S')
        //# prismatic: Should have ext='xyz' if/when multiple formats can have the same
        //# extension
        //F('prismatic', 'prismatic and computem XYZ-file', '1F')
        //F('py', 'Python file', '+F')
        //F('sys', 'qball sys file', '1F')
        //F('qbox', 'QBOX output file', '+F',
        //magic=b'*:simulation xmlns:')
        //F('res', 'SHELX format', '1S', ext='shelx')
        //F('rmc6f', 'RMCProfile', '1S', ext='rmc6f')
        //F('sdf', 'SDF format', '1F')
        //F('siesta-xv', 'Siesta .XV file', '1F',
        //glob='*.XV', module='siesta')
        //F('struct', 'WIEN2k structure file', '1S', module='wien2k')
        //F('struct_out', 'SIESTA STRUCT file', '1F', module='siesta')
        //F('traj', 'ASE trajectory', '+B', module='trajectory', ext='traj',
        //magic=[b'- of UlmASE-Trajectory', b'AFFormatASE-Trajectory'])
        //F('turbomole', 'TURBOMOLE coord file', '1F', glob='coord',
        //magic=b'$coord')
        //F('turbomole-gradient', 'TURBOMOLE gradient file', '+F',
        //module='turbomole', glob='gradient', magic=b'$grad')
        //F('v-sim', 'V_Sim ascii file', '1F', ext='ascii')
        //F('vasp', 'VASP POSCAR/CONTCAR', '1F',
        //ext='poscar', glob=['*POSCAR*', '*CONTCAR*', '*CENTCAR*'])
        //F('vasp-out', 'VASP OUTCAR file', '+F',
        //module='vasp', glob='*OUTCAR*')
        //F('vasp-xdatcar', 'VASP XDATCAR file', '+F',
        //module='vasp', glob='*XDATCAR*')
        //F('vasp-xml', 'VASP vasprun.xml file', '+F',
        //module='vasp', glob='*vasp*.xml')
        //F('vti', 'VTK XML Image Data', '1F', module='vtkxml')
        //F('vtu', 'VTK XML Unstructured Grid', '1F', module='vtkxml', ext='vtu')
        //F('wout', 'Wannier90 output', '1F', module='wannier90')
        //F('x3d', 'X3D', '1S')
        //F('xsd', 'Materials Studio file', '1F')
        //F('xsf', 'XCrySDen Structure File', '+F',
        //magic=[b'*\nANIMSTEPS', b'*\nCRYSTAL', b'*\nSLAB', b'*\nPOLYMER',
        //b'*\nMOLECULE', b'*\nATOMS'])
        //F('xtd', 'Materials Studio file', '+F')
        //# xyz: No `ext='xyz'` in the definition below.
        //#      The .xyz files are handled by the extxyz module by default.
        //F('xyz', 'XYZ-file', '+F')
    }
}
val extension2format = mutableMapOf<String, IOFormat>()

/** Return ioformat object or raise appropriate error. */
fun getIoformat(name: String): IOFormat {
    if (name !in ioFormats)
        throw UnknownFileTypeError (name)
    val fmt = ioFormats[name]!!
    // Make sure module is importable, since this could also raise an error.
//    fmt.module
    return fmt
}

// We define all the IO formats below.  Each IO format has a code,
// such as '1F', which defines some of the format's properties:
//
// 1=single atoms object
// +=multiple atoms objects
// F=accepts a file-descriptor
// S=needs a file-name str
// B=like F, but opens in binary mode


/**
 *  Read Atoms object(s) from file.
 *  filename: str or file
 *         Name of the file to read from or a file descriptor.
 *     index: int, slice or str
 *         The last configuration will be returned by default.  Examples:
 *
 *             * ``index=0``: first configuration
 *             * ``index=-2``: second to last
 *             * ``index=':'`` or ``index=slice(None)``: all
 *             * ``index='-3:'`` or ``index=slice(-3, None)``: three last
 *             * ``index='::2'`` or ``index=slice(0, None, 2)``: even
 *             * ``index='1::2'`` or ``index=slice(1, None, 2)``: odd
 *     format: str
 *         Used to specify the file-format.  If not given, the
 *         file-format will be guessed by the *filetype* function.
 *     parallel: bool
 *         Default is to read on master and broadcast to slaves.  Use
 *         parallel=False to read on all slaves.
 *     do_not_split_by_at_sign: bool
 *         If False (default) ``filename`` is splitted by at sign ``@``
 *
 *     Many formats allow on open file-like object to be passed instead
 *     of ``filename``. In this case the format cannot be auto-detected,
 *     so the ``format`` argument should be explicitly given.
 */
fun read(file: File,
         format: String? = null,
         parallel: Boolean = true,
         doNotSplitByAtSign: Boolean = false): Atoms {

    //    filename, index = parse_filename(filename, index, do_not_split_by_at_sign)
    //    if index is None:
    //    index = -1
    val format = format ?: TODO() // filetype(filename, read=isinstance(filename, str))

    if (format == "espresso-out")
        return readEspressoOut(file)
    TODO()
}