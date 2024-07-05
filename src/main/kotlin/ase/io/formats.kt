package ase.io

import ase.Atoms
import ase.calculators.lammps.Prism
import java.io.File
import kotlin.reflect.KFunction3

class IOFormat(val name: String,
               val description: String,
               val code: String,
               val moduleName: String,
               val encoding: String? = null) {
    init {
        check(code.length == 2 && code[0] in "+1" && code[1] in "BFS")
    }

    // (To be set by define_io_format())
    var extensions = emptyList<String>()
    var globs = emptyList<String>()
    var magic = emptyList<String>()
    var magicRegex: String? = null

    val canWrite: Boolean
        get() = writeFunc != null

    /** Whether this format is for a single Atoms object. */
    val single: Boolean
        get() = code[0] == '1'

    val writeFunc: KFunction3<File, Any, Map<String, Any>, Unit>? = when (name) {
        "lammps-data" -> ::writeLammpsData
        else -> null
    }

    fun write(file: File, images: Any, kwargs: Map<String, Any> = emptyMap()) = writeFunc!!.invoke(file, images, kwargs)

    val acceptsFd: Boolean
        get() = code[1] != 'S'
}

/** These will be filled at run-time. */
val ioFormats: Map<String, IOFormat> = buildMap {

    operator fun String.invoke(desc: String, code: String,
        //*,
                               module: String = replace('-', '_'), ext: Any? = null,
                               glob: Any? = null, magic: Any? = null, encoding: String? = null,
                               magicRegex: String? = null, external: Boolean = false): IOFormat {

        //        if not external :
        //        module = 'ase.io.' + module

        @Suppress("UNCHECKED_CAST")
        fun Any?.normalizePatterns() = when (this) {
            null -> emptyList()
            is String -> listOf(this)
            else -> this as List<String>
        }

        val fmt = IOFormat(this, desc, code, module, encoding).apply {
            extensions = ext.normalizePatterns()
            globs = glob.normalizePatterns()
            this.magic = magic.normalizePatterns()
        }
        if (magicRegex != null)
            fmt.magicRegex = magicRegex

        for (ext in fmt.extensions) {
            check(ext !in extension2format) { "extension \"$ext\" already registered" }
            extension2format[ext] = fmt
        }

        put(this, fmt)
        return fmt
    }

    /*F('abinit-gsr', 'ABINIT GSR file', '1S',
      module='abinit', glob='*o_GSR.nc')
    F('abinit-in', 'ABINIT input file', '1F',
      module='abinit', magic=b'*znucl *')
    F('abinit-out', 'ABINIT output file', '1F',
      module='abinit', magic=b'*.Version * of ABINIT')
    F('aims', 'FHI-aims geometry file', '1S', ext='in')
    F('aims-output', 'FHI-aims output', '+S',
      module='aims', magic=b'*Invoking FHI-aims ...')
    F('bundletrajectory', 'ASE bundle trajectory', '+S')
    F('castep-castep', 'CASTEP output file', '+F',
      module='castep', ext='castep')
    F('castep-cell', 'CASTEP geom file', '1F',
      module='castep', ext='cell')
    F('castep-geom', 'CASTEP trajectory file', '+F',
      module='castep', ext='geom')
    F('castep-md', 'CASTEP molecular dynamics file', '+F',
      module='castep', ext='md')
    F('castep-phonon', 'CASTEP phonon file', '1F',
      module='castep', ext='phonon')
    F('cfg', 'AtomEye configuration', '1F')
    F('cif', 'CIF-file', '+B', ext='cif')
    F('cmdft', 'CMDFT-file', '1F', glob='*I_info')
    F('cjson', 'Chemical json file', '1F', ext='cjson')
    F('cp2k-dcd', 'CP2K DCD file', '+B',
      module='cp2k', ext='dcd')
    F('cp2k-restart', 'CP2K restart file', '1F',
      module='cp2k', ext='restart')
    F('crystal', 'Crystal fort.34 format', '1F',
      ext=['f34', '34'], glob=['f34', '34'])
    F('cube', 'CUBE file', '1F', ext='cube')
    F('dacapo-text', 'Dacapo text output', '1F',
      module='dacapo', magic=b'*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n')
    F('db', 'ASE SQLite database file', '+S')
    F('dftb', 'DftbPlus input file', '1S', magic=b'Geometry')
    F('dlp4', 'DL_POLY_4 CONFIG file', '1F',
      module='dlp4', ext='config', glob=['*CONFIG*'])
    F('dlp-history', 'DL_POLY HISTORY file', '+F',
      module='dlp4', glob='HISTORY')
    F('dmol-arc', 'DMol3 arc file', '+S',
      module='dmol', ext='arc')
    F('dmol-car', 'DMol3 structure file', '1S',
      module='dmol', ext='car')
    F('dmol-incoor', 'DMol3 structure file', '1S',
      module='dmol')
    F('elk', 'ELK atoms definition from GEOMETRY.OUT', '1F',
      glob=['GEOMETRY.OUT'])
    F('elk-in', 'ELK input file', '1F', module='elk')
    F('eon', 'EON CON file', '+F',
      ext='con')
    F('eps', 'Encapsulated Postscript', '1S')
    F('espresso-in', 'Quantum espresso in file', '1F',
      module='espresso', ext='pwi', magic=[b'*\n&system', b'*\n&SYSTEM'])
    F('espresso-out', 'Quantum espresso out file', '+F',
      module='espresso', ext=['pwo', 'out'], magic=b'*Program PWSCF')
    F('exciting', 'exciting input', '1F', module='exciting', glob='input.xml')
    F('exciting', 'exciting output', '1F', module='exciting', glob='INFO.out')
    F('extxyz', 'Extended XYZ file', '+F', ext='xyz')
    F('findsym', 'FINDSYM-format', '+F')
    F('gamess-us-out', 'GAMESS-US output file', '1F',
      module='gamess_us', magic=b'*GAMESS')
    F('gamess-us-in', 'GAMESS-US input file', '1F',
      module='gamess_us')
    F('gamess-us-punch', 'GAMESS-US punchcard file', '1F',
      module='gamess_us', magic=b' $DATA', ext='dat')
    F('gaussian-in', 'Gaussian com (input) file', '1F',
      module='gaussian', ext=['com', 'gjf'])
    F('gaussian-out', 'Gaussian output file', '+F',
      module='gaussian', ext='log', magic=b'*Entering Gaussian System')
    F('acemolecule-out', 'ACE output file', '1S',
      module='acemolecule')
    F('acemolecule-input', 'ACE input file', '1S',
      module='acemolecule')
    F('gen', 'DFTBPlus GEN format', '1F')
    F('gif', 'Graphics interchange format', '+S',
      module='animation')
    F('gpaw-out', 'GPAW text output', '+F',
      magic=b'*  ___ ___ ___ _ _ _')
    F('gpumd', 'GPUMD input file', '1F', glob='xyz.in')
    F('gpw', 'GPAW restart-file', '1S',
      magic=[b'- of UlmGPAW', b'AFFormatGPAW'])
    F('gromacs', 'Gromacs coordinates', '1F',
      ext='gro')
    F('gromos', 'Gromos96 geometry file', '1F', ext='g96')
    F('html', 'X3DOM HTML', '1F', module='x3d')
    F('json', 'ASE JSON database file', '+F', ext='json', module='db')
    F('jsv', 'JSV file format', '1F')
    F('lammps-dump-text', 'LAMMPS text dump file', '+F',
      module='lammpsrun', magic_regex=b'.*?^ITEM: TIMESTEP$')
    F('lammps-dump-binary', 'LAMMPS binary dump file', '+B',
      module='lammpsrun')*/
    "lammps-data"("LAMMPS data file", "1F", module = "lammpsdata", encoding = "ascii")
    /*F('magres', 'MAGRES ab initio NMR data file', '1F')
    F('mol', 'MDL Molfile', '1F')
    F('mp4', 'MP4 animation', '+S',
      module='animation')
    F('mustem', 'muSTEM xtl file', '1F',
      ext='xtl')
    F('mysql', 'ASE MySQL database file', '+S',
      module='db')
    F('netcdftrajectory', 'AMBER NetCDF trajectory file', '+S',
      magic=b'CDF')
    F('nomad-json', 'JSON from Nomad archive', '+F',
      ext='nomad-json')
    F('nwchem-in', 'NWChem input file', '1F',
      module='nwchem', ext='nwi')
    F('nwchem-out', 'NWChem output file', '+F',
      module='nwchem', ext='nwo',
      magic=b'*Northwest Computational Chemistry Package')
    F('octopus-in', 'Octopus input file', '1F',
      module='octopus', glob='inp')
    F('onetep-out', 'ONETEP output file', '+F',
      module='onetep',
      magic=b'*Linear-Scaling Ab Initio Total Energy Program*')
    F('onetep-in', 'ONETEP input file', '1F',
      module='onetep',
      magic=[b'*lock species ',
      b'*LOCK SPECIES ',
      b'*--- INPUT FILE ---*'])
    F('proteindatabank', 'Protein Data Bank', '+F',
      ext='pdb')
    F('png', 'Portable Network Graphics', '1B')
    F('postgresql', 'ASE PostgreSQL database file', '+S', module='db')
    F('pov', 'Persistance of Vision', '1S')
    # prismatic: Should have ext='xyz' if/when multiple formats can have the same
    # extension
    F('prismatic', 'prismatic and computem XYZ-file', '1F')
    F('py', 'Python file', '+F')
    F('sys', 'qball sys file', '1F')
    F('qbox', 'QBOX output file', '+F',
      magic=b'*:simulation xmlns:')
    F('res', 'SHELX format', '1S', ext='shelx')
    F('rmc6f', 'RMCProfile', '1S', ext='rmc6f')
    F('sdf', 'SDF format', '1F')
    F('siesta-xv', 'Siesta .XV file', '1F',
      glob='*.XV', module='siesta')
    F('struct', 'WIEN2k structure file', '1S', module='wien2k')
    F('struct_out', 'SIESTA STRUCT file', '1F', module='siesta')
    F('traj', 'ASE trajectory', '+B', module='trajectory', ext='traj',
      magic=[b'- of UlmASE-Trajectory', b'AFFormatASE-Trajectory'])
    F('turbomole', 'TURBOMOLE coord file', '1F', glob='coord',
      magic=b'$coord')
    F('turbomole-gradient', 'TURBOMOLE gradient file', '+F',
      module='turbomole', glob='gradient', magic=b'$grad')
    F('v-sim', 'V_Sim ascii file', '1F', ext='ascii')
    F('vasp', 'VASP POSCAR/CONTCAR', '1F',
      ext='poscar', glob=['*POSCAR*', '*CONTCAR*', '*CENTCAR*'])
    F('vasp-out', 'VASP OUTCAR file', '+F',
      module='vasp', glob='*OUTCAR*')
    F('vasp-xdatcar', 'VASP XDATCAR file', '+F',
      module='vasp', glob='*XDATCAR*')
    F('vasp-xml', 'VASP vasprun.xml file', '+F',
      module='vasp', glob='*vasp*.xml')
    F('vti', 'VTK XML Image Data', '1F', module='vtkxml')
    F('vtu', 'VTK XML Unstructured Grid', '1F', module='vtkxml', ext='vtu')
    F('wout', 'Wannier90 output', '1F', module='wannier90')
    F('x3d', 'X3D', '1S')
    F('xsd', 'Materials Studio file', '1F')
    F('xsf', 'XCrySDen Structure File', '+F',
      magic=[b'*\nANIMSTEPS', b'*\nCRYSTAL', b'*\nSLAB', b'*\nPOLYMER',
      b'*\nMOLECULE', b'*\nATOMS'])
    F('xtd', 'Materials Studio file', '+F')
    # xyz: No `ext='xyz'` in the definition below.
    #      The .xyz files are handled by the extxyz module by default.
    F('xyz', 'XYZ-file', '+F')*/
}
val extension2format = mutableMapOf<String, IOFormat>()

/** Return ioformat object or raise appropriate error. */
fun getIoFormat(name: String): IOFormat {
    val fmt = ioFormats[name] ?: error("invalid name $name")
    // Make sure module is importable, since this could also raise an error.
    //    fmt.module
    return fmt
}

/**
 *     Wrapper around builtin `open` that will guess compression of a file
 *     from the filename and open it for reading or writing as if it were
 *     a standard file.
 *
 *     Implemented for ``gz``(gzip), ``bz2``(bzip2) and ``xz``(lzma).
 *
 *     Supported modes are:
 *        * 'r', 'rt', 'w', 'wt' for text mode read and write.
 *        * 'rb, 'wb' for binary read and write.
 *
 *     Parameters
 *     ==========
 *     filename: str
 *         Path to the file to open, including any extensions that indicate
 *         the compression used.
 *     mode: str
 *         Mode to open the file, same as for builtin ``open``, e.g 'r', 'w'.
 *
 *     Returns
 *     =======
 *     fd: file
 *         File-like object open with the specified mode.
 */
//fun openWithCompression(file: File/*, mode: str = 'r'*/) -> IO:
//"""
//
//    """
//
//# Compressed formats sometimes default to binary, so force text mode.
//if mode == 'r':
//mode = 'rt'
//elif mode == 'w':
//mode = 'wt'
//elif mode == 'a':
//mode = 'at'
//
//root, compression = get_compression(filename)
//
//if compression == 'gz':
//return gzip.open(filename, mode=mode)  # type: ignore[return-value]
//elif compression == 'bz2':
//return bz2.open(filename, mode=mode)
//elif compression == 'xz':
//return lzma.open(filename, mode)
//else:
//# Either None or unknown string
//return open(filename, mode)

/**
 * Write Atoms object(s) to file.
 *
 *     filename: str or file
 *         Name of the file to write to or a file descriptor.  The name '-'
 *         means standard output.
 *     images: Atoms object or list of Atoms objects
 *         A single Atoms object or a list of Atoms objects.
 *     format: str
 *         Used to specify the file-format.  If not given, the
 *         file-format will be taken from suffix of the filename.
 *     parallel: bool
 *         Default is to write on master only.  Use parallel=False to write
 *         from all slaves.
 *     append: bool
 *         Default is to open files in 'w' or 'wb' mode, overwriting
 *         existing files.  In some cases opening the file in 'a' or 'ab'
 *         mode (appending) is useful,
 *         e.g. writing trajectories or saving multiple Atoms objects in one file.
 *         WARNING: If the file format does not support multiple entries without
 *         additional keywords/headers, files created using 'append=True'
 *         might not be readable by any program! They will nevertheless be
 *         written without error message.
 *
 *     The use of additional keywords is format specific. write() may
 *     return an object after writing certain formats, but this behaviour
 *     may change in the future.
 */
fun write(filename: File,
          images: Any, //Union[Atoms, Sequence[Atoms]]
          format: String? = null,
          parallel: Boolean = true,
          append: Boolean = false,
          kwargs: Map<String, Any> = emptyMap()) {

    //    if isinstance(filename, PurePath):
    //    filename = str(filename)
    //
    //    if isinstance(filename, str):
    //    fd = None
    //    if filename == '-':
    //    fd = sys.stdout
    //    filename = None  # type: ignore[assignment]
    //    elif format is None:
    //    format = filetype(filename, read=False)
    //    assert isinstance(format, str)
    //    else:
    //    fd = filename  # type: ignore[assignment]
    //    if format is None:
    //    try:
    //        format = filetype(filename, read=False)
    //        assert isinstance(format, str)
    //        except UnknownFileTypeError:
    //        format = None
    //        filename = None  # type: ignore[assignment]

    val format = format ?: "json"
    val io = getIoFormat(format)

    return write(filename, null, format, io, images, parallel, append, kwargs)
}

private fun write(file: File, fd: Any?, format: String, io: IOFormat, images: Any,
                  parallel: Boolean? = null, append: Boolean = false, kwargs: Map<String, Any> = emptyMap()) {

    var images = images
    if (images is Atoms)
        images = listOf(images)

    if (io.single) {
        if ((images as List<Atoms>).size > 1)
            error("$format-format can only store 1 Atoms object.")
        images = images[0]
    }

    check(io.canWrite) { "Can't write to $format-format" }

    // Special case for json-format:
    if (format == "json" && ((images as List<Atoms>).size > 1 || append))
        TODO()
    //        if filename is not None:
    //            return io.write(filename, images, append=append, **kwargs)
    //    raise ValueError("Can't write more than one image to file-descriptor "
    //    'using json-format.')

    return when {
        io.acceptsFd -> {
            val openNew = fd == null
            //        try :
            //            if openNew:
            //            mode = 'wb' if io.isbinary else 'w'
            //            if append:
            //            mode = mode.replace('w', 'a')
            //            fd = open_with_compression(filename, mode)
            //            # XXX remember to re -enable compressed open
            //            # fd = io.open(filename, mode)
            return io.write(file, images, kwargs)
            //            finally:
            //            if openNew and fd is not None :
            //            fd.close()
        }
        else -> {
            TODO()
            //            if fd is not None:
            //            raise ValueError("Can't write {}-format to file-descriptor"
            //                .format(format))
            //            if io.can_append:
            //            return io.write(filename, images, append=append, **kwargs)
            //            elif append:
            //            raise ValueError("Cannot append to {}-format, write-function "
            //            "does not support the append keyword."
            //                .format(format))
            //            else:
            //            return io.write(filename, images, **kwargs)
        }
    }
}

