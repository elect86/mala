package lammps

import lammps.library_h.Style
import lammps.library_h.Type
import java.io.File
import java.lang.foreign.Arena
import java.lang.foreign.MemorySegment
import java.lang.foreign.ValueLayout
import java.util.ArrayList


/**
 * Create an instance of the LAMMPS Python class.
 *
 *   .. _mpi4py_docs: https://mpi4py.readthedocs.io/
 *
 *   This is a Python wrapper class that exposes the LAMMPS C-library
 *   interface to Python.  It either requires that LAMMPS has been compiled
 *   as shared library which is then dynamically loaded via the ctypes
 *   Python module or that this module called from a Python function that
 *   is called from a Python interpreter embedded into a LAMMPS executable,
 *   for example through the :doc:`python invoke <python>` command.
 *   When the class is instantiated it calls the :cpp:func:`lammps_open`
 *   function of the LAMMPS C-library interface, which in
 *   turn will create an instance of the :cpp:class:`LAMMPS <LAMMPS_NS::LAMMPS>`
 *   C++ class.  The handle to this C++ class is stored internally
 *   and automatically passed to the calls to the C library interface.
 *
 *   :param name: "machine" name of the shared LAMMPS library ("mpi" loads ``liblammps_mpi.so``, "" loads ``liblammps.so``)
 *   :type  name: string
 *   :param cmdargs: list of command line arguments to be passed to the :cpp:func:`lammps_open` function.  The executable name is automatically added.
 *   :type  cmdargs: list
 *   :param ptr: pointer to a LAMMPS C++ class instance when called from an embedded Python interpreter.  None means load symbols from shared library.
 *   :type  ptr: pointer
 *   :param comm: MPI communicator (as provided by `mpi4py <mpi4py_docs_>`_). ``None`` means use ``MPI_COMM_WORLD`` implicitly.
 *   :type  comm: MPI_Comm
 */
class Lammps(name: String = "", cmdArgs: ArrayList<String> = arrayListOf()) {

    // -------------------------------------------------------------------------
    // create an instance of LAMMPS

    val lmp: MemorySegment = Arena.ofConfined().use { offHeap ->
        cmdArgs.add(0, "lammps")
        // 4. Allocate a region of off-heap memory to store four pointers
        val pointers = offHeap.allocate(ValueLayout.ADDRESS, cmdArgs.size.toLong())
        // 5. Copy the strings from on-heap to off-heap
        for (i in cmdArgs.indices) {
            val cString = offHeap.allocateFrom(cmdArgs[i])
            pointers.setAtIndex(ValueLayout.ADDRESS, i.toLong(), cString)
        }
        library_h.lammps_open_no_mpi(cmdArgs.size, pointers, MemorySegment.NULL)
    }

    /**
     *     Read LAMMPS commands from a file.
     *
     *     This is a wrapper around the :cpp:func:`lammps_file` function of the C-library interface.
     *     It will open the file with the name/path `file` and process the LAMMPS commands line by line until
     *     the end. The function will return when the end of the file is reached.
     *
     *     :param path: Name of the file/path with LAMMPS commands
     *     :type path:  string
     */
    fun file(path: File) {
        Arena.ofConfined().use { offHeap ->
            val pString = offHeap.allocateFrom("/home/elect/PycharmProjects/mala/mala/descriptors/in.bgrid.python")
            library_h.lammps_file(lmp, pString)
        }
    }

    /**
     * Retrieve data from a LAMMPS compute
     *
     *     This is a wrapper around the :cpp:func:`lammps_extract_compute`
     *     function of the C-library interface.
     *     This function returns ``None`` if either the compute id is not
     *     recognized, or an invalid combination of :ref:`cstyle <py_style_constants>`
     *     and :ref:`ctype <py_type_constants>` constants is used. The
     *     names and functionality of the constants are the same as for
     *     the corresponding C-library function.  For requests to return
     *     a scalar or a size, the value is returned, otherwise a pointer.
     *
     *     :param cid: compute ID
     *     :type cid:  string
     *     :param cstyle: style of the data retrieve (global, atom, or local), see :ref:`py_style_constants`
     *     :type cstyle:  int
     *     :param ctype: type or size of the returned data (scalar, vector, or array), see :ref:`py_type_constants`
     *     :type ctype:  int
     *     :return: requested data as scalar, pointer to 1d or 2d double array, or None
     *     :rtype: c_double, ctypes.POINTER(c_double), ctypes.POINTER(ctypes.POINTER(c_double)), or NoneType
     */
    fun extractCompute(cId: String, cStyle: Style, cType: Type) = when(cType) {
        Type.scalar -> TODO()
        Type.vector -> TODO()
        Type.array -> Arena.ofConfined().use {
            library_h.lammps_extract_compute(lmp, it.allocateFrom(cId), cStyle.ordinal, cType.ordinal)
        }
        Type.sizeCols -> TODO()
        else -> TODO()
    }
}