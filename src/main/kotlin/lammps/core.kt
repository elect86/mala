package lammps

import jdk.jfr.MemoryAddress
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

    val address: MemorySegment = Arena.ofConfined().use { offHeap ->
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
}