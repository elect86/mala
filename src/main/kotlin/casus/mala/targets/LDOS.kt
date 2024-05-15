package casus.mala.targets

import casus.mala.common.ParametersInterface

/** Postprocessing / parsing functions for the local density of states. */
class LDOS(paramsInterface: ParametersInterface): Target(paramsInterface) {

    override val dataName = "LDOS"
    override val feaureSize = parameters.ldosGridsize

//    @classmethod
//    def from_numpy_file(cls, params, path, units="1/(eV*A^3)"):
//    """
//        Create an LDOS calculator from a numpy array saved in a file.
//
//        Parameters
//        ----------
//        params : mala.common.parameters.Parameters
//            Parameters used to create this LDOS object.
//
//        path : string
//            Path to file that is being read.
//
//        units : string
//            Units the LDOS is saved in.
//
//        Returns
//        -------
//        ldos_calculator : mala.targets.ldos.LDOS
//            LDOS calculator object.
//        """
//    return_ldos_object = LDOS(params)
//    return_ldos_object.read_from_numpy_file(path, units=units)
//    return return_ldos_object
//
//    @classmethod
//    def from_numpy_array(cls, params, array, units="1/(eV*A^3)"):
//    """
//        Create an LDOS calculator from a numpy array in memory.
//
//        By using this function rather then setting the local_density_of_states
//        object directly, proper unit coversion is ensured.
//
//        Parameters
//        ----------
//        params : mala.common.parameters.Parameters
//            Parameters used to create this LDOS object.
//
//        array : numpy.ndarray
//            Path to file that is being read.
//
//        units : string
//            Units the LDOS is saved in.
//
//        Returns
//        -------
//        ldos_calculator : mala.targets.ldos.LDOS
//            LDOS calculator object.
//        """
//    return_ldos_object = LDOS(params)
//    return_ldos_object.read_from_array(array, units=units)
//    return return_ldos_object
//
//    @classmethod
//    def from_cube_file(
//    cls, params, path_name_scheme, units="1/(eV*A^3)", use_memmap=None
//    ):
//    """
//        Create an LDOS calculator from multiple cube files.
//
//        The files have to be located in the same directory.
//
//        Parameters
//        ----------
//        params : mala.common.parameters.Parameters
//            Parameters used to create this LDOS object.
//
//        path_name_scheme : string
//            Naming scheme for the LDOS .cube files. Every asterisk will be
//            replaced with an appropriate number for the LDOS files. Before
//            the file name, please make sure to include the proper file path.
//
//        units : string
//            Units the LDOS is saved in.
//
//        use_memmap : string
//            If not None, a memory mapped file with this name will be used to
//            gather the LDOS.
//            If run in MPI parallel mode, such a file MUST be provided.
//        """
//    return_ldos_object = LDOS(params)
//    return_ldos_object.read_from_cube(
//    path_name_scheme, units=units, use_memmap=use_memmap
//    )
//    return return_ldos_object
//
//    @classmethod
//    def from_xsf_file(
//    cls, params, path_name_scheme, units="1/(eV*A^3)", use_memmap=None
//    ):
//    """
//        Create an LDOS calculator from multiple xsf files.
//
//        The files have to be located in the same directory.
//
//        Parameters
//        ----------
//        params : mala.common.parameters.Parameters
//            Parameters used to create this LDOS object.
//
//        path_name_scheme : string
//            Naming scheme for the LDOS .xsf files. Every asterisk will be
//            replaced with an appropriate number for the LDOS files. Before
//            the file name, please make sure to include the proper file path.
//
//        units : string
//            Units the LDOS is saved in.
//
//        use_memmap : string
//            If not None, a memory mapped file with this name will be used to
//            gather the LDOS.
//            If run in MPI parallel mode, such a file MUST be provided.
//        """
//    return_ldos_object = LDOS(params)
//    return_ldos_object.read_from_xsf(
//    path_name_scheme, units=units, use_memmap=use_memmap
//    )
//    return return_ldos_object
//
//    @classmethod
//    def from_openpmd_file(cls, params, path):
//    """
//        Create an LDOS calculator from an OpenPMD file.
//
//        Supports all OpenPMD supported file endings.
//
//        Parameters
//        ----------
//        params : mala.common.parameters.Parameters
//            Parameters used to create this LDOS object.
//
//        path : string
//            Path to OpenPMD file.
//
//        Returns
//        -------
//        ldos_calculator : mala.targets.ldos.LDOS
//            LDOS calculator object.
//        """
//    return_ldos_object = LDOS(params)
//    return_ldos_object.read_from_openpmd_file(path)
//    return return_ldos_object
}