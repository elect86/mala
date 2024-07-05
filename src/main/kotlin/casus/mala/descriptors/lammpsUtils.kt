package casus.mala.descriptors

/**
 *     Add a dicitionary of LAMMPS arguments in a command line argument string.
 *
 *     Parameters
 *     ----------
 *     cmdargs : list
 *         Command line argument string. Will be mutated by this function.
 *
 *     argdict : dict
 *         Dictionary to be added to LAMMPS command line argument string.
 *
 *     Returns
 *     -------
 *     cmdargs : list
 *         New command line argument string.
 */
fun setCmdlineVars(cmdargs: List<String>, argDict: Map<String, Any>): List<String> =
    cmdargs + argDict.flatMap { (k, v) -> listOf("-var", k, v.toString()) }