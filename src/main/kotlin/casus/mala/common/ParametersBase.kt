package casus.mala.common

abstract class ParametersBase {
    val gpu = false
    val horovod = false
    val mpi = false
    val device = "cpu"
    val openpmdConfiguration: List<String> = emptyList()
    val openpmdGranularity = 1
    val lammps = true
}

sealed interface ParametersInterface