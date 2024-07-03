package casus.mala.dataHandling

import casus.mala.common.Parameters
import casus.mala.descriptors.Descriptor
import casus.mala.targets.Target

/** Base class for all data handling (loading, shuffling, etc.). */
abstract class DataHandlerBase(parameters: Parameters,
                               targetCalculator: Target? = null,
                               descriptorCalculator: Descriptor? = null) {

    val parameters: Parameters.Data = parameters.data
    val useHorovod = parameters.useHorovod

    // Calculators used to parse data from compatible files.
    val targetCalculator = targetCalculator ?: Target.from(parameters)
    val descriptorCalculator = descriptorCalculator ?: Descriptor.from(parameters)

    // Dimensionalities of data.
    var inputDimension = 0L
    var outputDimension = 0L
    var nrSnapshots = 0

    /** Reset the entire data pipeline.
     *
     *  Useful when doing multiple investigations in the same python file. */
    open fun clearData() = parameters.snapshotDirectoriesList.clear()

    operator fun plusAssign(snapshot: Snapshot) {
        parameters.snapshotDirectoriesList += snapshot
    }

    //#############################
    // Private methods
    //#############################

    // Loading data
    //#####################

    // just a comfortable redirect for subclasses to not have to pass the needed argument for super methods
    open fun checkSnapshots() = checkSnapshots(null)

    /** Check the snapshots for consistency. */
    fun checkSnapshots(comm: String?) {
        nrSnapshots = parameters.snapshotDirectoriesList.size

        // Read the snapshots using a memorymap to see if there is consistency.
        var firstsnapshot = true
        for (snapshot in parameters.snapshotDirectoriesList) {
            //###################
            // Descriptors.
            //###################

            println(
                "Checking descriptor file " + snapshot.inputNpyFile + " at " + snapshot.inputNpyDirectory,
                //            min_verbosity = 1,
            )
            var tmpDimension = when (snapshot.snapshotType) {
                Snapshot.Type.numpy -> descriptorCalculator.readDimensionsFromNumpyFile(snapshot.inputNpyFile)
                Snapshot.Type.openpmd -> TODO()
                //                descriptorCalculator.read_dimensions_from_openpmd_file(
                //                os.path.join(
                //                    snapshot.input_npy_directory, snapshot.input_npy_file
                //                ),
                //                comm = comm,
                //            )
                else -> error("Unknown snapshot file type.")
            }

            // get the snapshot feature dimension - call it input dimension
            // for flexible grid sizes only this need be consistent
            val tmpInputDimension = tmpDimension.lastDimension
            val tmpGridDim = tmpDimension[0, 3]!!
            snapshot.gridDimensions = tmpGridDim
            snapshot.gridSize = snapshot.gridDimensions.prod()
            if (firstsnapshot)
                inputDimension = tmpInputDimension
            else
                check(inputDimension == tmpInputDimension) { "Invalid snapshot entered at " + snapshot.inputNpyFile }
            // ###################
            //  Targets.
            // ###################

            println(
                "Checking targets file " + snapshot.outputNpyFile + " at " + snapshot.outputNpyDirectory,
                //                min_verbosity = 1,
            )
            tmpDimension = when (snapshot.snapshotType) {
                Snapshot.Type.numpy -> targetCalculator.readDimensionsFromNumpyFile(snapshot.outputNpyFile)
                else -> TODO() /*targetCalculator.read_dimensions_from_openpmd_file(
                                os.path.join(
                                    snapshot.output_npy_directory,
                                    snapshot.output_npy_file,
                                ),
                                comm = comm,
                            )*/
            }
            // The first snapshot determines the data size to be used .
            // We need to make sure that snapshot size is consistent .
            val tmpOutputDimension = tmpDimension.lastDimension
            if (firstsnapshot)
                outputDimension = tmpOutputDimension
            else check(outputDimension == tmpOutputDimension) { "Invalid snapshot entered at " + snapshot.outputNpyFile }

            check(tmpDimension[0, 3].prod() == snapshot.gridSize) { "Inconsistent snapshot data provided." }

            if (firstsnapshot)
                firstsnapshot = false
        }
    }
}