package casus.mala.dataHandling

import ai.djl.ndarray.types.Shape
import java.io.File
import java.nio.file.Path

class Snapshot(
    // Inputs

    /** File with saved numpy input array. */
    val inputNpyFile: Path,
    /** File with saved numpy output array. */
    val outputNpyFile: Path,

    /** "Function" of the snapshot in the MALA workflow.
     *  Replaces the old approach of MALA to have a separate list.
     *  Default is None. */
    val snapshotFunction: Function,
    var inputUnits: String? = null,
    var outputUnits: String = "1/(eV*A^3)",
    /** File with the output of the original snapshot calculation.
     *  This is only needed when testing multiple snapshots. */
    val calculationOutput: File? = null,
    val snapshotType: Type = Type.openpmd) {

    // All the dimensionalities of the snapshot.
    lateinit var gridDimensions: Shape
    var gridSize = -1L
//    self.input_dimension = None
//    self.output_dimension = None

    /** Directory containing input_npy_directory. */
    val inputNpyDirectory: Path
        get() = inputNpyFile.parent

    /** Directory containing output_npy_file. */
    val outputNpyDirectory: Path
        get() = outputNpyFile.parent

    enum class Type { numpy, openpmd }
    enum class Function {
        /** This snapshot will be a testing snapshot. */
        testing,
        /** This snapshot will be a training snapshot. */
        training,
        /** This snapshot will be a validation snapshot. */
        validation
    }
}