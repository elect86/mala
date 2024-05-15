package casus.mala.dataHandling

import java.io.File

class Snapshot(
    // Inputs
    val inputNpyFile: File,
    val inputNpyDirectory: File,
    val outputNpyFile: File,
    val outputNpyDirectory: File,
    val snapshotFunction: String,
    val inputUnits: String = "",
    val outputUnits: String = "",
    val calculationOutput: String = "",
    val snapshotType: String = "openpmd") {
}