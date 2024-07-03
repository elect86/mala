package casus.mala.dataHandling

import java.io.File
import java.nio.file.Path

enum class Units { `1|(eV*A^3)`, `1|(Ry*Bohr^3)` }

operator fun File.div(string: String) = resolve(string)