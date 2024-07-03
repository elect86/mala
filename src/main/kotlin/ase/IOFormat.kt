package ase

data class IOFormat(val name: String,
                    val description: String,
                    val code: String,
                    val moduleName: String,
                    val encoding: String? = null) {

    init {
        check(code.length == 2)
        check(code[0] in listOf('+', '1'))
        check(code[1] in listOf('B', 'F', 'S'))
    }

    // (To be set by define_io_format())
    val extensions = ArrayList<String>()
    val globs = ArrayList<String>()
    val magic = ArrayList<String>()
    var magicRegex: String? = null
}