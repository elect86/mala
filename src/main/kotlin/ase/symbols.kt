package ase

fun symbols2numbers(symbols: Any): IntArray {
    var symbols = symbols
    if (symbols is String)
        symbols = TODO()// string2symbols (symbols)
    val numbers = ArrayList<Int>()
    symbols as List<String>
    for(s in symbols) {
//        if isinstance(s, str):
        numbers += atomicNumbers[s]!!
//        else:
//        numbers.append(int(s))
    }
    return numbers.toIntArray()
}