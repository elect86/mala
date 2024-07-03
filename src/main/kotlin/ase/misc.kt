package ase

import org.ejml.simple.SimpleMatrix
import kotlin.math.abs
import kotlin.math.floor

fun Array<FloatArray>.any(axis: Int): BooleanArray {
    check(size == 3 && this[0].size == 3 && this[1].size == 3 && this[2].size == 3)
    return when (axis) {
        0 -> BooleanArray(3) { x -> this[x].any { it != 0f } }
        1 -> BooleanArray(3) { this[0][it] != 0f && this[1][it] != 0f && this[2][it] != 0f }
        else -> error("invalid axis $axis")
    }
}

operator fun BooleanArray.not() = BooleanArray(size) { !this[it] }
fun BooleanArray.nonZero() = buildList { for (i in indices) if (!this@nonZero[i]) add(i) }

fun MutableList<FloatArray>.transpose(): Array<FloatArray> {
    val x = size
    val y = this[0].size
    check(y == 3)
    return Array(y) { y -> FloatArray(x) { x -> this[x][y] } }
}

fun Array<FloatArray>.transpose(): Array<FloatArray> {
    val x = size
    val y = this[0].size
    check(y == 3)
    return Array(y) { y -> FloatArray(x) { x -> this[x][y] } }
}

val SimpleMatrix.array: Array<FloatArray>
    get() = Array(numRows) { r -> FloatArray(numCols) { get(r, it).toFloat() } }

fun FloatArray.isClose(other: FloatArray, atol: Float) = BooleanArray(size) { abs(this[it] - other[it]) <= atol }

operator fun Array<IntArray>.times(ints: IntArray): Array<IntArray> {
    check(this[0].size == ints.size)
    return Array(size) { r ->
        IntArray(3) { c -> this[r][c] * ints[c] }
    }
}

fun Array<FloatArray>.flatNonZero(): List<Int> {
    val res = arrayListOf<Int>()
    for (r in indices)
        for (c in this[r].indices)
            if (this[r][c] != 0f)
                res += r * this[r].size + c
    return res
}

operator fun List<Int>.rem(int: Int) = map { it % int }

operator fun Array<FloatArray>.times(ints: IntArray) = Array(size) { x -> FloatArray(this[0].size) { this[x][it] * ints[it] } }
fun Array<FloatArray>.floor() = Array(size) { x -> FloatArray(this[0].size) { floor(this[x][it]) } }
fun Array<FloatArray>.toInt() = Array(size) { x -> IntArray(this[0].size) { this[x][it].toInt() } }
fun Array<IntArray>.zerosLike() = Array(size) { IntArray(this[0].size) }
fun IntArray.argSort(): IntArray = indices.sortedBy { this[it] }.toIntArray()
fun IntArray.bincount(): IntArray = IntArray(size) { count { i -> i == it } }
operator fun IntArray.get(mask: BooleanArray) = filterIndexed { idx, _ -> mask[idx] }.toIntArray()
//operator fun IntArray.times(int: Int) = IntArray(size) {}

val selector: (IntArray) -> Int = { it[0] * 10_000 + it[1] * 100 + it[2] }

operator fun List<IntArray>.times(floats: Array<FloatArray>) = Array(size) { r ->
    FloatArray(floats[0].size) { f ->
        var sum = 0f
        var i = 0
        this[r].sum { it * floats[i++][f] }
//        for (element in this[r])
//            sum += element * floats[i++][it]
//        sum
    }
}

inline fun IntArray.sum(selector: (Int) -> Float): Float {
    var sum = 0f
    for (element in this)
        sum += selector(element)
    return sum
}

operator fun FloatArray.plus(floats: Array<FloatArray>) = Array(floats.size) { r -> FloatArray(size) { this[it] + floats[r][it] } }
operator fun FloatArray.minus(floats: FloatArray) = FloatArray(size) { this[it] - floats[it] }
infix fun FloatArray.cross(other: FloatArray) = floatArrayOf(this[1] * other[2] - this[2] * other[1], this[2] * other[0] - this[0] * other[2], this[0] * other[1] - this[1] * other[0])

fun ArrayList<Float>.argWhere(filter: (Float) -> Boolean): List<Float> = buildList {
    val rec = this@argWhere
    for (i in rec.indices)
        if (filter(rec[i]))
            this += rec[i]
}

fun Int.factorial(): Long {
    var result = 1L
    for (factor in 2..this)
        result *= factor.toLong()
    return result
}

val now
    get() = System.currentTimeMillis() / 1_000f

fun FloatArray.cdist(xb: Array<FloatArray>, metric: String ="euclidean"/*, *, out=None, **kwargs*/) {
    val mA = 1
    val mB = xb.size
    val n = 3
}