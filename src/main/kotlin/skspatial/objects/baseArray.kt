package skspatial.objects

/** Private base class for spatial objects based on a single NumPy array. */
abstract class BaseArray(val array: Any) {
    init {
        if (array is Array<*>) {
            check(array.isNotEmpty()) { "The array must not be empty." }
            if (array[0] is FloatArray)
                check((array as Array<FloatArray>).none { a -> a.none { it.isInfinite() } }) { "The values must all be finite." }
        } else {
            if (array is FloatArray) {
                check(array.isNotEmpty()) { "The array must not be empty." }
                check(array.none { it.isInfinite() }) { "The values must all be finite." }
            }
        }
    }
}

/** Private base class for spatial objects based on a single 2D NumPy array. */
abstract class BaseArray1D(array: FloatArray) : BaseArray(array)

/** Private base class for spatial objects based on a single 2D NumPy array. */
abstract class BaseArray2D(array: Array<FloatArray>) : BaseArray(array) {


}