package skspatial.objects

import ase.io.dot
import ase.io.norm
import ase.minus

/**
 *     A vector implemented as a 1D array.
 *
 *     The array is a subclass of :class:`numpy.ndarray`.
 *
 *     Parameters
 *     ----------
 *     array : array_like
 *         Input array.
 *
 *     Attributes
 *     ----------
 *     dimension : int
 *         Dimension of the vector.
 *
 *     Raises
 *     ------
 *     ValueError
 *         If the array is empty, the values are not finite,
 *         or the dimension is not one.
 *
 *     Examples
 *     --------
 *     >>> from skspatial.objects import Vector
 *
 *     >>> vector = Vector([1, 2, 3])
 *
 *     >>> vector.dimension
 *     3
 *
 *     The object inherits methods from :class:`numpy.ndarray`.
 *
 *     >>> vector.mean()
 *     array(2.)
 *
 *     >>> Vector([])
 *     Traceback (most recent call last):
 *     ...
 *     ValueError: The array must not be empty.
 *
 *     >>> import numpy as np
 *
 *     >>> Vector([1, 2, np.nan])
 *     Traceback (most recent call last):
 *     ...
 *     ValueError: The values must all be finite.
 *
 *     >>> Vector([[1, 2], [3, 4]])
 *     Traceback (most recent call last):
 *     ...
 *     ValueError: The array must be 1D.
 */
class Vector(array: FloatArray) : BaseArray1D(array) {

    /**
     *         Return the unit vector in the same direction as the vector.
     *
     *         A unit vector is a vector with a magnitude of one.
     *
     *         Returns
     *         -------
     *         Vector
     *             Unit vector.
     *
     *         Raises
     *         ------
     *         ValueError
     *             If the magnitude of the vector is zero.
     *
     *         Examples
     *         --------
     *         >>> from skspatial.objects import Vector
     *
     *         >>> Vector([1, 0]).unit()
     *         Vector([1., 0.])
     *
     *         >>> Vector([-20, 0]).unit()
     *         Vector([-1.,  0.])
     *
     *         >>> Vector([1, 1]).unit().round(3)
     *         Vector([0.707, 0.707])
     *
     *         >>> Vector([1, 1, 1]).unit().round(3)
     *         Vector([0.577, 0.577, 0.577])
     *
     *         >>> Vector([0, 0]).unit()
     *         Traceback (most recent call last):
     *         ...
     *         ValueError: The magnitude must not be zero.
     */
    fun unit(): Vector {
        val magnitude = (array as FloatArray).norm()

        if (magnitude == 0f)
            error("The magnitude must not be zero.")

        return Vector(FloatArray(3) { array[it] / magnitude })
    }

    /**
     *         Return the scalar projection of an other vector onto the vector.
     *
     *         Parameters
     *         ----------
     *         other : array_like
     *             Other vector.
     *
     *         Returns
     *         -------
     *         np.float64
     *             Scalar projection.
     *
     *         Examples
     *         --------
     *         >>> from skspatial.objects import Vector
     *
     *         >>> Vector([0, 1]).scalar_projection([2, 1])
     *         1.0
     *
     *         >>> Vector([-1, -1]).scalar_projection([1, 0]).round(3)
     *         -0.707
     *
     *         >>> Vector([0, 100]).scalar_projection([9, 5])
     *         5.0
     *
     *         >>> Vector([5, 0]).scalar_projection([-10, 3])
     *         -10.0
     */
    infix fun scalarProjection(other: Vector): Float = (unit().array as FloatArray).dot(other.array as FloatArray)

    companion object {
        /**
         *         Instantiate a vector from point A to point B.
         *
         *         Parameters
         *         ----------
         *         point_a, point_b : array_like
         *             Points defining the vector.
         *
         *         Returns
         *         -------
         *         Vector
         *             Vector from point A to point B.
         *
         *         Examples
         *         --------
         *         >>> from skspatial.objects import Vector
         *
         *         >>> Vector.from_points([0, 0], [1, 0])
         *         Vector([1, 0])
         *
         *         >>> Vector.from_points([5, 2], [-2, 8])
         *         Vector([-7,  6])
         *
         *         >>> Vector.from_points([3, 1, 1], [7, 7, 0])
         *         Vector([ 4,  6, -1])
         */
        fun fromPoints(pointA: Point, pointB: Point): Vector = fromPoints(pointA.array as FloatArray, pointB.array as FloatArray)
        fun fromPoints(pointA: FloatArray, pointB: FloatArray): Vector {
            val arrayVectorAB = pointB - pointA
            return Vector(arrayVectorAB)
        }
    }
}