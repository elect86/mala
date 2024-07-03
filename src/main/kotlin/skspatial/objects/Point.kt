package skspatial.objects

/**
 *     A point in space implemented as a 1D array.
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
 *         Dimension of the point.
 *
 *     Raises
 *     ------
 *     ValueError
 *         If the array is empty, the values are not finite,
 *         or the dimension is not one.
 *
 *     Examples
 *     --------
 *     >>> from skspatial.objects import Point
 *
 *     >>> point = Point([1, 2, 3])
 *
 *     >>> point.dimension
 *     3
 *
 *     The object inherits methods from :class:`numpy.ndarray`.
 *
 *     >>> point.mean()
 *     array(2.)
 *
 *     >>> Point([])
 *     Traceback (most recent call last):
 *     ...
 *     ValueError: The array must not be empty.
 *
 *     >>> import numpy as np
 *
 *     >>> Point([1, 2, np.nan])
 *     Traceback (most recent call last):
 *     ...
 *     ValueError: The values must all be finite.
 *
 *     >>> Point([[1, 2], [3, 4]])
 *     Traceback (most recent call last):
 *     ...
 *     ValueError: The array must be 1D.
 */
class Point(array: FloatArray) : BaseArray1D(array)