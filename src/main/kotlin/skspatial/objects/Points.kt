package skspatial.objects

/**
 *     Multiple points in space implemented as a 2D array.
 *
 *     The array is a subclass of :class:`numpy.ndarray`.
 *     Each row in the array represents a point in space.
 *
 *     Parameters
 *     ----------
 *     points : array_like
 *         (N, D) array representing N points with dimension D.
 *
 *     Raises
 *     ------
 *     ValueError
 *         If the array is empty, the values are not finite,
 *         or the dimension is not two.
 *
 *     Examples
 *     --------
 *     >>> from skspatial.objects import Points
 *
 *     >>> points = Points([[1, 2, 0], [5, 4, 3]])
 *
 *     >>> points
 *     Points([[1, 2, 0],
 *             [5, 4, 3]])
 *
 *     >>> points.dimension
 *     3
 *
 *     The object inherits methods from :class:`numpy.ndarray`.
 *
 *     >>> points.mean(axis=0)
 *     array([3. , 3. , 1.5])
 *
 *     >>> Points([])
 *     Traceback (most recent call last):
 *     ...
 *     ValueError: The array must not be empty.
 *
 *     >>> import numpy as np
 *
 *     >>> Points([[1, 2], [1, np.nan]])
 *     Traceback (most recent call last):
 *     ...
 *     ValueError: The values must all be finite.
 *
 *     >>> Points([1, 2, 3])
 *     Traceback (most recent call last):
 *     ...
 *     ValueError: The array must be 2D.
 */
class Points(array: Array<FloatArray>): BaseArray(array) {
}