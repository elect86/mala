package skspatial.objects

import ase.cross
import kotlin.math.abs

/**
 *     A plane in space.
 *
 *     The plane is defined by a point and a normal vector.
 *
 *     Parameters
 *     ----------
 *     point : array_like
 *         Point on the plane.
 *     direction : array_like
 *         Normal vector of the plane.
 *     kwargs : dict, optional
 *         Additional keywords passed to :meth:`Vector.is_zero`.
 *         This method is used to ensure that the normal vector is not the zero vector.
 *
 *     Attributes
 *     ----------
 *     point : Point
 *         Point on the plane.
 *     normal : Vector
 *         Unit normal vector.
 *     vector : Vector
 *         Same as the normal.
 *     dimension : int
 *         Dimension of the plane.
 *
 *     Raises
 *     ------
 *     ValueError
 *         If the point and vector have different dimensions.
 *         If the vector is all zeros.
 *
 *     Examples
 *     --------
 *     >>> from skspatial.objects import Plane
 *
 *     >>> plane = Plane(point=[0, 0, 0], normal=[0, 0, 5])
 *
 *     >>> plane
 *     Plane(point=Point([0, 0, 0]), normal=Vector([0, 0, 5]))
 *
 *     >>> plane.normal
 *     Vector([0, 0, 5])
 *
 *     The normal can also be accessed with the ``vector`` attribute.
 *
 *     >>> plane.vector
 *     Vector([0, 0, 5])
 *
 *     The plane dimension is the dimension of the point and vector.
 *
 *     >>> plane.dimension
 *     3
 *
 *     >>> Plane([0, 0], [1, 0, 0])
 *     Traceback (most recent call last):
 *     ...
 *     ValueError: The point and vector must have the same dimension.
 *
 *     >>> Plane([1, 1], [0, 0])
 *     Traceback (most recent call last):
 *     ...
 *     ValueError: The vector must not be the zero vector.
 */
class Plane(val point: Point, val normal: Vector) : BaseLinePlane() {

    /**
     *         Return the signed distance from a point to the plane.
     *
     *         Parameters
     *         ----------
     *         point : array_like
     *             Input point.
     *
     *         Returns
     *         -------
     *         np.float64
     *             Signed distance from the point to the plane.
     *
     *         References
     *         ----------
     *         http://mathworld.wolfram.com/Point-PlaneDistance.html
     *
     *         Examples
     *         --------
     *         >>> from skspatial.objects import Plane
     *
     *         >>> plane = Plane([0, 0, 0], [0, 0, 1])
     *
     *         >>> plane.distance_point_signed([5, 2, 0])
     *         0.0
     *
     *         >>> plane.distance_point_signed([5, 2, 1])
     *         1.0
     *
     *         >>> plane.distance_point([5, 2, -4])
     *         4.0
     *         >>> plane.distance_point_signed([5, 2, -4])
     *         -4.0
     */
    fun distancePointSigned(point: Point): Float {
        val vectorToPoint = Vector.fromPoints(this.point, point)
        return normal scalarProjection vectorToPoint
    }

    /**
     *         Return the distance from a point to the plane.
     *
     *         Parameters
     *         ----------
     *         point : array_like
     *             Input point.
     *
     *         Returns
     *         -------
     *         np.float64
     *             Distance from the point to the plane.
     *
     *         References
     *         ----------
     *         http://mathworld.wolfram.com/Point-PlaneDistance.html
     *
     *         Examples
     *         --------
     *         >>> from skspatial.objects import Plane
     *
     *         >>> plane = Plane([0, 0, 0], [0, 0, 1])
     *
     *         >>> plane.distance_point([5, 2, 0])
     *         0.0
     *
     *         >>> plane.distance_point([5, 2, 1])
     *         1.0
     *
     *         >>> plane.distance_point([5, 2, -4])
     *         4.0
     */
    infix fun distancePoint(point: Point): Float = abs(distancePointSigned(point))

    companion object {

        /**
         *         Instantiate a plane from a point and two vectors.
         *
         *         The two vectors span the plane.
         *
         *         Parameters
         *         ----------
         *         point : array_like
         *             Point on the plane.
         *         vector_a, vector_b : array_like
         *             Input vectors.
         *         kwargs : dict, optional
         *             Additional keywords passed to :meth:`Vector.is_parallel`.
         *
         *         Returns
         *         -------
         *         Plane
         *             Plane containing input point and spanned by the two input vectors.
         *
         *         Raises
         *         ------
         *         ValueError
         *             If the vectors are parallel.
         *
         *         Examples
         *         --------
         *         >>> from skspatial.objects import Plane
         *
         *         >>> Plane.from_vectors([0, 0], [1, 0], [0, 1])
         *         Plane(point=Point([0, 0, 0]), normal=Vector([0, 0, 1]))
         *
         *         >>> Plane.from_vectors([0, 0], [1, 0], [2, 0])
         *         Traceback (most recent call last):
         *         ...
         *         ValueError: The vectors must not be parallel.
         */
        fun fromVectors(point: FloatArray, vectorA: FloatArray, vectorB: FloatArray, kwargs: Map<String, Any> = emptyMap()): Plane {
//            val vA = Vector(vectorA)

//            if vA.is_parallel(vector_b, ** kwargs):
//            raise ValueError ("The vectors must not be parallel.")

            // The cross product returns a 3D vector.
            val vectorNormal = vectorA cross vectorB

            // Convert the point to 3D so that it matches the vector dimension.
            val p = Point(point)//.set_dimension(3)

            return Plane(p, Vector(vectorNormal))
        }

        /**
         *         Instantiate a plane from three points.
         *
         *         The three points lie on the plane.
         *
         *         Parameters
         *         ----------
         *         point_a, point_b, point_c: array_like
         *             Three points defining the plane.
         *         kwargs: dict, optional
         *             Additional keywords passed to :meth:`Points.are_collinear`.
         *
         *         Returns
         *         -------
         *         Plane
         *             Plane containing the three input points.
         *
         *         Raises
         *         ------
         *         ValueError
         *             If the points are collinear.
         *
         *         Examples
         *         --------
         *         >>> from skspatial.objects import Plane
         *
         *         >>> Plane.from_points([0, 0], [1, 0], [3, 3])
         *         Plane(point=Point([0, 0, 0]), normal=Vector([0, 0, 3]))
         *
         *         The order of the points affects the direction of the normal vector.
         *
         *         >>> Plane.from_points([0, 0], [3, 3], [1, 0])
         *         Plane(point=Point([0, 0, 0]), normal=Vector([ 0,  0, -3]))
         *
         *         >>> Plane.from_points([0, 0], [0, 1], [0, 3])
         *         Traceback (most recent call last):
         *         ...
         *         ValueError: The points must not be collinear.
         */
        fun fromPoints(pointA: FloatArray, pointB: FloatArray, pointC: FloatArray, kwargs: Map<String, Any> = emptyMap()): Plane {
//            if Points([pointA, pointB, pointC]).are_collinear(** kwargs):
//            raise ValueError ("The points must not be collinear.")

            val vectorAB = Vector.fromPoints(pointA, pointB)
            val vectorAC = Vector.fromPoints(pointA, pointC)

            return fromVectors(pointA, vectorAB.array as FloatArray, vectorAC.array as FloatArray)
        }
    }
}