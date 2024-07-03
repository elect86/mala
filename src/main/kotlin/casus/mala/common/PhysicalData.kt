package casus.mala.common

import ai.djl.ndarray.NDArray
import ai.djl.ndarray.NDList
import ai.djl.ndarray.NDManager
import ai.djl.ndarray.types.Shape
import casus.mala.dataHandling.Units
import java.nio.file.Path
import kotlin.io.path.readBytes

/** Base class for physical data.
 *
 *  Implements general framework to read and write such data to and from files. */
abstract class PhysicalData(val paramsInterface: ParametersInterface) {

    var gridDimensions = IntArray(3)

    /** Get a string that describes the data (for e.g. metadata). */
    abstract val dataName: String

    val np = NDManager.newBaseManager() // TODO, remember to close it

    //#############################
    // Read functions.
    //   Write functions for now not implemented at this level
    //   because there is no need to.
    //#############################

    /** Read the data from a numpy file.
     *
     *  Parameters
     *  ----------
     *  @param path: Path to the numpy file.
     *
     *  @param units: Units the data is saved in.
     *
     *  @param array: If not None, the array to save the data into. The array has to be 4-dimensional.
     *
     *  @return: If array is None, a numpy array containing the data.
     *  Elsewise, None, as the data will be saved into the provided array. */
    fun readFromNumpyFile(path: Path, units: String? = null, array: NDArray? = null, reshape: Boolean = false): NDArray? =
        when (array) {
            null -> TODO()
            //        loaded_array = np.load(path)[:, :, :, self._feature_mask() :]
            //        self._process_loaded_array(loaded_array, units = units)
            //        return loaded_array
            else -> {
                if (reshape) {
                    val arrayDims = array.shape
                    array[":, :"] = np.load(path)[0][":, :, :, $featureMask :"].reshape(arrayDims)
                } else
                    array[":, :, :, :"] = np.load(path)[":, :, :, $featureMask :"]
                processLoadedArray(array, units = units)
                null
            }
        }

    /** Read only the dimensions from a numpy file.
     *
     * Parameters
     * ----------
     * path : string
     * Path to the numpy file.
     *
     * read_dtype : bool
     * If True, the dtype is read alongside the dimensions. */
    fun readDimensionsFromNumpyFile(path: Path, readDtype: Boolean = false): Shape =
        NDManager.newBaseManager().use {
            val loadedArray = NDList.decode(it, path.readBytes())[0]!!
            when {
                readDtype -> TODO() // return (self._process_loaded_dimensions(np.shape(loadedArray)), loadedArray.dtype, )
                else -> loadedArray.shape
            }
        }

    //#############################
    // Class-specific reshaping, processing, etc. of data.
    //    Has to be implemented by the classes themselves. E.g. descriptors may
    //    need to cut xyz-coordinates, LDOS/density may need unit conversion.
    //#############################

    protected abstract fun processLoadedArray(array: NDArray, units: String? = null)

    //    @abstractmethod
    //    def _process_loaded_dimensions(self, array_dimensions):
    //    pass

    //    @abstractmethod
    //    def _set_feature_size_from_array(self, array):
    //    pass

    //    def _process_additional_metadata(self, additional_metadata):
    //    pass

    abstract val featureMask: Int
}