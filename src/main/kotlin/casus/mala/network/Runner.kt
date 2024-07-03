package casus.mala.network

import casus.mala.common.ParametersNetwork
import casus.mala.common.Parameters
import casus.mala.common.ParametersRunning
import casus.mala.dataHandling.DataHandler

/**
 *  Parent class for all classes that in some sense "run" the network.
 *
 * That can be training, benchmarking, inference, etc.
 */
open class Runner
/**
 * Parameters
 * ----------
 * @param params: Parameters used to create this Runner object.
 *
 * @param network: Network which is being run.
 *
 * @data: DataHandler holding the data for the run. */
constructor(val parametersFull: Parameters,
            val network: Network,
            val data: DataHandler) {

    val parameters: ParametersRunning = parametersFull.running

}