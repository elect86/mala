package casus.mala.network

import ai.djl.engine.Engine.Companion.getEngine
import casus.mala.common.Parameters
import javax.xml.transform.Transformer

/**
 * Central network class for this framework, based on pytorch.nn.Module.
 *
 * The correct type of neural network will automatically be instantiated
 * by this class if possible. You can also instantiate the desired
 * network directly by calling upon the subclass.
 *
 * Parameters
 * ----------
 * @param parans: Parameters used to create this neural network. */
open class Network(params: Parameters) {

    // copy the network params from the input parameter object
    val useHorovod = params.useHorovod
    val miniBatchSize = params.running.miniBatchSize
    val params = params.network

    init {
        // if the user has planted a seed (for comparibility purposes) we should use it.
        params.manualSeed?.let {
            getEngine("PyTorch").setRandomSeed(it)
            //            torch.cuda.manual_seed(params.manual_seed)
        }
    }
    //    # initialize the parent class
    //    super(Network, self).__init__()

    //    # Mappings for parsing of the activation layers.
    //    self.activation_mappings = {
    //        "Sigmoid": nn.Sigmoid,
    //        "ReLU": nn.ReLU,
    //        "LeakyReLU": nn.LeakyReLU,
    //        "Tanh": nn.Tanh,
    //    }

    // initialize the layers
    val numberOfLayers = this.params.layerSizes.size - 1

    //    val lossFunc: Parameters.LossFunction
    init {
        // initialize the loss function
        if (this.params.lossFunction == Parameters.LossFunction.mse)
            TODO() // lossFunc = functional.mse_loss
        else
            error("Unsupported loss function.")
    }

    companion object {

        infix fun from(params: Parameters): Network =

        // Check if we're accessing through base class.
            // If not, we need to return the correct
            when (params.network.nnType) {
                Parameters.Nn.`feed-forward` -> FeedForwardNet(params)
                Parameters.Nn.transformer -> TODO() //Transformer(params)
                Parameters.Nn.lstm -> TODO() // LSTM (params)
                Parameters.Nn.gru -> TODO() //GRU(params)
            }
    }
}