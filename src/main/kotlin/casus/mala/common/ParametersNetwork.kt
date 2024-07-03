package casus.mala.common

/** Parameters necessary for constructing a neural network. */
class ParametersNetwork {
    /** Type of the neural network that will be used. */
    var nnType = Parameters.Nn.`feed-forward`

    /**  A list of integers detailing the sizes of the layer of the neural network. Please note that the input layer
     *  is included therein. Default: [10,10,10] */
    var layerSizes = IntArray(3) { 10 }

    /** A list of strings detailing the activation functions to be used by the neural network. If the dimension of
     * [layerActivations] is smaller than the dimension of layer_sizes-1, than the first entry is used for all layers. */
    var layerActivations = listOf(Parameters.Activation.Sigmoid)

    /** Loss function for the neural network */
    var lossFunction = Parameters.LossFunction.mse

    /** for LSTM/Gru + Transformer
     *  Number of hidden layers to be used in lstm or gru or transformer nets */
    var numHiddenLayers = 1

    /** for LSTM/Gru
     *  If True hidden and cell state is assigned to zeros for LSTM Network. false will keep the hidden state active */
    var noHiddenState = false

    /** Sets lstm network size based on bidirectional or just one direction */
    var bidirection = false

    /** for transformer net
     *  Dropout rate for transformer net
     *  0.0 ≤ dropout ≤ 1.0 */
    var dropout = 0.1f

    /** Number of heads to be used in Multi head attention network. This should be a divisor of input dimension */
    var numHeads = 10
}