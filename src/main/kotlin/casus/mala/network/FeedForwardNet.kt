package casus.mala.network

import casus.mala.common.Parameters

/** Initialize this network as a feed-forward network. */

class FeedForwardNet(params: Parameters) : Network(params) {

    //    self.layers = nn.ModuleList()

    // If we have only one entry in the activation list, we use it for the entire list.
    // We should NOT modify the list itself. This would break the hyperparameter algorithms.
    var useOnlyOneActivationType = false

    init {
        if (this.params.layerActivations.size == 1)
            useOnlyOneActivationType = true
        else if (this.params.layerActivations.size < numberOfLayers)
            error("Not enough activation layers provided.")
        else if (this.params.layerActivations.size > numberOfLayers)
            println("Too many activation layers provided. The last " + (this.params.layerActivations.size - numberOfLayers) +
                    "activation function(s) will be ignored."/*min_verbosity = 1,*/)

        // Add the layers.
        // As this is a feedforward layer we always add linear layers, and then an activation function
//        for (i in 0 ..< numberOfLayers) {
//            self.layers.append(
//                (
//                        nn.Linear(
//                            self.params.layer_sizes[i],
//                            self.params.layer_sizes[i + 1],
//                        )
//                )
//            )
//            try :
//                if use_only_one_activation_type:
//                self.layers.append(
//                    self.activation_mappings[
//                        self.params.layer_activations[0]
//                    ]()
//                )
//                else:
//                self.layers.append(
//                    self.activation_mappings[
//                        self.params.layer_activations[i]
//                    ]()
//                )
//                except KeyError :
//                raise Exception ("Invalid activation type seleceted.")
//            }
//            // Once everything is done, we can move the Network on the target device.
//            self.to(self.params._configuration["device"])
    }
}