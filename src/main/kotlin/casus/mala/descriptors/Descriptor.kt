package casus.mala.descriptors

import casus.mala.common.Parameters
import casus.mala.common.PhysicalData

/** Base class for all descriptors available in MALA.
 *
 *  Descriptors encode the atomic fingerprint of a DFT calculation. */
abstract class Descriptor(parameters: Parameters): PhysicalData(parameters) {

    companion object {
        infix fun from(params: Parameters): PhysicalData {
            val descriptors = when(params.descriptors.type) {
                Parameters.Descriptors.Type.Bispectrum -> Bispectrum(params)
                Parameters.Descriptors.Type.AtomicDensity -> AtomicDensity(params)
                Parameters.Descriptors.Type.MinterpyDescriptors -> TODO("")
            }
            return descriptors
        }
    }
}