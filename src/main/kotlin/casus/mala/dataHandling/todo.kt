package casus.mala.dataHandling

import ai.djl.ndarray.types.Shape

fun Shape.prod(): Long {
    var result = 1L
    for (d in 0..<dimension())
        result *= get(d)
    return result
}