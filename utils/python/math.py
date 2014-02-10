import numpy as np


def shannon_entropy(p):
    """Calculates shannon entropy in bits.

    **Parameters**

    p : np.array
        array of probabilities

    **Returns**

    shannon entropy in bits
    """
    return -np.sum(p * np.log2(p))
