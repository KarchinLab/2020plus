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


def max_shannon_entropy(n):
    """Returns max possible entropy given "n" mutations.

    The maximum possible entropy is the entropy of the
    uniform distribution. The uniform distribution has
    entropy equal to log(n) (I will use base 2).

    **Parameters**

    n : int
        total mutation counts

    **Returns**

    max possible shannon entropy in bits
    """
    if n <= 0:
        return 0.
    return float(np.log2(n))
