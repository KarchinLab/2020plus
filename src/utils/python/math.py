import numpy as np


def shannon_entropy(p):
    """Calculates shannon entropy in bits.

    Parameters
    ----------
    p : np.array
        array of probabilities

    Returns
    -------
    shannon entropy in bits
    """
    return -np.sum(p * np.log2(p))


def max_shannon_entropy(n):
    """Returns max possible entropy given "n" mutations.

    The maximum possible entropy is the entropy of the
    uniform distribution. The uniform distribution has
    entropy equal to log(n) (I will use base 2).

    Parameters
    ----------
    n : int
        total mutation counts

    Returns
    -------
    max possible shannon entropy in bits
    """
    if n <= 0:
        return 0.
    return float(np.log2(n))


def kl_divergence(p, q):
    """Compute the Kullback-Leibler (KL) divergence for discrete distributions.

    Parameters
    ----------
    p : np.array
        "Ideal"/"true" Probability distribution
    q : np.array
        Approximation of probability distribution p

    Returns
    -------
    kl : float
        KL divergence of approximating p with the distribution q
    """
    # make sure numpy arrays are floats
    p = p.astype(float)
    q = q.astype(float)

    # compute kl divergence
    kl = np.sum(np.where(p!=0, p*np.log2(p/q), 0))
    return kl


def js_divergence(p, q):
    """Compute the Jensen-Shannon Divergence between two discrete distributions.

    Parameters
    ----------
    p : np.array
        probability mass array (sums to 1)
    q : np.array
        probability mass array (sums to 1)

    Returns
    -------
    js_div : float
        js divergence between the two distrubtions
    """
    m = .5 * (p+q)
    js_div = .5*kl_divergence(p, m) + .5*kl_divergence(q, m)
    return js_div


def js_distance(p, q):
    """Compute the Jensen-Shannon distance between two discrete distributions.

    NOTE: JS divergence is not a metric but the sqrt of JS divergence is a
    metric and is called the JS distance.

    Parameters
    ----------
    p : np.array
        probability mass array (sums to 1)
    q : np.array
        probability mass array (sums to 1)

    Returns
    -------
    js_dist : float
        Jensen-Shannon distance between two discrete distributions
    """
    js_dist = np.sqrt(js_divergence(p, q))
    return js_dist
