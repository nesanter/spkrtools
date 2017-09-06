import sympy as y
from . import hz

def freqz(sos, j, dB = False):
    """
    Calculate frequency response of SOSs

    Args:
        sos: tuple of transforms in SOS form
        j: dependent variable (typically hz / nyquist)
        dB: return only magnitude in decibels (default False)

    Returns:
        response dependent on j
    """
    resp = 1
    for section in sos:
        z0 = section[0] * y.E ** (-1 * y.I * y.pi * j)
        z1 = section[1] * y.E ** (-2 * y.I * y.pi * j)
        z2 = section[2] * y.E ** (-3 * y.I * y.pi * j)
        p0 = section[3] * y.E ** (-1 * y.I * y.pi * j)
        p1 = section[4] * y.E ** (-2 * y.I * y.pi * j)
        p2 = section[5] * y.E ** (-3 * y.I * y.pi * j)

        z = z0 + z1 + z2
        p = p0 + p1 + p2

        resp *= (z / p)

    if dB:
        return 20 / y.log(10) * y.log(abs(resp))
    else:
        return resp

def lt(f0, q0, f1, q1, sr):
    """
    Calculate a Linkwitz Transform

    Args:
        f0: initial f3
        q0: initial Q
        f1: transformed f3
        q1: transformed Q
        sr: sample rate

    Returns:
        transform in SOS form
    """
    w0 = 2 * y.pi * f0
    w1 = 2 * y.pi * f1
    gn = 2 * sr#(2 * pi) / tan(pi / sr)

    a0 = w0 ** 2 + gn * w0 / q0 + gn ** 2
    a1 = 2 * (w0 ** 2 - gn ** 2)
    a2 = w0 ** 2 - gn * w0 / q0 + gn ** 2

    b0 = w1 ** 2 + gn * w1 / q1 + gn ** 2
    b1 = 2 * (w1 ** 2 - gn ** 2)
    b2 = w1 ** 2 - gn * w1 / q1 + gn ** 2

    return (a0 / b0, a1 / b0, a2 / b0, y.Integer(1), b1 / b0, b2 / b0)

def hp(fc, q, sr):
    """
    Calculate a 2nd order highpass filter

    Args:
        fc: filter cutoff (f3)
        q: filter Q
        sr: sample rate

    Returns:
        transform in SOS form
    """
    w0 = 2 * y.pi * fc / sr
    alpha = y.sin(w0) / (2 * q)
    cw0 = y.cos(w0)

    a0 = (1 + cw0) / 2
    a1 = -(1 + cw0)
    a2 = a0
    b0 = (1 + alpha)
    b1 = (-2 * cw0)
    b2 = (1 - alpha)

    return (a0 / b0, a1 / b0, a2 / b0, y.Integer(1), b1 / b0, b2 / b0)

def lp(fc, q, sr):
    """
    Calculate a 2nd order lowpass filter

    Args:
        fc: filter cutoff (f3)
        q: filter Q
        sr: sample rate

    Returns:
        transform in SOS form
    """
    w0 = 2 * y.pi * fc / sr
    alpha = y.sin(w0) / (2 * q)
    cw0 = y.cos(w0)

    a0 = (1 - cw0) / 2
    a1 = (1 - cw0)
    a2 = a0
    b0 = (1 + alpha)
    b1 = (-2 * cw0)
    b2 = (1 - alpha)

    return (a0 / b0, a1 / b0, a2 / b0, y.Integer(1), b1 / b0, b2 / b0)

def to_filterstr(coefs, gain = 1, leader = None, i = (0,1,2,4,5)):
    s = ''
    if leader:
        s += leader + ' '

    for c in i:
        s += float(coefs[i]).hex() + ' '

    s += str(gain)

    return s

def approx_gain(sos, domain=range(20, 20000)):
    f = y.lambdify(hz, abs(freqz(sos, hz / 24000)))
    g = max((f(n) for n in domain))
    return g

