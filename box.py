import sympy as y
from . import hz
from .ts import *

Vc = y.symbols('vc')

# derived quantities
Qtc = Qts * y.sqrt(1 + Vas / Vc)
fc = fs * y.sqrt(1 + Vas / Vc)
f3 = fc * y.sqrt((Qtc ** -2 - 2 + y.sqrt((Qtc ** -2 - 2) ** 2 + 4)) / 2)

# dbSPL for Xmax at given frequency
Ro_c = 1.18 / 345.0 # constant
Vd = sD * y.pi / 4 * Xmax / 1e7
Xspl = 112 + 10 * y.log(4 * y.pi ** 3 * Ro_c * Vd ** 2 * hz ** 4) / y.log(10)

# high frequency
hf_f3 = 1 / y.sqrt(Mmd * (1/(1/(Le / (Bl ** 2)) + (1/Cms))))
#hf_q = (1 / (1 / Rms + 1 / Re)) * y.sqrt(1/((1 / Le) + (1/Cms)) / (Bl**2) * Mms)

hf_q = 1 / ((1 / Rms + 1 / Re) * hf_f3)

