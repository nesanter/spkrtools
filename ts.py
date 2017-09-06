import sympy as y
from . import hz

# driver parameters
sD, Qts, Qms, Qes, Re, Le, Vas, fs, Xmax = y.symbols('sD Qts Qms Qes Re Le Vas fs Xmax') #cm2, unitless, l, hz, mm

# T/S derived quantities
Cms = Vas / (1.42 * sD ** 2) # 1.42 is adiabatic bulk modulus of air (fun fact)
Mms = ((2 * y.pi * fs) ** -2) / Cms
Rms = (1 / Qms) * y.sqrt(Mms / Cms)
Bl = y.sqrt(2 * y.pi * fs * Mms * Re / Qes)
Mmd = Mms - (2.67 * (y.sqrt(1e-2 * sD) / y.pi) ** 3 * 1.225e-3)

