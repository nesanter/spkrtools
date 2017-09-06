import sympy as y
from . import ts, box, filters
from . import hz

class System:
    """
    Container for sealed box single driver system

    Attributes:
        values: dictionary of initial T/S values
    """

    def __init__(self, sD, fs, Qts, Qes, Qms, Re, Le, Vas, Xmax, Vc):
        """
        Initialize new System using given T/S values
        """
        self.values = {
                ts.sD : sD,
                ts.fs : fs,
                ts.Qts : Qts,
                ts.Qes : Qes,
                ts.Qms : Qms,
                ts.Vas : Vas,
                ts.Xmax : Xmax,
                ts.Re : Re,
                ts.Le : Le,
                box.Vc : Vc
            }

    def f3(self, N=True):
        """
        Solve for system f3

        Args:
            N: return as reduced value (default True)

        Returns:
            f3
        """
        r = box.f3.subs(self.values)
        if N:
            return y.N(r)
        else:
            return r

    def fc(self, N = True):
        """
        Solve for system fc

        Args:
            N: return as reduced value (default True)

        Returns:
            fc
        """
        r = box.fc.subs(self.values)
        if N:
            return y.N(r)
        else:
            return r

    def Qtc(self, N = True):
        """
        Solve for system Qtc

        Args:
            N: return as reduced value (default True)

        Returns:
            Qtc
        """
        r = box.Qtc.subs(self.values)
        if N:
            return y.N(r)
        else:
            return r

    def lt(self, f1, q1, sr = 48000):
        """
        Calculate a Linkwitz Transform for the system

        Args:
            f1: adjusted f3
            q1: adjusted Qtc

        Returns:
            Tuple of SOS coefficients
        """
        return filters.lt(self.f3(N=False), self.Qtc(N=False), f1, q1, sr)

    def dB(self, lt = None, sr = 48000):
        """
        Calculate the normalized low-frequency acoustic response,
        optionally including a transform

        Args:
            lt: additional transform in SOS form (default None)
            sr: transform sample rate in Hz (default 48000)

        Returns:
            Response in decibels SPL dependent on 'hz'
        """
        hp = filters.hp(self.f3(N=False), self.Qtc(N=False), sr)
        if lt:
            return filters.freqz((hp, lt), hz / (sr / 2), dB=True)
        else:
            return filters.freqz((hp,), hz / (sr / 2), dB=True)

    def dB_at(self, athz, lt = None, sr = 48000, N = True):
        """
        Solve dB() for given frequency

        Args:
            athz: frequency to solve at in hz
            lt: as per dB()
            sr: as per dB()
            N: return as reduced value (default True)

        Returns:
            solution in decibels SPL
        """
        r = self.dB(lt, sr).subs({hz:athz})
        if N:
            return y.N(r)
        else:
            return r

    def f_xlimited(self, dBSPL, digits=10):
        """
        Solve for excursion-limited minimum frequency at given dbSPL

        Note that there is no option to include a transform,
        as SymPy fails to find a solution for that case

        Args:
            dbSPL: given dbSPL

        Returns:
            frequency in Hz
        """
        return tuple(y.N(y.solveset(
            box.Xspl.subs(self.values) - dBSPL, hz, domain=y.S.Reals) &
            y.Interval.Ropen(0, y.S.Infinity), digits))

    def lt_spl_xlimited(self, lt, athz = 20, N=True):
        """
        Solve for excursion-limited SPL at given frequency
        with given transform

        This is designed to test safe maximum SPL for when
        a Linkwitz Transform is applied

        Args:
            lt: the transform in SOS form
            athz: minimum test frequency in Hz (default 20)
            N: return reduced value (default True)

        Returns:
            maximum dbSPL
        """
        xbase = box.Xspl.subs(self.values).subs({hz : athz})
        adjust = y.N(filters.freqz((lt,), hz / 24000, dB=True).subs({hz:athz}))

        r = (xbase - adjust)
        if N:
            return y.N(r)
        else:
            return r

    def datasheet(self):
        values = self.values
        print("---- T/S -------")
        print("sD ", values[ts.sD], "cm2")
        print("fs ", values[ts.fs], "Hz")
        print("Qes", values[ts.Qes])
        print("Qms", values[ts.Qms])
        print("Qts", values[ts.Qts])
        print("Vas", values[ts.Vas], "l")
        print("Re ", values[ts.Re], "ohm")
        print("Le ", values[ts.Le], "mH")
        print("Xmx", values[ts.Xmax], "mm")
        print()
        print("---- Box -------")
        print("Sealed")
        print("Vc ", values[box.Vc], "l")
        print()
        print("---- Derived ---")
        print("Cms", y.N(ts.Cms.subs(values) * 1e3, 4), "mm/N")
        print("Rms", y.N(ts.Rms.subs(values), 4), "Ns/m")
        print("Mmd", y.N(ts.Mmd.subs(values) * 1e3, 4), "g")
        print("Bl ", y.N(ts.Bl.subs(values), 4), "Tm")
        print()
        print("---- System ----")
        print("f3 ", y.N(self.f3(N=False),3), "Hz")
        print("fc ", y.N(self.fc(N=False),3), "Hz")
        print("Qtc", y.N(self.Qtc(N=False),3))
        fX7 = self.f_xlimited(70, digits=3)
        fX8 = self.f_xlimited(80, digits=3)
        fX9 = self.f_xlimited(90, digits=3)
        if len(fX7) == 0:
            print("fX7", "(none)", "Hz @ 70dBSPL")
        else:
            print("fX7", fX7[0], "Hz @ 70dBSPL")
        if len(fX8) == 0:
            print("fX8", "(none)", "Hz @ 80dBSPL")
        else:
            print("fX8", fX8[0], "Hz @ 80dBSPL")
        if len(fX9) == 0:
            print("fX9", "(none)", "Hz @ 90dBSPL")
        else:
            print("fX9", fX9[0], "Hz @ 90dBSPL")
        print()
        print("---- Linkwitz --")
        print("Butterworth, 18Hz")
        print("SPL", y.N(self.lt_spl_xlimited(lt=self.lt(18, y.sqrt(0.5)), athz=20), 3), "dBSPL @ 20Hz)")
        print("Butterworth, 33Hz")
        print("SPL", y.N(self.lt_spl_xlimited(lt=self.lt(33, y.sqrt(0.5)), athz=20), 3), "dBSPL @ 20Hz)")
        print("SPL", y.N(self.lt_spl_xlimited(lt=self.lt(33, y.sqrt(0.5)), athz=33), 3), "dBSPL @ 33Hz)")
        print("Q = 0.5, 55Hz")
        print("SPL", y.N(self.lt_spl_xlimited(lt=self.lt(55, 0.5), athz=20), 3), "dBSPL @ 20Hz)")
        print("SPL", y.N(self.lt_spl_xlimited(lt=self.lt(55, 0.5), athz=33), 3), "dBSPL @ 33Hz)")
        print("SPL", y.N(self.lt_spl_xlimited(lt=self.lt(55, 0.5), athz=55), 3), "dBSPL @ 55Hz)")
        print("Butterworth, fs")
        fs = values[ts.fs]
        print("SPL", y.N(self.lt_spl_xlimited(lt=self.lt(fs, y.sqrt(0.5)), athz=20), 3), "dBSPL @ 20Hz)")
        print("SPL", y.N(self.lt_spl_xlimited(lt=self.lt(fs, y.sqrt(0.5)), athz=33), 3), "dBSPL @ 33Hz)")
        print("SPL", y.N(self.lt_spl_xlimited(lt=self.lt(fs, y.sqrt(0.5)), athz=55), 3), "dBSPL @ 55Hz)")
        print("SPL", y.N(self.lt_spl_xlimited(lt=self.lt(fs, y.sqrt(0.5)), athz=fs), 3), "dBSPL @ fs)")

class Filter:
    """
    Container class for SOS filter chain

    Attributes:
        sos: list-of-tuples (Nx6) coefficients
        sr: sample rate
    """

    def __init__(self, *sos, sr = 48000):
        self.sos = sos
        self.sr = sr

    def freqz(self, dB = True):
        return filters.freqz(self.sos, dB)

#    def __str__(self):
#        s = ''
#        for section in sos:
            
