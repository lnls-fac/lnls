import matplotlib.pyplot as _plt
import _optics

def plot_twiss(twiss):

    s = _optics.gettwiss(twiss,'s')
    betax = _optics.gettwiss(twiss,'betax')
    betay = _optics.gettwiss(twiss,'betay')
    etax  = _optics.gettwiss(twiss,'etax')

    slimit = (s[-1]/10.0)
    sel = (s <= slimit)
    _plt.plot(s[sel],betax[sel])
    _plt.plot(s[sel],betay[sel])
    _plt.plot(s[sel],100*etax[sel])
    _plt.xlabel('pos [m]'), _plt.ylabel('betax, betay [m] - etay [cm]')
    _plt.grid()
    _plt.xlim((0,slimit))
    _plt.show()
