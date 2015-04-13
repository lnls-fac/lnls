#!/usr/bin/env python3

import lnls as _lnls
import os as _os
import phase1_parameters
import _input
import _plot
import _optics

twiss_fname = _os.path.join(_lnls.folder_data,'sirius','si', 'beam_dynamics','oficial','v07','c05','multi.cod.tune.coup','trackcpp','rms01','twiss.txt')

phase1_parameters = phase1_parameters.IBSParameters()

elements, twiss, parameters = _input.read_twiss_file(twiss_fname, phase1_parameters)

_optics.calc_rf_acceptance(parameters)
_optics.calc_emittances(parameters)
print(parameters)



# _bo.calc_rf_acceptance(parameters)
# print(100 * parameters.rf_acceptance)
# print(parameters)
# _plot.plot_twiss(twiss)
