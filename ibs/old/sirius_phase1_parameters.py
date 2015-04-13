import mathphys as _mp

_harmonic_number = 864

parameters = {

    # -- beam --
    'En':  3.0e+9,                   # energy [eV]
    'I0': (100.0/_harmonic_number),  # total current per bunch [mA]

    # -- lattice --
    'C':   518.396,  # lattice circumference [m]

    # -- vaccum chamber pressure --
    'Pmed':    1.0,  # average pressure in the machine [nTorr]

    # -- RF cavities (6 Dampy Cavities) --
    'Vrf':     2.7e+6,           # total RF voltage [volts]
    'hh':      _harmonic_number, # hamonic number
    'Rs0':     19.8e+06,
    'Q0':      27000,
    'betac':   3.83,
    'Detune0': 100e+03,

    # -- harmonic cavity (SLS type SC HHC) --
    'mharm':        3,
    'Q0n':          2e+08,
    'Rsn':          176.8e+08,
    'DetuneHC':     77.7e3,

    'ap':           0.0001739520528,

    # # -- delta radiation integrals due to IDs --
    # 'I1':  0.087350856985953,
    # 'I2':  0.400501181860251,
    # 'I3':  0.029528249837988,
    # 'I4':  -0.132157247653490,
    # 'I5':  0.000010929804409,

    # 'I1':  0.0,
    # 'I2':  0.200481368,
    # 'I3':  0.036505767,
    # 'I4':  5.89612E-07,
    # 'I5':  1.2398E-06,
    # 'I6':  0.0,

    'I1':  0.0,
    'I2':  0.0,
    'I3':  0.0,
    'I4':  0.0,
    'I5':  0.0,
    'I6':  0.0,

    # -- transverse acceptances --
    'Ax':  7.579,
    'Ay':  2.348, # With in-vac undulators - 4.5 mm

    # -- betatron coupling --
    'k_beta':      0.01, # coupling in transverse motion
    'k_dw':        0.00, # coupling through dispersion waves

    # -- physical constants --
    'Cgamma':    _mp.constants.rad_cgamma, # 0.0000884628,
    'Cq':        _mp.constants.Cq, # 3.83194e-13,
    'Ccoul':     _mp.constants.elementary_charge, # 1.602176487e-19,
    'r0':        _mp.constants.electron_radius, # 2.8179402894e-15,
    'me':        _mp.units.joule_2_eV(_mp.constants.electron_rest_energy)/1e9, # 0.00051099891,
    'cluz':      _mp.constants.light_speed, # 299792458,

}
