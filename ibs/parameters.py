class _IBSParameters:

    def __str__(self):
        str = '-- electron beam'
        str += '\n   beam energy: {0:.1e} eV'.format(self.beam_energy)
        str += '\n   beam current per bunch: {0:.5f} mA'.format(self.beam_current_per_bunch)
        str += '\n   beam current total: {0} mA'.format(self.beam_current_per_bunch * self.harmonic_number)
        str += '\n   number of electrons per bunch: {0:.3e}'.format(self.nr_electrons_per_bunch)
        str += '\n-- lattice'
        str += '\n   ring circumference: {0} m'.format(self.circumference)
        str += '\n   mean gas pressure: {0} nTorr'.format(self.mean_gas_pressure)
        str += '\n   horizontal acceptance: {0} mm.mrad'.format(self.ax)
        str += '\n   vertical acceptance: {0} mm.mrad'.format(self.ay)
        str += '\n-- rf cavities'
        str += '\n   harmonic number: {0}'.format(self.harmonic_number)
        str += '\n   total rf voltage: {0:.2e} MV'.format(self.Vrf)
        str += '\n   Rs0: {0:.2e}'.format(self.Rs0)
        str += '\n   Q0: {0:.2e}'.format(self.Q0)
        str += '\n   betac: {0}'.format(self.betac)
        str += '\n   detune0: {0}'.format(self.detune0)
        str += '\n   mharm: {0}'.format(self.mharm)
        str += '\n   Q0n: {0}'.format(self.Q0n)
        str += '\n   Rsn: {0}'.format(self.Rsn)
        str += '\n   detuneHC: {0}'.format(self.detuneHC)
        str += '\n-- ids'
        str += '\n   I1 contribution from ids: {0:+.6e} m'.format(self.ids_i1)
        str += '\n   I2 contribution from ids: {0:+.6e} 1/m'.format(self.ids_i2)
        str += '\n   I3 contribution from ids: {0:+.6e} 1/m²'.format(self.ids_i3)
        str += '\n   I4 contribution from ids: {0:+.6e} 1/m'.format(self.ids_i4)
        str += '\n   I5 contribution from ids: {0:+.6e} 1/m'.format(self.ids_i5)
        str += '\n   I6 contribution from ids: {0:+.6e} 1/m'.format(self.ids_i6)
        str += '\n-- optics'
        str += '\n   compaction factor: {0:+.3e}'.format(self.mcf)
        str += '\n   betatron coupling: {0} %'.format(self.k_beta * 100)
        str += '\n   dispersion coupling: {0} %'.format(self.k_dw * 100)
        str += '\n   I1 from lattice: {0:+.6e} m'.format(self.latt_i1)
        str += '\n   I2 from lattice: {0:+.6e} 1/m'.format(self.latt_i2)
        str += '\n   I3 from lattice: {0:+.6e} 1/m²'.format(self.latt_i3)
        str += '\n   I4 from lattice: {0:+.6e} 1/m'.format(self.latt_i4)
        str += '\n   I5 from lattice: {0:+.6e} 1/m'.format(self.latt_i5)
        str += '\n   I6 from lattice: {0:+.6e} 1/m'.format(self.latt_i6)
        str += '\n   energy loss per turn: {0} keV'.format(self.U0)
        str += '\n   overvoltage q: {0}'.format(self.q)
        str += '\n   rf energy acceptance: {0} %'.format(100*self.rf_acceptance)
        str += '\n-- natural equilibrium parameters'
        str += '\n   natural energy spread: {0} %'.format(100*self.natural_sigmae)
        #str += '\n   natural bunch length: {0} mm'.format(1000*self.natural_sigmal)





        return str

    #
    # self.beam_energy  = 3.0e9   #[eV]
    # self.beam_current_per_bunch = 100.0 / self.harmonic_number # [mA]
    # self.circumference = 518.396
    # self.mean_gas_pressure = 1.0 # [nTorr]
    # self.Vrf = 2.7e6 # total RF voltage [volts]
    # self.Rs0 = 19.8e6
    # self.Q0 = 27000
    # self.betac = 3.83
    # self.detune0 = 100e3
    # self.mharm = 3
    # self.Q0n = 2e8
    # self.Rsn = 176.8e8
    # self.detuneHC = 77.7e3
    # self.ap = 0.0001739520528
    #
    # self.ids_i1 = 0.0
    # self.ids_i2 = 0.200481368
    # self.ids_i3 = 0.036505767
    # self.ids_i4 = 5.89612e-7
    # self.ids_i5 = 1.2398e-6
    # self.ids_i6 = 0.0
    #
    # self.ax = 7.579 # [mrad.mm]
    # self.ay = 2.348 # [mrad.mm]  (with in-vac undulators 4.5 mm)
    # self.k_beta = 0.01
    # self.k_dw = 0.00
