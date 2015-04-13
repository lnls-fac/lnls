import parameters as _parameters

class IBSParameters(_parameters._IBSParameters):

    def __init__(self):

        self.harmonic_number = 864

        self.beam_energy  = 3.0   #[GeV]
        self.beam_current_per_bunch = 100.0 / self.harmonic_number # [mA]
        self.circumference = 518.396
        self.mean_gas_pressure = 1.0 # [nTorr]

        self.Vrf = 2.7 # total RF voltage [MV]
        self.Rs0 = 19.8e6
        self.Q0 = 27000
        self.betac = 3.83
        self.detune0 = 100e3
        self.hcavities = False
        self.mharm = 3
        self.Q0n = 2e8
        self.Rsn = 176.8e8
        self.detuneHC = 77.7e3
        self.mcf = 0.0001739520528

        self.ids_i1 = 0.0
        self.ids_i2 = 0.200481368
        self.ids_i3 = 0.036505767
        self.ids_i4 = 5.89612e-7
        self.ids_i5 = 1.2398e-6
        self.ids_i6 = 0.0

        self.ax = 7.579 # [mrad.mm]
        self.ay = 2.348 # [mrad.mm]  (with in-vac undulators 4.5 mm)
        self.k_beta = 0.01
        self.k_dw = 0.00
