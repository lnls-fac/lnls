#!/usr/bin/env python-sirius

"""."""

import os as _os
import sys as _sys
from siriuspy import util as _util
import pickle as _pic
import numpy as _np


class Templates:
    """."""

    rawfield = """
    # ==========================
    # fma_rawfield.py input file
    # Date: TIMESTAMP
    # Accelerator Physics LNLS
    # ==========================

    # --- Summary ---
    #
    # this is the input file for fac-fma-rawfield.py script
    # this script reads a fieldmap from a 3D magnet model, stores it
    # for latter analysis and prints and plots basic information on the
    # field map. It is used to quickly inspect the fieldmap


    # --- Input parameters ---

    # each analysis has an identity label used for naming output files

      config_label             'CONFIG_NAME'


    # the next parameter specifies the type of magnet to be analysed.
    # each type may have its own particular algorithms to be applied

      magnet_type              'dipole'


    # the full name of the file that contains the field map

      fmap_filename            'FIELDMAP_FNAME'

    # Runge-Kutta algorithm used for the integration of the eqs. of motion needs to know
    # what to do when trajectory reaches the fieldmap bounds. It will either extrapolate the field
    # along the longitudinal (z) direction or consider it to have vanished. This is controlled with
    # the parameter below. Bear in mind that the calculation of extrapolation coefficients is very
    # time-consuming currently. As for the transverse directions (x and y), the RK algorithm will
    # generate exceptions.

      fmap_extrapolation_flag  False

    """

    trajectory = """
    # ==========================
    # fma_rawfield.py input file
    # Date: TIMESTAMP
    # Accelerator Physics LNLS
    # ==========================

    # --- Summary ---
    #
    # This is the input file for trajectory calculation based on a given
    # fieldmap which is performed with the script 'fac-fma-trajectory.py'
    # A controllable fixed-size Runge-Kutta algorithm is used to integrate
    # the equations of motion of a single electron in the presence of
    # the magnetic field as defined in the fieldmap.
    #
    # The implemented equations of motion are not approximated. Provided
    # a sufficiently fine RK step is chosen, this scripts may be used to
    # accurately obtain the trajectory of the electron with arbitrary energy
    #
    # Runge-Kutta algorithm used for the integration of the eqs. of motion needs to know
    # what to do when trajectory reaches the fieldmap bounds. It will either extrapolate the field
    # along the longitudinal (z) direction or consider it to have vanished.
    # As for the transverse directions (x and y), the RK algorithm will
    # generate exceptions.


    # --- Input parameters ---

    # each analysis has an identity label used for naming output files

      config_label             	'CONFIG_NAME'


    # beam energy

      beam_energy                     BEAM_ENERGY     # [GeV]


    # A trajectory can also be read from file. This is useful when the fieldmap of
    # 3D models with errors are being analysed. In this case we want to use as reference
    # trajectory a trajectory that was calculated from the 3D model without errors and
    # saved to file. If parameter 'traj_load_filename' is set to 'None' then a new
    # reference trajectory with be calculated with RK on the given fieldmap.

      traj_load_filename              None


    # If parameter 'traj_is_reference_traj' is set to True the algorithm will rescale the
    # fieldmap so that the total trajectory deflection will exactly match the nominal deflection

      traj_is_reference_traj          False
      model_nominal_angle             7.2                           # [deg] model nominal deflection angle of the magnet

    # There is the option to restrain the trajectory to the midplane (y = 0 mm) of the magnet

      traj_force_midplane_flag        True


    # There is the option to serach for a initial rx position that will result in a trajectory that
    # is centered in the good-field region of the magnet (around rx == 0)

      traj_center_sagitta_flag        False


    # The RK algorithm always integrates the trajectory from the center of the magnet (z = s = 0 mm)
    # The limits of the RK integration may be specified in various ways:
    # If only 'traj_rk_s_step' is given then the algorithm will integrate until
    # the z coordinate of the particle reaches the fieldmap bound.

      traj_init_rx                    RX_INIT                       # [mm] init_rx at z = s = 0 mm (center of magnet)
      traj_rk_s_step                  S_STEP                        # [mm]
      traj_rk_length                  None                          # [mm]
      traj_rk_nrpts                   None


    # whether to save trajectory to an ASCII file

      traj_save                       True
    """

    multipoles = """
    # ==========================
    # fma_rawfield.py input file
    # Date: TIMESTAMP
    # Accelerator Physics LNLS
    # ==========================

    # --- Summary ---
    #
    # this is the input file for the 'fac-fma-multipoles.py' script
    # this script calculates the multipoles around the reference trajectory.


    # --- Input parameters ---

    # each analysis has an identity label used for naming output files

      config_label                      'CONFIG_NAME'


    # the multipoles (m1,m2,...) to be calculated are defined by a list of position x exponents (n1,n2,...):
    # By = m1 * x^n1 + m2 * x^n2 + ...

      multipoles_normal_field_fitting_monomials      (0,1,2,3,4,5,6)                 # monomials to be included in the polynomial fit of multipoles
      multipoles_skew_field_fitting_monomials        ()

    # grid of perpendicular points around each point of the reference trajectory for the polynomial fit of By and Bx

      multipoles_perpendicular_grid     np.linspace(-12,12,65)          # grid of points on perpendicular line to ref trajectory [mm]

    # after multipole coeffs are calculated, their normalized strengths at perp. position r0 are calculated (as defined in tracy)

      multipoles_r0                     17.5                             # [mm] horizontal position at which polynomial fields are calculated relative to the principal multipole
      normalization_monomial            0
      normalization_is_skew             False

    # integrated residual field (converted to kick angle) calculated from fitted multipoles and
    # from integrated fieldmap are compared. The parameter below lists the monomials which are
    # supposed to define the main field. The rest makes up for the residual field

      normal_multipoles_main_monomials         (0,1,2)
      skew_multipoles_main_monomials           ()
    """

    model = """
    # ==========================
    # fma_rawfield.py input file
    # Date: 2018-07-20
    # Accelerator Physics LNLS
    # ==========================

    # --- Summary ---
    #
    # This script integrates fitted multipoles at each segment of the hard-edge model

    # --- Input parameters ---

    # each analysis has an identity label used for naming output files

      config_label                      'CONFIG_NAME'


    # list with lengths of model segments

       model_segmentation               (196, 192, 182, 10, 10, 13, 17, 20, 30, 50)
    """

    help = """
    NAME
           hallprobe.py - routines and libs to process fieldmap analysis

    SYNOPSIS
           hallprobe.py [CMD] [ARGS]

    DESCRIPTION
           The command options are

           help
           summary  [x0-9p1013mm|...] [current_label] [positive|negative]
           new-energy [x0-9p1013mm|...]
    """


# text templates
_rawfield = Templates.rawfield
_trajectory = Templates.trajectory
_multipoles = Templates.multipoles
_model = Templates.model
_help = Templates.help


class FMapAnalysis:
    """."""

    def __init__(self,
                 magnet,
                 curlabel,
                 path_analysis,
                 path_fmap,
                 beam_energy,
                 rx_init,
                 s_step):
        """."""
        self.magnet = magnet
        self.curlabel = curlabel
        self.path_analysis = path_analysis
        self.path_fmap = path_fmap
        self.beam_energy = beam_energy
        self.rx_init = rx_init
        self.s_step = s_step
        if not _os.path.exists(path_analysis):
            _os.makedirs(path_analysis)

    def files_create(self):
        """."""
        self._create_rawfield_in()
        self._create_trajectory_in()
        self._create_multipoles_in()
        self._create_model_in()

    def files_read_results_rawfield(self):
        """."""
        fname = self.path_analysis + '/rawfield.out'
        with open(fname, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if 'main_coil_current' in line:
                self.results_current = float(line.split()[1])

    def files_read_results_trajectory(self):
        """."""
        fname = self.path_analysis + '/trajectory.out'
        with open(fname, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if 'horizontal_deflection_angle' in line:
                self.results_angle = float(line.split()[1])
            elif 'rx position of reference' in line:
                self.results_rx_ref = float(line.split()[5])
            elif 'initial rx position' in line:
                self.results_rx_init = float(line.split()[5])

    def files_read_results_multipoles(self):
        """."""
        fname = self.path_analysis + '/multipoles.out'
        with open(fname, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if 'n=00' in line:
                self.results_intn0 = float(line.split()[2])
            elif 'n=01' in line:
                self.results_intn1 = float(line.split()[2])
            elif 'n=02' in line:
                self.results_intn2 = float(line.split()[2])

    def cmd_clean(self):
        """."""
        path = self.path_analysis
        cmd = 'cd ' + path + '; fac-fma-analysis.py clean'
        _os.system(cmd)

    def cmd_rawfield(self):
        """."""
        path = self.path_analysis
        cmd = 'cd ' + path + '; fac-fma-analysis.py rawfield'
        _os.system(cmd)

    def cmd_trajectory(self):
        """."""
        path = self.path_analysis
        cmd = 'cd ' + path + '; fac-fma-analysis.py trajectory'
        _os.system(cmd)

    def cmd_multipoles(self):
        """."""
        path = self.path_analysis
        cmd = 'cd ' + path + '; fac-fma-analysis.py multipoles'
        _os.system(cmd)

    def cmd_model(self):
        """."""
        path = self.path_analysis
        cmd = 'cd ' + path + '; fac-fma-analysis.py model'
        _os.system(cmd)

    def cmd_run(self):
        """."""
        path = self.path_analysis
        cmd = 'cd ' + path + '; fac-fma-analysis.py run'
        _os.system(cmd)

    def analysis_trajectory(self):
        """."""
        self.cmd_clean()
        self.files_create()
        self.cmd_rawfield()
        self.cmd_trajectory()
        self.files_read_results_rawfield()
        self.files_read_results_trajectory()

    def analysis_multipoles(self):
        """."""
        self.cmd_clean()
        self.files_create()
        self.cmd_rawfield()
        self.cmd_trajectory()
        self.cmd_multipoles()
        self.files_read_results_rawfield()
        self.files_read_results_trajectory()
        self.files_read_results_multipoles()

    def analysis_model(self):
        """."""
        self.cmd_clean()
        self.files_create()
        self.cmd_rawfield()
        self.cmd_trajectory()
        self.cmd_multipoles()
        self.cmd_model()
        self.files_read_results_rawfield()
        self.files_read_results_trajectory()
        self.files_read_results_multipoles()

    def analysis_concatenate_output_files(self):
        """."""
        path = self.path_analysis
        cmd = 'cd ' + path + '; cat rawfield.out >> analysis.txt'
        _os.system(cmd)
        cmd = 'cd ' + path + '; cat trajectory.out >> analysis.txt'
        _os.system(cmd)
        cmd = 'cd ' + path + '; cat multipoles.out >> analysis.txt'
        _os.system(cmd)
        cmd = 'cd ' + path + '; cat model.out >> analysis.txt'
        _os.system(cmd)

    # --- auxilliary method ---

    def _create_rawfield_in(self):
        text = self._process_text(_rawfield)
        fname = self.path_analysis + '/rawfield.in'
        with open(fname, 'w') as f:
            f.write(text)

    def _create_trajectory_in(self):
        text = self._process_text(_trajectory)
        fname = self.path_analysis + '/trajectory.in'
        with open(fname, 'w') as f:
            f.write(text)

    def _create_multipoles_in(self):
        text = self._process_text(_multipoles)
        fname = self.path_analysis + '/multipoles.in'
        with open(fname, 'w') as f:
            f.write(text)

    def _create_model_in(self):
        text = self._process_text(_model)
        fname = self.path_analysis + '/model.in'
        with open(fname, 'w') as f:
            f.write(text)

    def _process_text(self, text):
        t = text
        t = t.replace('TIMESTAMP', _util.get_timestamp())
        t = t.replace('MAGNET', self.magnet)
        t = t.replace('CONFIG_NAME', self.magnet+'-'+self.curlabel)
        t = t.replace('FIELDMAP_FNAME', self.path_fmap)
        t = t.replace('BEAM_ENERGY', str(self.beam_energy))
        t = t.replace('RX_INIT', str(self.rx_init))
        t = t.replace('S_STEP', str(self.s_step))
        return t


class FMapAnalysisProduction:
    """."""

    _path_analysis = (
        '/home/fac_files/lnls-ima/bo-dipoles/model-09'
        '/analysis/hallprobe/production/')

    _path_measurements = (
        '/home/fac_files/lnls-ima/bo-dipoles/model-09'
        '/measurement/magnetic/hallprobe/production/')

    def __init__(self,
                 magnet,
                 curlabel,
                 fmap_fname,
                 beam_energy,
                 s_step):
        """."""
        self.magnet = magnet
        self.curlabel = curlabel
        self.beam_energy = beam_energy
        self.path_analysis = FMapAnalysisProduction._path_analysis + \
            self.magnet + '/' + 'M1' + '/' + self.curlabel + '/'
        self.path_measurement = FMapAnalysisProduction._path_measurements + \
            self.magnet + '/' + 'M1' + '/'
        self.path_fmap = self.path_measurement + fmap_fname
        self.s_step = abs(s_step)
        self.analysis_neg = FMapAnalysis(self.magnet,
                                         self.curlabel,
                                         self.path_analysis + 'z-negative',
                                         self.path_fmap,
                                         self.beam_energy,
                                         -1.0*self.s_step)
        self.analysis_pos = FMapAnalysis(self.magnet,
                                         self.curlabel,
                                         self.path_analysis + 'z-positive',
                                         self.path_fmap,
                                         self.beam_energy,
                                         +1.0*self.s_step)

    def files_create(self):
        """."""
        self.analysis_pos.files_create()
        self.analysis_neg.files_create()

    def files_clean(self):
        """."""
        self.analysis_pos.files_clean()
        self.analysis_neg.files_clean()

    def run(self):
        """."""
        self.analysis_pos.files_run()
        self.analysis_neg.files_run()


def get_summary(dst, cur, pos):
    """."""
    data = (
        'Booster Dipoles Integrated Principal Multipoles\n'
        '================================================\n'
        '\n'
        'As calculated in {0:s}-half Runge-Kutta trajectory,\n'
        'defined by measured fieldmap with magnet excitated with current of'
        ' {1:s},\n'
        'corresponding to nominal particle energy of 3 GeV.\n'
        '\n').format(pos, cur)

    fmt = '{0:^10s} | {1:^11s} |  {2:^11s} | {3:^11s} | {4:^11s} |\n'
    data += fmt.format('Dipole', 'Angle [°]', 'Dint [T.m]', 'Gint [T]',
                       'Sint [T/m]')
    data += fmt.format('', '', '', '', '')

    fmt = ('{0:^10s} | {1:^+11.5f} |  {2:^+11.5f} | {3:^+11.5f} |'
           ' {4:^+11.5f} |\n')
    for i in range(4, 58):
        name = 'bd-{0:03d}'.format(i)
        fi = dst + '/' + name + '/M1/{0:s}/z-{1:s}/'.format(cur, pos)
        fname = './' + fi + cur + '.pkl'
        # print(fname)
        with open(fname, 'rb') as fi:
            config = _pic.load(fi)
        si = 1 if pos == 'positive' else -1
        multi = config.multipoles.normal_multipoles_integral
        theta = si*_np.arctan(config.traj.px[-1]/config.traj.pz[-1])*180/_np.pi
        data += fmt.format(name, theta, multi[0], multi[1], multi[2])

    print(data)
    with open(dst + '/README-{0:s}-Z{1:s}.md'.format(cur, pos), 'w') as fi:
        fi.write(data)


def get_new_energy(dst, cur):
    """."""
    theta = _np.zeros(54)
    for i in range(4, 58):
        for pos in ['positive', 'negative']:
            fi = dst + '/bd-{0:03d}/M1/' + cur + '/z-{1:s}/'.format(i, pos)
            with open(fi + cur + '.pkl', 'rb') as fi:
                config = _pic.load(fi)
            si = 1 if pos == 'positive' else -1
            theta[i-4] += \
                si*_np.arctan(config.traj.px[-1]/config.traj.pz[-1])*180/_np.pi

    print(theta)
    print('New Energy = {0:7.5f} GeV'.format(theta.mean()/-7.2 * 3))


def print_help():
    """."""
    print(_help)


def main():
    """."""
    n = len(_sys.argv)
    if n > 1 and 'summary' in _sys.argv[1]:
        dst = _sys.argv[2] if n > 2 else 'x0-9p1013mm'
        cur = _sys.argv[3] if n > 3 else '0991p63A'
        pos = _sys.argv[4] if n > 4 else 'positive'
        get_summary(dst, cur, pos)
    elif n > 1 and 'energy' in _sys.argv[1]:
        dst = _sys.argv[2] if n > 2 else 'x0-9p1013mm'
        cur = _sys.argv[3] if n > 3 else '0991p63A'
        get_new_energy(dst, cur)
    elif n > 1 and 'help' in _sys.argv[1]:
        print_help()
    else:
        print_help()


if __name__ == '__main__':
    main()