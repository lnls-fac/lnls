import os as _os


folder_root = _os.environ.get('FACROOT', '/home/fac_files')

folder_db = _os.path.join(folder_root, 'siriusdb')
folder_fac_code = _os.path.join(folder_root, 'lnls-fac')
folder_sirius_code = _os.path.join(folder_root, 'lnls-sirius')
folder_data = _os.path.join(folder_root, 'data')
folder_mml  = _os.path.join(folder_fac_code, 'MatlabMiddleLayer', 'Release')
folder_mml_sirius     = _os.path.join(folder_mml, 'lnls', 'fac_scripts', 'sirius')
folder_mml_sirius_ids = _os.path.join(folder_mml_sirius, 'insertion_devices')

folder_excitation_curves = _os.path.join(folder_sirius_code, 'control-system-constants/magnets/excitation-data')
folder_pulse_curves = _os.path.join(folder_sirius_code, 'control-system-constants/magnets/pulse-curve-data')
