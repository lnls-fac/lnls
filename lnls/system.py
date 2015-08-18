import os as _os


folder_root = _os.environ.get('FACROOT', '/home/fac_files')

folder_db = _os.path.join(folder_root, 'db')
folder_code = _os.path.join(folder_root, 'code')
folder_data = _os.path.join(folder_root, 'data')
folder_mml  = _os.path.join(folder_code, 'MatlabMiddleLayer', 'Release')
folder_mml_sirius     = _os.path.join(folder_mml, 'lnls', 'fac_scripts', 'sirius')
folder_mml_sirius_ids = _os.path.join(folder_mml_sirius, 'insertion_devices')

folder_excitation_curves = _os.path.join(folder_db, 'excitation_curves')
