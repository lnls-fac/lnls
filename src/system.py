import os as _os

folder_root = '/home/fac_files'

folder_data = _os.path.join(folder_root,'data')
folder_code = _os.path.join(folder_root,'code')
folder_mml  = _os.path.join(folder_code,'MatlabMiddleLayer','Release')
folder_mml_sirius     = _os.path.join(folder_mml, 'lnls', 'fac_scripts', 'sirius')
folder_mml_sirius_ids = _os.path.join(folder_mml_sirius, 'insertion_devices')
