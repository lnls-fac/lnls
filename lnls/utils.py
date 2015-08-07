import os as _os
import re as _re

def get_data_paths_list(top_folder=None, recursive=True, regexp=None):
    """List all files under 'top_folder' whose names obey regexp"""
    if top_folder is None: top_folder = os.getcwd()
    all_paths = []
    try:
        local_paths = [_os.path.join(top_folder, file) for file in _os.listdir(top_folder)]
    except PermissionError:
        local_paths = []
    regexp = _re.compile(regexp)
    for path in local_paths:
        if _os.path.isdir(path) and not _os.path.islink(path):
            if recursive: all_paths.extend(get_data_paths_list(path, True, regexp))
        else:
            #print(path)
            if regexp.match(_os.path.basename(path)): all_paths.append(path)
    return all_paths
