import os as _os
import re as _re

def files_get_matches(folder=None, recursive=True, strs_in=None, strs_out=None):
    if folder is None:
        folder = os.getcwd()
    if strs_in is None:
        strs_in = ('.dat','BOB_')
    elif isinstance(strs_in, str):
        strs_in = (strs_in,)
    elif isinstance(strs_out, str):
        strs_out = (strs_out,)
    files = []
    local_files = _os.listdir(folder)
    for fname in local_files:
        path = _os.path.join(folder, fname)
        if _os.path.isdir(path):
            if recursive:
                files.extend(files_get_matches(path, recursive, strs_in, strs_out))
        else:
            strs_in_flag  = [1 if token in path else 0 for token in strs_in]
            if strs_out:
                strs_out_flag = [1 if token not in path else 0 for token in strs_out]
            else:
                strs_out_flag = [1,]
            if all(strs_in_flag) and (not strs_out or all(strs_out_flag)):
                files.append(path)

    # --- sort list of filenames ---
    fname = _os.path.basename(files[0])
    if fname[0].isdigit():
        files = sorted(files)
    else:
        files = sorted(files, key=lambda v:v[-17:])

    return files
