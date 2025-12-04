"""Util module."""

import gzip as _gzip
import os as _os
import pickle as _pickle


def files_get_matches(folder=None, recursive=True, strs_in=None, strs_out=None):
    if folder is None:
        folder = _os.getcwd()
    if strs_in is None:
        strs_in = ('.dat', 'BOB_')
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
            strs_in_flag = [1 if token in path else 0 for token in strs_in]
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
        files = sorted(files, key=lambda v: v[-17:])

    return files


def save_pickle(filename, **kwargs):
    """Save variables in gzip compressed pickle format as a dictionary.
    INPUTS:
      filename : path to and name (with or without extension) of the file to
      save
      kwargs   : variables to save:

    Examples:
    >>> a = dict({'casa':[1,2], 'bla':3.4})
    >>> b = ['fla',3.42,True]
    >>>save_pickle('teste',a=a,b=b)
    >>>save_pickle('teste2',var1=a,var2=b)
    """
    if not filename.endswith('.pickle'):
        filename += '.pickle'

    with _gzip.open(filename, 'wb') as fi:
        _pickle.dump(kwargs, fi, _pickle.HIGHEST_PROTOCOL)


def load_pickle(filename):
    """Load gzip compressed files in pickle format

    Input:
      filename : path to and name (with or without extension) of the file to
    load

    Output:
        dictionary with the variables of the file.

    Examples:
        >>> a = dict({'casa':[1,2], 'bla':3.4})
        >>> b = ['fla',3.42,True]
        >>>save_pickle('teste',a=a,b=b)
        >>>save_pickle('teste2',var1=a,var2=b)
        >>>vars = load_pickle('teste2.pickle')
        >>> vars['var1']
        {'casa':[1,2], 'bla':3.4}
        >>>vars['var2']
        ['fla',3.42,True]
        >>> load_pickle('teste')
        {'a':{'casa':[1,2], 'bla':3.4},'b':['fla',3.42,True]}
    """
    if not filename.endswith('.pickle'):
        filename += '.pickle'

    with _gzip.open(filename, 'rb') as f1:
        data = _pickle.load(f1)
    return data
