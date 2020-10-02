import os
import sys
import logging

logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
from typing import Union


def check_file_exists(filename: str = None):
    """Checks if file exists in path"""
    if filename is None:
        return

    if not os.path.isfile(filename):
        raise ValueError("{} file does not exist.".format(filename))
    return filename


def ensure_folder_exists(filename):
    if os.path.isdir(filename):
        pass
    elif os.path.isfile(filename):
        ensure_folder_exists(os.path.dirname(filename))
    else:
        head, tail = os.path.split(filename)
        if head and not os.path.isdir(head):
            ensure_folder_exists(head)
        if tail:
            os.mkdir(filename)


def setup_output_directory(directory: Union[str, None]):
    """
    Creates output directory from given
    name.
    :param Union[str, None] directory: Name of
        the output dir. If `None`, 'out_VETA'
        will be created.
    """
    out_dir = directory if directory is not None else os.path.join(os.getcwd(), "out_VETA")
    assert not os.path.isdir(out_dir), "Output directory {} exists. " \
                                       "Remove it or set a new one.".format(out_dir)
    return out_dir


def print_clinvar_levels():
    """
    Prints the clinvar levels available
    to filter the input variants for all
    the analysis when the reference dataset
    is the clinvar database.
    """
    print("clinvar ---- Whole clinvar with Pathogenic and Benign assignments" + "\n" +
          "clinvar_l ---- Same as clinvar but with likely assignments" + "\n" +
          "1s --- clinvar with 1 star or more" + "\n" +
          "2s --- 2 stars or more" + "\n" +
          "3s --- 3 stars or more" + "\n" +
          "4s --- 4 stars" + "\n" +
          "1s_l ---- 1 star or more with likely assignments" + "\n" +
          "2s_l ---- 2 stars or more with likely assignments" + "\n" +
          "3s_l ---- 3 stars or more with likely assignments" + "\n" +
          "4s_l ---- 4 stars with likely assignments")
    exit(1)
