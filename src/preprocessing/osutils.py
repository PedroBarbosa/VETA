import logging
import os
import tempfile
from itertools import zip_longest
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
    # assert not os.path.isdir(out_dir), "Output directory {} exists. " \
    #                                    "Remove it or set a new one.".format(out_dir)
    return out_dir


def split_file_in_chunks(file: str,
                         header: list,
                         n_variants: int,
                         chunk_size: int = 2000):
    """
    Receives
    :param str file: Input VCF to split without the header
    :param list header: VCF header as an iterable
    :param int n_variants: Number of variants in the input `file`
    :param int chunk_size: Number of variants per each new smaller file. Default: `2000`

    :return list: List with tmp file names created
    """
    # Output Files
    outlist = []

    # Ceil division
    n_chunks = -(-n_variants // chunk_size)
    logging.info("More than 5000 variants found ({}). "
                 "Spliting file in {} chunks of {} variants".format(n_variants, n_chunks, chunk_size))

    # Decode header
    _header = [x.decode('utf-8') for x in header]

    # Process variants
    with open(file) as f:
        for i, g in enumerate(grouper(n_chunks, f, fillvalue=''), 1):
            with tempfile.NamedTemporaryFile('w', delete=False) as tmp:
                tmp.writelines(_header)
                tmp.writelines(g)
            outlist.append(tmp.name)

    return outlist


def grouper(n, file_iterable, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    args = [iter(file_iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def print_clinvar_levels():
    """
    Prints the clinvar levels available
    to filter the input variants for all
    the analysis when the reference dataset
    is the clinvar database.
    """
    print("0s ---- Whole clinvar with Pathogenic and Benign assignments" + "\n" +
          "1s --- clinvar with 1 star or more" + "\n" +
          "2s --- 2 stars or more" + "\n" +
          "3s --- 3 stars or more" + "\n" +
          "4s --- 4 stars" + "\n" +
          "0s_l ---- 0 star or more with likely assignments" + "\n" +
          "1s_l ---- 1 star or more with likely assignments" + "\n" +
          "2s_l ---- 2 stars or more with likely assignments" + "\n" +
          "3s_l ---- 3 stars or more with likely assignments" + "\n" +
          "4s_l ---- 4 stars with likely assignments")
    exit(1)
