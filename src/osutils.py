import os
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')

def check_file_exists(filename):
    if not os.path.exists(filename):
        logging.info("{} file does not exist.".format(filename))
        exit(1)


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