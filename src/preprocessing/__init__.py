import os.path
from .utils import *
from .clinvar import *
from .vcf import *
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')


