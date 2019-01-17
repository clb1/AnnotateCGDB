import sys
from collections import namedtuple
import lxml.etree
import requests
from math import factorial
from scipy.stats import rankdata
import numpy as np
from operator import itemgetter
import tempfile
import commonCGDB
import time
import random

from primer3.bindings import designPrimers, setP3Globals

import Primer3RegionSpecification, Signature, IsoformSignatureGenerator

# frac_duplexed is fraction of primer that are calculated to be duplexed to its exact complementary template at the reaction temperature
# In addition to thermodynamic considerations, want to take into account primer 3' "goodness". A 3' is "better" when it has
# a balane of A+T vs C+G and when it has G or C as its 3' nucleotide
PrimerSingleton = namedtuple("PrimerSingleton", ["template_5p_pos", "seq", "len", "Tm", "frac_duplexed", "thermo_penalty", "aux_3p_penalty"] )

