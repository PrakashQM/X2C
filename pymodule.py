import psi4
import re
import os
import inputparser
import math
import warnings
from driver import *
from wrappers import *
from molutil import *
import p4util
from p4xcpt import *


def run_x2c(name, **kwargs):
    r"""Function encoding sequence of PSI module and plugin calls so that
    x2c can be called via :py:func:`~driver.energy`. For post-scf plugins.

    >>> energy('x2c')

    """
    lowername = name.lower()
    kwargs = p4util.kwargs_lower(kwargs)

    # Your plugin's psi4 run sequence goes here
    psi4.plugin('x2c.so')
    #scf_helper(name, **kwargs)

# Integration with driver routines
procedures['energy']['x2c'] = run_x2c

