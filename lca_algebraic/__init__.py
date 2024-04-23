#
# This file defines several utility functions above brightway2 to be used by notebooks
#

from .activity import *
from .base_utils import *
from .database import *
from .interpolation import *
from .io import *
from .lca import *
from .log import *
from .methods import *
from .params import *
from .settings import Settings
from .stats import *
from .units import *

# Backport the test of pypardiso from bw2calc to emit a warning
try:
    from pypardiso import factorized, spsolve
except ImportError:
    from scipy.sparse.linalg import factorized, spsolve
    from scipy.sparse.linalg._dsolve import linsolve

    if not linsolve.useUmfpack:
        logger.warn(
            """
        Did not findPypardisio or Umfpack. Matrix computation may be very slow.

        If you are on an Intel architecture, please install pypardiso as explained in the docs :
        https://docs.brightway.dev/en/latest/content/installation/index.html

        > pip install pypardiso
        or
        > conda install pypardiso
        """
        )


# Global print options
np.set_printoptions(threshold=30)
pd.options.display.float_format = "{:,g}".format
