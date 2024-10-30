import time

import brightway2 as bw
import bw2data
import bw2io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import SALib
import seaborn as sns
from IPython.display import HTML, display
from scipy.stats import binned_statistic
from sympy import *
from tabulate import tabulate

# +
# Custom utils defined for inter-acv
import lca_algebraic as agb

# from lca_algebraic.params import FixedParamMode
# from lca_algebraic.stats import *
# from lca_algebraic.stats import _generate_random_params, _compute_stochastics
# -


# Larger space in notebook for large graphs
display(HTML("<style>.container { width:70% !important; }</style>"))

# Some options for pandas formatting
pd.options.display.max_rows = 200
pd.options.display.float_format = "{:,.3g}".format
