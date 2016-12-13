from __future__ import absolute_import, division, print_function
from argparse import ArgumentParser

try:
    import IPython
except ImportError:
    IPython = None
try:
    import Gnuplot
except ImportError:
    Gnuplot = None
try:
    import matplotlib
except ImportError:
    matplotlib = None


__all__ = [
    'parser_common',
    'parser_plot_common',
    ]

# common command-line parser
# ------------------------------------------------------------------------------
parser_common = ArgumentParser(add_help=False)
group = parser_common.add_mutually_exclusive_group()
group.add_argument('-v', '--verbose', action='count')
group.add_argument('-q', '--quiet', action='store_true')
parser_common.add_argument('--log-level', choices=['notset', 'debug', 'info', 'warning', 'error', 'critical'])
# ------------------------------------------------------------------------------


# common command-line parser for plotting scripts
# ------------------------------------------------------------------------------
parser_plot_common = ArgumentParser(add_help=False)
# parser_plot_common.add_argument('-o', '--output')
group = parser_plot_common.add_argument_group(title='plot arguments')
group.add_argument('-i', '--interact', action='store_true')
if Gnuplot is not None:
    group.add_argument('--interact-gnuplot', action='store_true')
if IPython is not None:
    group.add_argument('--interact-ipython', action='store_true')
if matplotlib is not None:
    group.add_argument('--matplotlib', action='store_true')
    group.add_argument('--interact-matplotlib', action='store_true')
    # group.add_argument('--mpl-style', type=mpl_style_arg)
# ------------------------------------------------------------------------------
