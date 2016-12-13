from __future__ import absolute_import, division, print_function
import os
import re
from argparse import ArgumentTypeError

from ..regex import regex_int


__all__ = [
    'extant_file',
    'extant_path',
    'positive_int_arg',
    'nonnegative_int_arg',
    'band_pair_arg',
    'band_range_arg',
    'spin_arg',
    ]


def extant_file(f):
    if not os.path.isfile(f):
        raise ArgumentTypeError('the file "%s" does not exist' % f)
    return f


def extant_path(s):
    if not os.path.exists(s):
        raise ArgumentTypeError('the path "%s" does not exist' % s)
    return s


def positive_int_arg(s):
    try:
        i = int(s)
        if i < 1:
            raise ArgumentTypeError('Integer must be positive')
        return i
    except ArgumentTypeError:
        raise
    except:
        raise ArgumentTypeError('Must be a positive integer')


def nonnegative_int_arg(s):
    try:
        i = int(s)
        if i < 0:
            raise ArgumentTypeError('Integer must be nonnegative')
        return i
    except ArgumentTypeError:
        raise
    except:
        raise ArgumentTypeError('Must be a nonnegative integer')


def band_pair_arg(s):
    if not re.compile('{0},{0}'.format(regex_int.pattern)):
        raise ArgumentTypeError('Band pairs must be in the form ibnd1,ibnd2')
    try:
        ibnd1, ibnd2 = map(int, s.split(','))
        if ibnd1 < 0 or ibnd2 < 0:
            raise ArgumentTypeError('Band indices must be nonnegative')
        return ibnd1, ibnd2
    except ArgumentTypeError:
        raise
    except:
        raise ArgumentTypeError('Enter band pairs in the form ibnd1,ibnd2')


def band_range_arg(s):
    if not re.compile('{0}-{0}'.format(regex_int.pattern)):
        raise ArgumentTypeError('Band range must be in the form ibnd1-ibnd2')
    try:
        ibnd1, ibnd2 = map(int, s.split('-'))
        if ibnd1 < 1 or ibnd2 < 1:
            raise ArgumentTypeError('Band indices must be >= 1')
        return ibnd1, ibnd2
    except ArgumentTypeError:
        raise
    except:
        raise ArgumentTypeError('Enter band pairs in the form ibnd1-ibnd2')


def spin_arg(s):
    try:
        ispn = int(s)
        if ispn < 0 or ispn > 1:
            raise ArgumentTypeError('Spin index should be 0 or 1')
        return ispn
    except ArgumentTypeError:
        raise
    except:
        raise ArgumentTypeError('Enter spin index (0 or 1)')
