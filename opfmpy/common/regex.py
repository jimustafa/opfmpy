"""Common regular expressions"""
from __future__ import absolute_import, division, print_function
import re


pattern_float = r'[-+]?[0-9]*\.?[0-9]+'

#: integers
regex_int = re.compile(r'[+-]?[0-9]+')

#: floats
regex_float = re.compile(pattern_float)

#: floats with exponent
regex_float_exp = re.compile(r'[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?')

#: floats with exponent
regex_float_exp_fortran = re.compile(r'[-+]?[0-9]*\.?[0-9]+([eEdD][-+]?[0-9]+)?')

#: boolean
regex_bool = re.compile(r'(1|0|t|f|on|off|true|false)', re.IGNORECASE)

#: Fortran logical
regex_bool_fortran = re.compile(r'(\.(t|f|true|false)\.)', re.IGNORECASE)

#: 3-tuple
regex_tuple3 = re.compile(r'[(]\s*(' + pattern_float + r')\s*,\s*(' + pattern_float + r')\s*,\s*(' + pattern_float + r')\s*[)]')
