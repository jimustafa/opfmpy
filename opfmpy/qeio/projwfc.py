from __future__ import absolute_import, division, print_function
import re
from xml.etree import cElementTree as ET

import numpy as np


pattern_fr = re.compile(
    r'state\s+#\s*(?P<state_idx>\d+):\s+'
    r'atom\s+(?P<atom_idx>\d+)\s+'
    r'[(](?P<atom_symbol>.+)[)],\s+'
    r'wfc\s+(?P<wfn_idx>\d+)\s+'
    r'[(]j=(?P<j>.+)\s+l=(?P<l>.+)\s+m_j=(?P<mj>.+)[)]'
    )

pattern_sr = re.compile(
    r'state\s+#\s*(?P<state_idx>\d+):\s+'
    r'atom\s+(?P<atom_idx>\d+)\s+'
    r'[(](?P<atom_symbol>.+)[)],\s+'
    r'wfc\s+(?P<wfn_idx>\d+)\s+'
    r'[(]l=(?P<l>.+)\s+m=(?P<m>.+)[)]'
    )

orbitals = {
    0: {
        1: 's',
    },
    1: {
        1: 'pz',
        2: 'px',
        3: 'py',
    },
    2: {
        1: 'dz2',
        2: 'dxz',
        3: 'dyz',
        4: 'dx2-y2',
        5: 'dxy',
    },
}


def read_atomic_proj(fname):
    xml_data = ET.parse(fname)

    nkpts = int(xml_data.find('HEADER/NUMBER_OF_K-POINTS').text)
    nbnds = int(xml_data.find('HEADER/NUMBER_OF_BANDS').text)
    nspns = int(xml_data.find('HEADER/NUMBER_OF_SPIN_COMPONENTS').text)
    nproj = int(xml_data.find('HEADER/NUMBER_OF_ATOMIC_WFC').text)

    atomic_proj = np.zeros((nkpts, nbnds, nproj), dtype=complex)

    for kpoint in xml_data.findall('PROJECTIONS/*'):
        ikpt = int(re.compile(r'K[-]POINT[.](?P<ikpt>\d+)').match(kpoint.tag).group('ikpt')) - 1
        for proj in kpoint.findall('./*'):
            iproj = int(re.compile(r'ATMWFC[.](?P<iproj>\d+)').match(proj.tag).group('iproj')) - 1
            atomic_proj[ikpt, :, iproj] = np.fromstring(proj.text.replace(',', ' '), sep='\n').view(complex)

    return atomic_proj


def read_atomic_states(fname):
    with open(fname, 'r') as f:
        contents = f.read()

    states_raw = pattern_sr.findall(contents)
    if not states_raw:
        states_raw = pattern_fr.findall(contents)
        if not states_raw:
            raise Exception
        else:
            states = [m.groupdict() for m in pattern_fr.finditer(contents)]
            map(lambda state: state.update({'j': float(state['j'])}), states)
            map(lambda state: state.update({'l': int(state['l'])}), states)
            map(lambda state: state.update({'mj': float(state['mj'])}), states)
            map(lambda state: state.update({'orbital': 'j=%f,l=%d,mj=%f' % (state['j'], state['l'], state['mj'])}), states)
    else:
        states = [m.groupdict() for m in pattern_sr.finditer(contents)]
        map(lambda state: state.update({'l': int(state['l'])}), states)
        map(lambda state: state.update({'m': int(state['m'])}), states)
        map(lambda state: state.update({'orbital': orbitals[state['l']][state['m']]}), states)

    map(lambda state: state.update({'state_idx': int(state['state_idx'])-1}), states)
    map(lambda state: state.update({'atom_idx': int(state['atom_idx'])-1}), states)
    map(lambda state: state.update({'atom_symbol': state['atom_symbol'].strip()}), states)
    map(lambda state: state.update({'wfn_idx': int(state['wfn_idx'])-1}), states)

    return states
