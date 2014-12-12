#
# Maps.py
#
# Utilities for dealing with event rate maps in analysis
#
# author: Tim Arlen   <tca3@psu.edu>
#

import os
import numpy as np
from pisa.analysis.stats.LLHStatistics import get_random_map

def get_pseudo_data_fmap(template_maker,fiducial_params,seed=None,chan='both',
                         no_osc=False):
    '''
    Creates a true template from fiducial_params, then uses Poisson statistics
    to vary the expected counts per bin to create a pseudo data set.
    If seed is provided, the random state is seeded with seed before the map is
    created.
    '''
    
    true_template = template_maker.get_template(fiducial_params) if not no_osc else template_maker.get_template_no_osc(fiducial_params)
    true_fmap = flatten_map(true_template,chan=chan)
    fmap = get_random_map(true_fmap, seed=seed)

    return fmap

def flatten_map(template,chan='both'):
    '''
    Takes a final level true (expected) template of trck/cscd, and returns a
    single flattened map of trck appended to cscd, with all zero bins
    removed.

    IMPORTANT: returns a SINGLE flattened map of:
      if chan == 'both': both trck/casc combined
      if chan == 'trck': trck channel only
      if chan == 'cscd': cscd channel only

    '''
    cscd = template['cscd']['map'].flatten()
    trck = template['trck']['map'].flatten()

    if chan == 'both':
        fmap = np.append(cscd,trck)
        fmap = np.array(fmap)[np.nonzero(fmap)]
    elif chan == 'trck':
        fmap = np.array(trck)[np.nonzero(trck)]
    elif chan == 'cscd':
        fmap = np.array(cscd)[np.nonzero(cscd)]

    return fmap

def get_seed():
    '''
    Returns a random seed from /dev/urandom that can be used to seed the random
    state, e.g. for the poisson random variates.
    '''
    return int(os.urandom(8).encode('hex'),16)

