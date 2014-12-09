#! /usr/bin/env python
#
# RunOscillationAnalysis.py
#
# Runs the oscillation analysis, very similar to the LLROptimizer
# analysis, but different enough that it will be it's own script
#
# author: Tim Arlen - tca3@psu.edu
#
# date:   02-July-2014
#

import numpy as np
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from pisa.utils.log import logging, profile, physics, set_verbosity
from pisa.utils.jsons import from_json,to_json
from pisa.analysis.llr.LLHAnalysis import find_max_llh_bfgs, find_llh_no_osc, find_llh_mc_true
from pisa.analysis.stats.Maps import get_pseudo_data_fmap, get_seed
from pisa.analysis.TemplateMaker import TemplateMaker
from pisa.utils.params import get_values, select_hierarchy

parser = ArgumentParser(description='''Runs the LLR optimizer-based analysis varying a number of systematic parameters
defined in settings.json file and saves the likelihood values for all
combination of hierarchies.''',
                        formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-t','--template_settings',type=str,
                    metavar='JSONFILE', required = True,
                    help='''Settings related to the template generation and systematics.''')
parser.add_argument('-m','--minimizer_settings',type=str,
                    metavar='JSONFILE', required = True,
                    help='''Settings related to the optimizer used in the LLR analysis.''')
parser.add_argument('-n','--ntrials',type=int, default = 1,
                    help="Number of trials to run")
parser.add_argument('-s','--save-steps',action='store_true',default=False,
                    dest='save_steps',
                    help="Save all steps the optimizer takes.")
parser.add_argument('--save_pseudo_data',action='store_true',default=False,
                    help='Saves pseudo data in output file.')
parser.add_argument('--fit_no_osc',action='store_true',default=False,
                    help='Fits null hypothesis (no oscillations) for the given pseudo data set.')
parser.add_argument('-o','--outfile',type=str,default='llh_data.json',metavar='JSONFILE',
                    help="Output filename.")
parser.add_argument('-v', '--verbose', action='count', default=None,
                    help='set verbosity level')
args = parser.parse_args()

set_verbosity(args.verbose)

#Read in the settings
template_settings = from_json(args.template_settings)
minimizer_settings  = from_json(args.minimizer_settings)

#Workaround for old scipy versions
import scipy
if scipy.__version__ < '0.12.0':
    logging.warn('Detected scipy version %s < 0.12.0'%scipy.__version__)
    if 'maxiter' in minimizer_settings:
      logging.warn('Optimizer settings for \"maxiter\" will be ignored')
      minimizer_settings.pop('maxiter')

#Get the parameters
params = template_settings['params']

#store results from all the trials
trials = []

template_maker = TemplateMaker(get_values(params),**template_settings['binning'])

for itrial in xrange(1,args.ntrials+1):
    profile.info("start trial %d"%itrial)
    logging.info(">"*10 + "Running trial: %05d"%itrial + "<"*10)


    # //////////////////////////////////////////////////////////////////////
    # For each trial, generate two pseudo-data experiemnts (one for each
    # hierarchy), and for each find the best matching template in each of the
    # hierarchy hypothesis.
    # //////////////////////////////////////////////////////////////////////
    results = {}
    for data_tag, data_normal in [('data_NMH',True),('data_IMH',False)]:

        results[data_tag] = {}
        # 0) get a random seed and store with the data
        #results[data_tag]['seed'] = get_seed()
        #print "  seed used: ",results[data_tag]['seed']
        # 1) get a pseudo data fmap from fiducial model (best fit vals of params).
        fmap = get_pseudo_data_fmap(template_maker,
                                    get_values(select_hierarchy(params,
                                                                normal_hierarchy=data_normal)))
                                    #seed=results[data_tag]['seed'])

        if args.save_pseudo_data: results[data_tag]['pseudo_data'] = fmap

        if args.fit_no_osc:
            logging.info("Finding llh for no oscillations...")
            llh_no_osc = find_llh_no_osc(fmap,template_maker,params)
            results[data_tag]['fit_no_osc'] = {}
            results[data_tag]['fit_no_osc']['llh'] = llh_no_osc

        # Fit to true osc params:
        logging.info("Finding llh for true injected params...")
        llh_mc_true = find_llh_mc_true(fmap,template_maker,params,normal_hierarchy=data_normal)
        results[data_tag]['mc_true'] = {}
        results[data_tag]['mc_true']['llh'] = llh_mc_true

        # 2) find max llh (and best fit free params) from matching pseudo data
        #    to templates.
        for hypo_tag, hypo_normal in [('hypo_NMH',True),('hypo_IMH',False)]:

            physics.info("Finding best fit for %s under %s assumption"%(data_tag,hypo_tag))
            profile.info("start optimizer")
            llh_data = find_max_llh_bfgs(fmap,template_maker,params,
                                        minimizer_settings,args.save_steps,normal_hierarchy=hypo_normal)
            profile.info("stop optimizer")

            #Store the LLH data
            results[data_tag][hypo_tag] = llh_data


    #Store this trial
    trials += [results]
    profile.info("stop trial %d"%itrial)

#Assemble output dict
output = {'trials' : trials,
          'template_settings' : template_settings,
          'minimizer_settings' : minimizer_settings}
#And write to file
to_json(output,args.outfile)
