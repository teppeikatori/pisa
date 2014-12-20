#! /usr/bin/env python
#
# RunOscillationAnalysis.py
#
# Runs the oscillation analysis, very similar to the LLROptimizer
# NMH analysis, but different enough that it will be it's own script
#
# author: Tim Arlen - tca3@psu.edu
#
# date:   11-December-2014
#

import numpy as np
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys

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
parser.add_argument('-c','--chan',type=str,default='both',
                    help='which channel to use in the fit. For now, can be one of [trck,cscd,both]')
parser.add_argument('--save_pseudo_data',action='store_true',default=False,
                    help='Saves pseudo data in output file.')
parser.add_argument('--fit_null',action='store_true',default=False,
                    help='Fits null hypothesis (no oscillations) for the given pseudo data set, to allow computation of significance of numu appearance.')
parser.add_argument('--fit_mc_true',action='store_true',default=False,
                    help='Fits mc true values for given pseudo data set to find chi2-like distribution.')
parser.add_argument('-o','--outfile',type=str,default='llh_data.json',metavar='JSONFILE',
                    help="Output filename.")
parser.add_argument('-v', '--verbose', action='count', default=None,
                    help='set verbosity level')
args = parser.parse_args()

cmd_str = ""
for arg in sys.argv: cmd_str+=arg+" "

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
    # For each trial, generate pseudo-data experiments. If we are fitting
    # for the NULL hypothesis (no oscillations) then generate data sets
    # for 'no_osc', 'nmh', and 'imh'. Otherwise, will only generate pseudo
    # data for 'nmh' and 'imh', since oscillations are assumed.
    #
    # Pseudo data sets are then fit to templates assuming these alternative
    # hypotheses.
    # //////////////////////////////////////////////////////////////////////

    data_types = [('data_NMH',True),('data_IMH',False)]
    if args.fit_null: data_types.append(('data_NULL',None))

    results = {}
    for data_tag, data_normal in data_types:
        results[data_tag] = {}

        if data_normal is None:
            logging.info("Defining pseudo data from template assuming NO OSCILLATIONS...")
            fmap = get_pseudo_data_fmap(template_maker,get_values(params),chan=args.chan,no_osc=True)
        else:
            logging.info("Defining pseudo data from template assuming oscillations...")
            params_hierarchy = select_hierarchy(params,normal_hierarchy=data_normal)
            fmap = get_pseudo_data_fmap(template_maker,get_values(params_hierarchy),chan=args.chan)

        if args.save_pseudo_data: results[data_tag]['pseudo_data'] = fmap

        if args.fit_mc_true:
            # Find LLH of MC True-NEED TO IMPLEMENT NUISANCE PARAM FIT HERE
            logging.info("Finding llh for true injected params...")
            llh_mc_true = find_llh_mc_true(fmap,template_maker,params,normal_hierarchy=data_normal)
            results[data_tag]['mc_true'] = llh_mc_true

        hypo_types = [('hypo_NMH',True),('hypo_IMH',False)]
        if args.fit_null: hypo_types.append(('hypo_NULL',None))
        for hypo_tag, hypo_normal in hypo_types:

            physics.info("Finding best fit for %s under %s assumption"%(data_tag,hypo_tag))
            profile.info("start optimizer")
            # HOW do we handle 'NULL' case here? Well, I think we can
            # give normal_hierarchy=None, and if that happens,
            # find_max_llh_bfgs will use
            # template_maker.get_template_no_osc() instead??
            if hypo_normal == None:
                no_osc = True
                hypo_normal = True
            else:
                no_osc = False
            llh_data = find_max_llh_bfgs(fmap,template_maker,params,
                                         minimizer_settings,args.save_steps,
                                         normal_hierarchy=hypo_normal,no_osc=no_osc,
                                         chan=args.chan)
            profile.info("stop optimizer")

            #Store the LLH data
            results[data_tag][hypo_tag] = llh_data


    #Store this trial
    trials += [results]
    profile.info("stop trial %d"%itrial)

#Assemble output dict
output = {'trials' : trials,
          'template_settings' : template_settings,
          'minimizer_settings' : minimizer_settings,
          'command': cmd_str}
#And write to file
to_json(output,args.outfile)
