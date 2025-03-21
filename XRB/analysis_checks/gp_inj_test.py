#!/usr/bin/env python

import argparse
import csky as cy
import numpy as np
import pandas as pd
import sys, os
import matplotlib.pyplot as plt
import histlite as  hl



def load_data(dataset='combo'):
    repo = cy.selections.Repository()
    cspec = cy.selections.DNNCascadeDataSpecs.DNNC_10yr                             
    psspec = cy.selections.PSDataSpecs.ps_v4[3:]
    ana = cy.get_analysis(repo, 'version-004-p02', psspec, 'version-001-p01', cspec)
    #ana_cascade = cy.get_analysis(repo, 'version-001-p01', cspec)
    return ana

def run_gp_test_stacking(ana, n_trials = 1, ns=1697.85, seed=0, cpus =12):

    def get_gp_conf(ana, template_str, ): 
        template_repo = cy.selections.Repository(
        local_root='/data/ana/analyses/NuSources/2021_DNNCascade_analyses')
        if template_str == 'pi0':

            # get default gamma
            gamma = 2.7

            template = template_repo.get_template('Fermi-LAT_pi0_map')
            gp_conf = {
                'template': template,
                'flux': cy.hyp.PowerLawFlux(gamma),
                'randomize': ['ra'],
                'fitter_args': dict(gamma=gamma),
                'sigsub': True,
                'update_bg' : True,
                'fast_weight' : False,
                'extra_keep' : ['energy' , 'mjd'],
                }
            return gp_conf                                                

    def get_stacking_config(ana, src_gamma, fix_gamma, thresh, lag, weight):
        source_file  = '/home/ssclafani/XRB_Analysis/XRB/sources/lc_sources_reselected_IC86.hdf'
        sources = pd.read_hdf(source_file)
        ras = []
        decs = []
        lcs = []
        flux_weight = []
        for name in sources.name_disp:
            print(name)
            source = sources.loc[sources['name_disp'] == name]
            bins = np.array(source.lc_bins_10)[0]
            fluxes = np.array(source.lc_values_10)[0]
            index  = np.where((bins>=ana.mjd_min)& (bins<=ana.mjd_max))[0]
            if len(index) > 2:
                ras.append(source.ra_deg)
                decs.append(source.dec_deg)
                flux_weight.append(source.mean_rate)
                if index[0] != 0:
                        bins_in_data = np.append(np.append(ana.mjd_min-7.,bins[index]),ana.mjd_max+7.)
                        flux_in_data = np.append(fluxes[index[0]-1],fluxes[index])
                else:
                        bins_in_data = np.append(bins[index],ana.mjd_max+7.)
                        flux_in_data = fluxes[index] 
                lc = hl.Hist(bins_in_data, flux_in_data)
                lcs.append(lc)
            else:
                print('Not Including: ' + name)
        print(len(decs), len(ras), len(lcs))
        if weight == 'equal':
            src_weight = 1./(len(decs)) * np.ones_like(decs)
        elif weight == 'flux':
            src_weight = np.concatenate(flux_weight / np.sum(flux_weight))
        src = cy.utils.Sources(dec = np.concatenate(decs), ra = np.concatenate(ras), weight=src_weight, deg=True)

        print(f'inj lag: {lag}')
        print(f'inj thresh: {thresh}')                                                                 
        if fix_gamma:
            fitter_dict = {}
            for i in range(len(ras)):
                num_str = '0000'[0:4 - len(str(i))] + str(i)
                fitter_dict['lag_' + num_str] = 0.0
                fitter_dict['thresh_' + num_str] = 0.0
                fitter_dict['gamma'] = float(fix_gamma)
            fitter_dict['_fmin_method'] = 'minuit'
            print('Fixed gamma: {}'.format(fix_gamma))
            conf = dict(
                    src = src,
                    kind = 'binned',                
                    time='lc',
                    lcs=lcs,
                    range_lag=(-7.,7.),
                    sig='lc',
                    concat_evs=True,
                    fitter_args = fitter_dict,
                    extra_keep = ['energy'],
                    use_grl=True,
                    n_seeds_thresh=20,
                    sigsub = True,
                    n_seeds_lag = 20,
                    update_bg = True,
             )
        else:
            fitter_dict = {}
            for i in range(len(ras)):
                num_str = '0000'[0:4 - len(str(i))] + str(i)
                fitter_dict['lag_' + num_str] = 0.0         
                fitter_dict['thresh_' + num_str] = 0.0
            fitter_dict['_fmin_method'] = 'minuit'
            print('fitting gamma')
            conf = dict(
                    src = src,
                    kind = 'binned',                
                    time='lc',
                    lcs=lcs,
                    range_lag=(-7.,7.),
                    sig='lc',
                    #concat_evs=True,
                    extra_keep = ['energy'],
                    n_seeds_lag = 20,
                    update_bg = True,
                    extended = False,
                    sigsub = True,
                    fitter_args = fitter_dict,
                    )                                                                                                                   
        inj_conf=dict(lcs=lcs, threshs=np.zeros_like(lcs), lags=np.zeros_like(lcs), flux = cy.hyp.PowerLawFlux(gamma=src_gamma))

        return conf, inj_conf


    gp_conf = get_gp_conf(ana, 'pi0')
    gp_tr = cy.get_trial_runner(ana=ana, conf=gp_conf)


    stacking_conf, stacking_inj_conf = get_stacking_config(
                                    ana, 
                                    src_gamma = 2.0, 
                                    fix_gamma = None, 
                                    thresh = 0, 
                                    lag = 0, 
                                    weight = 'flux')
    stacking_tr  = cy.get_trial_runner(ana=ana, conf=stacking_conf, mp_cpus = cpus)


    stacking_tr.sig_injs[0]  = gp_tr.sig_injs[0]
    stacking_tr.sig_injs[1]  = gp_tr.sig_injs[1]    

    trials = stacking_tr.get_many_fits(n_trials=n_trials, n_sig=ns, seed=seed, mp_cpus = cpus)
    
    dir = '/data/user/ssclafani/data/analyses/XRB_baseline_v0.3/gp_tests'
    np.save(f'{dir}/stacking_trials_{n_trials:08d}_ns_{ns}_seed_{seed:08d}', trials.as_array)
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--n_sig', type = float, default = 0)
    parser.add_argument('--n_trials', type = int, default = 1)
    parser.add_argument('--seed', type = int, default=0)
    parser.add_argument('--cpus', type = int, default=1)
    args = parser.parse_args()

    ana = load_data(dataset='combo')
    trials = run_gp_test_stacking(ana, args.n_trials, args.n_sig, args.seed, args.cpus)


