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
    cspec = cy.selections.DNNCascadeDataSpecs.DNNC_12yr                             
    psspec = cy.selections.PSDataSpecs.ps_v4_15yr[3:]
    ana = cy.get_analysis(repo, 'version-004-p03', psspec, 'version-001-p02', cspec)
    return ana

def run_gp_test_stacking(ana, n_trials = 1, ns=0, seed=0, cpus =12, gp_inject=False):

    def get_stacking_config(ana, src_gamma, fix_gamma, thresh, lag, weight, inject_gp):
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
            if inject_gp:
                inj_conf=dict(lcs=lcs,
                            threshs=np.zeros_like(lcs),
                            lags=np.zeros_like(lcs),
                            flux = cy.hyp.PowerLawFlux(gamma=src_gamma),
                            gp_filenames = ['/data/user/ssclafani/GP_injected_trials/trial_tracks_IC86.npy',
                                  '/data/user/ssclafani/GP_injected_trials/trial_cascades_IC86.npy'])
            else:
                inj_conf=dict(lcs=lcs, threshs=np.zeros_like(lcs), lags=np.zeros_like(lcs), flux = cy.hyp.PowerLawFlux(gamma=src_gamma))
        return conf, inj_conf

    if gp_inject:
        print('Injecting gp into bkg')
        inject_gp = True
    else:
        print('Not Injecting gp into bkg')
        inject_gp = False 

    stacking_conf, stacking_inj_conf = get_stacking_config(
                                    ana, 
                                    src_gamma = 2.0, 
                                    fix_gamma = None, 
                                    thresh = 0, 
                                    lag = 0, 
                                    weight = 'flux',
                                    inject_gp = inject_gp)
    stacking_tr  = cy.get_trial_runner(ana=ana, conf=stacking_conf, inj_conf = stacking_inj_conf, mp_cpus = cpus) 
    #stacking_tr.sig_injs[0]  = gp_tr.sig_injs[0]
    #stacking_tr.sig_injs[1]  = gp_tr.sig_injs[1]    

    trials = stacking_tr.get_many_fits(n_trials=n_trials, n_sig=ns, seed=seed, mp_cpus = cpus)
    print(np.median(trials.ts)) 
    if inject_gp:
        dir = cy.utils.ensure_dir('/data/user/ssclafani/data/analyses/XRB_baseline_v0.5/gp_tests/gp_inj/')
        np.save(f'{dir}/stacking_trials_{n_trials:08d}_ns_{ns}_seed_{seed:08d}', trials.as_array)
    else:
        dir = cy.utils.ensure_dir('/data/user/ssclafani/data/analyses/XRB_baseline_v0.5/gp_tests/no_g/')
        np.save(f'{dir}/stacking_trials_{n_trials:08d}_ns_{ns}_seed_{seed:08d}', trials.as_array)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--n_sig', type = float, default = 2012.25)
    parser.add_argument('--n_trials', type = int, default = 1)
    parser.add_argument('--seed', type = int, default=0)
    parser.add_argument('--cpus', type = int, default=1)
    parser.add_argument('--gp_inject', action = 'store_true', default=False)
    args = parser.parse_args()

    ana = load_data(dataset='combo')
    trials = run_gp_test_stacking(ana, args.n_trials, args.n_sig, args.seed, args.cpus, args.gp_inject)


