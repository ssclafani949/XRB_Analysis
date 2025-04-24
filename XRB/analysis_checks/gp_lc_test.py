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

def run_gp_lc_test(ana, name, n_trials = 1, ns=0, seed=0, cpus =12, gp_inject=False):

    def get_lc_config(ana, name, src_gamma, fix_gamma, cutoff_GeV, lag, thresh, gp_inject = False):
        source_file  = '/home/ssclafani/XRB_Analysis/XRB/sources/lc_sources_reselected_IC86.hdf'
        sources = pd.read_hdf(source_file)
        source = sources.loc[sources['name_disp'] == name]
        print(name)
        print(source.dec_deg)
        src = cy.utils.Sources(dec = source.dec_deg, ra = source.ra_deg, deg=True)

        print('Source RA: {} DEC: {}'.format(np.degrees(src.ra[0]), np.degrees(src.dec[0])))

        bins = np.array(source.lc_bins_10)[0]
        fluxes = np.array(source.lc_values_10)[0]
        index        = np.where((bins>=ana.mjd_min)& (bins<=ana.mjd_max))[0]
        print (bins[index])
        print (ana.mjd_min)
        print (ana.mjd_max)
        if len(index) > 2:                                                                    
            if index[0] != 0:
                bins_in_data = np.append(np.append(ana.mjd_min-7.,bins[index]),ana.mjd_max+7.)
                flux_in_data = np.append(fluxes[index[0]-1],fluxes[index])
            else:
                bins_in_data = np.append(bins[index],ana.mjd_max+7.)
                flux_in_data = fluxes[index] 
        else:
            print('No bins in data')
            print(index)
        print(f'lag: {lag}')
        print(f'thresh: {thresh}')
        #dir = cy.utils.ensure_dir ('{}/lc/{}'.format (base_dir, name))
        lc = hl.Hist(bins_in_data, flux_in_data)
        if fix_gamma:                                                               
            print('Fixed gamma: {}'.format(fix_gamma))
            conf = dict(
                src = src,
                kind = 'binned',                
                time='lc',
                lcs=lc,
                range_lag=(-7.,7.),
                #range_thresh=(0.,.5),
                sig='lc',
                concat_evs=True,
                extra_keep = ['energy'],
                #concat_evs=False,
                use_grl=True,
                n_seeds_thresh=20,
                n_seeds_lag = 20,
                update_bg = True,
                sigsub = True,
                fitter_args = dict(_fmin_method='minuit', gamma = float(fix_gamma)),
            )
        else:                                                                       
            print('fitting gamma')
            conf = dict(
                src = src,
                kind = 'binned',                
                time='lc',
                lcs=lc,
                range_lag=(-7.,7.),
                sig='lc',
                concat_evs=True,
                extra_keep = ['energy'],
                n_seeds_thresh=20,
                n_seeds_lag = 20,
                update_bg = True,
                fitter_args = dict(_fmin_method='minuit'),
                sigsub = True,
            )
        if gp_inject:
            inj_conf=dict(lcs=lc, 
                        threshs=thresh, 
                        lags=lag, 
                        flux = cy.hyp.PowerLawFlux(gamma=src_gamma),
                        gp_filenames = ['/data/user/ssclafani/GP_injected_trials/trial_tracks_IC86.npy',
                                        '/data/user/ssclafani/GP_injected_trials/trial_cascades_IC86.npy'])
        else:
            inj_conf=dict(lcs=lc, threshs=thresh, lags=lag, flux = cy.hyp.PowerLawFlux(gamma=src_gamma))

        del sources, lc, 
        return conf, inj_conf                                                                       

    if gp_inject:
        print('Injecting gp into bkg')
        inject_gp = True
    else:
        print('Not Injecting gp into bkg')
        inject_gp = False 

    lc_conf, lc_inj_conf = get_lc_config(
                                    ana, 
                                    name,
                                    src_gamma = 2.0, 
                                    fix_gamma = None, 
                                    thresh = 0, 
                                    lag = 0, 
                                    cutoff_GeV = np.inf,
                                    gp_inject = inject_gp)

    tr  = cy.get_trial_runner(ana=ana, conf=lc_conf, inj_conf = lc_inj_conf, mp_cpus = cpus) 

    trials = tr.get_many_fits(n_trials=n_trials, n_sig=ns, seed=seed, mp_cpus = cpus)
    print(np.median(trials.ts)) 
    if inject_gp:
        dir = cy.utils.ensure_dir('/data/user/ssclafani/data/analyses/XRB_baseline_v0.5/gp_tests/lc/{name}/gp/')
        np.save(f'{dir}/stacking_trials_{n_trials:08d}_ns_{ns}_seed_{seed:08d}', trials.as_array)
    else:
        dir = cy.utils.ensure_dir('/data/user/ssclafani/data/analyses/XRB_baseline_v0.5/gp_tests/lc/{name}/no_gp/')
        np.save(f'{dir}/stacking_trials_{n_trials:08d}_ns_{ns}_seed_{seed:08d}', trials.as_array)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--name', type = str, default = 'Cyg_X_dash_1')
    parser.add_argument('--n_sig', type = float, default = 0)
    parser.add_argument('--n_trials', type = int, default = 1)
    parser.add_argument('--seed', type = int, default=0)
    parser.add_argument('--cpus', type = int, default=1)
    parser.add_argument('--gp_inject', action = 'store_true', default=False)
    args = parser.parse_args()

    ana = load_data(dataset='combo')
    trials = run_gp_lc_test(ana, args.name, args.n_trials, args.n_sig, args.seed, args.cpus, args.gp_inject)


