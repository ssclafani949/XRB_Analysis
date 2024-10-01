#!/usr/bin/env python

from __future__ import print_function
import csky as cy
from csky import coord
import numpy as np
import pickle
import pandas as pd
import datetime
from submitter import Submitter
import histlite as hl
now = datetime.datetime.now
import matplotlib.pyplot as plt
import click, sys, os, time
flush = sys.stdout.flush


#repo = cy.selections.Repository(local_root='/data/i3store/users/analyses')
#ana_dir = cy.utils.ensure_dir('/data/i3store/users/ssclafani/XRB/analyses')
#base_dir = cy.utils.ensure_dir('/data/i3store/users/ssclafani/XRB_stacking_ss/')
repo = cy.selections.Repository()
ana_dir = cy.utils.ensure_dir('/data/user/ssclafani/XRB/analyses')
base_dir = cy.utils.ensure_dir('/data/user/ssclafani/XRB_stacking_ss_test/')
source_file  = '/home/ssclafani/XRB_Analysis/XRB/sources/lc_sources_reselected.hdf'
job_basedir = '/scratch/ssclafani/logs/' 


class State (object):
    def __init__ (self, ana_name, ana_dir, save,  base_dir):
        self.ana_name, self.ana_dir, self.save = ana_name, ana_dir, save
        self.base_dir = base_dir
        #self._ana = None

    @property
    def ana (self):
        print(self.ana_name)
        if self.ana_name == 'combo':
            cspec = cy.selections.DNNCascadeDataSpecs.DNNC_10yr
            psspec = cy.selections.PSDataSpecs.ps_v4
            ana = cy.get_analysis(repo, 'version-004-p02', psspec, 'version-001-p01', cspec,
                        dir=self.ana_dir)
            if self.save:
                cy.utils.ensure_dir (self.ana_dir)
                ana.save (self.ana_dir)
            ana.name = self.ana_name
            self._ana = ana
        elif self.ana_name == 'dnnc':
            cspec = cy.selections.DNNCascadeDataSpecs.DNNC_10yr
            ana = cy.get_analysis(repo, 'version-001-p01', cspec)
            if self.save:
                cy.utils.ensure_dir (self.ana_dir)
                ana.save (self.ana_dir)
            ana.name = self.ana_name
            self._ana = ana
        elif self.ana_name == 'gfu':
            spec = cy.selections.GFUDataSpecs.gfu_IC86_2018
            ana = cy.get_analysis(repo, 'version-002-p05', spec)
            if self.save:
                cy.utils.ensure_dir (self.ana_dir)
                ana.save (self.ana_dir)
            ana.name = self.ana_name
            self._ana = ana
        elif self.ana_name == 'tracks':
            psspec = cy.selections.PSDataSpecs.ps_v4
            ana = cy.get_analysis(repo, 'version-004-p02', psspec)
            if self.save:
                cy.utils.ensure_dir (self.ana_dir)
                ana.save (self.ana_dir)
            ana.name = self.ana_name
            self._ana = ana
        else:
            print('Ana Name {} not valid'.format(self.ana_name))
        return self._ana

    @property
    def state_args (self):
        return '--ana {} --ana-dir {} --base-dir {}'.format (
            self.ana_name, self.ana_dir, self.base_dir)

pass_state = click.make_pass_decorator (State)

@click.group (invoke_without_command=True, chain=True)
@click.option ('-a', '--ana', 'ana_name', default='combo', help='Dataset title')
@click.option ('--ana-dir', default=ana_dir, type=click.Path ())
@click.option ('--save/--nosave', default=False)
@click.option ('--base-dir', default=base_dir,
               type=click.Path (file_okay=False, writable=True))
@click.pass_context
def cli (ctx, ana_name, ana_dir, save, base_dir):
    ctx.obj = State.state = State (ana_name, ana_dir, save, base_dir)


@cli.result_callback ()
def report_timing (result, **kw):
    exe_t1 = now ()
    print ('c7: end at {} .'.format (exe_t1))
    print ('c7: {} elapsed.'.format (exe_t1 - exe_t0))

@cli.command ()
@pass_state
def setup_ana (state):
    state.ana

@cli.command()
@click.option('--n-trials', default=1000, type=int)
@click.option ('-n', '--n-sig', default=0, type=float)
@click.option ('--poisson/--nopoisson', default=True)
@click.option ('--fix_gamma', default=None)
@click.option ('--src_gamma', default=2.0, type=float)
@click.option ('--thresh', default=0.0, type=float)
@click.option ('--lag', default=0.0, type=float)
@click.option ('--name', default='Cyg_X_dash_1', type=str)
@click.option ('--seed', default=None, type=int)
@click.option ('--cutoff', default=np.inf, type=float, help='exponential cutoff energy, (TeV)')
@click.option ('--cpus', default=1, type=int)
@pass_state
def do_lc_trials ( state, name, n_trials, fix_gamma, src_gamma, thresh, lag, n_sig, poisson, seed, cutoff,  cpus, logging=True):
    ana = state.ana
    if seed is None:
        seed = int (time.time () % 2**32)
    random = cy.utils.get_random (seed)
   
    cutoff_GeV = cutoff * 1e3 
    sources = pd.read_hdf(source_file)
    source = sources.loc[sources['name_disp'] == name]
    print(source.dec_deg)
    print(source.dec_deg)
    src = cy.utils.Sources(dec = source.dec_deg, ra = source.ra_deg, deg=True)

    #def source(name):
    #    source_list  = np.genfromtxt("/home/qliu/binaries/analysis/xb_list/list/source_list.txt",dtype=str,delimiter=',')
    #    names        = [source_list[i][1] for i in range(len(source_list))]                                              
    #    names_display = [name.replace(' ', '_') for name in names] 
    #    index        = names.index(name)
    #    return       cy.utils.Sources(dec=float(source_list[index][5]),ra=float(source_list[index][6]),deg=True)
    
    #src          = source(name) 
    print('Source RA: {} DEC: {}'.format(np.degrees(src.ra[0]), np.degrees(src.dec[0])))

    #lightcurve   = np.genfromtxt("/data/user/qliu/binaries/SwiftBatLC/lc/watchdog_list/outburst/bin_outburst_fit_{}".format(name.replace(' ','').replace('+','p')))
    #bins         = lightcurve[:,0]
    #fluxes       = lightcurve[:,1]
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
    print(seed)
    print(f'lag: {lag}')
    print(f'thresh: {thresh}')
    dir = cy.utils.ensure_dir ('{}/lc/{}'.format (state.base_dir, name))
    lc = hl.Hist(bins_in_data, flux_in_data)
    if fix_gamma:
        print('Fixed gamma: {}'.format(fix_gamma))
        conf = dict(
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
                #inj_class=cy.inj.PointSourceBinnedTimeInjector,
                #inj_kw=dict( lcs=lc, threshs=thresh,lags=lag, gamma=src_gamma),
        )
    else:
        print('fitting gamma')
        conf = dict(
                kind = 'binned',                
                time='lc',
                lcs=lc,
                range_lag=(-7.,7.),
                #range_thresh=(0.,.5),
                sig='lc',
                concat_evs=True,
                extra_keep = ['energy'],
                #concat_evs=False,
                #use_grl=True,
                n_seeds_thresh=20,
                #n_seeds_thresh=40,
                n_seeds_lag = 20,
                update_bg = True,
                sigsub = True,
                #fitter_args=dict(_fmin_method='minuit'),
                #fitter_args= dict(lag=lag, gamma=src_gamma),
                #sig_kw = dict(gamma=src_gamma),
                #cut_n_sigma=3,
                #inj_class=cy.inj.PointSourceBinnedTimeInjector,
                #inj_kw=dict( lcs=lc, threshs=thresh,lags=lag, gamma=src_gamma),
        )
    inj_conf=dict(lcs=lc, threshs=thresh, lags=lag, flux = cy.hyp.PowerLawFlux(gamma=src_gamma))
    tr = cy.get_trial_runner(conf=conf,  inj_conf = inj_conf, ana=ana, src=src, flux=cy.hyp.PowerLawFlux(src_gamma, energy_cutoff = cutoff_GeV), mp_cpus=cpus, dir=dir)
    cy.describe(tr)
    t0 = now ()
    print ('Beginning trials at {} ...'.format (t0))
    flush ()
    trials = tr.get_many_fits (
        n_trials, n_sig=n_sig, poisson=poisson, seed=seed, logging=logging)
    t1 = now ()
    print ('Finished trials at {} ...'.format (t1))
    print (trials if n_sig else cy.dists.Chi2TSD (trials))
    print (t1 - t0, 'elapsed.')
    flush ()
    if n_sig:
        if fix_gamma:
            out_dir = cy.utils.ensure_dir (
                    '{}/{}/lc/{}/trials/fix_gamma_{}/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/nsig/{}'.format (
                        state.base_dir, state.ana_name, name, fix_gamma, src_gamma, thresh, lag, cutoff,
                          n_sig))
        else:
            out_dir = cy.utils.ensure_dir (
                '{}/{}/lc/{}/trials/fit_gamma/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/nsig/{}'.format (
                    state.base_dir, state.ana_name, name,  src_gamma, thresh, lag, cutoff,
                      n_sig))
    else:        
        if fix_gamma:
            out_dir = cy.utils.ensure_dir ('{}/{}/lc/{}/trials/fix_gamma_{}/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/bg'.format (
                state.base_dir, state.ana_name, name, fix_gamma, src_gamma, thresh, lag, cutoff))
 
        else:
            out_dir = cy.utils.ensure_dir ('{}/{}/lc/{}/trials/fit_gamma/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/bg'.format (
                state.base_dir, state.ana_name, name,  src_gamma, thresh, lag, cutoff))

    out_file = '{}/trials_{:07d}__seed_{:010d}.npy'.format (
        out_dir, n_trials, seed)
    print ('-> {}'.format (out_file))
    np.save (out_file, trials.as_array)


@cli.command()
@click.option('--n-trials', default=1000, type=int)
@click.option ('-n', '--n-sig', default=0, type=float)
@click.option ('--poisson/--nopoisson', default=True)
@click.option ('--fix_gamma', default=None)
@click.option ('--src_gamma', default=2.0, type=float)
@click.option ('--thresh', default=0.0, type=float)
@click.option ('--lag', default=0.0, type=float)
@click.option ('--seed', default=None, type=int)
@click.option ('--cutoff', default=np.inf, type=float, help='exponential cutoff energy, (TeV)')
@click.option ('--cpus', default=1, type=int)
@pass_state
def do_stacking_trials ( state, n_trials, fix_gamma, src_gamma, thresh, lag, n_sig, poisson, seed, cutoff, cpus, logging=True):
    ana = state.ana
    if seed is None:
        seed = int (time.time () % 2**32)
    random = cy.utils.get_random (seed)
   
    cutoff_GeV = cutoff * 1e3 
    sources = pd.read_hdf(source_file)
    ras = []
    decs = []
    lcs = []
    for name in sources.name_disp:
        print(name)
        source = sources.loc[sources['name_disp'] == name]
        bins = np.array(source.lc_bins_10)[0]
        fluxes = np.array(source.lc_values_10)[0]
        index  = np.where((bins>=ana.mjd_min)& (bins<=ana.mjd_max))[0]
        if len(index) > 2:
            ras.append(source.ra_deg)
            decs.append(source.dec_deg)
            if index[0] != 0:
                    bins_in_data = np.append(np.append(ana.mjd_min-7.,bins[index]),ana.mjd_max+7.)
                    flux_in_data = np.append(fluxes[index[0]-1],fluxes[index])
            else:
                    bins_in_data = np.append(bins[index],ana.mjd_max+7.)
                    flux_in_data = fluxes[index] 
            lc = hl.Hist(bins_in_data, flux_in_data)
            lcs.append(lc)
    print(len(decs), len(ras), len(lcs))
    src = cy.utils.Sources(dec = np.concatenate(decs), ra = np.concatenate(ras), deg=True)
    print(seed)
    print(f'inj lag: {lag}')
    print(f'inj thresh: {thresh}')
    dir = cy.utils.ensure_dir ('{}/lc/{}'.format (state.base_dir, name))

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
                kind = 'binned',                
                time='lc',
                lcs=lcs,
                #range_lag=(-7.,7.),
                #range_thresh=(0.,.5),
                sig='lc',
                concat_evs=True,
                fitter_args = fitter_dict,
                extra_keep = ['energy'],
                #concat_evs=False,
                use_grl=True,
                n_seeds_thresh=20,
                sigsub = True,
                n_seeds_lag = 20,
                update_bg = True,
                inj_class=cy.inj.PointSourceBinnedTimeInjector,
                inj_kw=dict(src,  lcs=lc, threshs=thresh, lags=lag, gamma=src_gamma),
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
                #range_thresh=(0.,.5),
                sig='lc',
                #concat_evs=True,
                extra_keep = ['energy'],
                #concat_evs=False,
                #use_grl=True,
                #n_seeds_thresh=60,
                n_seeds_lag = 20,
                update_bg = True,
                sigsub = True,
                #gammas = np.arange(0,4.01,.0125),
                fitter_args = fitter_dict,
                #fitter_args=dict(_fmin_method='minuit'),
                #fitter_args= dict(lag=lag, gamma=src_gamma),
                #sig_kw = dict(gamma=src_gamma),
                #cut_n_sigma=3,
                #inj_class=cy.inj.PointSourceBinnedTimeInjector,
                #inj_kw=dict(src, lcs=lcs, gamma=src_gamma),
        )
    inj_conf=dict(lcs=lcs, threshs=np.zeros_like(lcs), lags=np.zeros_like(lcs), flux = cy.hyp.PowerLawFlux(gamma=src_gamma))
    tr = cy.get_trial_runner(conf=conf,  inj_conf=inj_conf, ana=ana, #flux=cy.hyp.PowerLawFlux(src_gamma, energy_cutoff = cutoff_GeV), 
                                mp_cpus=cpus, dir=dir)
    t0 = now ()
    #cy.describe(tr)
    print ('Beginning trials at {} ...'.format (t0))
    flush ()
    trials = tr.get_many_fits (
        n_trials, n_sig=n_sig, poisson=poisson, seed=seed, logging=logging)
    print(f'TS: {trials.ts}')
    print(f'ns: {trials.ns}')
    print(trials.seed)
    if 'gamma' in trials.keys():
        print(f'Gamma: {trials.gamma}')
    t1 = now ()
    print ('Finished trials at {} ...'.format (t1))
    print (trials if n_sig else cy.dists.Chi2TSD (trials))
    print(np.median(trials['ns']))
    print (t1 - t0, 'elapsed.')
    flush ()
    if n_sig:
        if fix_gamma:
            out_dir = cy.utils.ensure_dir (
                    '{}/{}/lc/stacking/trials/fix_gamma_{}/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/nsig/{}'.format (
                        state.base_dir, state.ana_name, fix_gamma, src_gamma, thresh, lag, cutoff,
                          n_sig))
        else:
            out_dir = cy.utils.ensure_dir (
                '{}/{}/lc/stacking/trials/fit_gamma/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/nsig/{}'.format (
                    state.base_dir, state.ana_name,  src_gamma, thresh, lag, cutoff,
                      n_sig))
    else:        
        if fix_gamma:
            out_dir = cy.utils.ensure_dir ('{}/{}/lc/stacking/trials/fix_gamma_{}/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/bg'.format (
                state.base_dir, state.ana_name, fix_gamma, src_gamma, thresh, lag, cutoff))
 
        else:
            out_dir = cy.utils.ensure_dir ('{}/{}/lc/stacking/trials/fit_gamma/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/bg'.format (
                state.base_dir, state.ana_name, src_gamma, thresh, lag, cutoff))

    out_file = '{}/trials_{:07d}__seed_{:010d}.npy'.format (
        out_dir, n_trials, seed)
    print ('-> {}'.format (out_file))
    np.save (out_file, trials.as_array)


@cli.command ()
@click.option ('--fit/--nofit', default=True)
@click.option ('--hist/--nohist', default=False)
@click.option ('-n', '--n', default=0, type=int)
@pass_state
def collect_lc_trials (state, fit, hist, n):
    print(state.ana_name)    
    sources = pd.read_hdf(source_file)
    names = sources.name_disp
    #names = names[:50]
    #names = np.r_[np.array(names[:10]), ['GX_339_dash_4', 'GX_1_plus_4']]
    #names = ['Cyg_X_dash_1'] 
    kw = {}
    if hist:
        TSD = cy.dists.BinnedTSD
        suffix = '_hist'
        kw['keep_trials'] = False
    elif fit:
        TSD = cy.dists.Chi2TSD
        suffix = '_chi2'
    else:
        TSD = cy.dists.TSD
        suffix = ''
    combined_outfile = '{}/{}/lc/TSD{}.dict'.format (
        state.base_dir, state.ana_name,  suffix)
    bg = {}
    bgs = {}                                                                    
    for name in names:    
        outfile = '{}/{}/lc/{}/TSD{}.dict'.format (
            state.base_dir, state.ana_name, name, suffix)
        key = name 
        print ('\r{} ...'.format (key))
        #flush ()
        print('{}/{}/lc/{}'.format(state.base_dir, state.ana_name, key))
        #if (not bg['dec'][dec_deg]) or bg['dec'][dec_deg].n_total < n:
        #print('{}/dec/{}/'.format(bg_dir, key))
        #print(cy.bk.get_all ('{}/dec/{}/'.format (bg_dir, key), '*.npy', merge=np.concatenate))
        post_convert = (lambda x: TSD (cy.utils.Arrays (x), **kw))
        bgs['{}'.format(key)] = cy.bk.get_all('{}/{}/lc/{}/trials/'.format(state.base_dir, state.ana_name, key),'trials*.npy', 
                          merge=np.concatenate, post_convert=post_convert)
        #bg['name'] = bgs             
        #print ('\rDone.' + 20 * ' ')
        #flush ()
        #print ('->', outfile)
        #with open (outfile, 'wb') as f:
        #     pickle.dump (bg, f, -1)
    print(bgs.keys())
    bg['name'] = bgs             
    print ('\rDone.' + 20 * ' ')
    flush ()
    print ('->', combined_outfile)
    with open (combined_outfile, 'wb') as f:
         pickle.dump (bg, f, -1)

@cli.command ()
@click.option ('--fit/--nofit', default=True)
@click.option ('--hist/--nohist', default=False)
@click.option ('-n', '--n', default=0, type=int)
@pass_state
def collect_stacking_trials (state, fit, hist, n):
    kw = {}
    if hist:
        TSD = cy.dists.BinnedTSD
        suffix = '_hist'
        kw['keep_trials'] = False
    elif fit:
        TSD = cy.dists.Chi2TSD
        suffix = '_chi2'
    else:
        TSD = cy.dists.TSD
        suffix = ''
    outfile = '{}/{}/lc/stacking/TSD{}.dict'.format (
        state.base_dir, state.ana_name,  suffix)
    bg = {}
    bgs = {}                                            
    print('{}/{}/lc/'.format(state.base_dir, state.ana_name))
    post_convert = (lambda x: TSD (cy.utils.Arrays (x), **kw))
    bg = cy.bk.get_all('{}/{}/lc/stacking/trials/'.format(state.base_dir, state.ana_name),'trials*.npy', 
                      merge=np.concatenate, post_convert=post_convert)
    print ('\rDone.' + 20 * ' ')
    flush ()
    print ('->', outfile)
    with open (outfile, 'wb') as f:
         pickle.dump (bg, f, -1)


@cli.command ()
@click.option ('--n-trials', default=10000, type=int)
@click.option ('--n-jobs', default=10, type=int)
@click.option ('-n', '--n-sig', 'n_sigs', multiple=True, default=[0], type=float)
@click.option ('--src_gamma',  default=2.0, type=float)
@click.option ('--fix_gamma',  default=None, type=float)
@click.option ('--poisson/--nopoisson', default=True)
@click.option ('--cutoff', default=np.inf, type=float, help='exponential cutoff energy, (TeV)')
@click.option ('--dry/--nodry', default=False)
@click.option ('--seed', default=0)
@pass_state
def submit_do_lc_trials (
        state, n_trials, n_jobs, n_sigs, src_gamma, fix_gamma, poisson, cutoff, dry, seed):
    ana_name = state.ana_name
    T = time.time ()
    #job_basedir = '/scratch/ssclafani/logs/' 
    poisson_str = 'poisson' if poisson else 'nopoisson'
    job_dir = '{}/{}/lc_trials/T_{:17.6f}'.format (
        job_basedir, ana_name,  T)
    sub = Submitter (job_dir=job_dir, memory=9, max_jobs=1000)
    #env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'
    commands, labels = [], []
    this_script = os.path.abspath (__file__)
    #lags = np.arange(-6, 6.1, 3)
    lags = [0, 6] 
    threshs = [0]
    src_gammas = [src_gamma]
    fix_gammas = [fix_gamma] 
    
    sources = pd.read_hdf(source_file)
    names = sources.name_disp[:5]
    #names = ['4U_2206_plus_54', 'Cyg_X_dash_1', 'Cen_X_dash_3', 'Cyg_X_dash_3', 'Sgr_Astar', 'Vela_X_dash_1']
    #names = np.r_[np.array(names[5:])]
    #names = ['XTE_J1855_dash_026']
    #names = ['Cyg_X_dash_1']
    if fix_gamma:
        for name in names:
            for src_gamma in src_gammas:
                for fix_gamma in fix_gammas:
                    for thresh in threshs:
                        for lag in lags:
                            for n_sig in n_sigs:
                                for i in range (n_jobs):
                                    s = i + seed
                                    #fmt = '{} --ana-dir \'\' --base-dir={} do-ps-trials --dec={:+08.3f} --n-trials={}'
                                    fmt = ' {} {} do-lc-trials --name {} --fix_gamma={} --src_gamma={} --cutoff {} --n-trials={}' \
                                            ' --n-sig={} --thresh={} --lag={}' \
                                            ' --{} --seed={}'
                                    command = fmt.format ( this_script, state.state_args, name, fix_gamma, src_gamma, cutoff, n_trials,
                                                          n_sig, thresh, lag, poisson_str, s)
                                    fmt = 'xrb{}_fix_gamma_{}__trials_{:07d}__n_sig_{:08.3f}__' \
                                            'src_gamma_{:.3f}_thresh_{}_lag_{}_cutoff_{}seed_{:04d}'
                                    label = fmt.format (name, fix_gamma, n_trials, n_sig, src_gamma,
                                                        thresh, lag, cutoff,  s)
                                    #print(label)
                                    commands.append (command)
                                    labels.append (label)
    else:
        for name in names:
            for src_gamma in src_gammas:
                for thresh in threshs:
                    for lag in lags:
                        for n_sig in n_sigs:
                            for i in range (n_jobs):
                                s = i + seed
                                #fmt = '{} --ana-dir \'\' --base-dir={} do-ps-trials --dec={:+08.3f} --n-trials={}'
                                fmt = '{} {} do-lc-trials --name {}  --src_gamma={} --n-trials={}' \
                                        ' --n-sig={} --thresh={} --lag={}' \
                                        ' --{} --seed={}'
                                #command = fmt.format (env_shell, this_script, state.state_args, name, src_gamma, n_trials,
                                #                      n_sig, thresh, lag, poisson_str, s)
                                command = fmt.format ( this_script, state.state_args, name, src_gamma, n_trials,
                                                      n_sig, thresh, lag, poisson_str, s)
                                fmt = 'xrb{}_fit_gamma__trials_{:07d}__n_sig_{:08.3f}__' \
                                        'src_gamma_{:.3f}_thresh_{}_lag_{}_seed_{:04d}'
                                label = fmt.format (name,  n_trials, n_sig, src_gamma,
                                                    thresh, lag,  s)
                                print(label)
                                commands.append (command)
                                labels.append (label)

    sub.dry = dry
    sub.submit_npx4 (commands, labels)

@cli.command ()
@click.option ('--n-trials', default=10000, type=int)
@click.option ('--n-jobs', default=10, type=int)
@click.option ('-n', '--n-sig', 'n_sigs', multiple=True, default=[0], type=float)
@click.option ('--src_gamma',  default=2.0, type=float)
@click.option ('--fix_gamma',  default=None, type=float)
@click.option ('--poisson/--nopoisson', default=True)
@click.option ('--cutoff', default=np.inf, type=float, help='exponential cutoff energy, (TeV)')
@click.option ('--dry/--nodry', default=False)
@click.option ('--seed', default=0)
@pass_state
def submit_do_stacking_trials (
        state, n_trials, n_jobs, n_sigs, src_gamma, fix_gamma, poisson, cutoff, dry, seed):
    ana_name = state.ana_name
    T = time.time ()
    #job_basedir = '/scratch/ssclafani/logs/' 
    poisson_str = 'poisson' if poisson else 'nopoisson'
    job_dir = '{}/{}/lc_trials/T_{:17.6f}'.format (
        job_basedir, ana_name,  T)
    sub = Submitter (job_dir=job_dir, ncpu=2, memory=13, max_jobs=1000)
    #env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'
    commands, labels = [], []
    this_script = os.path.abspath (__file__)
    
    if fix_gamma:
        for n_sig in n_sigs:
            for i in range (n_jobs):
                s = i + seed
                fmt = ' {} {} do-stacking-trials --fix_gamma={} --src_gamma={} --cutoff {} --n-trials={}' \
                        ' --n-sig={}' \
                        ' --{} --seed={}'
                command = fmt.format ( this_script, state.state_args, fix_gamma, src_gamma, cutoff, n_trials,
                                      n_sig, poisson_str, s)
                fmt = 'xrb_fix_gamma_{}__trials_{:07d}__n_sig_{:08.3f}__' \
                        'src_gamma_{:.3f}__cutoff_{}seed_{:04d}'
                label = fmt.format ( fix_gamma, n_trials, n_sig, src_gamma,
                                    cutoff,  s)
                commands.append (command)
                labels.append (label)
    else:
        for n_sig in n_sigs:
            for i in range (n_jobs):
                s = i + seed
                fmt = '{} {} do-stacking-trials  --src_gamma={} --n-trials={}' \
                        ' --n-sig={}' \
                        ' --{} --seed={}'
                command = fmt.format ( this_script, state.state_args,  src_gamma, n_trials,
                                      n_sig,  poisson_str, s)
                fmt = 'xrb_stacking_fit_gamma__trials_{:07d}__n_sig_{:08.3f}__' \
                        'src_gamma_{:.3f}_seed_{:04d}'
                label = fmt.format (n_trials, n_sig, src_gamma, s)
                print(label)
                commands.append (command)
                labels.append (label)

    sub.dry = dry
    sub.submit_npx4 (commands, labels)




if __name__ == '__main__':
    exe_t0 = now ()
    print ('c7: start at {} .'.format (exe_t0))
    cli ()
