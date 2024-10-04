#!/usr/bin/env python

import csky as cy
import numpy as np
import pickle
import pandas as pd
import datetime
from submitter import Submitter
import histlite as hl
now = datetime.datetime.now
import click, sys, os, time
flush = sys.stdout.flush
import config as cg
import socket

hostname = socket.gethostname()

repo, ana_dir, base_dir, job_basedir, source_file = cg.repo, cg.ana_dir, cg.base_dir, cg.job_basedir, cg.source_file
#repo = cy.selections.Repository(local_root='/data/i3store/users/analyses')
#ana_dir = cy.utils.ensure_dir('/data/i3store/users/ssclafani/XRB/analyses')
#base_dir = cy.utils.ensure_dir('/data/i3store/users/ssclafani/XRB_stacking_ss/')
#repo = cy.selections.Repository()
#ana_dir = cy.utils.ensure_dir('/data/user/ssclafani/XRB/analyses')
#base_dir = cy.utils.ensure_dir('/data/user/ssclafani/XRB_stacking_ss_test/')
#source_file  = '/home/ssclafani/XRB_Analysis/XRB/sources/lc_sources_reselected.hdf'
#job_basedir = '/scratch/ssclafani/logs/' 


class State (object):
    def __init__ (self, ana_name, ana_dir, save,  base_dir):
        self.ana_name, self.ana_dir, self.save = ana_name, ana_dir, save
        self.base_dir = base_dir
        #self._ana = None

    @property
    def ana (self):
        print(self.ana_name)
        if self.ana_name == 'combo':
            print(repo.local_root)
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
@click.option('--n-trials', default=100, type=int)
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
    print(seed)
    conf, inj_conf  =  cg.get_ps_config(
                                ana,
                                name, 
                                src_gamma, 
                                fix_gamma, 
                                cutoff_GeV, 
                                lag, 
                                thresh
                                ) 
    tr = cy.get_trial_runner(
                        conf=conf,
                        inj_conf = inj_conf, 
                        ana=ana, 
                        flux=cy.hyp.PowerLawFlux(src_gamma, energy_cutoff = cutoff_GeV), 
                        mp_cpus=cpus, 
                        dir=dir
                        )

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
    conf, inj_conf = cg.get_stacking_config(
                                    ana,
                                    src_gamma, 
                                    fix_gamma, 
                                    thresh, 
                                    lag
                                    )
    tr = cy.get_trial_runner(
                            conf=conf,  
                            inj_conf=inj_conf, 
                            ana=ana,
                            flux=cy.hyp.PowerLawFlux(src_gamma, energy_cutoff = cutoff_GeV), 
                            dir=dir)
    t0 = now ()
    print ('Beginning trials at {} ...'.format (t0))
    flush ()
    trials = tr.get_many_fits (
        n_trials, n_sig=n_sig, poisson=poisson, seed=seed, logging=logging, mp_cpus = cpus)
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
    sources = pd.read_hdf(source_file)
    names = sources.name_disp
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
    poisson_str = 'poisson' if poisson else 'nopoisson'
    job_dir = '{}/{}/lc_trials/T_{:17.6f}'.format (
        job_basedir, ana_name,  T)
    sub = Submitter (job_dir=job_dir, memory=9, max_jobs=1000, config = submit_cfg_file)
    #env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'
    commands, labels = [], []
    this_script = os.path.abspath (__file__)
    #lags = np.arange(-6, 6.1, 3)
    lags = [0] 
    threshs = [0]
    src_gammas = [src_gamma]
    fix_gammas = [fix_gamma] 
    
    sources = pd.read_hdf(cg.source_file)
    names = sources.name_disp
    #names = sources.name_disp[:5]
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
    if 'condor' in hostname:
        sub.submit_condor00 (commands, labels)
    else:
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
    if 'condor' in hostname:
        sub.submit_condor00 (commands, labels)
    else:
        sub.submit_npx4 (commands, labels)

@cli.command()
@click.option ('--fix_gamma', default=None)
@click.option ('--src_gamma', default=2.0, type=float)
@click.option ('--nsigma', default=None, type=float)
@click.option ('--thresh', default=0.0, type=float)
@click.option ('--lag', default=0.0, type=float)
@click.option ('--cutoff', default=np.inf, type=float, help='exponential cutoff energy, (TeV)')
@click.option ('--cpus', default=1, type=int)
@click.option ('--save/--nosave', default=False)
@pass_state
def find_stacking_nsig ( state,fix_gamma, src_gamma, nsigma, thresh, lag, cutoff, cpus, logging=True, save=False):
    ana = state.ana
    cutoff_GeV = 1e3 * cutoff
    sigfile = '{}/{}/lc/stacking/TSD_chi2.dict'.format (state.base_dir, state.ana_name)
    sig = np.load (sigfile, allow_pickle=True)
    if fix_gamma:
        sig_trials = cy.bk.get_best(sig, 'fix_gamma_{}'.format(fit_gamma),  'src_gamma_{}'.format(src_gamma),
                                                'thresh_{}'.format(thres), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff), 'nsig')
        b = cy.bk.get_best(sig, 'fix_gamma_{}'.format(fit_gamma),  'src_gamma_{}'.format(src_gamma),
                                                    'thresh_{}'.format(thresh), 'lag_{}'.format(lag),  'cutoff_{}'.format(cutoff), 'bg')
    else:
        print('fit gamma')
        sig_trials = cy.bk.get_best(sig, 'fit_gamma',  'src_gamma_{}'.format(src_gamma),
                                                    'thresh_{}'.format(thresh), 'lag_{}'.format(lag),  'cutoff_{}'.format(cutoff), 'nsig')
    
        b = cy.bk.get_best(sig,  'fit_gamma',  'src_gamma_{}'.format(src_gamma),
                                                    'thresh_{}'.format(thresh), 'lag_{}'.format(lag),  'cutoff_{}'.format(cutoff), 'bg')
    if logging:
        print(b)
    conf, inj_conf = cg.get_stacking_config(
                                    ana,
                                    src_gamma, 
                                    fix_gamma, 
                                    thresh, 
                                    lag
                                    )
    tr = cy.get_trial_runner(
                            conf=conf,  
                            inj_conf=inj_conf, 
                            ana=ana,
                            flux=cy.hyp.PowerLawFlux(src_gamma, energy_cutoff = cutoff_GeV), 
                            dir=dir)
        # determine ts threshold
    if nsigma !=None:
        beta = 0.5
        print('sigma = {}'.format(nsigma))
        ts = b.isf_nsigma(nsigma)
    else:
        print('Getting sensitivity')
        beta = 0.9
        ts = b.median()
    
    # include background trials in calculation
    trials = {0: b}
    trials.update(sig_trials)
    for key in trials.keys():
        trials[key] = trials[key].trials
    # get number of signal events
    # (arguments prevent additional trials from being run)
    result = tr.find_n_sig(ts, beta, max_batch_size=0, logging=logging, trials=trials)
    
    f = tr.to_E2dNdE(result['n_sig'], E0=1, unit=1e3)
    if logging:
        print(ts, beta, result['n_sig'], f)
    print(f) 
    print(result['n_sig'])
    if save:
        if nsigma:
            np.save(f'{state.base_dir}/{state.ana_name}/lc/stacking/trials/fit_gamma/src_gamma_{src_gamma}/{nsigma}sigma_dp.npy')
        else:
            np.save(f'{state.base_dir}/{state.ana_name}/lc/stacking/trials/fit_gamma/src_gamma_{src_gamma}/sens.npy')


@cli.command()
@click.option ('--fix_gamma', default=None)
@click.option ('--src_gamma', default=2.0, type=float)
@click.option ('--nsigma', default=None, type=float)
@click.option ('--thresh', default=0.0, type=float)
@click.option ('--lag', default=0.0, type=float)
@click.option ('--cutoff', default=np.inf, type=float, help='exponential cutoff energy, (TeV)')
@click.option ('--cpus', default=1, type=int)
@click.option ('--save/--nosave', default=False)
@pass_state
def plot_stacking_bias ( state,fix_gamma, src_gamma, nsigma, thresh, lag, cutoff, cpus, logging=True, save=False):
    import matplotlib.pyplot as plt

    cutoff_GeV = 1e3 * cutoff
    sigfile = '{}/{}/lc/stacking/TSD_chi2.dict'.format (state.base_dir, state.ana_name)
    sig = np.load (sigfile, allow_pickle=True)
    if fix_gamma:
        sig_trials = cy.bk.get_best(sig, 'fix_gamma_{}'.format(fit_gamma),  'src_gamma_{}'.format(src_gamma),
                                                'thresh_{}'.format(thres), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff), 'nsig')
        b = cy.bk.get_best(sig, 'fix_gamma_{}'.format(fit_gamma),  'src_gamma_{}'.format(src_gamma),
                                                    'thresh_{}'.format(thresh), 'lag_{}'.format(lag),  'cutoff_{}'.format(cutoff), 'bg')
    else:
        print('fit gamma')
        sig_trials = cy.bk.get_best(sig, 'fit_gamma',  'src_gamma_{}'.format(src_gamma),
                                                    'thresh_{}'.format(thresh), 'lag_{}'.format(lag),  'cutoff_{}'.format(cutoff), 'nsig')
    
        b = cy.bk.get_best(sig,  'fit_gamma',  'src_gamma_{}'.format(src_gamma),
                                                    'thresh_{}'.format(thresh), 'lag_{}'.format(lag),  'cutoff_{}'.format(cutoff), 'bg')
    if logging:
        print(b)
    
    trials = {0: b}
    trials.update(sig_trials)
    n_sigs = [n for n in sorted(trials.keys())]
    
    for key in sorted(trials.keys()):
        trials[key] = trials[key].trials
    trials = [trials[k] for k in sorted(trials.keys())]
    
    for (n_sig, t) in zip(sorted(n_sigs), trials):
        t['ntrue'] = np.repeat(n_sig, len(t))
    
    allt = cy.utils.Arrays.concatenate(trials)
    
    if logging:
        for ntrue in sorted(np.unique(allt['ntrue'])):
            mask = allt['ntrue'] == ntrue
            median = np.median(allt[mask].ns)
            gamma = np.median(allt[mask].gamma)
            print(f'Inj: {ntrue:3.2f}   | Median: {median:.2f}   | gamma {gamma:.2f}')
 
    if fix_gamma:
        fig, ax = plt.subplots(1, 1, figsize=(4,3))
    else:
        fig, axs = plt.subplots(1, 2, figsize=(6,3))
    
    dns = np.mean(np.diff(n_sigs))
    ns_bins = np.r_[n_sigs - 0.5*dns, n_sigs[-1] + 0.5*dns]
    expect_kw = dict(color='C0', ls='--', lw=1, zorder=-10)
    if not fix_gamma:
        expect_gamma = src_gamma 
        ax = axs[0]
        ax.set_xlim(0,500)
        ax.set_ylim(0,500)
    else:
        ax.set_xlim(0,200)
        ax.set_ylim(0,200)
    h = hl.hist((allt.ntrue, allt.ns), bins=(ns_bins, 70))
    hl.plot1d(ax, h.contain_project(1),errorbands=True, drawstyle='default')
    ax.set_xlabel(r'n$_{inj}$')
    ax.set_ylabel(r'n$_s$')
    ax.set_aspect('equal')
    
    lim = ns_bins[[0, -1]]
    ax.set_xlim(ax.set_ylim(lim))
    ax.plot(lim, lim, **expect_kw)
    ax.set_aspect('equal')


    if not fix_gamma:
        ax = axs[1]
        h = hl.hist((allt.ntrue, allt.gamma), bins=(ns_bins, 70))
        hl.plot1d(ax, h.contain_project(1),errorbands=True, drawstyle='default')
        ax.axhline(expect_gamma, **expect_kw)
        ax.set_xlim(axs[0].get_xlim())
        ax.set_xlabel(r'n$_{inj}$')
        ax.set_ylabel(r'$\gamma$')
        #ax.set_aspect('equal')

    plt.tight_layout()
    if save:
        np.save(cy.utils.ensure_dir(f'{state.base_dir}/{state.ana_name}/lc/stacking/plots/fit_gamma/src_gamma_{src_gamma}/bias.png'))
                

if __name__ == '__main__':
    exe_t0 = now ()
    print ('c7: start at {} .'.format (exe_t0))
    cli ()
