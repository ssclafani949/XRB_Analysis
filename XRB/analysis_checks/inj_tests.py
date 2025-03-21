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
submit_cfg_file = cg.submit_cfg_file
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
            psspec = cy.selections.PSDataSpecs.ps_v4[3:]
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
@click.option ('--inj_name', default='Cyg_X_dash_1', type=str)
@click.option ('--test_name', default='IGR_J00370_plus_6122', type=str)
@click.option ('--seed', default=None, type=int)
@click.option ('--cutoff', default=np.inf, type=float, help='exponential cutoff energy, (TeV)')
@click.option ('--cpus', default=1, type=int)
@pass_state
def make_inj_trials ( state, inj_name,test_name, n_trials, 
    fix_gamma, src_gamma, thresh, lag, n_sig, poisson, seed, cutoff,  cpus, logging=True):
    ana = state.ana
    if seed is None:
        seed = int (time.time () % 2**32)
    random = cy.utils.get_random (seed)
   
    cutoff_GeV = cutoff * 1e3 
    print(seed)
    test_conf, _  =  cg.get_ps_config(
                                ana,
                                test_name, 
                                src_gamma, 
                                fix_gamma, 
                                cutoff_GeV, 
                                lag, 
                                thresh
                                ) 
    inj_conf, _  =  cg.get_ps_config(
                                ana,
                                inj_name, 
                                src_gamma, 
                                fix_gamma, 
                                cutoff_GeV, 
                                lag, 
                                thresh
                                ) 
    tr = cy.get_trial_runner(
                        conf=test_conf,
                        inj_conf = inj_conf, 
                        ana=ana, 
                        flux=cy.hyp.PowerLawFlux(src_gamma, energy_cutoff = cutoff_GeV), 
                        mp_cpus=cpus, 
                        dir=dir
                        )
    trials = tr.get_many_fits(n_trials, n_sig, logging=logging)
    t0 = now ()
    print ('Beginning trials at {} ...'.format (t0))
    flush ()
    if n_sig:
        if fix_gamma:
            out_dir = cy.utils.ensure_dir (
                '{}/{}/lc/inj_name/{}/test_name/{}/inj_trials/fix_gamma_{}/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/nsig/{}'.format (
                         state.base_dir, state.ana_name, inj_name, test_name,  fix_gamma, src_gamma, thresh, lag, cutoff,
                           n_sig))
        else:
             out_dir = cy.utils.ensure_dir (
                 '{}/{}/lc/inj_name/{}/test_name/{}/inj_trials/fit_gamma/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/nsig/{}'.format (
                     state.base_dir, state.ana_name, inj_name, test_name,  src_gamma, thresh, lag, cutoff, n_sig))
    else:        
        if fix_gamma:
             out_dir = cy.utils.ensure_dir (
                '{}/{}/lc/inj_name/{}/test_name/{}/inj_trials/fix_gamma_{}/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/bg'.format (
                 state.base_dir, state.ana_name, inj_name, test_name, fix_gamma, src_gamma, thresh, lag, cutoff))
  
        else:
             out_dir = cy.utils.ensure_dir (
                '{}/{}/lc/inj_name/{}/test_name/{}/inj_trials/fit_gamma/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/bg'.format (
                state.base_dir, state.ana_name, inj_name, test_name,  src_gamma, thresh, lag, cutoff))                                         
 
    out_file = '{}/trials_{:07d}__seed_{:010d}.npy'.format (
         out_dir, n_trials, seed)
    print ('-> {}'.format (out_file))
    np.save (out_file, trials.as_array)

@cli.command()
@click.option('--n-trials', default=100, type=int)
@click.option ('-n', '--n-sig', default=0, type=float)
@click.option ('--poisson/--nopoisson', default=True)
@click.option ('--fix_gamma', default=None)
@click.option ('--src_gamma', default=2.0, type=float)
@click.option ('--thresh', default=0.0, type=float)
@click.option ('--lag', default=0.0, type=float)
@click.option ('--inj_name_1', default='AX_J1749_dot_1_dash_2639', type=str)
@click.option ('--inj_name_2', default='Swift_J1745_dot_1_dash_2624', type=str)
@click.option ('--test_name', default='Swift_J1745_dot_1_dash_2624', type=str)
@click.option ('--seed', default=None, type=int)
@click.option ('--cutoff', default=np.inf, type=float, help='exponential cutoff energy, (TeV)')
@click.option ('--cpus', default=1, type=int)
@pass_state
def make_inj_trials_stacking ( state, inj_name_1, inj_name_2 ,test_name, n_trials, 
    fix_gamma, src_gamma, thresh, lag, n_sig, poisson, seed, cutoff,  cpus, logging=True):
    ana = state.ana
    if seed is None:
        seed = int (time.time () % 2**32)
    random = cy.utils.get_random (seed)
   
    cutoff_GeV = cutoff * 1e3 
    print(seed)
    test_conf, _  =  cg.get_ps_config(
                                ana,
                                test_name, 
                                src_gamma, 
                                fix_gamma, 
                                cutoff_GeV, 
                                lag, 
                                thresh
                                ) 
    inj_conf, _  =  cg.get_stacking_config_AB_test(
                                ana,
                                src_gamma,
                                fix_gamma,
                                inj_name_1, 
                                inj_name_2, 
                                thresh,
                                lag,
                                weight = 'equal'
                                ) 
    tr = cy.get_trial_runner(
                        conf=test_conf,
                        inj_conf = inj_conf, 
                        ana=ana, 
                        flux=cy.hyp.PowerLawFlux(src_gamma, energy_cutoff = cutoff_GeV), 
                        mp_cpus=cpus, 
                        dir=dir
                        )
    trials = tr.get_many_fits(n_trials, n_sig, logging=logging)
    t0 = now ()
    print ('Beginning trials at {} ...'.format (t0))
    flush ()
    if n_sig:
        out_dir = cy.utils.ensure_dir (
                 '{}/{}/lc/inj_name/{}_{}/test_name/{}/inj_trials/fit_gamma/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/nsig/{}'.format (
                     state.base_dir, state.ana_name, inj_name_1, inj_name_2, test_name,  src_gamma, thresh, lag, cutoff, n_sig))
    else:        
        out_dir = cy.utils.ensure_dir (
                '{}/{}/lc/inj_name/{}_{}/test_name/{}/inj_trials/fit_gamma/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/bg'.format (
                state.base_dir, state.ana_name, inj_name_1, inj_name_2, test_name,  src_gamma, thresh, lag, cutoff))                                         
 
    out_file = '{}/trials_{:07d}__seed_{:010d}.npy'.format (
         out_dir, n_trials, seed)
    print ('-> {}'.format (out_file))
    np.save (out_file, trials.as_array)

@cli.command()
@click.option('--n-trials', default=100, type=int)
@click.option ('-n', '--n-sig', default=0, type=float)
@click.option ('--poisson/--nopoisson', default=True)
@click.option ('--fix_gamma', default=None)
@click.option ('--src_gamma', default=2.0, type=float)
@click.option ('--thresh', default=0.0, type=float)
@click.option ('--lag', default=0.0, type=float)
@click.option ('--weight', default='equal', type=str)
@click.option ('--seed', default=None, type=int)
@click.option ('--cpus', default=1, type=int)
@click.option ('--template_str', default="pi0", type=str)
@pass_state
def inj_gp_trials_stacking ( state,  n_trials, fix_gamma, src_gamma, thresh, lag, weight,
         n_sig, poisson, seed,  cpus, template_str, logging=True):
    ana = state.ana
    if seed is None:
        seed = int (time.time () % 2**32)
    random = cy.utils.get_random (seed)
    cutoff = cutoff_GeV =  np.inf
    stacking_conf, stacking_inj_conf  = cg.get_stacking_config(    
                                ana,
                                src_gamma, 
                                fix_gamma, 
                                thresh, 
                                lag,
                                weight
                                )
    gp_conf = cg.get_gp_conf(
                    ana,
                    template_str,
                    )

    tr = cy.get_trial_runner(
                        #conf = stacking_conf,
                        conf = gp_conf, 
                        ana=ana, 
                        mp_cpus=cpus, 
                        dir=dir
                        )
    print(gp_conf)
    print(stacking_conf)
    trials = tr.get_many_fits(n_trials, n_sig, logging=logging)
    t0 = now ()
    print ('Beginning trials at {} ...'.format (t0))
    flush ()
    if n_sig:
        if fix_gamma:
            out_dir = cy.utils.ensure_dir (
                    '{}/{}/lc/gp_inj/{}/stacking/trials/fix_gamma_{}/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/weight_{}/nsig/{}'.format (     
                        state.base_dir, state.ana_name, template_str, fix_gamma, src_gamma, thresh, lag, cutoff,
                          n_sig))
    else:
            out_dir = cy.utils.ensure_dir (
                    '{}/{}/lc/gp_inj/{}/stacking/trials/fix_gamma_{}/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/weight_{}/bg'.format (     
                        state.base_dir, state.ana_name, template_str, fix_gamma, src_gamma, thresh, lag, cutoff,
                          ))
 
    out_file = '{}/trials_{:07d}__seed_{:010d}.npy'.format (
         out_dir, n_trials, seed)
    print ('-> {}'.format (out_file))
    np.save (out_file, trials.as_array)



@cli.command ()
@click.option ('--n-trials', default=10000, type=int)
@click.option ('--n-jobs', default=10, type=int)
@click.option ('-n', '--n-sig', 'n_sigs', multiple=True, default=[0], type=float)
@click.option ('--src_gamma',  default=2.0, type=float)
@click.option ('--fix_gamma',  default=None, type=float)
@click.option ('--inj_name', default='XTE_J1710_dash_281', type=str)
@click.option ('--poisson/--nopoisson', default=True)
@click.option ('--cutoff', default=np.inf, type=float, help='exponential cutoff energy, (TeV)')
@click.option ('--dry/--nodry', default=False)
@click.option ('--seed', default=0)
@pass_state
def submit_do_inj_trials (
        state, inj_name, n_trials, n_jobs, n_sigs, src_gamma, fix_gamma, poisson, cutoff, dry, seed):
    ana_name = state.ana_name
    T = time.time ()
    poisson_str = 'poisson' if poisson else 'nopoisson'
    job_dir = '{}/{}/lc_trials/T_{:17.6f}'.format (
        job_basedir, ana_name,  T)
    sub = Submitter (job_dir=job_dir, memory=9, max_jobs=1000, config = submit_cfg_file)
    #env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'
    commands, labels = [], []
    this_script = os.path.abspath (__file__)
    lags = [0] 
    threshs = [0]
    src_gammas = [src_gamma]
    fix_gammas = [fix_gamma] 
    
    sources = pd.read_hdf(cg.source_file)
    names = sources.name_disp
    for test_name in names:
        for src_gamma in src_gammas:
            for thresh in threshs:
                for lag in lags:
                    for n_sig in n_sigs:
                        for i in range (n_jobs):
                            s = i + seed
                            #fmt = '{} --ana-dir \'\' --base-dir={} do-ps-trials --dec={:+08.3f} --n-trials={}'
                            fmt = '{} {} make-inj-trials --inj_name={} --test_name={}  --src_gamma={} --n-trials={}' \
                                    ' --n-sig={} --thresh={} --lag={}' \
                                    ' --{} --seed={}'
                            command = fmt.format ( this_script, state.state_args, inj_name, test_name, src_gamma, n_trials,
                                                  n_sig, thresh, lag, poisson_str, s)
                            fmt = 'xrb{}_{}_fit_gamma__trials_{:07d}__n_sig_{:08.3f}' \
                                    'src_gamma_{:.3f}_thresh_{}_lag_{}_seed_{:04d}'
                            label = fmt.format (inj_name, test_name,  n_trials, n_sig, src_gamma,
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
def submit_do_AB_test_trials (
        state, n_trials, n_jobs, n_sigs, src_gamma, fix_gamma, poisson, cutoff, dry, seed):
    ana_name = state.ana_name
    T = time.time ()
    poisson_str = 'poisson' if poisson else 'nopoisson'
    job_dir = '{}/{}/AB_test_lc_trials/T_{:17.6f}'.format (
        job_basedir, ana_name,  T)
    sub = Submitter (job_dir=job_dir, memory=9, max_jobs=1000, config = submit_cfg_file)
    #env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'
    commands, labels = [], []
    this_script = os.path.abspath (__file__)
    lags = [0] 
    threshs = [0]
    src_gammas = [src_gamma]
    fix_gammas = [fix_gamma] 
    
    for src_gamma in src_gammas:
        for thresh in threshs:
            for lag in lags:
                for n_sig in n_sigs:
                    for i in range (n_jobs):
                        s = i + seed
                        fmt = '{} {} make-inj-trials-stacking  --src_gamma={} --n-trials={}' \
                                ' --n-sig={} --thresh={} --lag={}' \
                                ' --{} --seed={}'
                        command = fmt.format ( this_script, state.state_args, src_gamma, n_trials,
                                              n_sig, thresh, lag, poisson_str, s)
                        fmt = 'xrb_fit_gamma__trials_{:07d}__n_sig_{:08.3f}' \
                                'src_gamma_{:.3f}_thresh_{}_lag_{}_seed_{:04d}'
                        label = fmt.format (n_trials, n_sig, src_gamma,
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
@click.option ('--inj_name', default='XTE_J1710_dash_281', type=str)
@click.option ('--fit/--nofit', default=True)
@click.option ('--hist/--nohist', default=False)
@pass_state
def collect_inj_trials(state, inj_name, fit, hist):
    sources = pd.read_hdf(cg.source_file)
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
    #combined_outfile = '{}/{}/lc/TSD{}.dict'.format (
    #    state.base_dir, state.ana_name,  suffix)
    bg = {}
    bgs = {}                                                                    
    key = inj_name 
    
    outfile = '{}/{}/lc/inj_name/{}/TSD{}.dict'.format (
        state.base_dir, state.ana_name, inj_name, suffix)
    print ('\r{} ...'.format (key))
    print('{}/{}/lc/{}'.format(state.base_dir, state.ana_name, key))
    post_convert = (lambda x: TSD (cy.utils.Arrays (x), **kw))
    bgs['{}'.format(key)] = cy.bk.get_all('{}/{}/lc/inj_name/{}'.format(state.base_dir,
                            state.ana_name, key),'trials*.npy', 
                            merge=np.concatenate, post_convert=post_convert)
    print(bgs.keys())
    bg['inj_name'] = bgs             
    print ('\rDone.' + 20 * ' ')
    flush ()
    print ('->', outfile)
    with open (outfile, 'wb') as f:
         pickle.dump (bg, f, -1)

@cli.command ()
@click.option ('--inj_name_1', default='AX_J1749_dot_1_dash_2639', type=str)
@click.option ('--inj_name_2', default='Swift_J1745_dot_1_dash_2624', type=str)
@click.option ('--test_name', default='Swift_J1745_dot_1_dash_2624', type=str)
@click.option ('--fit/--nofit', default=True)
@click.option ('--hist/--nohist', default=False)
@pass_state
def collect_stacking_trials(state, inj_name_1, inj_name_2, test_name, fit, hist):
    sources = pd.read_hdf(cg.source_file)
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
    #combined_outfile = '{}/{}/lc/TSD{}.dict'.format (
    #    state.base_dir, state.ana_name,  suffix)
    bg = {}
    bgs = {}                                                                    
    key = test_name
    outfile = '{}/{}/lc/inj_name/{}_{}/test_name/{}/TSD{}.dict'.format (
        state.base_dir, state.ana_name, inj_name_1, inj_name_2, test_name, suffix)
    post_convert = (lambda x: TSD (cy.utils.Arrays (x), **kw))
    bgs['{}'.format(key)] = cy.bk.get_all('{}/{}/lc/inj_name/{}_{}/test_name/{}/'.format(state.base_dir,
                            state.ana_name,inj_name_1, inj_name_2, test_name),'trials*.npy', 
                            merge=np.concatenate, post_convert=post_convert)
    print(bgs.keys())
    bg['test_name'] = bgs             
    print ('\rDone.' + 20 * ' ')
    flush ()
    print ('->', outfile)
    with open (outfile, 'wb') as f:
         pickle.dump (bg, f, -1)

@cli.command()
@click.option ('--name', default=None)
@click.option ('--fix_gamma', default=None)
@click.option ('--src_gamma', default=2.0, type=float)
@click.option ('--nsigma', default=None, type=float)
@click.option ('--thresh', default=0.0, type=float)
@click.option ('--lag', default=0.0, type=float)
@click.option ('--cutoff', default=np.inf, type=float, help='exponential cutoff energy, (TeV)')
@click.option ('--cpus', default=1, type=int)
@click.option ('--save/--nosave', default=False)
@pass_state
def find_lc_nsig(
    state, name, fix_gamma, src_gamma, nsigma, thresh, lag, cutoff, cpus,
    logging=True, save= False):
    ana = state.ana
    cutoff_GeV = 1e3 * cutoff
    sigfile = '{}/{}/lc/TSD_chi2.dict'.format (state.base_dir, state.ana_name)

    if not name:
        sens = {}
        sources = pd.read_hdf(cg.source_file)
        for name in sources.name_disp:
            sig = np.load (sigfile, allow_pickle=True)
            print(name)
            print(sig['name'][name])
            if fix_gamma:
                if sig['name'][name]:
                    sig_trials = cy.bk.get_best(sig, 'name', name, 'fix_gamma_{}'.format(fix_gamma), 'src_gamma_{}'.format(src_gamma),
                                            'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff), 'nsig')
                    try:
                        b = cy.bk.get_best(sig, 'name', name, 'fix_gamma_{}'.format(fix_gamma), 'src_gamma_{}'.format(src_gamma),
                                            'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff), 'bg')
                    except:
                        print( 0,0)
                else:
                    print( 0, 0)
            else:
                print('fit gamma')
                try: 
                    sig_trials = cy.bk.get_best(sig, 'name', name, 'fit_gamma', 'src_gamma_{}'.format(src_gamma), 
                                'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff), 'nsig')
                    b = cy.bk.get_best(sig, 'name', name, 'fit_gamma',  'src_gamma_{}'.format(src_gamma),
                                    'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff_GeV), 'bg')
                    if logging:
                        print(b)
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
                    if nsigma !=None:
                        beta = 0.5
                        print('sigma = {}'.format(nsigma))
                        ts = b.isf_nsigma(nsigma)
                    else:
                        print('Getting sensitivity')
                        beta = 0.9
                        ts = b.median()
                    #print(ts)

                    # include background trials in calculation
                    trials = {0: b}
                    trials.update(sig_trials)
                    for key in trials.keys():
                        trials[key] = trials[key].trials
                    # get number of signal events
                    # (arguments prevent additional trials from being run)
                    
                    result = tr.find_n_sig(ts, beta, max_batch_size=0, logging=logging, trials=trials)
                    flux = tr.to_E2dNdE(result['n_sig'], E0=1, unit=1e3)
                    #flux = tr.to_dNdE(result['n_sig'], E0=1, unit=1e3)
                    # return flux
                    if logging:
                        print(ts, beta, result['n_sig'], flux)
                    sens[name] = flux
                    result['flux_E2dNdE_1TeV'] = flux
                    if save:
                        if nsigma:
                            np.save(f'{state.base_dir}/{state.ana_name}/lc/{name}/trials/fit_gamma/src_gamma_{src_gamma}/thresh_{thresh}/lag_{lag}/{nsigma}sigma_dp.npy', result)
                        else:
                            np.save(f'{state.base_dir}/{state.ana_name}/lc/{name}/trials/fit_gamma/src_gamma_{src_gamma}/thresh_{thresh}/lag_{lag}/sens.npy', result)
                except:
                    pass 
        if save:
            np.save(f'{state.base_dir}/{state.ana_name}/lc/sens_lag_{lag}_thresh_{thresh}.npy', sens)

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
def find_stacking_nsig(
    state, fix_gamma, src_gamma, nsigma, thresh, lag, cutoff, cpus,
    logging=True, save= False):
    ana = state.ana
    cutoff_GeV = 1e3 * cutoff
    inj_1_name = 'AX_J1749_dot_1_dash_2639'
    inj_2_name = 'Swift_J1745_dot_1_dash_2624'
    test_name = 'Swift_J1745_dot_1_dash_2624'
    sigfile = '{}/{}/lc/inj_name/{}_{}/test_name/{}/TSD_chi2.dict'.format (state.base_dir,
        state.ana_name, inj_1_name, inj_2_name, test_name)

    sig = np.load (sigfile, allow_pickle=True)
    print('fit gamma')
    try: 
        sig_trials = cy.bk.get_best(sig, 'test_name', test_name, 'inj_trials', 'fit_gamma', 'src_gamma_{}'.format(src_gamma), 
                    'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff), 'nsig')
        b = cy.bk.get_best(sig, 'test_name', test_name, 'inj_trials', 'fit_gamma', 'src_gamma_{}'.format(src_gamma), 
                        'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff_GeV), 'bg')
        if logging:
            print(b)
        conf, inj_conf  =  cg.get_ps_config(
                                    ana,
                                    test_name, 
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
        if nsigma !=None:
            beta = 0.5
            print('sigma = {}'.format(nsigma))
            ts = b.isf_nsigma(nsigma)
        else:
            print('Getting sensitivity')
            beta = 0.9
            ts = b.median()
        #print(ts)

        # include background trials in calculation
        trials = {0: b}
        trials.update(sig_trials)
        for key in trials.keys():
            trials[key] = trials[key].trials
        # get number of signal events
        # (arguments prevent additional trials from being run)
        
        result = tr.find_n_sig(ts, beta, max_batch_size=0, logging=logging, trials=trials)
        flux = tr.to_E2dNdE(result['n_sig'], E0=1, unit=1e3)
        #flux = tr.to_dNdE(result['n_sig'], E0=1, unit=1e3)
        # return flux
        if logging:
            print(ts, beta, result['n_sig'], flux)
        sens[name] = flux
        result['flux_E2dNdE_1TeV'] = flux
    except:
        pass 
    if save:
        np.save(f'{state.base_dir}/{state.ana_name}/inj_name/{inj_name_1}_{inj_name_2}/test_name/{test_name}/sens_lag_{lag}_thresh_{thresh}.npy', sens)


                

if __name__ == '__main__':
    exe_t0 = now ()
    print ('c7: start at {} .'.format (exe_t0))
    cli ()
