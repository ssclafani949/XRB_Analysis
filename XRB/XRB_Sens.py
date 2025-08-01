#!/usr/bin/env python
# This script is used for sensitivity analysis in XRB (X-ray binaries) using the csky framework.

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
overlap_file = cg.overlap_file
template_file = cg.template_file

class State (object): 
    """
    A class representing the state of an analysis with universal parameters.

    This class encapsulates the configuration and methods necessary to manage 
    different types of analyses, including the ability to save results and 
    apply specific techniques such as downsampling.

    Attributes:
    - ana_name (str): Name of the analysis.
    - ana_dir (str): Directory where analysis data is stored.
    - save (bool): Option to save analysis objects.
    - base_dir (str): Base directory for the project.
    - gp_downsample (bool): Option to enable the Galactic Plane "downsample" technique.
    - template_file (str): Path to the template file for the analysis (necessary for Downsample PDFs).
    """

    def __init__ (self, ana_name, ana_dir, save, base_dir, gp_downsample):
        """
        Initializes the State object with the given parameters.

        Parameters:
        - ana_name (str): Name of the analysis.
        - ana_dir (str): Directory where analysis data is stored.
        - save (bool): Option to save analysis objects.
        - base_dir (str): Base directory for the project.
        - gp_downsample (bool): Option to enable the Galactic Plane "downsample" technique.
        """
        self.ana_name, self.ana_dir, self.save = ana_name, ana_dir, save
        self.base_dir = base_dir
        self.gp_downsample = gp_downsample 
        self.template_file = template_file

    @property
    def ana (self):
        """
        Retrieves the analysis object based on the specified analysis name.

        Depending on the value of `ana_name`, this method fetches the appropriate 
        analysis data specifications, processes the data, and optionally saves 
        the results to the specified directory.

        Returns:
        - The analysis object corresponding to the specified analysis name.
        """
        print(self.ana_name)
        if self.ana_name == 'combo':
            print(repo.local_root)
            cspec = cy.selections.DNNCascadeDataSpecs.DNNC_12yr
            psspec = cy.selections.PSDataSpecs.ps_v4_15yr[3:]
            ana = cy.get_analysis(repo, 'version-004-p03', psspec, 'version-001-p02', cspec)
            ind_list = np.load(overlap_file)
            ana[0].data = cy.utils.Arrays(
                ana[0].data.as_dataframe.drop(
                    labels=ind_list, 
                    errors='ignore',  # Ignore events that are not present since we are just using IC86
                ))
            if self.save:
                cy.utils.ensure_dir(self.ana_dir)
                ana.save(self.ana_dir)
            ana.name = self.ana_name
            self._ana = ana
        elif self.ana_name == 'dnnc':
            cspec = cy.selections.DNNCascadeDataSpecs.DNNC_10yr
            ana = cy.get_analysis(repo, 'version-001-p01', cspec)
            if self.save:
                cy.utils.ensure_dir(self.ana_dir)
                ana.save(self.ana_dir)
            ana.name = self.ana_name
            self._ana = ana
        elif self.ana_name == 'gfu':
            spec = cy.selections.GFUDataSpecs.gfu_IC86_2018
            ana = cy.get_analysis(repo, 'version-002-p05', spec)
            if self.save:
                cy.utils.ensure_dir(self.ana_dir)
                ana.save(self.ana_dir)
            ana.name = self.ana_name
            self._ana = ana
        elif self.ana_name == 'tracks':
            psspec = cy.selections.PSDataSpecs.ps_v4
            ana = cy.get_analysis(repo, 'version-004-p02', psspec)
            if self.save:
                cy.utils.ensure_dir(self.ana_dir)
                ana.save(self.ana_dir)
            ana.name = self.ana_name
            self._ana = ana
        else:
            print('Ana Name {} not valid'.format(self.ana_name))

        if self.gp_downsample:
            print('Building Downsample PDFs')
            for a in ana:
                cy.utils.get_gp_prob(a, 
                    norm=2.18e-18, 
                    gamma=2.7, 
                    template_file=self.template_file,
                    bins=[np.linspace(2,8,50), np.linspace(-1,1,50)])

        return self._ana

    @property
    def state_args (self):
        """
        Constructs command-line arguments for the analysis state.

        Returns:
        - str: A formatted string containing the command-line arguments 
               for the analysis state.
        """
        return '--ana {} --ana-dir {} --base-dir {} --gp_downsample {}'.format (
            self.ana_name, self.ana_dir, self.base_dir, self.gp_downsample)

pass_state = click.make_pass_decorator (State)

@click.group (invoke_without_command=True, chain=True)
@click.option ('-a', '--ana', 'ana_name', default='combo', help='Dataset title')
@click.option ('--ana-dir', default=ana_dir, type=click.Path ())
@click.option ('--save/--nosave', default=True)
@click.option ('--base-dir', default=base_dir,
               type=click.Path (file_okay=False, writable=True))
@click.option ('--gp_downsample/--nogp_downsample', default=True)
@click.pass_context
def cli (ctx, ana_name, ana_dir, save, base_dir, gp_downsample):
    ctx.obj = State.state = State (ana_name, ana_dir, save, base_dir, gp_downsample)


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
@click.option ('--save_trials/--nosave_trials', default=False)
@click.option ('--inject_gp/--noinject_gp',  default=False)
@pass_state
def do_lc_trials(state, name, n_trials, fix_gamma, src_gamma, thresh, lag, n_sig, poisson, seed,
                  cutoff, cpus, save_trials, inject_gp, logging=True):
    """
    Perform individual trials for for XRB sources with lightcurve pdf.

    Parameters:
    - state: analysis state variables 
    - name: The name of the source.
    - n_trials: The number of trials to run.
    - fix_gamma: Fixed gamma value if provided (False).
    - src_gamma: Source gamma value.
    - thresh: Threshold value for analysis.
    - lag: Lag time for analysis.
    - n_sig: Number of signal events to inject.
    - poisson: Boolean flag if injecting via Poisson statistics.
    - seed: Random seed for reproducibility.
    - cutoff: Exponential cutoff energy in TeV.
    - cpus: Number of CPUs to use for parallel processing.
    - save_trials: Boolean flag to save trial results.
    - inject_gp: Boolean flag to enable Galactic Plane in background.
    - logging: Boolean flag to enable logging.
    """
    
    ana = state.ana
    if inject_gp:
        gp_inj_str = 'gp_inj'
    else:
        gp_inj_str = 'no_gpinj'
    
    # Set the random seed if not provided
    if seed is None:
        seed = int(time.time() % 2**32)
    
    random = cy.utils.get_random(seed)
   
    cutoff_GeV = cutoff * 1e3 
    print(seed)
    
    # Get the configuration for the analysis
    conf, inj_conf = cg.get_ps_config(
        ana,
        name, 
        src_gamma, 
        fix_gamma, 
        cutoff_GeV, 
        lag, 
        thresh,
        inject_gp=inject_gp
    ) 
    
    # Initialize the trial runner
    tr = cy.get_trial_runner(
        conf=conf,
        inj_conf=inj_conf, 
        ana=ana, 
        flux=cy.hyp.PowerLawFlux(src_gamma, energy_cutoff=cutoff_GeV), 
        mp_cpus=cpus, 
        dir=dir
    )

    # Start timing the trials
    t0 = now()
    print('Beginning trials at {} ...'.format(t0))
    flush()
    
    # Run the trials and get the results
    trials = tr.get_many_fits(
        n_trials, n_sig=n_sig, poisson=poisson, seed=seed, logging=logging
    )
    
    # End timing the trials
    t1 = now()
    print('Finished trials at {} ...'.format(t1))
    print(trials if n_sig else cy.dists.Chi2TSD(trials))
    print(t1 - t0, 'elapsed.')
    flush()
    
    # Determine the output directory based on parameters
    if n_sig:
        if fix_gamma:
            out_dir = cy.utils.ensure_dir(
                '{}/{}/lc/{}/trials/fix_gamma_{}/{}/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/nsig/{}'.format(
                    state.base_dir, state.ana_name, name, fix_gamma, gp_inj_str, src_gamma, thresh, lag, cutoff,
                    n_sig))
        else:
            out_dir = cy.utils.ensure_dir(
                '{}/{}/lc/{}/trials/fit_gamma/{}/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/nsig/{}'.format(
                    state.base_dir, state.ana_name, name, gp_inj_str, src_gamma, thresh, lag, cutoff,
                    n_sig))
    else:        
        if fix_gamma:
            out_dir = cy.utils.ensure_dir(
                '{}/{}/lc/{}/trials/fix_gamma_{}/{}/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/bg'.format(
                    state.base_dir, state.ana_name, name, fix_gamma, gp_inj_str, src_gamma, thresh, lag, cutoff))
        else:
            out_dir = cy.utils.ensure_dir(
                '{}/{}/lc/{}/trials/fit_gamma/{}/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/bg'.format(
                    state.base_dir, state.ana_name, name, gp_inj_str, src_gamma, thresh, lag, cutoff))

    # Save the trial results to a file
    out_file = '{}/trials_{:07d}__seed_{:010d}.npy'.format(
        out_dir, n_trials, seed)
    print('-> {}'.format(out_file))
    np.save(out_file, trials.as_array)


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
@click.option ('--weight',  default='equal')
@click.option ('--inject_gp/--noinject_gp',  default=False)
@pass_state
def do_stacking_trials(state, n_trials, fix_gamma, src_gamma, thresh, lag, n_sig, poisson, seed, cutoff, cpus, weight, inject_gp, logging=True):
    """
    Perform stacking trials for XRB analysis.

    Parameters:
    - state: The current state of the analysis.
    - n_trials (int): Number of trials to perform (default: 1000).
    - fix_gamma (float or None): Fixed gamma value for the analysis (default: None).
    - src_gamma (float): Source gamma value (default: 2.0).
    - thresh (float): Threshold value for the analysis (default: 0.0).
    - lag (float): Lag time for the analysis (default: 0.0).
    - n_sig (float): Number of signal events to inject (default: 0).
    - poisson (bool): Whether to use Poisson statistics (default: True).
    - seed (int or None): Random seed for reproducibility (default: None).
    - cutoff (float): Exponential cutoff energy in TeV (default: np.inf).
    - cpus (int): Number of CPUs to use for parallel processing (default: 1).
    - weight (str): Weighting scheme for the trials (default: 'equal').
    - inject_gp (bool): Whether to inject Gaussian processes (default: False).
    - logging (bool): Whether to enable logging (default: True).

    This function initializes the analysis, sets up the trial runner, and performs the stacking trials,
    saving the results to a specified directory based on the parameters provided.
    """
    
    # Initialize the analysis state
    ana = state.ana
    if seed is None:
        seed = int(time.time() % 2**32)
    
    random = cy.utils.get_random(seed)
    
    # Determine the injection string based on the inject_gp option
    gp_inj_str = 'gp_inj' if inject_gp else 'no_gpinj'
    
    # Convert cutoff energy from TeV to GeV
    cutoff_GeV = cutoff * 1e3 

    # Get stacking configuration
    conf, inj_conf = cg.get_stacking_config(
        ana,
        src_gamma, 
        fix_gamma, 
        thresh, 
        lag,
        weight,
        inject_gp=inject_gp 
    )
    
    # Initialize the trial runner
    tr = cy.get_trial_runner(
        conf=conf,  
        inj_conf=inj_conf, 
        ana=ana,
        flux=cy.hyp.PowerLawFlux(src_gamma, energy_cutoff=cutoff_GeV), 
        dir=dir
    )
    
    # Start timing the trials
    t0 = now()
    print('Beginning trials at {} ...'.format(t0))
    flush()
    
    # Perform the trials
    trials = tr.get_many_fits(
        n_trials, n_sig=n_sig, poisson=poisson, seed=seed, logging=logging, mp_cpus=cpus
    )
    
    # Output results
    print(f'TS: {trials.ts}')
    print(f'ns: {trials.ns}')
    if 'gamma' in trials.keys():
        print(f'Gamma: {trials.gamma}')
    
    # End timing the trials
    t1 = now()
    print('Finished trials at {} ...'.format(t1))
    print(trials if n_sig else cy.dists.Chi2TSD(trials))
    print(np.median(trials['ns']))
    print(t1 - t0, 'elapsed.')
    flush()
    
    # Determine output directory based on parameters
    if n_sig:
        if fix_gamma:
            out_dir = cy.utils.ensure_dir(
                '{}/{}/lc/stacking/trials/fix_gamma_{}/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/weight_{}/{}/nsig/{}/'.format(
                    state.base_dir, state.ana_name, fix_gamma, src_gamma, thresh, lag, cutoff, gp_inj_str, n_sig
                )
            )
        else:
            out_dir = cy.utils.ensure_dir(
                '{}/{}/lc/stacking/trials/fit_gamma/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/weight_{}/{}/nsig/{}'.format(
                    state.base_dir, state.ana_name, src_gamma, thresh, lag, cutoff, weight, gp_inj_str, n_sig
                )
            )
    else:        
        if fix_gamma:
            out_dir = cy.utils.ensure_dir(
                '{}/{}/lc/stacking//trials/fix_gamma_{}/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/weight_{}/{}{}//bg'.format(
                    state.base_dir, state.ana_name, fix_gamma, src_gamma, thresh, lag, cutoff, weight, gp_inj_str
                )
            )
        else:
            out_dir = cy.utils.ensure_dir(
                '{}/{}/lc/stacking/trials/fit_gamma/src_gamma_{}/thresh_{}/lag_{}/cutoff_{}/weight_{}/{}/bg'.format(
                    state.base_dir, state.ana_name, src_gamma, thresh, lag, cutoff, weight, gp_inj_str
                )
            )
    
    # Save the results if specified
    if state.save:
        out_file = '{}/trials_{:07d}__seed_{:010d}.npy'.format(
            out_dir, n_trials, seed
        )
        print('-> {}'.format(out_file))
        np.save(out_file, trials.as_array)


@cli.command()
@click.option('--fit/--nofit', default=True, help="Enable fitting of the light curves.")
@click.option('--hist/--nohist', default=False, help="Enable histogram generation.")
@click.option('-n', '--n', default=0, type=int, help="Minimum number of trials required.")
@pass_state
def collect_lc_trials(state, fit, hist, n):
    """
    Collect light curve trials for sources defined in the source file.

    This function reads source data from an HDF file, processes light curve trials
    based on the specified options (fit, hist), and saves the results to a dictionary
    file. The trials can be processed as either binned histograms or chi-squared fits.

    Parameters:
    state (State): The state object containing base directory and analysis name.
    fit (bool): Flag to indicate if fitting should be performed.
    hist (bool): Flag to indicate if histogram generation should be performed.
    n (int): Minimum number of trials required for processing.

    Returns:
    None
    """
    
    # Read source data from the specified HDF file
    sources = pd.read_hdf(cg.source_file)
    names = sources.name_disp
    kw = {}
    
    # Determine the type of TSD to use based on the options provided
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
    
    # Define the output file path for combined results
    combined_outfile = '{}/{}/lc/TSD{}.dict'.format(state.base_dir, state.ana_name, suffix)
    bg = {}
    bgs = {}
    
    # Iterate over each source name to collect trials
    for name in names:
        outfile = '{}/{}/lc/{}/TSD{}.dict'.format(state.base_dir, state.ana_name, name, suffix)
        key = name 
        print('\r{} ...'.format(key))
        
        # Print the current processing path
        print('{}/{}/lc/{}'.format(state.base_dir, state.ana_name, key))
        
        # Define a post-conversion function for the trials
        post_convert = (lambda x: TSD(cy.utils.Arrays(x), **kw))
        
        # Collect all trials for the current source
        bgs['{}'.format(key)] = cy.bk.get_all('{}/{}/lc/{}/trials/'.format(state.base_dir, state.ana_name, key), 'trials*.npy', 
                          merge=np.concatenate, post_convert=post_convert)
    
    # Print the keys of the collected background trials
    print(bgs.keys())
    bg['name'] = bgs             
    print('\rDone.' + 20 * ' ')
    flush()
    
    # Print the path of the combined output file
    print('->', combined_outfile)
    
    # Save the collected background trials to the combined output file
    with open(combined_outfile, 'wb') as f:
         pickle.dump(bg, f, -1)

@cli.command ()
@click.option ('--fit/--nofit', default=True)
@click.option ('--hist/--nohist', default=False)
@click.option ('-n', '--n', default=0, type=int)
@pass_state
def collect_stacking_trials (state, fit, hist, n):
"""
    Collect stacking trials for sources defined in the source file.

    This function reads source data from an HDF file, processes light curve trials
    based on the specified options (fit, hist), and saves the results to a dictionary
    file. The trials can be processed as either binned histograms or chi-squared fits.

    Parameters:
    state (State): The state object containing base directory and analysis name.
    fit (bool): Flag to indicate if fitting should be performed.
    hist (bool): Flag to indicate if histogram generation should be performed.
    n (int): Minimum number of trials required for processing.

    Returns:
    None
    """
    kw = {}
    # Determine the type of TSD based on the user's choice of histogram or fitting
    if hist:
        TSD = cy.dists.BinnedTSD  # Use BinnedTSD for histogram generation
        suffix = '_hist'  # Suffix for output file
        kw['keep_trials'] = False  # Do not keep trials for histogram
    elif fit:
        TSD = cy.dists.Chi2TSD  # Use Chi2TSD for fitting
        suffix = '_chi2'  # Suffix for output file
    else:
        TSD = cy.dists.TSD  # Default TSD
        suffix = ''  # No suffix for output file

    # Define the output file path for saving results
    outfile = '{}/{}/lc/stacking/TSD{}.dict'.format(
        state.base_dir, state.ana_name, suffix)
    
    bg = {}  # Background dictionary to store results
    bgs = {}  # Additional background dictionary (currently unused)
    
    print('{}/{}/lc/'.format(state.base_dir, state.ana_name))  # Print the directory path
    
    # Function to convert data using the selected TSD
    post_convert = (lambda x: TSD(cy.utils.Arrays(x), **kw))
    
    # Retrieve all trials and apply the post-conversion function
    bg = cy.bk.get_all('{}/{}/lc/stacking/trials/'.format(state.base_dir, state.ana_name), 'trials*.npy', 
                       merge=np.concatenate, post_convert=post_convert)
    
    print('\rDone.' + 20 * ' ')  # Indicate completion
    flush()  # Flush the output buffer
    print('->', outfile)  # Print the output file path
    
    # Save the background results to a file
    with open(outfile, 'wb') as f:
        pickle.dump(bg, f, -1)  # Use pickle to serialize the background data


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
@click.option ('--inject_gp/--noinject_gp',  default=False)
@pass_state
def submit_do_lc_trials (
        state, n_trials, n_jobs, n_sigs, src_gamma, fix_gamma, poisson, cutoff, dry, seed, inject_gp):
    ana_name = state.ana_name
    T = time.time ()
    downsample_str = 'gp_downsample' if state.gp_downsample else 'nogp_downsample'
    inj_gp_str = 'inject_gp' if inject_gp else 'noinject_gp'
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
    names = sources.name_disp[30:]
    #names = ['Swift_J1745_dot_1_dash_2624']
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
                                            ' --{} --{} --seed={} --save'
                                    command = fmt.format ( this_script, downsample_str, name, fix_gamma, src_gamma, cutoff, n_trials,
                                                          n_sig, thresh, lag, poisson_str, inj_gp_str, s)
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
                                fmt = '{} --{} do-lc-trials --name {}  --src_gamma={} --n-trials={}' \
                                        ' --n-sig={} --thresh={} --lag={}' \
                                        ' --{} --{} --seed={}'
                                #command = fmt.format (env_shell, this_script, state.state_args, name, src_gamma, n_trials,
                                #                      n_sig, thresh, lag, poisson_str, s)
                                command = fmt.format ( this_script, downsample_str,  name, src_gamma, n_trials,
                                                      n_sig, thresh, lag, poisson_str, inj_gp_str, s)
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
@click.option ('--weight', default='equal')
@click.option ('--inject_gp/--noinject_gp',  default=False)
@pass_state
def submit_do_stacking_trials (
        state, n_trials, n_jobs, n_sigs, src_gamma, fix_gamma, poisson, cutoff, dry, seed, weight, inject_gp):
    ana_name = state.ana_name
    gp_downsample = state.gp_downsample
    inj_gp = inject_gp
    downsample_str = 'gp_downsample' if state.gp_downsample else 'nogp_downsample'
    inj_gp_str = 'inject_gp' if inject_gp else 'noinject_gp'
    T = time.time ()
    #job_basedir = '/scratch/ssclafani/logs/' 
    poisson_str = 'poisson' if poisson else 'nopoisson'
    job_dir = '{}/{}/lc_trials/T_{:17.6f}'.format (
        job_basedir, ana_name,  T)
    sub = Submitter (job_dir=job_dir, ncpu=1, memory=13, max_jobs=1000, config = submit_cfg_file)
    #env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'
    commands, labels = [], []
    this_script = os.path.abspath (__file__)
    if fix_gamma:
        for n_sig in n_sigs:
            for i in range (n_jobs):
                s = i + seed
                fmt = ' {} {} do-stacking-trials --fix_gamma={} --src_gamma={} --cutoff {} --n-trials={}' \
                        ' --n-sig={} --weight={}' \
                        ' --{} --seed={}'
                command = fmt.format ( this_script, state.state_args, fix_gamma, src_gamma, cutoff, n_trials,
                                      n_sig, weight, poisson_str, s)
                fmt = 'xrb_fix_gamma_{}__trials_{:07d}__n_sig_{:08.3f}__' \
                        'src_gamma_{:.3f}__{}_cutoff_{}seed_{:04d}'
                label = fmt.format ( fix_gamma, n_trials, n_sig, src_gamma, weight,
                                    cutoff,  s)
                commands.append (command)
                labels.append (label)
    else:
        for n_sig in n_sigs:
            for i in range (n_jobs):
                s = i + seed
                fmt = '{} --{} do-stacking-trials  --src_gamma={} --n-trials={}' \
                        ' --n-sig={} --weight={}' \
                        ' --{}  --{} --seed={}'
                command = fmt.format ( this_script, downsample_str, src_gamma, n_trials,
                                      n_sig, weight, poisson_str, inj_gp_str,  s)
                fmt = 'xrb_stacking_fit_gamma__trials_{:07d}__n_sig_{:08.3f}__{}' \
                        'src_gamma_{:.3f}_seed_{:04d}'
                label = fmt.format (n_trials, n_sig, weight, src_gamma, s)
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
@click.option ('--weight', default = 'equal')
@click.option ('--save/--nosave', default=False)
@pass_state
def find_stacking_nsig(state, fix_gamma, src_gamma, nsigma, thresh, lag, cutoff, cpus, weight, logging=True, save=False):
    """
    Find the number of signal events in a stacking analysis.

    This function performs a stacking analysis to determine the number of signal events 
    based on the provided parameters. It utilizes background trials and calculates the 
    sensitivity of the signal based on the specified thresholds and gamma values.

    Parameters:
    - state: The current state of the analysis.
    - fix_gamma (float or None): The fixed gamma value for the analysis. If None, the function will fit gamma.
    - src_gamma (float): The source gamma value for the analysis.
    - nsigma (float or None): The number of sigma for the sensitivity calculation. If None, sensitivity is calculated.
    - thresh (float): The threshold value for the analysis.
    - lag (float): The lag value for the analysis.
    - cutoff (float): The exponential cutoff energy in TeV.
    - cpus (int): The number of CPUs to use for the analysis.
    - weight (str): The weighting method to use ('equal' or other).
    - logging (bool): Whether to log the output.
    - save (bool): Whether to save the results.

    Returns:
    - None: The function prints the results and saves them if specified.
    """
    ana = state.ana
    cutoff_GeV = 1e3 * cutoff
    sigfile = '{}/{}/lc/stacking/TSD_chi2.dict'.format(state.base_dir, state.ana_name)
    sig = np.load(sigfile, allow_pickle=True)
    print(sig)
    if fix_gamma:
        sig_trials = cy.bk.get_best(sig, 'fix_gamma_{}'.format(fit_gamma), 'src_gamma_{}'.format(src_gamma),
                                     'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff), 'weight_{}'.format(weight), 'nsig')
        b = cy.bk.get_best(sig, 'fix_gamma_{}'.format(fit_gamma), 'src_gamma_{}'.format(src_gamma),
                           'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff), 'weight_{}'.format(weight), 'bg')
    else:
        print('fit gamma')
        sig_trials = cy.bk.get_best(sig, 'fit_gamma', 'src_gamma_{}'.format(src_gamma),
                                     'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff), 'weight_{}'.format(weight), 'no_gpinj', 'nsig')

        b = cy.bk.get_best(sig, 'fit_gamma', 'src_gamma_{}'.format(src_gamma),
                           'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff), 'weight_{}'.format(weight), 'gp_inj', 'bg')
    if logging:
        print(b)
    conf, inj_conf = cg.get_stacking_config(
        ana,
        src_gamma,
        fix_gamma,
        thresh,
        lag,
        weight)
    tr = cy.get_trial_runner(
        conf=conf,
        inj_conf=inj_conf,
        ana=ana,
        flux=cy.hyp.PowerLawFlux(src_gamma, energy_cutoff=cutoff_GeV),
        dir=dir)
    # determine ts threshold
    if nsigma is not None:
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
    result['flux_E2dNdE_1TeV'] = f
    if save:
        cy.utils.ensure_dir(f'{state.base_dir}/{state.ana_name}/lc/stacking/trials/fit_gamma/src_gamma_{src_gamma}/weight/{weight}/')
        if nsigma:
            np.save(f'{state.base_dir}/{state.ana_name}/lc/stacking/trials/fit_gamma/src_gamma_{src_gamma}/weight/{weight}/{nsigma}sigma_dp.npy', result)
        else:
            np.save(f'{state.base_dir}/{state.ana_name}/lc/stacking/trials/fit_gamma/src_gamma_{src_gamma}/weight/{weight}/sens.npy', result)

@cli.command()
@click.option('--name', default=None)
@click.option('--fix_gamma', default=None)
@click.option('--src_gamma', default=2.0, type=float)
@click.option('--nsigma', default=None, type=float)
@click.option('--thresh', default=0.0, type=float)
@click.option('--lag', default=0.0, type=float)
@click.option('--cutoff', default=np.inf, type=float, help='exponential cutoff energy, (TeV)')
@click.option('--cpus', default=1, type=int)
@click.option('--save/--nosave', default=False)
@pass_state
def find_lc_nsig(
    state, name, fix_gamma, src_gamma, nsigma, thresh, lag, cutoff, cpus, logging=True, save=False, sens=True):
    """
    Find the sensitivity of the light curve based on the provided parameters.

    Parameters:
    - state: The current state containing analysis information.
    - name: The name of the source to analyze. If None, all sources will be analyzed.
    - fix_gamma: Fixed gamma value for the analysis.
    - src_gamma: Source gamma value (default is 2.0).
    - nsigma: Number of sigma for sensitivity calculation.
    - thresh: Threshold value for the analysis (default is 0.0).
    - lag: Lag time for the analysis (default is 0.0).
    - cutoff: Exponential cutoff energy in TeV (default is infinity).
    - cpus: Number of CPUs to use for the analysis (default is 1).
    - logging: Boolean flag to enable logging (default is True).
    - save: Boolean flag to save the results (default is False).
    - sens: Boolean flag to indicate if sensitivity should be calculated (default is True).

    This function loads signal and bkg tirals,
    and calculates the sensitivity or discovery potential for the specified sources.
    Results can be saved to disk if the save flag is set to True.
    Iterates over all sources if name is not provided.
    """
    
    ana = state.ana
    cutoff_GeV = 1e3 * cutoff
    sigfile = '{}/{}/lc/TSD_chi2.dict'.format(state.base_dir, state.ana_name)

    if name:
        name_list = [name]
    else:
        sens = {}
        sources = pd.read_hdf(cg.source_file)
        name_list = sources.name_disp
        print(name_list) 

    # Iterate over each source name (just one if name is provided)
    for name in name_list:
        sig = np.load(sigfile, allow_pickle=True)
        print(name)
        print(sig['name'][name])
        
        if fix_gamma:
            if sig['name'][name]:
                sig_trials = cy.bk.get_best(sig, 'name', name, 'fix_gamma_{}'.format(fix_gamma), 'src_gamma_{}'.format(src_gamma),
                                        'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff), 'no_gpinj', 'nsig')
                try:
                    b = cy.bk.get_best(sig, 'name', name, 'fix_gamma_{}'.format(fix_gamma), 'src_gamma_{}'.format(src_gamma),
                                        'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff), 'gp_inj', 'bg')
                except:
                    print(0, 0)
            else:
                print(0, 0)
        else:
            try:
                sig_trials = cy.bk.get_best(sig, 'name', name, 'fit_gamma', 'no_gpinj', 'src_gamma_{}'.format(src_gamma), 
                            'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff), 'nsig')
                b = cy.bk.get_best(sig, 'name', name, 'fit_gamma', 'gp_inj',  'src_gamma_2.0'.format(src_gamma),
                                'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff_GeV),  'bg')
                if logging:
                    print(b)
                
                # Get configuration for sensitivity calculation
                conf, inj_conf = cg.get_ps_config(
                                            ana,
                                            name, 
                                            src_gamma, 
                                            fix_gamma, 
                                            cutoff_GeV, 
                                            lag, 
                                            thresh
                                            ) 
                
                # Initialize trial runner
                tr = cy.get_trial_runner(
                                    conf=conf,
                                    inj_conf=inj_conf, 
                                    ana=ana, 
                                    flux=cy.hyp.PowerLawFlux(src_gamma, energy_cutoff=cutoff_GeV), 
                                    mp_cpus=cpus, 
                                    dir=dir
                                    )
                
                if nsigma is not None:
                    beta = 0.5
                    print('sigma = {}'.format(nsigma))
                    ts = b.isf_nsigma(nsigma)
                else:
                    print('Getting sensitivity')
                    beta = 0.9
                    ts = b.median()

                # Include background trials in calculation
                trials = {0: b}
                trials.update(sig_trials)
                for key in trials.keys():
                    trials[key] = trials[key].trials
                
                # Get number of signal events
                result = tr.find_n_sig(ts, beta, max_batch_size=0, logging=logging, trials=trials)
                
                # Convert to flux
                flux = tr.to_E2dNdE(result['n_sig'], E0=1, unit=1e3)
                
                if logging:
                    print(ts, beta, result['n_sig'], flux)
                
                sens[name] = flux
                result['flux_E2dNdE_1TeV'] = flux
                
                # Save results if the save flag is set
                if save:
                    if nsigma:
                        np.save(f'{state.base_dir}/{state.ana_name}/lc/{name}/trials/fit_gamma/no_gpinj/src_gamma_{src_gamma}/thresh_{thresh}/lag_{lag}/{nsigma}sigma_dp.npy', result)
                    else:
                        print('saving..')
                        np.save(f'{state.base_dir}/{state.ana_name}/lc/{name}/trials/fit_gamma/no_gpinj/src_gamma_{src_gamma}/thresh_{thresh}/lag_{lag}/sens.npy', result)
            except:
                print('skipping...')
    
    # Save total results
    if save and not name:
        if nsigma:
            np.save(f'{state.base_dir}/{state.ana_name}/lc/dp{nsigma}sig_lag_{lag}_thresh_{thresh}.npy', sens)
        else:
            np.save(f'{state.base_dir}/{state.ana_name}/lc/sens_lag_{lag}_thresh_{thresh}.npy', sens)


@cli.command()
@click.option ('--fix_gamma', default=None)
@click.option ('--src_gamma', default=2.0, type=float)
@click.option ('--thresh', default=0.0, type=float)
@click.option ('--lag', default=0.0, type=float)
@click.option ('--cutoff', default=np.inf, type=float, help='exponential cutoff energy, (TeV)')
@click.option ('--cpus', default=1, type=int)
@click.option ('--weight', default = 'equal')
@click.option ('--save/--nosave', default=False)
@click.option ('--sens/--nosens', default=True)                                                                 
@pass_state
def plot_stacking_bias ( state,fix_gamma, src_gamma, thresh, lag, cutoff, cpus, weight, logging=True, save=False, sens=True):
    import matplotlib.pyplot as plt

    cutoff_GeV = 1e3 * cutoff
    sigfile = '{}/{}/lc/stacking/TSD_chi2.dict'.format (state.base_dir, state.ana_name)
    sig = np.load (sigfile, allow_pickle=True)
    if fix_gamma:
        sig_trials = cy.bk.get_best(sig, 'fix_gamma_{}'.format(fit_gamma),  'src_gamma_{}'.format(src_gamma),
                                                'thresh_{}'.format(thres), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff), 'weight_{}'.format(weight), 'nsig')
        b = cy.bk.get_best(sig, 'fix_gamma_{}'.format(fit_gamma),  'src_gamma_{}'.format(src_gamma),
                                                    'thresh_{}'.format(thresh), 'lag_{}'.format(lag),  'cutoff_{}'.format(cutoff), 'weight_{}'.format(weight), 'bg')
    else:
        print('fit gamma')
        sig_trials = cy.bk.get_best(sig, 'fit_gamma',  'src_gamma_{}'.format(src_gamma),
                                                    'thresh_{}'.format(thresh), 'lag_{}'.format(lag),  'cutoff_{}'.format(cutoff), 'weight_{}'.format(weight), 'no_gpinj', 'nsig')
    
        b = cy.bk.get_best(sig,  'fit_gamma',  'src_gamma_{}'.format(src_gamma),
                                                    'thresh_{}'.format(thresh), 'lag_{}'.format(lag),  'cutoff_{}'.format(cutoff), 'weight_{}'.format(weight), 'gp_inj', 'bg')
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
    if not fix_gamma:
        expect_gamma = src_gamma 
        ax = axs[0]
    else:
        ax = axs[0]
    
    dns = np.mean(np.diff(n_sigs))
    ns_bins = np.r_[n_sigs - 0.5*dns, n_sigs[-1] + 0.5*dns]
    expect_kw = dict(color='C0', ls='--', lw=1, zorder=-10)
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
    if sens:
        sens = np.load(f'{state.base_dir}/{state.ana_name}/lc/stacking/trials/fit_gamma/src_gamma_{src_gamma}/weight/{weight}/sens.npy', allow_pickle = True)[()]
        nsigma = 5.0
        dp = np.load(f'{state.base_dir}/{state.ana_name}/lc/stacking/trials/fit_gamma/src_gamma_{src_gamma}/weight/{weight}/{nsigma}sigma_dp.npy', allow_pickle = True)[()]
        try:
            for ax in axs:
                ax.axvline(sens['n_sig'], ls = '--', c='C0')
                ax.axvline(dp['n_sig'], ls = ':', c='C0')
        except(FileNotFoundError):
            pass 
    for ax in axs:
        if src_gamma == 2.0:
            ax.set_xlim(0,50)
            axs[0].set_ylim(0,50)
        else:
            ax.set_xlim(0,100)
            axs[0].set_ylim(0,100)
    plt.suptitle(f'Weight {weight} $\gamma$={src_gamma:.2}')
    plt.tight_layout()
    if save:
        save_dir = cy.utils.ensure_dir(f'{state.base_dir}/{state.ana_name}/lc/stacking/plots/fit_gamma/src_gamma_{src_gamma}/')
        print(f'Saving to ... {save_dir}')
        plt.savefig(f'{save_dir}/bias_g{src_gamma:.2f}_{weight}.png')
                
@cli.command()
@click.option ('--name', default=None)
@click.option ('--fix_gamma', default=None)
@click.option ('--src_gamma', default=2.0, type=float)
@click.option ('--thresh', default=0.0, type=float)
@click.option ('--lag', default=0.0, type=float)
@click.option ('--cutoff', default=np.inf, type=float, help='exponential cutoff energy, (TeV)')
@click.option ('--cpus', default=1, type=int)
@click.option ('--logging/--nologging', default=True)
@click.option ('--save/--nosave', default=False)
@click.option ('--sens/--nosens', default=False)
@pass_state
def plot_lc_bias (
    state, name, fix_gamma, src_gamma,  thresh, lag, cutoff, cpus,
    logging=True, save= False, sens=False):
    import matplotlib.pyplot as plt

    cutoff_GeV = 1e3 * cutoff
    sigfile = '{}/{}/lc/TSD_chi2.dict'.format (state.base_dir, state.ana_name)
    sig = np.load (sigfile, allow_pickle=True)

    def get_bias(name, fix_gamma, src_gamma, thresh, lag, cutoff, cpus, sens):
        if fix_gamma:
            sig_trials = cy.bk.get_best(sig, 'fix_gamma_{}'.format(fit_gamma),  'src_gamma_{}'.format(src_gamma),
                                                    'thresh_{}'.format(thres), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff), 'weight_{}'.format(weight), 'nsig')
            b = cy.bk.get_best(sig, 'fix_gamma_{}'.format(fit_gamma),  'src_gamma_{}'.format(src_gamma),
                                                        'thresh_{}'.format(thresh), 'lag_{}'.format(lag),  'cutoff_{}'.format(cutoff), 'weight_{}'.format(weight), 'bg')
        else:
            print('fit gamma')
 
            sig_trials = cy.bk.get_best(sig, 'name', name, 'fit_gamma', 'src_gamma_{}'.format(src_gamma), 
                        'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff), 'nsig')
            b = cy.bk.get_best(sig, 'name', name, 'fit_gamma',  'src_gamma_{}'.format(src_gamma),
                            'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff_GeV), 'bg')
        if logging:
            print(b)
        trials = {0: b}
        trials.update(sig_trials)
        n_sigs = [n for n in sorted(trials.keys())]
        
        for key in sorted(trials.keys()):
            trials[key] = trials[key].trials
            #trials['ntrue'] = trials[key]
        trials = [trials[k] for k in sorted(trials.keys())]
        print(trials)
        for (n_sig, t) in zip(sorted(n_sigs), trials):
            t['ntrue'] = np.repeat(n_sig, len(t))
        allt = cy.utils.Arrays.concatenate(trials)
            
        if logging:
            for ntrue in sorted(np.unique(allt['ntrue'])):
                mask = allt['ntrue'] == ntrue
                median = np.median(allt[mask].ns)
                gamma = np.median(allt[mask].gamma)
                l = np.median(allt[mask].lag)
                print(f'Inj: {ntrue:3.2f}   | Median: {median:.2f}   | gamma {gamma:.2f} | lag {lag:.2f}')


        fig, axs = plt.subplots(1, 4, figsize=(12,3))

        name_disp = name.replace('_plus_', '+')
        name_disp = name_disp.replace('_dash_', '-')
        name_disp = name_disp.replace('_', ' ')
        
        fig.suptitle(name_disp)
        print(n_sigs)
        print(allt)
        
        dns = np.mean(np.diff(n_sigs))
        ns_bins = np.r_[n_sigs - 0.5*dns, n_sigs[-1] + 0.5*dns]
        print(n_sigs)
        print(ns_bins)
        expect_kw = dict(color='C0', ls='--', lw=1, zorder=-10)
        expect_gamma = src_gamma
        expect_thresh = thresh
        expect_lag = lag
        
        ax = axs[0]
        h = hl.hist((allt.ntrue, allt.ns), bins=(ns_bins, 70))
        hl.plot1d(ax, h.contain_project(1),errorbands=True, drawstyle='default')
        ax.set_xlabel(r'n$_{inj}$')
        ax.set_ylabel(r'n$_s$')
        ax.set_aspect('equal')
        
        lim = ns_bins[[0, -1]]
        ax.set_xlim(ax.set_ylim(lim))
        ax.plot(lim, lim, **expect_kw)
        ax.set_aspect('equal')
        ax.set_xlim(0,50)
        ax.set_ylim(0,50)

        ax = axs[1]
        h = hl.hist((allt.ntrue, allt.gamma), bins=(ns_bins, 70))
        hl.plot1d(ax, h.contain_project(1),errorbands=True, drawstyle='default')
        ax.axhline(expect_gamma, **expect_kw)
        ax.set_xlim(axs[0].get_xlim())
        ax.set_xlabel(r'n$_{inj}$')
        ax.set_ylabel(r'$\gamma$')
        ax.set_xlim(0,50)
        #ax.set_aspect('equal')
        
        ax = axs[2]
        h = hl.hist((allt.ntrue, allt.lag), bins=(ns_bins, 70))
        hl.plot1d(ax, h.contain_project(1),errorbands=True, drawstyle='default')
        ax.axhline(expect_lag, **expect_kw)
        ax.set_xlim(axs[0].get_xlim())
        ax.set_xlabel(r'n$_{inj}$')
        ax.set_ylabel(r'lag')
        ax.set_xlim(0,50)
        #ax.set_aspect('equal')
        
        ax = axs[3]
        h = hl.hist((allt.ntrue, allt.thresh), bins=(ns_bins, 70))
        hl.plot1d(ax, h.contain_project(1),errorbands=True, drawstyle='default')
        ax.axhline(expect_thresh, **expect_kw)
        ax.set_xlim(axs[0].get_xlim())
        ax.set_xlabel(r'n$_{inj}$')
        ax.set_ylabel(r'threshold')
        ax.set_xlim(0,50)
        #ax.set_aspect('equal')
        #for ax in axs:
        #    ax.set_xlabel(r'$n_\text{inj}$')
        #    ax.grid()
        #axs[0].set_ylabel(r'$n_s$')
        #axs[1].set_ylabel(r'$\gamma$')
        if sens:
            sens_file = f'{state.base_dir}/{state.ana_name}/lc/{name}/trials/fit_gamma/src_gamma_{src_gamma}/thresh_{thresh}/lag_{lag}/sens.npy'
            try:
                result = np.load(sens_file, allow_pickle = True)[()]
                for ax in axs:
                    ax.axvline(result['n_sig'], ls = '--')
            except(FileNotFoundError):
                pass 
        plt.tight_layout()  
        if save:
            save_dir = cy.utils.ensure_dir(f'{state.base_dir}/{state.ana_name}/lc/plots/{name}/fit_gamma/src_gamma_{src_gamma}/')
            if logging:
                print(f'saving... {save_dir}')
            plt.savefig(f'{save_dir}/{name}_bias_g{src_gamma:.2f}_l{lag:.2f}_t{thresh}.png')

    if not name:
        sources = pd.read_hdf(cg.source_file)
        for n in sources.name_disp:
            print(n)
            try:
                get_bias(n, fix_gamma, src_gamma, thresh, lag, cutoff, cpus, sens)
            except(AttributeError):
                print('MISSING...')
                pass
    else:
        get_bias(name, fix_gamma, src_gamma, thresh, lag, cutoff, cpus, sens)

                

if __name__ == '__main__':
    exe_t0 = now ()
    print ('c7: start at {} .'.format (exe_t0))
    cli ()
