#!/usr/bin/env python
# Unblinding Script 

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
@click.option('--name', default = 'Cyg_X_dash_1')
@click.option('--n-trials', default=1, type=int)
@click.option ('--poisson/--nopoisson', default=False)
@click.option ('--fix_gamma', default=None)
@click.option ('--src_gamma', default=2.0, type=float)
@click.option ('--thresh', default=0.0, type=float)
@click.option ('--lag', default=0.0, type=float)
@click.option ('--seed', default=None, type=int)
@click.option ('--cutoff', default=np.inf, type=float, help='exponential cutoff energy, (TeV)')
@click.option ('--cpus', default=1, type=int)
@click.option ('--save_trials/--nosave_trials', default=True)
@click.option ('--inject_gp/--noinject_gp',  default=False)
@click.option ('--unblind/--nounblind',  default=False)
@pass_state
def unblind_lc_trials(state, name, n_trials, fix_gamma, src_gamma, thresh, lag, poisson, seed,
                  cutoff, cpus, save_trials, inject_gp, unblind, logging=True):
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
    if name == 'all':
        print('Unblinding ALL sources')
        sources = pd.read_hdf(cg.source_file)
        names = sources.name_disp
    else:
        names = [name]  
    cutoff_GeV = cutoff * 1e3 
    print(seed)
    for n in names:
 
        # Get the configuration for the analysis
        conf, inj_conf = cg.get_ps_config(
            ana,
            n, 
            src_gamma, 
            fix_gamma, 
            cutoff_GeV, 
            lag, 
            thresh,
            inject_gp=inject_gp
        ) 
        truth = unblind
        if truth:
            print('WARNING: UNBLINDING')   
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
        if truth:
            trials = tr.get_many_fits(
                n_trials, n_sig=0, poisson=False, seed=seed, logging=logging, TRUTH=True,
            )
            print(trials.as_table)
            out_dir = cy.utils.ensure_dir(
                '{}/{}/lc/{}/results/'.format(
                    state.base_dir, state.ana_name, n))
            out_file = '{}/result.npy'.format(
                out_dir)
            sigfile = '{}/{}/lc/TSD_chi2.dict'.format(state.base_dir, state.ana_name)
            sig = np.load(sigfile, allow_pickle=True)
            b = cy.bk.get_best(sig, 'name', n, 'fit_gamma', 'gp_inj',  'src_gamma_2.0'.format(src_gamma),        
                             'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff_GeV),  'bg')
            print(f'p-value calculation with {len(b.trials.ts)} trials')
            pval = np.mean(b.trials.ts > trials.ts)
            print(f'{name} p-value: {pval}')
            trials['pval'] = [pval]
            np.save(out_file, trials.as_array)
        else:    
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
@click.option('--n-trials', default=1, type=int)
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
@click.option ('--save_trials/--nosave_trials', default=True)
@click.option ('--unblind/--nounblind',  default=False)
@pass_state
def unblind_stacking_trials(state, n_trials, fix_gamma, src_gamma, thresh, lag, n_sig, poisson, seed, cutoff, cpus, weight, inject_gp, save_trials,
            unblind, logging=True):
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
    if unblind:
        print('Unblinding')
        trials = tr.get_many_fits(
            n_trials, n_sig=n_sig, poisson=poisson, seed=seed, logging=logging, mp_cpus=cpus, TRUTH=False)
        print(trials.as_table)
        sigfile = '{}/{}/lc/stacking/TSD_chi2.dict'.format(state.base_dir, state.ana_name)
        sig = np.load(sigfile, allow_pickle=True)
        b = cy.bk.get_best(sig, 'fit_gamma', 'src_gamma_{}'.format(src_gamma),
                   'thresh_{}'.format(thresh), 'lag_{}'.format(lag), 'cutoff_{}'.format(cutoff), 'weight_{}'.format(weight), 'gp_inj', 'bg')
        print(f'Calculating pvalue with {len(b.trials)} trials')
        pval = np.mean(b.trials.ts > trials.ts)
        print(f'{weight} p-value: {pval}')
        trials['pval'] = [pval]
        
        out_dir = cy.utils.ensure_dir(
                    '{}/{}/lc/stacking/trials/results/'.format(
                        state.base_dir, state.ana_name
                    )
                )
        out_file = f'{out_dir}/results_{weight}.npy'
        np.save(out_file, trials.as_array)

    else:
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

                

if __name__ == '__main__':
    exe_t0 = now ()
    print ('c7: start at {} .'.format (exe_t0))
    cli ()
