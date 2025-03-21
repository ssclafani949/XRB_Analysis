from submitter import Submitter
import os, time


job_basedir = '/scratch/ssclafani'
T = time.time ()
job_dir = '{}/gp_tests/T_{:17.6f}'.format( job_basedir,  T)
source_file  = '/home/ssclafani/XRB_Analysis/XRB/sources/lc_sources_reselected_IC86.hdf'
submit_cfg_file = 'XRB_Analysis/XRB/submitter_config_npx'
cpus = 4

sub = Submitter (job_dir = job_dir, memory=26, ncpu = cpus, max_jobs = 1000, config = submit_cfg_file)
commands, labels = [], []

trial_script = os.path.abspath('gp_inj_test.py')
print(trial_script)
seed = 0 
n_trials = 100
n_sigs = [0, 1697]
n_jobs = 10


for n_sig in n_sigs:
    for i in range (n_jobs):
        s = i + seed
        command = f'{trial_script} --n_sig {n_sig} --n_trials {n_trials} --seed {s} --cpus {cpus}'
        label = f'gp_test_n_sig{n_sig:06d}__{n_trials:08d}__seed{s:08d}'
        commands.append(command)
        labels.append(label)

sub.submit_npx4 (commands, labels) 
