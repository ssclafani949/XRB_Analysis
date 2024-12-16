# config.py
import os
import socket
import numpy as np
import csky as cy
import getpass
import pandas as pd
import histlite as hl

hostname = socket.gethostname()
username = getpass.getuser()
print('Running as User: {} on Hostname: {}'.format(username, hostname))
job_base = 'XRB_baseline_v0.2'
if 'condor00' in hostname or 'cobol' in hostname:
    submit_cfg_file = 'XRB_Analysis/XRB/submitter_config_umd'
    repo = cy.selections.Repository(
        local_root='/data/i3store/users/analyses')
    ana_dir = cy.utils.ensure_dir(
        '/data/i3store/users/{}/data/analyses'.format(username))
    base_dir = cy.utils.ensure_dir(
        '/data/i3store/users/{}/data/analyses/{}'.format(username, job_base))
    job_basedir = '/data/i3home/{}/submitter_logs'.format(username)
    source_file  = '/data/i3home/ssclafani/XRB_Analysis/XRB/sources/lc_sources_reselected.hdf'
elif 'gpu' in hostname:
    if os.path.exists( '/data/i3home/'):
        print('I am on UMD')
        repo = cy.selections.Repository(
            local_root='/data/i3store/users/analyses')
        ana_dir = cy.utils.ensure_dir(
            '/data/i3store/users/{}/data/analyses'.format(username))
        base_dir = cy.utils.ensure_dir(
            '/data/i3store/users/{}/data/analyses/{}'.format(username, job_base))
        job_basedir = '/data/i3home/{}/submitter_logs'.format(username)
        source_file  = '/data/i3home/ssclafani/XRB_Analysis/XRB/sources/lc_sources_reselected.hdf'
        submit_cfg_file = 'XRB_Analysis/XRB/submitter_config_umd'
    elif os.path.exists( '/data/user/'):
        print('I am on NPX')
        repo = cy.selections.Repository()
        ana_dir = cy.utils.ensure_dir('/data/user/{}/data/analyses'.format(username))
        base_dir = cy.utils.ensure_dir('/data/user/{}/data/analyses/{}'.format(username, job_base))
        ana_dir = '{}/ana'.format (base_dir)
        job_basedir = '/scratch/{}/'.format(username) 
        source_file  = '/home/ssclafani/XRB_Analysis/XRB/sources/lc_sources_reselected.hdf'
        submit_cfg_file = 'XRB_Sens/submitter_config_npx'
    else:
        print('Could not find direcotry')
else:
    repo = cy.selections.Repository()
    ana_dir = cy.utils.ensure_dir('/data/user/{}/data/analyses'.format(username))
    base_dir = cy.utils.ensure_dir('/data/user/{}/data/analyses/{}'.format(username, job_base))
    ana_dir = '{}/ana'.format (base_dir)
    job_basedir = '/scratch/{}/'.format(username) 
    source_file  = '/home/ssclafani/XRB_Analysis/XRB/sources/lc_sources_reselected.hdf'
    submit_cfg_file = 'XRB_Sens/submitter_config_npx'


# path at which source catalogs are located
catalog_dir = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'catalogs')

# Path to submit config file. This needs to be a relative path to $HOME
# Example content of this file:
#    eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.0.1/setup.sh`
#    source  ~/path/to/venv/bin/activate



def get_ps_config(ana, name, src_gamma, fix_gamma, cutoff_GeV, lag, thresh):
    
    sources = pd.read_hdf(source_file)
    source = sources.loc[sources['name_disp'] == name]
    print(source.dec_deg)
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
    dir = cy.utils.ensure_dir ('{}/lc/{}'.format (base_dir, name))
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
    inj_conf=dict(lcs=lc, threshs=thresh, lags=lag, flux = cy.hyp.PowerLawFlux(gamma=src_gamma))

    del sources, lc, 
    return conf, inj_conf

def get_stacking_config(ana, src_gamma, fix_gamma, thresh, lag, weight):
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
    print(len(decs), len(ras), len(lcs))
    if weight == 'equal':
        src_weight = 1./(len(decs)) * np.ones_like(decs)
    elif weight == 'flux':
        src_weight = np.concatenate(flux_weight / np.sum(flux_weight))
    src = cy.utils.Sources(dec = np.concatenate(decs), ra = np.concatenate(ras), deg=True)

    print(f'inj lag: {lag}')
    print(f'inj thresh: {thresh}')                                                                 
    dir = cy.utils.ensure_dir ('{}/lc/{}'.format (base_dir, name))
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
                concat_evs=True,
                extra_keep = ['energy'],
                n_seeds_lag = 20,
                update_bg = True,
                sigsub = True,
                fitter_args = fitter_dict,
                )                                                                                                                   
    inj_conf=dict(lcs=lcs, threshs=np.zeros_like(lcs), lags=np.zeros_like(lcs), flux = cy.hyp.PowerLawFlux(gamma=src_gamma))

    return conf, inj_conf


