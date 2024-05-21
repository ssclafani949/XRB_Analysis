Scripts to download and format lightcurves for XRB analysis
Sources are named in swift_select.txt and maxi_select.txt

Swift Sources:

Run download_lcs_swift.py this will download all swift lightcurves to the
directory fits_files.

Run bayseian_block.py with "--name all" to generate block files and plots for
many of the swift sources

Some sources have been renamed in the database, so fix_missing.py will generate
replacement names based on the source RA and DEC.  This will be saved to
replacement_names.py

Rerun bayesian_block.py with "--name missing" to run for the missing light
curves.

the lightcurves are saved in the directory swift_lc

Maxi Sources:
run download_lcs_maxi.py
This script should download from maxi and then create the lightcurves in one
shot

the .dat files from maxi are saved in the directory maxi_lc_raw.
the bayesian blocks curves are saved in the directory maxi_lc

