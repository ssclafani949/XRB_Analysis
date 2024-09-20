import requests
import numpy as np

import os
from bayesian_block import block 

#url = 'http://www.maxi.jaxa.jp/obs/agn_etc/data/'
url = 'http://maxi.riken.jp/star_data/'
#J0052-724/J0052-724_g_lc_1orb_all.dat'
cat = ['hm_NS','hm_un','lm_BH','lm_NS','lm_un']

source_path = 'list/'

def download_file(url, directory):
    local_filename = os.path.join(directory, url.split('/')[-1])
    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(local_filename, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        print(f"Downloaded {local_filename}")
    except requests.HTTPError as e:
        print(f"Failed to download {url}: {e}")

save_directory = 'maxi_lc_raw/'
subdir = 'maxi_lc/'

def download_and_block(run_names, save_directory):
    
    for j, k in enumerate(run_names):
        index = np.where(name2==k)[0]
        k1 = name1[index][0]
        print('{}/{}'.format(j, len(run_names)))
        print(k)
        print(k1)
        if k == 'IGR J13020-6359':
            k1 ='J1302-638'
        k1 = k1.replace(' ','')
        fullurl = url + k1 + '/' + k1 +'_g_lc_1orb_all.dat'
        print(fullurl)
        if not os.path.exists(save_directory + '/' + k1 +'_g_lc_1orb_all.dat'): 
            print('File does not exist')
            print('Downloading...')
            download_file(fullurl, save_directory)
        else:
            print('File found at:')
            print(save_directory + '/' + k1 +'_g_lc_1orb_all.dat')
            print('looking for file at ')
            print(subdir + 'bin_outburst_fit_{}'.format(k.replace(' ','')))
        if not os.path.exists(subdir +'bin_outburst_fit_{}'.format(k.replace(' ',''))):
            try: 
                block(lc_path = save_directory+'/' +k1, 
                    name = k.replace(' ',''),
                    subdir = subdir,
                    MAXI=True)
            except(FileNotFoundError):
                print('MISSING INPUT FILE')
        else:
            print('File already made') 


do_all= True 
if do_all:
    csv = np.genfromtxt(source_path+'maxi_cat.csv',dtype=str,delimiter=',')
    name1 = csv[:,3]
    name2 = csv[:,4]

    name1 = np.append(name1, ['J1758-338',  'J1734-304', 'J0632+058'])
    name2 = np.append(name2, ['4U 1755-338', 'Terzan1', 'HESS J0632+057'])
    run_names = name2
    download_and_block(run_names, save_directory + '/other/')        

else:
    for i in cat:
        source = np.genfromtxt(source_path+i+'_det_maxi.txt',dtype=str,delimiter=',')
        if len(source[0]) == 5:
            run_names = [source[j][1] for j in range(len(source))]
        else:
            run_names = [source[1]]
        csv = np.genfromtxt(source_path+'maxi_cat.csv',dtype=str,delimiter=',')
        name1 = csv[:,3]
        name2 = csv[:,4]
        if not os.path.exists(save_directory + i):
            os.makedirs(save_directory + i )
        download_and_block(run_names, save_directory + '/' + i)        
