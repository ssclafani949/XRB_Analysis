import requests
import numpy as np

import os
from bayesian_block import block 

url = 'http://www.maxi.jaxa.jp/obs/agn_etc/data/'
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
for i in cat:
    source = np.genfromtxt(source_path+i+'_det_maxi.txt',dtype=str,delimiter=',')
    if len(source[0]) == 5:
        name = [source[j][1] for j in range(len(source))]
    else:
        name = [source[1]]

    csv = np.genfromtxt(source_path+'maxi_cat.csv',dtype=str,delimiter=',')
    name1 = csv[:,3]
    name2 = csv[:,4]
    if not os.path.exists(save_directory + i):
        os.makedirs(save_directory + i )
    for k in name:
        index = np.where(name2==k)[0]
        k1 = name1[index][0]
        print(k)
        print(k1)
        if k == 'IGR J13020-6359':
            k1 ='J1302-638'
        k1 = k1.replace(' ','')
        fullurl = url + k1 + '/' + k1 +'_g_lc_1orb_all.dat'
        print(fullurl) 
        download_file(fullurl, save_directory + i)
        
        print(save_directory+'/' + i +'/'+ k1)
        print(k.replace(' ',''))
        block(lc_path = save_directory+i+'/' +k1, 
            name = k.replace('',''),
            subdir = 'maxi_lc/',
            MAXI=True)
