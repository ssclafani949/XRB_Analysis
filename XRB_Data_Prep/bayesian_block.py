#Generate bayesian blocks and the light curve plot. One light curve plot is generated and two block files: one allows negative numbers and another force all negative values to be zero.
import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import argparse
from astropy.stats import bayesian_blocks
from astropy.time import Time
from astropy.io import fits

p = argparse.ArgumentParser(description="maximum likelihood block method",
                            formatter_class=argparse.RawTextHelpFormatter)
p.add_argument("--lc_path",type=str,
                help="folder of the light curve file")
p.add_argument("--name",type=str,
                help="name of the source")
p.add_argument("--subdir",type=str,
                help="sub directory of output files")
p.add_argument("--MAXI",type=bool,default=False,
                help="swift lc or MAXI lc, default is swift")
args = p.parse_args()


def func(x,a):
        return a

def block(lc_path,name,subdir,MAXI=False):
        r""" make bayesian blocks. Output file is stored and light curve is plotted.
        Parameters
        ------------
        lc_path         : str
                                  folder of the lightcurve
        name            : str
                                  source name
        subdir          : str
                                  subdir of the source
        MAXI            : bool
                                  if true it is a MAXI light curve; if not it is a swift light cuve. 
        
        ------------
        
        """
        
        if not MAXI:
                #parsing input files
                light = fits.open(lc_path+name.replace(' ','').replace('_','').replace('+','p')+'.lc.fits')
                t = light[1].data['TIME']
                rate = light[1].data['RATE']
                lc_error = light[1].data['Error']
                flag = light[1].data['DATA_FLAG']
                index = np.where(flag==0)[0]
                t = t[index]
                rate = rate[index]
                lc_error = lc_error[index]
        
        elif MAXI:
                light     = np.genfromtxt(lc_path + '_g_lc_1orb_all.dat' )
                print (light)
                t        = light[:,0]
                rate     = light[:,-2]
                lc_error = light[:,-1]
        
        data = {}
        data['MJD']=np.empty(len(t), dtype=[('MJD',float),('rate',float), ('error', float)])
        data['MJD']['MJD'] = t 
        data['MJD']['rate'] = rate 
        data['MJD']['error'] = lc_error 

        # MJD to date 
        t_min = Time(min(data['MJD']['MJD']),format = 'mjd')
        t_min = t_min.datetime.year
        t_max = Time(max(data['MJD']['MJD']),format = 'mjd')
        t_max = t_max.datetime.year
        years = np.arange(t_min,t_max+1,1)
        years_mjd = [str(x)+'-01-01T00:00:00' for x in years]
        years_mjd = Time(years_mjd,format='isot',scale = 'utc')
        years_mjd = years_mjd.mjd
        
        prior  = 1.32+0.577*np.log10(len(data['MJD']['rate'])) #prior from 1207.5578 
        blocks = bayesian_blocks(t,rate,lc_error, fitness='measures',gamma = np.exp(-prior))  #standard
        

        z=[] # time of blocks   
        y=[] # value of blocks

        file = open('./{}/bin_fit_'.format(subdir)+name.replace(' ','').replace('_',''),'w')
        filep = open('./{}/bin_outburst_fit_'.format(subdir)+name.replace(' ','').replace('_',''),'w') #file of positive values

        #calculate value of each block and save files
        for j in range(len(blocks)-1):
                a = [x for x in data['MJD'] if blocks[j] <= x['MJD'] and blocks[j+1] >= x['MJD']]
                a = list(zip(*a))
                times = list(a[0])
                count = list(a[1])
                error_bar = list(a[2])
                file.write(str(blocks[j])+'\t')
                #fit constant within each block 
                if len(times) == 1:
                        fit   = a[1][0]
                else:
                        popt,pcov = curve_fit(func,times,count,p0=np.mean(count),sigma = error_bar)
                        fit = func(times,*popt)
                file.write(str(fit)+'\n')
                z.append(blocks[j])
                y.append(fit)
                filep.write(str(blocks[j])+'\t')
                fit = fit if fit >= 0. else 0.
                filep.write(str(fit)+'\t')


        file.write(str(blocks[j+1])+'\t')
        file.write(str(0.)+'\n')
        filep.write(str(blocks[j+1])+'\t')
        filep.write(str(0.)+'\n')
        z.append(blocks[j+1])
        y.append(0.)
        
        #plot light curves and bayesian block fit. 
        plt.clf()
        fig=plt.figure()
        plt.title(name.replace('p', '+').replace('m', '-'),fontsize=40,y = 1.08)
        
        if MAXI:
                label = 'MAXI 10-20 keV'
        else:   
                label = 'Swift-BAT 15-50 keV'
        plt.errorbar(t,rate,xerr=0,yerr=lc_error,fmt='.',color= 'violet',label=label,markersize=4,elinewidth=1.,alpha=0.7)
        y = np.array(y)

        plt.step(z,y,color='red',where='post',linewidth = 3,label = 'Bayesian Blocks')

        for i in range(len(years)):
                plt.axvline (x = years_mjd[i],color = 'black',linestyle='-.',linewidth=1)
                plt.text(years_mjd[i],max(y)*1.02,years[i],fontsize=20)

        plt.xlabel("MJD",fontsize=35)
        plt.xlim(min(data['MJD']['MJD']),max(data['MJD']['MJD']))
        plt.ylim(min(rate),max(rate))
        plt.xticks(fontsize=24)
        plt.ylabel(r"$\mathrm{count}\,\mathrm{cm}^{-2}s^{-1}$",fontsize=35)
        plt.yticks(fontsize=24)
        plt.legend(prop={'size':25},loc='best')
        plt.gca().get_legend().get_frame().set_color('w')
        plt.gca().get_legend().get_frame().set_edgecolor('k')
        fig.set_size_inches(32,11)
        plt.savefig('./{}/'.format(subdir)+name+'.png')
        plt.close()

if __name__ == "__main__":
    lc_path = args.lc_path
    name    = args.name
    subdir = args.subdir
    MAXI    = args.MAXI

    print( 'lc_path {}'.format(lc_path))
    print( 'source {}'.format(name))
    print( 'subdir {}'.format(subdir))

    if name == 'all':
        #do all 
        names = np.genfromtxt('swift_select.txt', dtype='str', delimiter = '/n')
        missing = []
        for n in names:
            n = n.replace(' ', '')
            try:
                block(lc_path,n, subdir, MAXI=MAXI)
            except:
                print('{} - not found'.format(n))
                missing.append(n)
        print(missing)
    elif name == 'missing': #just do the missing items from all
        names = np.load('replacement_names.npy')[:,0]
        missing = []
        for n in names:
            print(n)
            try:
                block(lc_path,n, subdir, MAXI=MAXI)
            except:
                print('{} - not found'.format(n))
                missing.append(n)
        print(missing)

    else:
        block(lc_path,name,subdir,MAXI=MAXI)
