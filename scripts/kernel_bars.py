# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 14:18:08 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
import glob as glob

def SeasonalCalcs(density,flux):
    """Function to calculate the climatological seasonal means (DJF,MAM,JJA,SON)
    of a variable.
    
    Args:
        variable_array (array): variable of which the climatological seasonal 
                                means are required
    
    Returns:
        variable_seasons (array): climatological seasonal means of the input 
                                  variable
    """
    # Calculate the mean of each trio of months for each year, missing out the 
    # first year as there is no data for December the year before the data starts
    MAM_d = pl.sum(density[:,2:5],axis=1) # March, April, May
    JJA_d = pl.sum(density[:,5:8],axis=1) # June, July, August
    SON_d = pl.sum(density[:,8:11],axis=1) # September, October, November
    DJF_d = pl.sum(density[:,1:3],axis=1) + pl.sum(density[:,-1],axis=0)# December, January, February
    
    MAM_d_s = pl.sum(MAM_d,axis=0)
    JJA_d_s = pl.sum(JJA_d,axis=0)
    SON_d_s = pl.sum(SON_d,axis=0)
    DJF_d_s = pl.sum(DJF_d,axis=0)
    
    MAM_f = pl.mean(flux[:,2:5],axis=1)
    JJA_f = pl.mean(flux[:,5:8],axis=1)
    SON_f = pl.mean(flux[:,8:11],axis=1)
    DJF_f = pl.zeros_like(SON_f)
    DJF_f[0] = pl.mean(flux[:,0],axis=0)
    DJF_f[1] = pl.mean(flux[:,-2],axis=0)
    DJF_f[2] = pl.mean(flux[:,-1],axis=0)
    
    # Calculate the climatological mean of each season:  
    MAM_f_mn = pl.mean(MAM_f,axis=0)
    JJA_f_mn = pl.mean(JJA_f,axis=0)
    SON_f_mn = pl.mean(SON_f,axis=0)
    DJF_f_mn = pl.mean(DJF_f,axis=0)
    
    # Stick all seasons in one array before outputting:
    density_sns = pl.array([DJF_d_s,MAM_d_s,JJA_d_s,SON_d_s])
    flux_sns = pl.array([DJF_f_mn,MAM_f_mn,JJA_f_mn,SON_f_mn])
    
    return density_sns,flux_sns

pl.close('all')
exec(open('/home/users/qx911590/np838619/PminusE_data/ERA_Int/functions.py').read())

ncfile = Dataset('/home/users/qx911590/np838619/Watershed/quv_annmean.nc','r')
eralat = ncfile.variables['lat'][:]
eralon = ncfile.variables['lon'][:]
ncfile.close()

clusdir = '/storage/shared/glusterfs/scenario/users/qx911590/traj/monthly_means/'
homedir = '/home/users/qx911590/np838619/Trajectory/kernel/'
sheds = ['amr','afr','eaa']
catch = 'all'
#bas1 = 'ind'; bas2 = 'pac'


filenames1 = pl.zeros([5,12],dtype='S99')
filenames2 = pl.zeros_like(filenames1)
filenames3 = pl.zeros_like(filenames1)

years = ['2010','2011','2012','2013','2014']

for yr in range(len(years)):
    flist = glob.glob(clusdir+years[yr]+'/'+sheds[0]+'/ordens*')#PrintFiles(clusdir+years[yr]+'/'+sheds[0]+'/','ordens')
    flist = pl.sort(flist)
    filenames1[yr] = flist
    
    flist = glob.glob(clusdir+years[yr]+'/'+sheds[1]+'/ordens*')
    flist = pl.sort(flist)
    filenames2[yr] = flist
    
    flist = glob.glob(clusdir+years[yr]+'/'+sheds[2]+'/ordens*')
    flist = pl.sort(flist)
    filenames3[yr] = flist

dens1 = pl.zeros([5,12,256,512]); dens2 = pl.zeros_like(dens1)
dens3 = pl.zeros_like(dens1)
flux1 = pl.zeros_like(dens1); flux2 = pl.zeros_like(dens1)
flux3 = pl.zeros_like(dens1)
flx_bas1A = pl.zeros_like(dens1); flx_bas1B = pl.zeros_like(dens1)
flx_bas2A = pl.zeros_like(dens1); flx_bas2B = pl.zeros_like(dens1)
flx_bas3A = pl.zeros_like(dens1); flx_bas3B = pl.zeros_like(dens1)

#for i in range(len(sheds)):
for yr in range(len(years)):
    for m in range(filenames1.shape[1]):
        ncfile = Dataset(filenames1[yr,m],'r')
        D = ncfile.variables['origin density '+sheds[0]+catch+' all'][:]
        F = ncfile.variables['flux density '+sheds[0]+catch+' all'][:]
        B1 = ncfile.variables['flux density '+sheds[0]+'pac'+' all'][:]
        B2 = ncfile.variables['flux density '+sheds[0]+'atl'+' all'][:]
        ncfile.close()
        dens1[yr,m] = D
        flux1[yr,m] = F
        flx_bas1A[yr,m] = B1; flx_bas1B[yr,m] = B2
        
        ncfile = Dataset(filenames2[yr,m],'r')
        D = ncfile.variables['origin density '+sheds[1]+catch+' all'][:]
        F = ncfile.variables['flux density '+sheds[1]+catch+' all'][:]
        B1 = ncfile.variables['flux density '+sheds[1]+'atl'+' all'][:]
        B2 = ncfile.variables['flux density '+sheds[1]+'ind'+' all'][:]
        ncfile.close()
        dens2[yr,m] = D
        flux2[yr,m] = F
        flx_bas2A[yr,m] = B1; flx_bas2B[yr,m] = B2
        
        ncfile = Dataset(filenames3[yr,m],'r')
        D = ncfile.variables['origin density '+sheds[2]+catch+' all'][:]
        F = ncfile.variables['flux density '+sheds[2]+catch+' all'][:]
        B1 = ncfile.variables['flux density '+sheds[2]+'ind'+' all'][:]
        B2 = ncfile.variables['flux density '+sheds[2]+'pac'+' all'][:]
        ncfile.close()
        dens3[yr,m] = D
        flux3[yr,m] = F
        flx_bas3A[yr,m] = B1; flx_bas3B[yr,m] = B2

dens1_sns, flx1_sns = SeasonalCalcs(dens1,flux1)
dens1_sns, bas1A_sns = SeasonalCalcs(dens1,flx_bas1A) # Pacific
dens1_sns, bas1B_sns = SeasonalCalcs(dens1,flx_bas1B) # Atlantic

dens2_sns, flx2_sns = SeasonalCalcs(dens2,flux2)
dens2_sns, bas2A_sns = SeasonalCalcs(dens2,flx_bas2A) # Atlantic
dens2_sns, bas2B_sns = SeasonalCalcs(dens2,flx_bas2B) # Indian

dens3_sns, flx3_sns = SeasonalCalcs(dens3,flux3)
dens3_sns, bas3A_sns = SeasonalCalcs(dens3,flx_bas3A) # Indian
dens3_sns, bas3B_sns = SeasonalCalcs(dens3,flx_bas3B) # Pacific


flx1_sums = pl.sum(flx1_sns,axis=(1,2))/(10**9)
bas1A_sums = pl.sum(bas1A_sns,axis=(1,2))/(10**9)
bas1B_sums = pl.sum(bas1B_sns,axis=(1,2))/(10**9)

flx2_sums = pl.sum(flx2_sns,axis=(1,2))/(10**9)
bas2A_sums = pl.sum(bas2A_sns,axis=(1,2))/(10**9)
bas2B_sums = pl.sum(bas2B_sns,axis=(1,2))/(10**9)

flx3_sums = pl.sum(flx3_sns,axis=(1,2))/(10**9)
bas3A_sums = pl.sum(bas3A_sns,axis=(1,2))/(10**9)
bas3B_sums = pl.sum(bas3B_sns,axis=(1,2))/(10**9)


fig, ax  = pl.subplots(3,1,figsize=(10,10))
names = ['DJF','MAM','JJA','SON']
pos = pl.arange(4); width=0.25

ax1 = pl.subplot(311)
ax1.bar(pos,bas1A_sums,width,color='g',label='Pacific',ec='k')
ax1.bar(pos+width,flx1_sums,width,color='k',label='total')
ax1.bar(pos+2*width,bas1B_sums,width,color='r',label='Atlantic',ec='k')
ax1.grid(axis='y',ls=':',color='k',lw=0.25); pl.ylim(-0.6,0.2); pl.yticks(fontsize=16)
ax1.set_yticklabels(['$-$0.6',' ','$-$0.4',' ','$-$0.2',' ',0.0,' ',0.2])
pl.xticks(pos+3*width/2); pl.ylabel('Sv',fontsize=22,labelpad=1)
ax1.xaxis.set_major_formatter(pl.NullFormatter())
pl.title('(a) Americas',fontsize=18)
ax1.legend(loc=3,ncol=2,columnspacing=0.5,fontsize=16,framealpha=1)

ax2 = pl.subplot(312)
ax2.bar(pos,bas2A_sums,width,color='r',label='Atlantic',ec='k')
ax2.bar(pos+width,flx2_sums,width,color='k',label='total')
ax2.bar(pos+2*width,bas2B_sums,width,color='b',label='Indian',ec='k')
ax2.grid(axis='y',ls=':',color='k'); pl.ylim(-0.4,0.2); pl.yticks(fontsize=16)
ax2.set_yticklabels(['$-$0.4',' ','$-$0.2',' ',0.0,' ' ,0.2])
pl.xticks(pos+3*width/2); pl.ylabel('Sv',fontsize=22,labelpad=1)
ax2.xaxis.set_major_formatter(pl.NullFormatter())
pl.title('(b) Africa',fontsize=18)
ax2.legend(loc=3,ncol=2,columnspacing=0.5,fontsize=16,framealpha=1)

ax3 = pl.subplot(313)
ax3.bar(pos,bas3A_sums,width,color='b',label='Indian',ec='k')
ax3.bar(pos+width,flx3_sums,width,color='k',label='total')
ax3.bar(pos+2*width,bas3B_sums,width,color='g',label='Pacific',ec='k')
ax3.grid(axis='y',ls=':',color='k'); pl.yticks(fontsize=16); pl.ylim(-0.4,0.8)
pl.xticks(pos+3*width/2); pl.ylabel('Sv',fontsize=22,labelpad=1)
ax3.set_xticklabels(names,fontsize=22)
pl.title('(c) South-East Asia',fontsize=18)
ax3.legend(loc=2,ncol=2,columnspacing=0.5,fontsize=16,framealpha=1)

pl.tight_layout()

#pl.savefig(homedir+'flxbars_sns.png')