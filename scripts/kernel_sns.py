# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 20:28:36 2017

@author: np838619
"""


from __future__ import division
import pylab as pl
from netCDF4 import Dataset
import cartopy
import cartopy.crs as ccrs
from matplotlib.colors import Normalize
from matplotlib import ticker
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
shed = 'amr'; catch = 'all'
bas1 = 'pac'; bas2 = 'atl'

filenames = pl.zeros([5,12],dtype='S99')

years = ['2010','2011','2012','2013','2014']

for yr in range(len(years)):
    flist = glob.glob(clusdir+years[yr]+'/'+shed+'/ordens*')#PrintFiles(clusdir+years[yr]+'/'+shed+'/','ordens')
    flist = pl.sort(flist)
    filenames[yr] = flist

density = pl.zeros([5,12,256,512])
flux = pl.zeros_like(density)
flx_bas1 = pl.zeros_like(density); flx_bas2 = pl.zeros_like(density)

for yr in range(len(years)):
    for m in range(filenames.shape[1]):
        ncfile = Dataset(filenames[yr,m],'r')
        D = ncfile.variables['origin density '+shed+catch+' all'][:]
        F = ncfile.variables['flux density '+shed+catch+' all'][:]
        B1 = ncfile.variables['flux density '+shed+bas1+' all'][:]
        B2 = ncfile.variables['flux density '+shed+bas2+' all'][:]
        ncfile.close()
        density[yr,m] = D
        flux[yr,m] = F
        flx_bas1[yr,m] = B1; flx_bas2[yr,m] = B2

dens_sns, flx_sns = SeasonalCalcs(density,flux)
dens_sns, fb1_sns = SeasonalCalcs(density,flx_bas1)
dens_sns, fb2_sns = SeasonalCalcs(density,flx_bas2)

levels =  [-100000, -85000, -70000, -55000, -40000, -25000, -10000,-1000,
          1000,10000, 25000, 40000, 55000, 70000, 85000, 100000]
#levels = [-500000,-400000,-300000,-200000,-100000,-50000,-20000,-10000,
#          10000,20000,50000,100000,200000,300000,400000,500000]
norm = pl.Normalize(-100000,100000,clip=False)
lons,lats = pl.meshgrid(eralon,eralat)
proj = ccrs.Robinson(central_longitude=30)

fig,ax = pl.subplots(2,2,figsize=(11,8))

ax1 = pl.subplot(221,projection=proj); ax1.coastlines()
cs = ax1.contourf(lons, lats, flx_sns[0],transform=ccrs.PlateCarree(),norm=norm,
                 extend='both',cmap='seismic',levels=levels)
ax1.set_title('(a) DJF',fontsize=18); ax1.gridlines()
#ax1.text(-1.96e7,-2.9e5,'0$^\circ$',fontsize=15)
#ax1.text(-2.01e7,-3.9e6,'-30$^\circ$',fontsize=15)
#ax1.text(-1.83e7,-6.7e6,'-60$^\circ$',fontsize=15)
#ax1.text(-1.97e7,3.0e6,'30$^\circ$',fontsize=15)
#ax1.text(-1.75e7,6.0e6,'60$^\circ$',fontsize=15)

ax2 = pl.subplot(222,projection=proj); ax2.coastlines()
cs = ax2.contourf(lons, lats, flx_sns[1],transform=ccrs.PlateCarree(),norm=norm,
                 extend='both',cmap='seismic',levels=levels)
ax2.set_title('(b) MAM',fontsize=18); ax2.gridlines()
#ax2.text(-1.96e7,-2.9e5,'0$^\circ$',fontsize=18)

ax3 = pl.subplot(223,projection=proj); ax3.coastlines()
cs = ax3.contourf(lons, lats, flx_sns[2],transform=ccrs.PlateCarree(),norm=norm,
                 extend='both',cmap='seismic',levels=levels)
ax3.set_title('(c) JJA',fontsize=18); ax3.gridlines()
#ax3.text(-1.96e7,-2.9e5,'0$^\circ$',fontsize=15)
#ax3.text(-2.01e7,-3.9e6,'-30$^\circ$',fontsize=15)
#ax3.text(-1.83e7,-6.7e6,'-60$^\circ$',fontsize=15)
#ax3.text(-1.97e7,3.0e6,'30$^\circ$',fontsize=15)
#ax3.text(-1.75e7,6.0e6,'60$^\circ$',fontsize=15)
#ax3.text(6.1e4,-1.01e7,'0$^\circ$',fontsize=15)
#ax3.text(-5.0e6,-1.01e7,'-60$^\circ$',fontsize=15)
#ax3.text(-1.0e7,-1.01e7,'-120$^\circ$',fontsize=15)
#ax3.text(2.7e6,-1.01e7,'60$^\circ$',fontsize=15)
#ax3.text(6.2e6,-1.01e7,'120$^\circ$',fontsize=15)

ax4 = pl.subplot(224,projection=proj); ax4.coastlines()
cs = ax4.contourf(lons, lats, flx_sns[3],transform=ccrs.PlateCarree(),norm=norm,
                 extend='both',cmap='seismic',levels=levels)
ax4.set_title('(d) SON',fontsize=18); ax4.gridlines()
#ax4.text(6.1e4,-1.01e7,'0$^\circ$',fontsize=15)
#ax4.text(-5.0e6,-1.01e7,'-60$^\circ$',fontsize=15)
#ax4.text(-1.0e7,-1.01e7,'-120$^\circ$',fontsize=15)
#ax4.text(2.7e6,-1.01e7,'60$^\circ$',fontsize=15)
#ax4.text(6.2e6,-1.01e7,'120$^\circ$',fontsize=15)

pl.subplots_adjust(left=0.05,right=0.95,hspace=-0.3,wspace=0.05,top=1.0)

f = pl.gcf()
colax = f.add_axes([0.11,0.12,0.76,0.05])
clb = pl.colorbar(cs, cax=colax,orientation='horizontal')
clb.set_ticks(levels)
ticklabs = pl.asarray(levels); ticklabs = ticklabs/1000
clb.set_ticklabels(ticklabs.astype('int'))
clb.ax.tick_params(labelsize=14)
clb.update_ticks()#; cb.ax.set_aspect(0.09)
clb.set_label('10$^3$ kg s$^{-1}$ steradian$^{-1}$',fontsize=20)

#pl.savefig(homedir+shed+'/'+shed+'_sns_maps_'+catch+'_all.png')

flx_sums = pl.sum(flx_sns,axis=(1,2))/(10**9)
fb1_sums = pl.sum(fb1_sns,axis=(1,2))/(10**9)
fb2_sums = pl.sum(fb2_sns,axis=(1,2))/(10**9)
colors = ['b','g','r','grey']
names = ['DJF','MAM','JJA','SON']
pos = pl.arange(4); width=0.25

fig,ax = pl.subplots()
#for i in range(4):
ax.bar(pos,fb1_sums,width,color='b',label='Indian')
ax.bar(pos+width,flx_sums,width,color='k',label='total')
ax.bar(pos+2*width,fb2_sums,width,color='g',label='Pacific')
pl.legend(loc=2,fontsize=18)
pl.ylim(-0.4,0.7); pl.grid(axis='y')
pl.ylabel('Sv',fontsize=22,labelpad=1)
pl.yticks(pl.linspace(-0.4,0.7,12),fontsize=16)
pl.xticks(pos+3*width/2)
ax.set_xticklabels(names,fontsize=22)
pl.subplots_adjust(left=0.14,right=0.91)
pl.show()
#pl.savefig('/home/np838619/Trajectory/kernel/'+shed+'/flx_sns_bars_'+catch+'_all.png')
