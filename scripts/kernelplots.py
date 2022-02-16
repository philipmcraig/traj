# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 10:41:16 2017

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
import cartopy
import cartopy.crs as ccrs
from matplotlib.colors import Normalize
from matplotlib import ticker
import matplotlib.gridspec as gridspec
import matplotlib.path as mpath

class MidpointNormalize(Normalize):
    """Found this on the internet. Use to centre colour bar at zero.
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        a, b = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return pl.ma.masked_array(pl.interp(value, a, b))

pl.close('all')
exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

clusdir = '/glusterfs/scenario/users/np838619/traj/monthly_means/clim/'
sheddir = '/home/np838619/Watershed/shed_defs/'
trajdir = '/home/np838619/Trajectory/kernel/'
#catchments = ['all','atl','ind','pac','arc','sou']
#categories = ['all','CAT I','CAT II']
shed = 'sop'; catch = 'all'; cat = 'CAT II'
sheds = ['amr','afr','eaa','ara','ari','arp','soa','soi','sop']
locs = ['Am','AfMe','EAA','ArA','ArI','ArP','SOA','SOI','SOP']
names = ['(a) Americas','(b) Africa','(c) South-East Asia',
         '(d) Arctic Atlantic','(e) Arctic Indian','(f) Arctic Pacific',
         '(g) Southern Atlantic','(h) Southern Indian','(i) Southern Pacific' ]
#rlspts = pl.genfromtxt('/home/np838619/Watershed/shed_defs/SOP_traj_release_new.txt',
#                       skip_header=5)

ncfile = Dataset('/home/np838619/Watershed/quv_annmean.nc','r')
eralat = ncfile.variables['lat'][:]
eralon = ncfile.variables['lon'][:]
ncfile.close()


#for shed in sheds:
filenames = PrintFiles(clusdir+'/'+shed+'/','ordens')
filenames = pl.sort(filenames)

#ntraj10 = pl.genfromtxt(clusdir+'2010/amr/NTRAJ/amr_all201001_all.csv')
#ntraj11 = pl.genfromtxt(clusdir+'2011/amr/NTRAJ/amr_all201101_all.csv')
#ntraj12 = pl.genfromtxt(clusdir+'2012/amr/NTRAJ/amr_all201201_all.csv')
#ntraj13 = pl.genfromtxt(clusdir+'2013/amr/NTRAJ/amr_all201301_all.csv')
#ntraj14 = pl.genfromtxt(clusdir+'2014/amr/NTRAJ/amr_all201401_all.csv')

#N = pl.sum(ntraj10) + pl.sum(ntraj11) + pl.sum(ntraj12) + pl.sum(ntraj13) + pl.sum(ntraj14)
#    for catch in catchments:
#        for cat in categories:
density = pl.zeros([len(sheds),len(filenames),eralat.size,eralon.size])
flux = pl.zeros_like(density)

for s in range(len(sheds)):
    filenames = PrintFiles(clusdir+'/'+sheds[s]+'/','ordens')
    filenames = pl.sort(filenames)
    for name in range(filenames.size):
        ncfile = Dataset(clusdir+sheds[s]+'/'+filenames[name],'r')
        #lon = ncfile.variables['lon'][:]
        #lat = ncfile.variables['lat'][:]
        D = ncfile.variables['origin density '+sheds[s]+catch+' '+cat][:]
        F = ncfile.variables['flux density '+sheds[s]+catch+' '+cat][:]
        ncfile.close()
        density[s,name] = D
        flux[s,name] = F

Ds = pl.nanmean(density,axis=1)
Fm = pl.nanmean(flux,axis=1)

#for i in range(fm.shape[0]):
#    for j in range(fm.shape[1]):
#        if fm[i,j] == 0.:
#            fm[i,j] = pl.float64('nan')

eralon[-1] = 360
cl = [-70,0,90,0,60,180,0,80,180]

#pl.figure(1)
fig, ax = pl.subplots(3,3,figsize=(15,11))
levels1 = [1,50,100,150,200,250,300,350]
levels2 = [1,25,50,75,100,125,150]
norm1 = pl.Normalize(1,350,clip=True); norm2 = pl.Normalize(1,150,clip=True)
lons,lats = pl.meshgrid(eralon,eralat)
for i in range(9):
    axx = pl.subplot(3,3,i+1,projection=ccrs.Robinson(central_longitude=cl[i]))
    axx.coastlines(color='grey'); axx.gridlines()
    if i in (0,1,2,6,7,8):
        norm = norm1; levels = levels1
    else:
        norm = norm2; levels = levels2
    cs = axx.contourf(lons, lats, Zero2Nan(Ds[i]),transform=ccrs.PlateCarree(),
                      norm=norm,levels=levels,extend='both',cmap='magma_r')
    rlspts = pl.genfromtxt(sheddir+locs[i]+'_traj_release_new.txt',skip_header=5)
    axx.plot(rlspts[:,0],rlspts[:,1],color='k',transform=ccrs.Geodetic(),lw=1.5)
    pl.title(names[i],fontsize=18)
    #if i == 1:
    cb = pl.colorbar(cs,orientation='horizontal',pad=0.03,fraction=0.08)
    cb.ax.tick_params(labelsize=14); cb.ax.set_aspect(0.09)
    if i == 7:
        cb.set_label('trajectories/steradian/month',fontsize=15,labelpad=4)
    
#f = pl.gcf()
#colax = f.add_axes([0.00,0.08,1.0,0.03])
#cb = pl.colorbar(cs,orientation='horizontal',cax=colax)#pad=0.05,fraction=0.10,
#cb.set_label('trajectories/steradian/month',fontsize=20)
#cb.ax.tick_params(labelsize=16); cb.ax.set_aspect(0.09)
#    #ax1.plot(rlspts[:,0],rlspts[:,1],color='k',lw=2,transform=ccrs.Geodetic())
#    #pl.title('Annual mean origin density '+shed+' '+catch+' '+cat)
pl.suptitle('CAT II trajectory origins',fontsize=18,y=0.995)
pl.tight_layout(); pl.subplots_adjust(hspace=0.05,bottom=0.06,top=0.99)
#pl.savefig(trajdir+'panel_'+catch+'_'+cat+'.png')

#cmap, clevs = get_eofColormap(Fm)
#pl.figure(2)
##levels=pl.linspace(-100000,100000,20)
##levels = [-100000,-90000,-80000,-70000,-60000,-50000,-40000,-25000,-10000,
#          #   -2000,2000,10000,25000,40000,50000,60000,70000,80000,90000,100000]
levels = [-100000, -85000, -70000, -55000, -40000, -25000, -10000,-1000,
          1000,10000, 25000, 40000, 55000, 70000, 85000, 100000]
norm = pl.Normalize(-100000,100000,clip=False)
sd = [3,4,5]
names = ['(a) Atlantic sector','(b) Indian sector','(c) Pacific sector']#,'(d) Southern Atlantic']
cl = [0,80,180]

#fig, ax = pl.subplots(2,2,figsize=(10,10))
pl.figure(figsize=(10,10))
gs = gridspec.GridSpec(2, 4)
ig = [gs[0,:2],gs[0,2:],gs[1,1:3]]#[0,0,1]; iy = [:2,2:,1:3]
for i in range(3):
    axx = pl.subplot(ig[i],projection=ccrs.Robinson(central_longitude=cl[i]))
    axx.coastlines()
#    if i in (1,2):
#        A = -1
#    else:
#        A = 1
    cs = axx.contourf(lons, lats, Fm[sd[i]],transform=ccrs.PlateCarree(),norm=norm,
                     extend='both',cmap='seismic',levels=levels)
    axx.gridlines()
    rlspts = pl.genfromtxt(sheddir+locs[sd[i]]+'_traj_release_new.txt',skip_header=5)
    axx.plot(rlspts[:,0],rlspts[:,1],color='k',transform=ccrs.Geodetic())
    pl.title(names[i],fontsize=18)

f = pl.gcf()
colax = f.add_axes([0.1,0.2,0.8,0.03])                   
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)
cb.set_ticks(levels)
ticklabs = pl.asarray(levels); ticklabs = ticklabs/1000
cb.set_ticklabels(ticklabs.astype('int'))
cb.ax.tick_params(labelsize=14)
cb.update_ticks()#; cb.ax.set_aspect(0.09)
cb.set_label('10$^3$ kg/s/steradian',fontsize=20)
#pl.title('Annual mean flux density '+shed+' '+catch+' '+cat)
pl.suptitle('Arctic ('+cat+' trajectory origins)',fontsize=18,y=0.82)
pl.tight_layout(); pl.subplots_adjust(hspace=-0.6)
#pl.savefig(trajdir+'flxdens_panel_arc_'+cat+'_new.png')

#pl.close()
