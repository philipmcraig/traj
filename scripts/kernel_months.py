# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 11:50:30 2017

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset

def OpenNCfile(clusdir,shed,yr,m,eralat,eralon):
    newnc = Dataset(clusdir+'monthly_means/'+yr+'/'+shed+'/ordens_'+shed+yr+m+'mm.nc','w')
    
    lat_dim = newnc.createDimension('lat',eralat.size)
    lat_in = newnc.createVariable('lat',pl.float64,('lat',))
    lat_in.units = 'degrees_north'
    lat_in.long_name = 'latitude'
    lat_in[:] = eralat
    
    lon_dim = newnc.createDimension('lon',eralon.size)
    lon_in = newnc.createVariable('lon',pl.float64,('lon',))
    lon_in.units = 'degrees_east'
    lon_in.long_name = 'longitude'
    lon_in[:] = eralon
    
    newnc.close()
    
def SavetoNC(densfield,fluxfield,shed,catch,cat,yr,m):
    # which category of trajectories? all, CAT 1 or CAT 2?
    if cat == 'all':
        field_ex = ' all'
    elif cat == '1':
        field_ex = ' CAT I'
    elif cat == '2':
        field_ex = ' CAT II'
    
    clusdir = '/net/glusterfs/scenario/users/np838619/traj/monthly_means/'
    ncfile = Dataset(clusdir+yr+'/'+shed+'/ordens_'+shed+yr+m+'mm.nc','a')
    
    D = ncfile.createVariable('origin density ' + shed + catch + field_ex,pl.float64,
                                                          ('lat','lon'))
    D.units = 'trajectories per steradian'
    D.standard_name = 'trajectory origin density ' + shed + catch + field_ex
    D[:,:] = densfield
    
    F = ncfile.createVariable('flux density ' + shed + catch + field_ex,pl.float64,
                                                              ('lat','lon'))
    F.units = 'kg/s per steradian'
    F.standard_name = 'moisture flux density ' + shed + catch + field_ex
    F[:,:] = fluxfield
    
    ncfile.close() 

clusdir = '/net/glusterfs/scenario/users/np838619/traj/'

years = ['2010','2011','2012','2013','2014']
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
sheds = ['amr','afr','eaa','ara','ari','arp','soa','soi','sop']

ncfile = Dataset('/home/np838619/Watershed/quv_annmean.nc','r')
eralat = ncfile.variables['lat'][:]
eralon = ncfile.variables['lon'][:]
ncfile.close()

for yr in years:
    for shed in sheds:
        for m in months:
            print '***PROCESSING ' + yr + ' ' + m + ' ' + shed + '***'
            ncfile = Dataset(clusdir+yr+'/'+shed+'/density/ordens_'+shed+yr+m+'.nc','r')
            den_all_all = ncfile.variables['origin density '+shed+'all all'][:]
            flx_all_all = ncfile.variables['flux density '+shed+'all all'][:]
            den_all_c1 = ncfile.variables['origin density '+shed+'all CAT I'][:]
            flx_all_c1 = ncfile.variables['flux density '+shed+'all CAT I'][:]
            den_all_c2 = ncfile.variables['origin density '+shed+'all CAT II'][:]
            flx_all_c2 = ncfile.variables['flux density '+shed+'all CAT II'][:]
            den_atl_all = ncfile.variables['origin density '+shed+'atl all'][:]
            flx_atl_all = ncfile.variables['flux density '+shed+'atl all'][:]
            den_atl_c1 = ncfile.variables['origin density '+shed+'atl CAT I'][:]
            flx_atl_c1 = ncfile.variables['flux density '+shed+'atl CAT I'][:]
            den_atl_c2 = ncfile.variables['origin density '+shed+'atl CAT II'][:]
            flx_atl_c2 = ncfile.variables['flux density '+shed+'atl CAT II'][:]
            den_ind_all = ncfile.variables['origin density '+shed+'ind all'][:]
            flx_ind_all = ncfile.variables['flux density '+shed+'ind all'][:]
            den_ind_c1 = ncfile.variables['origin density '+shed+'ind CAT I'][:]
            flx_ind_c1 = ncfile.variables['flux density '+shed+'ind CAT I'][:]
            den_ind_c2 = ncfile.variables['origin density '+shed+'ind CAT II'][:]
            flx_ind_c2 = ncfile.variables['flux density '+shed+'ind CAT II'][:]
            den_pac_all = ncfile.variables['origin density '+shed+'pac all'][:]
            flx_pac_all = ncfile.variables['flux density '+shed+'pac all'][:]
            den_pac_c1 = ncfile.variables['origin density '+shed+'pac CAT I'][:]
            flx_pac_c1 = ncfile.variables['flux density '+shed+'pac CAT I'][:]
            den_pac_c2 = ncfile.variables['origin density '+shed+'pac CAT II'][:]
            flx_pac_c2 = ncfile.variables['flux density '+shed+'pac CAT II'][:]
            den_arc_all = ncfile.variables['origin density '+shed+'arc all'][:]
            flx_arc_all = ncfile.variables['flux density '+shed+'arc all'][:]
            den_arc_c1 = ncfile.variables['origin density '+shed+'arc CAT I'][:]
            flx_arc_c1 = ncfile.variables['flux density '+shed+'arc CAT I'][:]
            den_arc_c2 = ncfile.variables['origin density '+shed+'arc CAT II'][:]
            flx_arc_c2 = ncfile.variables['flux density '+shed+'arc CAT II'][:]
            den_sou_all = ncfile.variables['origin density '+shed+'sou all'][:]
            flx_sou_all = ncfile.variables['flux density '+shed+'sou all'][:]
            den_sou_c1 = ncfile.variables['origin density '+shed+'sou CAT I'][:]
            flx_sou_c1 = ncfile.variables['flux density '+shed+'sou CAT I'][:]
            den_sou_c2 = ncfile.variables['origin density '+shed+'sou CAT II'][:]
            flx_sou_c2 = ncfile.variables['flux density '+shed+'sou CAT II'][:]
            ncfile.close()
            
            # sum trajectory density fields along time axis:
            den_all_all = pl.nansum(den_all_all,axis=0)
            den_all_c1 = pl.nansum(den_all_c1,axis=0)
            den_all_c2 = pl.nansum(den_all_c2,axis=0)
            den_atl_all = pl.nansum(den_atl_all,axis=0)
            den_atl_c1 = pl.nansum(den_atl_c1,axis=0)
            den_atl_c2 = pl.nansum(den_atl_c2,axis=0)
            den_ind_all = pl.nansum(den_ind_all,axis=0)
            den_ind_c1 = pl.nansum(den_ind_c1,axis=0)
            den_ind_c2 = pl.nansum(den_ind_c2,axis=0)
            den_pac_all = pl.nansum(den_pac_all,axis=0)
            den_pac_c1 = pl.nansum(den_pac_c1,axis=0)
            den_pac_c2 = pl.nansum(den_pac_c2,axis=0)
            den_arc_all = pl.nansum(den_arc_all,axis=0)
            den_arc_c1 = pl.nansum(den_arc_c1,axis=0)
            den_arc_c2 = pl.nansum(den_arc_c2,axis=0)
            den_sou_all = pl.nansum(den_sou_all,axis=0)
            den_sou_c1 = pl.nansum(den_sou_c1,axis=0)
            den_sou_c2 = pl.nansum(den_sou_c2,axis=0)
            
            # average flux density fields along time axis:
            flx_all_all = pl.nanmean(flx_all_all,axis=0)
            flx_all_c1 = pl.nanmean(flx_all_c1,axis=0)
            flx_all_c2 = pl.nanmean(flx_all_c2,axis=0)
            flx_atl_all = pl.nanmean(flx_atl_all,axis=0)
            flx_atl_c1 = pl.nanmean(flx_atl_c1,axis=0)
            flx_atl_c2 = pl.nanmean(flx_atl_c2,axis=0)
            flx_ind_all = pl.nanmean(flx_ind_all,axis=0)
            flx_ind_c1 = pl.nanmean(flx_ind_c1,axis=0)
            flx_ind_c2 = pl.nanmean(flx_ind_c2,axis=0)
            flx_pac_all = pl.nanmean(flx_pac_all,axis=0)
            flx_pac_c1 = pl.nanmean(flx_pac_c1,axis=0)
            flx_pac_c2 = pl.nanmean(flx_pac_c2,axis=0)
            flx_arc_all = pl.nanmean(flx_arc_all,axis=0)
            flx_arc_c1 = pl.nanmean(flx_arc_c1,axis=0)
            flx_arc_c2 = pl.nanmean(flx_arc_c2,axis=0)
            flx_sou_all = pl.nanmean(flx_sou_all,axis=0)
            flx_sou_c1 = pl.nanmean(flx_sou_c1,axis=0)
            flx_sou_c2 = pl.nanmean(flx_sou_c2,axis=0)
            
            # open an nc file:
            OpenNCfile(clusdir,shed,yr,m,eralat,eralon)
            SavetoNC(den_all_all,flx_all_all,shed,'all','all',yr,m)
            SavetoNC(den_all_c1,flx_all_c1,shed,'all','1',yr,m)
            SavetoNC(den_all_c2,flx_all_c2,shed,'all','2',yr,m)
            SavetoNC(den_atl_all,flx_atl_all,shed,'atl','all',yr,m)
            SavetoNC(den_atl_c1,flx_atl_c1,shed,'atl','1',yr,m)
            SavetoNC(den_atl_c2,flx_atl_c2,shed,'atl','2',yr,m)
            SavetoNC(den_ind_all,flx_ind_all,shed,'ind','all',yr,m)
            SavetoNC(den_ind_c1,flx_ind_c1,shed,'ind','1',yr,m)
            SavetoNC(den_ind_c2,flx_ind_c2,shed,'ind','2',yr,m)
            SavetoNC(den_pac_all,flx_pac_all,shed,'pac','all',yr,m)
            SavetoNC(den_pac_c1,flx_pac_c1,shed,'pac','1',yr,m)
            SavetoNC(den_pac_c2,flx_pac_c2,shed,'pac','2',yr,m)
            SavetoNC(den_arc_all,flx_arc_all,shed,'arc','all',yr,m)
            SavetoNC(den_arc_c1,flx_arc_c1,shed,'arc','1',yr,m)
            SavetoNC(den_arc_c2,flx_arc_c2,shed,'arc','2',yr,m)
            SavetoNC(den_sou_all,flx_sou_all,shed,'sou','all',yr,m)
            SavetoNC(den_sou_c1,flx_sou_c1,shed,'sou','1',yr,m)
            SavetoNC(den_sou_c2,flx_sou_c2,shed,'sou','2',yr,m)
            print '*** ' + yr + ' ' + m + ' ' + shed + ' WRITTEN TO FILE***'