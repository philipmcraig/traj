# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 16:30:31 2017

@author: np838619
"""

import pylab as pl
from netCDF4 import Dataset


def OpenNCfile(clusdir,shed,m,eralat,eralon):
    newnc = Dataset(clusdir+'clim/'+shed+'/ordens_'+shed+m+'mm.nc','w')
    
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
    
def SavetoNC(densfield,fluxfield,shed,catch,cat,m):
    # which category of trajectories? all, CAT 1 or CAT 2?
    if cat == 'all':
        field_ex = ' all'
    elif cat == '1':
        field_ex = ' CAT I'
    elif cat == '2':
        field_ex = ' CAT II'
    
    clusdir = '/net/glusterfs/scenario/users/np838619/traj/monthly_means/clim'
    ncfile = Dataset(clusdir+'/'+shed+'/ordens_'+shed+m+'mm.nc','a')
    
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

clusdir = '/net/glusterfs/scenario/users/np838619/traj/monthly_means/'

years = ['2010','2011','2012','2013','2014']
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
sheds = ['amr','afr','eaa','ara','ari','arp','soa','soi','sop']

ncfile = Dataset('/home/np838619/Watershed/quv_annmean.nc','r')
eralat = ncfile.variables['lat'][:]
eralon = ncfile.variables['lon'][:]
ncfile.close()

for shed in sheds:
    for m in months:
        print '***PROCESSING ' + shed + ' ' + m + '*****'
        den_all_a_yrs = pl.zeros([len(years),eralat.size,eralon.size])
        flx_all_a_yrs = pl.zeros_like(den_all_a_yrs)
        den_all_1_yrs = pl.zeros_like(den_all_a_yrs)
        flx_all_1_yrs = pl.zeros_like(den_all_a_yrs)
        den_all_2_yrs = pl.zeros_like(den_all_a_yrs)
        flx_all_2_yrs = pl.zeros_like(den_all_a_yrs)
        den_atl_a_yrs = pl.zeros_like(den_all_a_yrs)
        flx_atl_a_yrs = pl.zeros_like(den_all_a_yrs)
        den_atl_1_yrs = pl.zeros_like(den_all_a_yrs)
        flx_atl_1_yrs = pl.zeros_like(den_all_a_yrs)
        den_atl_2_yrs = pl.zeros_like(den_all_a_yrs)
        flx_atl_2_yrs = pl.zeros_like(den_all_a_yrs)
        den_ind_a_yrs = pl.zeros_like(den_all_a_yrs)
        flx_ind_a_yrs = pl.zeros_like(den_all_a_yrs)
        den_ind_1_yrs = pl.zeros_like(den_all_a_yrs)
        flx_ind_1_yrs = pl.zeros_like(den_all_a_yrs)
        den_ind_2_yrs = pl.zeros_like(den_all_a_yrs)
        flx_ind_2_yrs = pl.zeros_like(den_all_a_yrs)
        den_pac_a_yrs = pl.zeros_like(den_all_a_yrs)
        flx_pac_a_yrs = pl.zeros_like(den_all_a_yrs)
        den_pac_1_yrs = pl.zeros_like(den_all_a_yrs)
        flx_pac_1_yrs = pl.zeros_like(den_all_a_yrs)
        den_pac_2_yrs = pl.zeros_like(den_all_a_yrs)
        flx_pac_2_yrs = pl.zeros_like(den_all_a_yrs)
        den_arc_a_yrs = pl.zeros_like(den_all_a_yrs)
        flx_arc_a_yrs = pl.zeros_like(den_all_a_yrs)
        den_arc_1_yrs = pl.zeros_like(den_all_a_yrs)
        flx_arc_1_yrs = pl.zeros_like(den_all_a_yrs)
        den_arc_2_yrs = pl.zeros_like(den_all_a_yrs)
        flx_arc_2_yrs = pl.zeros_like(den_all_a_yrs)
        den_sou_a_yrs = pl.zeros_like(den_all_a_yrs)
        flx_sou_a_yrs = pl.zeros_like(den_all_a_yrs)
        den_sou_1_yrs = pl.zeros_like(den_all_a_yrs)
        flx_sou_1_yrs = pl.zeros_like(den_all_a_yrs)
        den_sou_2_yrs = pl.zeros_like(den_all_a_yrs)
        flx_sou_2_yrs = pl.zeros_like(den_all_a_yrs)
        for yr in range(len(years)):
            y = years[yr]
            ncfile = Dataset(clusdir+y+'/'+shed+'/ordens_'+shed+y+m+'mm.nc','r')
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
            
            den_all_a_yrs[yr] = den_all_all
            flx_all_a_yrs[yr] = flx_all_all
            den_all_1_yrs[yr] = den_all_c1
            flx_all_1_yrs[yr] = flx_all_c1
            den_all_2_yrs[yr] = den_all_c2
            flx_all_2_yrs[yr] = flx_all_c2
            den_atl_a_yrs[yr] = den_atl_all
            flx_atl_a_yrs[yr] = flx_atl_all
            den_atl_1_yrs[yr] = den_atl_c1
            flx_atl_1_yrs[yr] = flx_atl_c1
            den_atl_2_yrs[yr] = den_atl_c2
            flx_atl_2_yrs[yr] = flx_atl_c2
            den_ind_a_yrs[yr] = den_ind_all
            flx_ind_a_yrs[yr] = flx_ind_all
            den_ind_1_yrs[yr] = den_ind_c1
            flx_ind_1_yrs[yr] = flx_ind_c1
            den_ind_2_yrs[yr] = den_ind_c2
            flx_ind_2_yrs[yr] = flx_ind_c2
            den_pac_a_yrs[yr] = den_pac_all
            flx_pac_a_yrs[yr] = flx_pac_all
            den_pac_1_yrs[yr] = den_pac_c1
            flx_pac_1_yrs[yr] = flx_pac_c1
            den_pac_2_yrs[yr] = den_pac_c2
            flx_pac_2_yrs[yr] = flx_pac_c2
            den_arc_a_yrs[yr] = den_arc_all
            flx_arc_a_yrs[yr] = flx_arc_all
            den_arc_1_yrs[yr] = den_arc_c1
            flx_arc_1_yrs[yr] = flx_arc_c1
            den_arc_2_yrs[yr] = den_arc_c2
            flx_arc_2_yrs[yr] = flx_arc_c2
            den_sou_a_yrs[yr] = den_sou_all
            flx_sou_a_yrs[yr] = flx_sou_all
            den_sou_1_yrs[yr] = den_sou_c1
            flx_sou_1_yrs[yr] = flx_sou_c1
            den_sou_2_yrs[yr] = den_sou_c2
            flx_sou_2_yrs[yr] = flx_sou_c2
        
        den_all_a_tot = pl.nansum(den_all_a_yrs,axis=0)
        den_all_1_tot = pl.nansum(den_all_1_yrs,axis=0)
        den_all_2_tot = pl.nansum(den_all_2_yrs,axis=0)
        den_atl_a_tot = pl.nansum(den_atl_a_yrs,axis=0)
        den_atl_1_tot = pl.nansum(den_atl_1_yrs,axis=0)
        den_atl_2_tot = pl.nansum(den_atl_2_yrs,axis=0)
        den_ind_a_tot = pl.nansum(den_ind_a_yrs,axis=0)
        den_ind_1_tot = pl.nansum(den_ind_1_yrs,axis=0)
        den_ind_2_tot = pl.nansum(den_ind_2_yrs,axis=0)
        den_pac_a_tot = pl.nansum(den_pac_a_yrs,axis=0)
        den_pac_1_tot = pl.nansum(den_pac_1_yrs,axis=0)
        den_pac_2_tot = pl.nansum(den_pac_2_yrs,axis=0)
        den_arc_a_tot = pl.nansum(den_arc_a_yrs,axis=0)
        den_arc_1_tot = pl.nansum(den_arc_1_yrs,axis=0)
        den_arc_2_tot = pl.nansum(den_arc_2_yrs,axis=0)
        den_sou_a_tot = pl.nansum(den_sou_a_yrs,axis=0)
        den_sou_1_tot = pl.nansum(den_sou_1_yrs,axis=0)
        den_sou_2_tot = pl.nansum(den_sou_2_yrs,axis=0)
        
        flx_all_a_mn = pl.nanmean(flx_all_a_yrs,axis=0)
        flx_all_1_mn = pl.nanmean(flx_all_1_yrs,axis=0)
        flx_all_2_mn = pl.nanmean(flx_all_2_yrs,axis=0)
        flx_atl_a_mn = pl.nanmean(flx_atl_a_yrs,axis=0)
        flx_atl_1_mn = pl.nanmean(flx_atl_1_yrs,axis=0)
        flx_atl_2_mn = pl.nanmean(flx_atl_2_yrs,axis=0)
        flx_ind_a_mn = pl.nanmean(flx_ind_a_yrs,axis=0)
        flx_ind_1_mn = pl.nanmean(flx_ind_1_yrs,axis=0)
        flx_ind_2_mn = pl.nanmean(flx_ind_2_yrs,axis=0)
        flx_pac_a_mn = pl.nanmean(flx_pac_a_yrs,axis=0)
        flx_pac_1_mn = pl.nanmean(flx_pac_1_yrs,axis=0)
        flx_pac_2_mn = pl.nanmean(flx_pac_2_yrs,axis=0)
        flx_arc_a_mn = pl.nanmean(flx_arc_a_yrs,axis=0)
        flx_arc_1_mn = pl.nanmean(flx_arc_1_yrs,axis=0)
        flx_arc_2_mn = pl.nanmean(flx_arc_2_yrs,axis=0)
        flx_sou_a_mn = pl.nanmean(flx_sou_a_yrs,axis=0)
        flx_sou_1_mn = pl.nanmean(flx_sou_1_yrs,axis=0)
        flx_sou_2_mn = pl.nanmean(flx_sou_2_yrs,axis=0)
        
        OpenNCfile(clusdir,shed,m,eralat,eralon)
        SavetoNC(den_all_a_tot,flx_all_a_mn,shed,'all','all',m)
        SavetoNC(den_all_1_tot,flx_all_1_mn,shed,'all','1',m)
        SavetoNC(den_all_2_tot,flx_all_2_mn,shed,'all','2',m)
        SavetoNC(den_atl_a_tot,flx_atl_a_mn,shed,'atl','all',m)
        SavetoNC(den_atl_1_tot,flx_atl_1_mn,shed,'atl','1',m)
        SavetoNC(den_atl_2_tot,flx_atl_2_mn,shed,'atl','2',m)
        SavetoNC(den_ind_a_tot,flx_ind_a_mn,shed,'ind','all',m)
        SavetoNC(den_ind_1_tot,flx_ind_1_mn,shed,'ind','1',m)
        SavetoNC(den_ind_2_tot,flx_ind_2_mn,shed,'ind','2',m)
        SavetoNC(den_pac_a_tot,flx_pac_a_mn,shed,'pac','all',m)
        SavetoNC(den_pac_1_tot,flx_pac_1_mn,shed,'pac','1',m)
        SavetoNC(den_pac_2_tot,flx_pac_2_mn,shed,'pac','2',m)
        SavetoNC(den_arc_a_tot,flx_arc_a_mn,shed,'arc','all',m)
        SavetoNC(den_arc_1_tot,flx_arc_1_mn,shed,'arc','1',m)
        SavetoNC(den_arc_2_tot,flx_arc_2_mn,shed,'arc','2',m)
        SavetoNC(den_sou_a_tot,flx_sou_a_mn,shed,'sou','all',m)
        SavetoNC(den_sou_1_tot,flx_sou_1_mn,shed,'sou','1',m)
        SavetoNC(den_sou_2_tot,flx_sou_2_mn,shed,'sou','2',m)
        
        print '*** ' + shed + ' ' + m + ' WRITTEN TO FILE*****'