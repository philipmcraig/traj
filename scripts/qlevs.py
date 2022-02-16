# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 16:43:03 2016

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset

def EtaToP(a,b,psurf):
    """
    """
    p0 = 10**5
    A = pl.zeros([len(a)-1]); B = pl.zeros_like(A)
    for i in range(len(A)):
        A[i] = 0.5*(a[i+1]+a[i])
        B[i] = 0.5*(b[i+1]+b[i])
    
    P = A*p0 + B*psurf
    
    return P

def BilinInterp(relpt,lon,lat,flux):
    """Function to interpolate to trajectory release point from 4 nearest grid
    points. e.g. https://en.wikipedia.org/wiki/Bilinear_interpolation
    
    Args:
        relpt (array): longitude & latitude co-ordinates of release point
        lon (array): longitude array from ERA-Interim
        lat (array): latitude array from ERA-Interim
        flux (array): variable from ERA-Interm requiring interpolation
    
    Returns:
        F (float): variable interpolated to release point
    """
    # First find p1,p2,p3,p4: the points forming a rectangle around relpt
    # Start with empty arrays for p1,p2,p3,p4; p = (x,y,flux):
    p1 = pl.zeros([3]); p2 = pl.zeros([3]); p3 = pl.zeros([3]); p4 = pl.zeros([3])
    # if release point longitude co-ordinate greater than max. ERA-Int longitude
    if relpt[0] > lon[-1]:
        a = -1 # take max ERA-Interim longitude as nearest longitude index
    else:
        a = NearestIndex(lon,relpt[0]) # nearest longitude index
    b = NearestIndex(lat,relpt[1]) # nearest latitude index
    
    if lon[a] < relpt[0]: # nearest lon west of relpt
        p1[0] = lon[a]; p3[0] = lon[a];  p2[0] = lon[a+1]; p4[0] = lon[a+1]
    elif lon[a] > relpt[0]: # nearest lon east of relpt
        p2[0] = lon[a]; p4[0] = lon[a]; p1[0] = lon[a-1]; p3[0] = lon[a-1]
        
    # does not take 0 meridian into account yet

    
    if lat[b] < relpt[1]: # nearest lat south of relpt
        p1[1] = lat[b]; p2[1] = lat[b]; p3[1] = lat[b-1]; p4[1] = lat[b-1]
    elif lat[b] > relpt[1]: # nearest lat north of relpt
        p3[1] = lat[b]; p4[1] = lat[b]; p1[1] = lat[b+1]; p2[1] = lat[b+1]
    
    # values of flux at p1,p2,p3,p4:
    nrth_lat = pl.where(lat==p3[1]); sth_lat = pl.where(lat==p1[1])
    west_lon = pl.where(lon==p1[0]); east_lon = pl.where(lon==p2[0])
    p1[2] = flux[sth_lat[0][0],west_lon[0][0]]
    p2[2] = flux[sth_lat[0][0],east_lon[0][0]]
    p3[2] = flux[nrth_lat[0][0],west_lon[0][0]]
    p4[2] = flux[nrth_lat[0][0],east_lon[0][0]]
    
    # if release point longitude co-ordinate greater than max. ERA-Int longitude
    if relpt[0] > lon[-1]:
        # dx is 360 + min ERA-Int longitude minus max ERA-Int longitude
        dx = (360. + lon[0]) - lon[-1]
    else:
        dx = p2[0] - p1[0] # dx is difference between 2 nearest longitudes
    dy = p3[1] - p2[1] # dy is difference between 2 nearest latitudes
    
    # if release point longitude co-ordinate greater than max. ERA-Int longitude
    if relpt[0] > lon[-1]:
        # need 360 + p2 longitude
        f1 = (((360+p2[0])-relpt[0])/dx)*p1[2] + ((relpt[0]-p1[0])/dx)*p2[2]
        f2 = (((360+p2[0])-relpt[0])/dx)*p3[2] + ((relpt[0]-p1[0])/dx)*p4[2]
    else:
        f1 = ((p2[0]-relpt[0])/dx)*p1[2] + ((relpt[0]-p1[0])/dx)*p2[2]
        f2 = ((p2[0]-relpt[0])/dx)*p3[2] + ((relpt[0]-p1[0])/dx)*p4[2]
    
    F = ((p3[1]-relpt[1])/dy)*f1 + ((relpt[1]-p2[1])/dy)*f2
    
    return F

def Add360(rlspts):
    """Function to add 360° to longitude co-ordinates below zero.
    
    Args:
        rlspts (array): longitude & latitude co-ordinates (degrees) of trajectory
                        release points, some longitude co-ordinates may be negative
    
    Returns:
        rlspts (array): longitude & latitude co-ordinates (degrees) of trajectory
                        release points, all longitude co-ordinates are positive
    """
    for i in range(rlspts.shape[0]): # loop over release points
        if rlspts[i,0] < 0.: # check for negative co-ordinates
            rlspts[i,0] = rlspts[i,0] + 360. # add 360. if co-ordinate < 0
    
    return rlspts

pl.close('all')
exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())


trajdir = '/home/np838619/Trajectory/'
clusdir = '/glusterfs/scenario/users/np838619/tjnew/TESTS/'
#eradir = '/panfs/jasmin/era/era-in/netc/'
sheddir = '/home/np838619/Watershed/shed_defs/'
file_plev = 'ggap'; file_surf = 'ggas'
year = '2014'; mnth_name = 'jul'

mnth_no = '08'; day = '31'; time = '0000'
fname1 = file_plev + year + mnth_no + day + time
fname2 = file_surf + year + mnth_no + day + time

nc1 = Dataset(clusdir+fname2+'.nc','r')
psurf = nc1.variables['SP'][:]
nc1.close()

nc2 = Dataset(clusdir+fname1+'.nc','r')
lat = nc2.variables['latitude'][:]
lon = nc2.variables['longitude'][:]
pres = nc2.variables['p'][:]
q = nc2.variables['Q'][:]
temp = nc2.variables['T'][:]; PV = nc2.variables['PV'][:]
u = nc2.variables['U'][:]; v = nc2.variables['V'][:]
nc2.close()

# remove 1D axes:
pres = pl.squeeze(pres); q = pl.squeeze(q); psurf = pl.squeeze(psurf)
temp = pl.squeeze(temp); PV = pl.squeeze(PV); u = pl.squeeze(u); v = pl.squeeze(v)

rlspts = pl.genfromtxt(sheddir+'ArP_traj_release_new.txt',skip_header=5)
rlspts = Add360(rlspts)
#loc = (360.-108.98438,42.456036) #-113.90625 -108.98438 42.456036
locinds = (NearestIndex(lon,rlspts[0,0]),NearestIndex(lat,rlspts[0,1]))
NREL = rlspts.shape[0]

uv_tj = pl.genfromtxt(clusdir+'uv_out.txt')
u_tj = uv_tj[:,0]; u_tj = pl.reshape(u_tj,(17,NREL))
v_tj = uv_tj[:,1]; v_tj = pl.reshape(v_tj,(17,NREL))
pres_tj = pl.genfromtxt(clusdir+'p_out.txt')
pres_tj = pl.reshape(pres_tj,(17,NREL))
sp_tj = pl.genfromtxt(clusdir+'sp_out.txt')
q_tj = pl.genfromtxt(clusdir+'qout.txt')
q_tj = pl.reshape(q_tj,(17,NREL))
pv_tj = pl.genfromtxt(clusdir+'pvout.txt')
pv_tj = pl.reshape(pv_tj,(17,NREL))
temp_tj = pl.genfromtxt(clusdir+'tempout.txt')
temp_tj = pl.reshape(temp_tj,(17,NREL))

ut_mn = pl.mean(u_tj,axis=1); vt_mn = pl.mean(v_tj,axis=1)
qt_mn = pl.mean(q_tj,axis=1); pvt_mn = pl.mean(pv_tj,axis=1)
tt_mn = pl.mean(temp_tj,axis=1); prt_mn = pl.mean(pres_tj,axis=1)

u_era = pl.zeros([pres.size,NREL]); v_era = pl.zeros_like(u_era)
sp_era = pl.zeros([NREL]); q_era = pl.zeros_like(u_era)
pv_era = pl.zeros_like(u_era); temp_era = pl.zeros_like(u_era)
for pt in range(NREL):
    sp_era[pt] = BilinInterp(rlspts[pt],lon,lat,psurf)
    for lev in range(pres.size):
        u_era[lev,pt] = BilinInterp(rlspts[pt],lon,lat,u[lev])
        v_era[lev,pt] = BilinInterp(rlspts[pt],lon,lat,v[lev])
        q_era[lev,pt] = BilinInterp(rlspts[pt],lon,lat,q[lev])
        pv_era[lev,pt] = BilinInterp(rlspts[pt],lon,lat,PV[lev])
        temp_era[lev,pt] = BilinInterp(rlspts[pt],lon,lat,temp[lev])

ue_mn = pl.mean(u_era,axis=1); ve_mn = pl.mean(v_era,axis=1)
qe_mn = pl.mean(q_era,axis=1); pve_mn = pl.mean(pv_era,axis=1)
te_mn = pl.mean(temp_era,axis=1)

coeffs = pl.genfromtxt(trajdir + 'coeffs.txt')
a = coeffs[:,1]; b = coeffs[:,2]

pfull = EtaToP(a,b,sp_era[0])
trajlevs = pl.linspace(0.95,0.15,17)

#for i in range(NREL):
#    pl.close('all')
#    
#    pl.figure(1)
#    pl.plot(u_tj[:,i],pres_tj[:,i],label='$u_{traj}$',color='b',zorder=3,linewidth=2.)
#    pl.plot(u_era[:,i],pres,label='$u_{ERA}$',color='b',ls='--')
#    pl.plot(v_tj[:,i],pres_tj[:,i],label='$v_{traj}$',color='g',zorder=3,linewidth=2.)
#    pl.plot(v_era[:,i],pres,label='$v_{ERA}$',color='g',ls='--')
#    pl.xlabel('m/s',fontsize=26); pl.ylabel('hPa',fontsize=26)
#    pl.ylim(1100,0); pl.legend(loc=0)
#    pl.savefig(clusdir+'test3/attr_checks/uv_pts/uv_rel'+str(i+1)+'.png')
#    
#    pl.figure(2)
#    pl.plot(q_tj[:,i],pres_tj[:,i],label='$q_{traj}$',linewidth=2,color='b')
#    pl.plot(q_era[:,i],pres,label='$q_{ERA}$',linewidth=2.,color='b',ls='--')
#    pl.xlabel('kg/kg',fontsize=26); pl.ylabel('hPa',fontsize=26)
#    pl.ylim(1100,0); pl.legend(loc=0)
#    pl.savefig(clusdir+'test3/attr_checks/q_pts/q_rel'+str(i+1)+'.png')
#    
#    pl.figure(3)
#    pl.plot(pv_tj[:,i],pres_tj[:,i],label='$PV_{traj}$',linewidth=2,color='b')
#    pl.plot(pv_era[:,i]*(10**6),pres,label='$PV_{ERA}$',linewidth=2.,color='b',ls='--')
#    pl.xlabel('PVU',fontsize=26); pl.ylabel('hPa',fontsize=26)
#    pl.ylim(1100,100); pl.xlim(-2,2); pl.legend(loc=0)
#    pl.savefig(clusdir+'test3/attr_checks/pv_pts/pv_rel'+str(i+1)+'.png')
#    
#    pl.figure(4)
#    pl.plot(temp_tj[:,i],pres_tj[:,i],label='$T_{traj}$',linewidth=2,color='b')
#    pl.plot(temp_era[:,i],pres,label='$T_{ERA}$',linewidth=2.,color='b',ls='--')
#    pl.xlabel('K',fontsize=26); pl.ylabel('hPa',fontsize=26)
#    pl.ylim(1100,0); pl.legend(loc=0)
#    pl.savefig(clusdir+'test3/attr_checks/temp_pts/temp_rel'+str(i+1)+'.png')