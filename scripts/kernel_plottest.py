# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 08:54:52 2017

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
import cartopy
import cartopy.crs as ccrs

def SphereDot(trajor,x):
    """Dot product in spherical co-ordinates
    """
    phi1 = trajor[1]; phi2 = x[1]
    lam1 = trajor[0]; lam2 = x[0]
    dotprod = pl.cos(phi1)*pl.cos(phi2)*(pl.sin(lam1)*pl.sin(lam2) + 
                            pl.cos(lam1)*pl.cos(lam2)) + pl.sin(phi1)*pl.sin(phi2)
    
    return dotprod

def KernelCalc(newlam,newphi,point,a,R):
    """Kernel calculation.
    """
    theta_grid = pl.zeros([newlam.size,newphi.size])
    K = pl.zeros_like(theta_grid)
    
    for i in range(newlam.size):
        for j in range(newphi.size):
            x_grid = (newlam[i],newphi[j],a)
            theta_grid[i,j] = pl.arccos(SphereDot(point,x_grid))
            r = a*theta_grid[i,j] # value of r
            if r > R: # check if r in kernel
                r = 0.; K[i,j] = 0.
            else:
                K[i,j] = (2/pl.pi*(R**2))*(1-(r**2)/(R**2)) # kernel function
    
    return K

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

pl.close('all')
R = 200000 # kernel width
a = 6.371*10**6 # radius of earth
lat = 88.9; lon = 359 # lat, lon co-ordinates of a point
phi = pl.radians(lat); lam = pl.radians(lon) # convert point to radians
xi = pl.array([lon,lat,a]) # kernel location
#x = pl.array([48,44,a])

nc1 = Dataset('/net/glusterfs_essc/scenario/users/np838619/tjnew/TESTS/'+'ggas201408310000.nc','r')
eralon = nc1.variables['longitude'][:]
eralat = nc1.variables['latitude'][:]
nc1.close()

lonx = pl.zeros([eralon.size+56])
lonx[:28] = -1*pl.flipud(eralon[1:29]); lonx[28:540] = eralon;
lonx[540:] = eralon[:28] + 360
#lonx[:eralon.size] = -1*pl.flipud(eralon); lonx[eralon.size:2*eralon.size] = eralon
#lonx[2*eralon.size:] = eralon + 360
eraphi = pl.radians(eralat); eralam = pl.radians(lonx) # lat,lon in radians

lamloc = NearestIndex(lonx,xi[0]) # nearest longitude grid point
philoc = NearestIndex(eralat,xi[1]) # nearest latitude grid point
dlam = eralam[lamloc] - eralam[lamloc-1] # longitude spacing
dphi = eraphi[philoc-1] - eraphi[philoc] # latitude spacing

di = R/(2*dlam*a*pl.cos(phi)); dj = R/(2*dphi*a)

xi_rad = pl.zeros_like(xi)
xi_rad[:2] = pl.radians(xi[:2]); xi_rad[-1] = a

q = 3; s = 2
newlam = eralam#[lamloc-3*round(di):lamloc+3*round(di)+1]
newphi = eraphi[:philoc+4]#[philoc-round(dj)-q+1:philoc+round(dj)+q]

K = KernelCalc(newlam,newphi,xi_rad,a,R)
#pl.imshow(K)

D = pl.zeros([eralat.size,lonx.size]); D[:] = pl.float64('nan')
#D[philoc-4*round(dj)+1:philoc+4*round(dj),lamloc-3*round(di):lamloc+3*round(di)+1] = K.T
D[:philoc+4,:] = K.T

ax = pl.axes(projection=ccrs.NorthPolarStereo(central_longitude=20))
ax.coastlines()
lons,lats = pl.meshgrid(lonx,eralat)
#levels=pl.linspace(0.01,40,20)
norm = pl.Normalize(0.001,0.024,clip=False)
cb = ax.contourf(lons, lats, D/(10e11),transform=ccrs.PlateCarree(),norm=norm,
                 extend='max')
pl.colorbar(cb)
ax.set_extent([-180, 180, 90, 80], ccrs.PlateCarree())
pl.tight_layout()