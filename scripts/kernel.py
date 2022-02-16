# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 08:46:36 2017

@author: np838619
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap
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

def Haver(pt1,pt2):
    """Function to calculate distance on sphere using the Haversine formula
    e.g. https://en.wikipedia.org/wiki/Haversine_formula
    
    Args:
        pt1 (array): longitude, latitude co-ordinates of a trajectory release point
        pt2 (array): longitude, latitude co-ordinates of a trajectory release point
    
    Returns:
        d (float): distance between two points on the sphere in metres
    """
    # convert both point to radians:
    pt1 = pl.radians(pt1); pt2 = pl.radians(pt2)
    
    # get dphi and dlamda
    phi1 = pt1[1]; phi2 = pt2[1]; dphi = pl.absolute(phi2 - phi1)
    lam1 = pt1[0]; lam2 = pt2[0]; dlam = pl.absolute(lam2 - lam1)
    
    # calculate dsigma:
    dsig = 2*pl.arcsin(pl.sqrt((pl.sin(dphi/2)**2 + pl.cos(phi1)*pl.cos(phi2)*pl.sin(dlam/2)**2)))
    
    a = 6.37e6 # radius of Earth (m)
    
    d = a*dsig # calculate distance
    
    return d

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

pl.close('all')
R = 1000000 # kernel width
a = 6.371*10**6 # radius of earth
lat = 87; lon = 182 # lat, lon co-ordinates of a point
phi = pl.radians(lat); lam = pl.radians(lon) # convert point to radians
xi = pl.array([lon,lat,a]) # kernel location
#x = pl.array([48,44,a])
x1 = pl.array([55,-6,a])
x2 = pl.array([58,-9,a])
x3 = pl.array([63,0,a])
xi = pl.array([x1,x2,x3])

nc1 = Dataset('/net/glusterfs_essc/scenario/users/np838619/tjnew/TESTS/'+'ggas201408310000.nc','r')
eralon = nc1.variables['longitude'][:]
eralat = nc1.variables['latitude'][:]
nc1.close()

lonx = pl.zeros([3*eralon.size])
#lonx[:27] = -1*pl.flipud(eralon[1:28]); lonx[27:27+eralon.size] = eralon[:]
#lonx[27+eralon.size:] = eralon[1:28]+360
lonx[:eralon.size] = -1*pl.flipud(eralon); lonx[eralon.size:2*eralon.size] = eralon
lonx[2*eralon.size:] = eralon + 360

#latx = pl.zeros([eralat.size+2])
#latx[0] = 90.; latx[1:-1] = eralat; latx[-1] = -90.

eraphi = pl.radians(eralat); eralam = pl.radians(lonx) # lat,lon in radians
# Nearest grid points to point:
nearest = pl.zeros([xi.shape[0],2])
for pt in range(xi.shape[0]):
    nearest[pt,0] = NearestIndex(lonx,xi[pt,0])
    nearest[pt,1] = NearestIndex(eralat,xi[pt,1])
dlam = eralam[nearest[0,0]] - eralam[nearest[0,0]-1]
dphi = eraphi[nearest[0,1]-1] - eraphi[nearest[0,1]]
# number of grid boxes in kernel:
di = R/(2*dlam*a*pl.cos(phi)); dj = R/(2*dphi*a)

xi_rad = pl.radians(xi); xi_rad[:,2] = a # kernel location in radians
#x_rad = pl.array([pl.radians(x[0]),pl.radians(x[1]),a])
#theta = pl.arccos(SphereDot(xi_rad,xi))
#r = a*theta
#if r > R:
#    r = 0.

Ki = pl.zeros([xi.shape[0],eralat.size,lonx.size])

for pt in range(xi.shape[0]):
    lamloc = nearest[pt,0]; philoc = nearest[pt,1]
    # patch of ERA-Interim grid
    q = 10; s = 12 # number of extra points in lat,lon directions to plot entire kernel
    # longitude, latitude co-ordinates for patch in radians:
    if philoc < 18.:
        newlam = eralam
        newphi = eraphi[:18+philoc]
    elif philoc > 238.:
        newlam = eralam
        newphi = eraphi[238-(256-philoc):]
    else:    
        newlam = eralam[lamloc-round(di)-s:lamloc+round(di)+s]
        newphi = eraphi[philoc-round(dj)-q:philoc+round(dj)+q]

    theta_grid = pl.zeros([newlam.size,newphi.size]) # empty array for thetas
    K = pl.zeros_like(theta_grid) # empty array for values of kernel
    for i in range(newlam.size): # loop over lon points
        for j in range(newphi.size): # loop over lat points
            x_grid = (newlam[i],newphi[j],a) # grid point in radians
            theta_grid[i,j] = pl.arccos(SphereDot(xi_rad[pt],x_grid)) # theta
            r = a*theta_grid[i,j] # value of r
            if r > R: # check if r in kernel
                r = 0.
                K[i,j] = 0.
            else:
                K[i,j] = (2/pl.pi*(R**2))*(1-(r**2)/(R**2)) # kernel function
    Ki[pt,philoc-round(dj)-q:philoc+round(dj)+q,lamloc-round(di)-s:lamloc+round(di)+s] = K.T

Ks = pl.sum(Ki,axis=0)

# array of nan as nan plots white, same size as ERA-Interim grid:
#D = pl.zeros([eralat.size,lonx.size]); D[:] = pl.float64('nan')
## put kernel in right place:
#if philoc < 18.:
#    D[:18+philoc] = K.T
#elif philoc > 238.:
#    D[238-(256-philoc):] = K.T
#else:
#    D[philoc-round(dj)-q:philoc+round(dj)+q,lamloc-round(di)-s:lamloc+round(di)+s] = K.T

#m = Basemap(projection='robin',lon_0=0.)#,boundinglat=20)
#lons,lats = pl.meshgrid(lonx,eralat)
#X,Y = m(lons,lats)
#m.drawcoastlines()
#levels = pl.linspace(0,0.64,10)
#m.contourf(X,Y,D/(10e11),levels=levels); m.plot(0,90,latlon=True,marker='x')
#m.colorbar()
#m.drawmeridians(eralon); m.drawparallels(eralat)

#ax = pl.axes(projection=ccrs.NorthPolarStereo(central_longitude=0))
#ax.coastlines()
#lons,lats = pl.meshgrid(lonx,eralat)
#levels=pl.linspace(0,0.63,10)
#norm = pl.Normalize(0,0.63,clip=False)
#cb = ax.contourf(lons, lats, D/(10e11),transform=ccrs.PlateCarree(),norm=norm,
#                 extend='max',levels=levels)
#pl.colorbar(cb)
#ax.set_extent([-180, 180, 90, 60], ccrs.PlateCarree())
#pl.tight_layout()
#ax.plot(0,-90,color='r',marker='x',transform=ccrs.Geodetic())
#c = pl.matplotlib.patches.Circle((xi[0],xi[1]),radius=R,fill=False,lw=2.,color='r',transform=ccrs.Geodetic())
#c = pl.Circle((9000000,45),R)
#ax.add_artist(c)

#r = a*pl.radians(phi)

#empty = [] # list for indices of empty arrays
#for rls in range(len(atlcount)): # loop over no. of releases
#    if atlcount[rls].size == 0.: # if no traj have origin at this release
#        empty.append(rls) # ... add this release index to list
#
#if len(empty) == 0.: # if list has zero length ...
#    B = pl.vstack(atlcount)
#else: # otherwise ...
#    B = pl.delete(atlcount,empty) # . delete the empty arrays
#    B = pl.vstack(atlcount)
#
#origins = pl.zeros([B.shape[0],3])
#origins[:,:2] = B[:,3:]; origins[:,-1] = a
#
##K = (2/(pl.pi*R**2))*(1-(r**2)/(R**2))
#K = pl.zeros([origins.shape[0],eralat.shape[0],eralon.shape[0]])

#for k in range(origins.shape[0]):
#    xk = origins[k]
#    ni = NearestIndex(eralat,origins[k,0]); nj = NearestIndex(eralon,origins[k,1])
#    for i in range(ni-1,ni+2):
#        if nj == 511:
#            p = (nj-1,nj,0)
#            for j in range(len(p)):
#                if p[j] == 0:
#                    xk = pl.array([eralat[i],eralon[p[j]]+360,a])
#                    f = (origins[k,1],origins[k,0]); g = (xk[1],xk[0])
#                    r = Haver(f,g); K[k,i,j] = (2/(pl.pi*R**2))*(1-(r**2/R**2))
#                else:
#                    xk = pl.array([eralat[i],eralon[p[j]],a])
#                    f = (origins[k,1],origins[k,0]); g = (xk[1],xk[0])
#                    r = Haver(f,g); K[k,i,j] = (2/(pl.pi*R**2))*(1-(r**2/R**2))
#        else:
#            for j in range(nj-1,nj+2):
#                xk = pl.array([eralat[i],eralon[j],a])
#                f = (origins[k,1],origins[k,0]); g = (xk[1],xk[0])
#                r = Haver(f,g)
#                #r = pl.arccos(pl.dot(origins[k],xk)/a**2)
#                K[k,i,j] = (2/(pl.pi*R**2))*(1-(r**2/R**2))