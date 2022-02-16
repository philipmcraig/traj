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
import timeit

def SphereDot(trajor,x):
    """Dot product in spherical co-ordinates
    
    Args:
        trajor (array-like): longitude, latitude co-ordinates of trajectory
                             origin (kernel centre) in radians
        x (array-like): longitude, latitude co-ordinates of grid point
    
    Returns:
        dotprod (float): value of scalar product in spherical co-ordinates
                         between trajectory origin and grid point
    """
    phi1 = trajor[1]; phi2 = x[1] # latitude co-ordinates
    lam1 = trajor[0]; lam2 = x[0] # longitude co-ordinates
    # dot product in spherical co-ordinates:
    dotprod = pl.cos(phi1)*pl.cos(phi2)*(pl.sin(lam1)*pl.sin(lam2) + 
                            pl.cos(lam1)*pl.cos(lam2)) + pl.sin(phi1)*pl.sin(phi2)
    
    return dotprod

def HalfGrid2(lat):
    """Function to create latitude half-grid
    
    Args:
        lat (array): latitude array in radians
    
    Returns:
        lath (array): half-grid latitude array in radians
    """
    # set up empty array, size one more than lat
    lath = pl.zeros([lat.shape[0]+1])
    # 
    lath[0] = lat[0]+0.5*(lat[0]-lat[1]); lath[-1] = lat[-1] - 0.5*(lat[-2]-lat[-1])
    # loop over lat_half from index 1 to index -2:
    for i in range(1,lath.shape[0]-1): # loops over 256
        lath[i] = 0.5*(lat[i]+lat[i-1])
    
    return lath

def KernelCalc(newlam,newphi,point,a,R):
    """Kernel calculation.
    
    Args:
        newlam (array): subset of ERA-Interim longitude to calculate kernel
        newphi (array): subset of ERA-Interim latitude to calculate kernel
        point (array-like): longitude, latitude co-ordinates of trajectory
                            origin (kernel centre) in radians
        a (float): radius of Earth in metres
        R (float): radius of kernel in metres
    
    Returns:
        K (array): values of kernel function inside subset of ERA-Interim grid
    """
    # make empty arrays for theta & K:
    theta = pl.zeros([newlam.size,newphi.size])
    K = pl.zeros_like(theta)#; N = pl.zeros_like(K)
    
#    lat_half = HalfGrid2(newphi)
#    delta_lambda = newlam[1] - newlam[0]
    
    for i in range(newlam.size): # loop over longitude subset
        for j in range(newphi.size): # loop over latitude subset
            #dmu = pl.sin(lat_half[j]) - pl.sin(lat_half[j+1])
            #delta_phi = lat_half[j] - lat_half[j+1]
            x_grid = (newlam[i],newphi[j],a) # point from subset of ERA-I grid
            # calculate theta using dot product in spherical co-ordinates:
            theta[i,j] = pl.arccos(SphereDot(point,x_grid))
            r = a*theta[i,j] # value of r
            if r > R: # check if r in kernel
                r = 0.; K[i,j] = 0. # set r & K as zero if outside kernel
            else: # otherwise calculate value of kernel
                N = 1/(2*pl.pi*(1-(2*a/R)*pl.sin(R/a) + (2*a**2/R**2)*(1-pl.cos(R/a))))
                K[i,j] = N*(1-(r**2)/(R**2))#*dphi*delta_lambda # kernel function
    
    return K

def WriteNetCDF(field,lon,lat,trajfile):
    """
    """
    # open nc file, need to input a string to function for filename
    
    #lat,lon dimensions
    #create lat,lon variables
    # lat,lon units
    # lat,lon = arrays
    
    #time dimensions, create time variable, time units
    
    # create field variable, field units, = array
    
    # close nc file

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

pl.close('all')
R = 200000 # kernel radius
a = 6.371*10**6 # radius of earth

NREL = rlslabs.shape[0]
dl = pl.zeros([NREL]) # empty arrays for distance weights
for pt in range(1,NREL-1): # exclude first & last release points
    # use Haversine formula to calculate dl
    dl[pt] = 0.5*(Haversine(rlslabs[pt-1,:2],rlslabs[pt,:2]) + 
                    Haversine(rlslabs[pt,:2],rlslabs[pt+1,:2]))
# use first 2 points & last 2 points for first & last dl
dl[0] = Haversine(rlslabs[0,:2],rlslabs[1,:2])
dl[-1] = Haversine(rlslabs[-2,:2],rlslabs[-1,:2])

sfx = sf3*dl#; sfx[:] = 1

# ERA-Interim longitude and latitude:
nc1 = Dataset('/net/glusterfs_essc/scenario/users/np838619/tjnew/TESTS/'+'ggas201408310000.nc','r')
eralon = nc1.variables['longitude'][:]
eralat = nc1.variables['latitude'][:]
nc1.close()

# To deal with points near zero longitude, extend the longitude array
# extended longitude array ranges from -320 to 380
lonx = pl.zeros([eralon.size+56])
lonx[:28] = -1*pl.flipud(eralon[1:29]); lonx[28:540] = eralon;
lonx[540:] = eralon[:28] + 360
eraphi = pl.radians(eralat); eralam = pl.radians(lonx) # lat,lon in radians

start_time = timeit.default_timer()
t = 0

CAT = 'all'
catchcount = newlabs # origin catchment
NTRAJ = 0. # number of trajectories for normalization, start with zero
Kz = pl.zeros([len(catchcount),256,512]) # empty Kernel array for all releases
F = pl.zeros([len(catchcount),eraphi.size,eralam.size])
for rls in range(0,len(catchcount)): # loop over trajectory releases
    #NTRAJ = NTRAJ + len(catchcount[rls])
    if catchcount[rls].size == 0.: # if no origins in release...
        continue # ... skip to next iteration in loop
    # extract origin point from catchcount array, flip so lon-lat:
    points = catchcount[rls][:,2:]; points = pl.fliplr(points)

    # empty array for Kernel for this release for all points
    Ki = pl.zeros([points.shape[0],eralat.size,lonx.size])
    cross = sfx[rls].flatten()#; H = pl.zeros_like(sfx[0])#; f2 = pl.zeros_like(cross)

    for pt in range(points.shape[0]): # loop over points
        if CAT == 'all':
            if catchcount[rls][pt,0] == -1:
                continue
            elif catchcount[rls][pt,0] == 0.:
                continue
            else:
                pass
        elif catchcount[rls][pt,0] != int(CAT):
            continue
#        elif catchcount[rls][pt,0] == -1.:
#            continue
        NTRAJ = NTRAJ + 1
        trajno = pt#catchcount[rls][pt,0]#; tj.append(trajno)
        H = pl.zeros_like(F[rls])
        
        flx = (1/9.81)*cross[trajno]#; f2[trajno] = (1/9.81)*cross[trajno]
        lamloc = NearestIndex(lonx,points[pt,0]) # nearest longitude grid point
        philoc = NearestIndex(eralat,points[pt,1]) # nearest latitude grid point
        dlam = eralam[lamloc] - eralam[lamloc-1] # longitude spacing
        dphi = eraphi[philoc-1] - eraphi[philoc] # latitude spacing
#        H[philoc,lamloc] = flx
#        F[rls] = F[rls] + H

        point_rad = pl.zeros([3])
        point_rad[:2] = pl.radians(points[pt]); point_rad[-1] = a
        
        phi = point_rad[1]
        # number of grid boxes in kernel:
        di = R/(2*dlam*a*pl.cos(phi)); dj = R/(2*dphi*a)
        
        # patch of ERA-Interim grid
        if philoc < 3.: # if close to north pole
            newlam = eralam # use entire longitude array
            newphi = eraphi[:4+philoc] # use first 4 + philoc latitude points
            K = KernelCalc(newlam,newphi,point_rad,a,R) # calculate Kernel
            Ki[pt,:4+philoc] = K.T # stick into Ki
            H[:4+philoc] = flx*K.T
        elif philoc > 253.: # if close to south pole
            newlam = eralam # use entire longitude array
            newphi = eraphi[philoc-2:] # use last philoc-2 latitude points 
            K = KernelCalc(newlam,newphi,point_rad,a,R) # calculate Kernel
            Ki[pt,philoc-2:] = K.T # stick into Ki
            H[philoc-2:] = flx*K.T
        else: # otherwise...
            # extract a subset of the ERA-Interim grid around the kernel centre
            # 3*di & d*dj give enough points to fit kernel inside grid subset
            newlam = eralam[lamloc-3*round(di):lamloc+3*round(di)+1]
            newphi = eraphi[philoc-4*round(dj)-1:philoc+4*round(dj)]
            K = KernelCalc(newlam,newphi,point_rad,a,R) # calculate Kernel
            # stick into Ki:
            Ki[pt,philoc-4*round(dj)-1:philoc+4*round(dj),lamloc-3*round(di):lamloc+3*round(di)+1] = K.T
            H[philoc-4*round(dj)-1:philoc+4*round(dj),lamloc-3*round(di):lamloc+3*round(di)+1] = flx*K.T
    
        F[rls] = F[rls] + H
    Ks = pl.nansum(Ki,axis=0) # sum over all points for this trajectory release
    # make it for the ERA-Interim grid:
    Kx = pl.zeros([256,512])
    Kx[:,:28] = Ks[:,540:]; Kx[:,484:] = Ks[:,:28]; Kx = Kx + Ks[:,28:540]
    Kz[rls] = Kx # stick into Kz

Fx = pl.zeros([len(catchcount),256,512])
Fx[:,:,:28] = F[:,:,540:]; Fx[:,:,484:] = F[:,:,:28]; Fx = Fx + F[:,:,28:540]
#Fx = Fx/(10**9)
Ktot = pl.nansum(Kz,axis=0) # sum over all trajectory releases

elapsed = timeit.default_timer() - start_time
print elapsed


lon2 = eralon; lon2[-1] = 360 # change final longitude point to avoid plotting discontinuity
#m = Basemap(projection='robin',lon_0=0.)#,boundinglat=20)
#lons,lats = pl.meshgrid(lon2,eralat)
#X,Y = m(lons,lats)
#m.drawcoastlines()
##levels = pl.linspace(0.01,37,20)
##norm = pl.Normalize(0.01,37,clip=False)
#cb=m.contourf(X,Y,Ktot/(10e11))
#m.colorbar(cb)
#m.drawmeridians(eralon); m.drawparallels(eralat)

eralam = pl.radians(eralon)
lat_half = HalfGrid(eraphi)
nlon = eralon.shape[0] # number of longitude points
delta_lambda = (2*pl.pi)/nlon
D = pl.zeros([eraphi.size,eralam.size]); T = pl.zeros([rls+1,eraphi.size,eralam.size])
# loop over latitude and longitude
for i in range(lat_half.size-1): # loops over 256
    #latpair = (lat_half[i],lat_half[i+1])
    dmu = pl.sin(lat_half[i]) - pl.sin(lat_half[i+1])
    delta_phi = lat_half[i] - lat_half[i+1] 
    for j in range(eralon.size): # loops over 512
        D[i,j] = Ktot[i,j]*dmu*delta_lambda#/NTRAJ

for r in range(rls+1):
    for i in range(lat_half.size-1): # loops over 256
    #latpair = (lat_half[i],lat_half[i+1])
        dmu = pl.sin(lat_half[i]) - pl.sin(lat_half[i+1])
        for j in range(eralon.size): # loops over 512
            T[r,i,j] = Fx[r,i,j]*dmu*delta_lambda

print pl.sum(D)/NTRAJ
Tmn = pl.nanmean(T[:rls+1],axis=0)


pl.figure(1)
ax = pl.axes(projection=ccrs.Robinson(central_longitude=0))
ax.coastlines()
lons,lats = pl.meshgrid(lon2,eralat)
levels=pl.linspace(0.00001,0.00035,10)
norm = pl.Normalize(0.00001,0.00035,clip=False)
cb = ax.contourf(lons, lats, D/NTRAJ,transform=ccrs.PlateCarree(),norm=norm,
                 levels=levels,extend='max')
pl.colorbar(cb)
pl.tight_layout()

#B = (NTRAJ/(pl.sum(D)))*pl.nansum(Kz*Fx,axis=0)

pl.figure(2)
ax2 = pl.axes(projection=ccrs.Robinson(central_longitude=0))
ax2.coastlines()
lons,lats = pl.meshgrid(lon2,eralat)
v1 = -1.0e-4; v2 = 1.0e-4
levels=pl.linspace(v1,v2,20)
norm = pl.Normalize(v1,v2,clip=False)#MidpointNormalize(midpoint=0)#
cb = ax2.contourf(lons, lats, Tmn/(10**9),transform=ccrs.PlateCarree(),norm=norm,
                 extend='both',cmap='seismic',levels=levels)
pl.colorbar(cb)
ax2.plot(rlslabs[:,0],rlslabs[:,1],transform=ccrs.Geodetic(),color='k')
#ax.set_extent([-180, 180, 90, 45], ccrs.PlateCarree())
pl.tight_layout()
#ax.plot(0,-90,color='r',marker='x',transform=ccrs.Geodetic())