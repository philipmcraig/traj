# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 20:42:49 2017

@author: np838619
"""

R = 200000 # kernel radius
a = 6.371*10**6 # radius of earth

NRLS = newlabs.shape[0] # number of trajectory releases
NREL = rlslabs.shape[0] # number of release points
dl = pl.zeros([NREL]) # empty arrays for distance weights
for pt in range(1,NREL-1): # exclude first & last release points
    # use Haversine formula to calculate dl
    dl[pt] = 0.5*(Haversine(rlslabs[pt-1,:2],rlslabs[pt,:2]) + 
                    Haversine(rlslabs[pt,:2],rlslabs[pt+1,:2]))
# use first 2 points & last 2 points for first & last dl
dl[0] = Haversine(rlslabs[0,:2],rlslabs[1,:2])
dl[-1] = Haversine(rlslabs[-2,:2],rlslabs[-1,:2])

sfx = sf3*dl # vertical cross-section of moisture flux, weighted by dl

# To deal with points near zero longitude, extend the longitude array
# extended longitude array ranges from -320 to 380
lonx = pl.zeros([eralon.size+56])
lonx[:28] = -1*pl.flipud(eralon[1:29]); lonx[28:540] = eralon;
lonx[540:] = eralon[:28] + 360

eraphi = pl.radians(eralat); eralam = pl.radians(lonx) # lat,lon in radians

# empty Kernel array for all releases:
Kz = pl.zeros([NRLS,eralat.size,eralon.size])
F = pl.zeros([NRLS,eraphi.size,eralam.size])
NTRAJ = 0
catchcount = newlabs
for rls in range(NRLS): # loop over trajectory releases
    if catchcount[rls].size == 0.: # if no origins in release...
        continue # ... skip to next iteration in loop
    # extract origin point from catchcount array, flip so lon-lat:
    points = catchcount[rls][:,2:]; points = pl.fliplr(points)

    # empty array for Kernel for this release for all points
    Ki = pl.zeros([points.shape[0],eralat.size,lonx.size])
    cross = sfx[rls].flatten() # flatten cross-section so 1D

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
    
        NTRAJ = NTRAJ + 1
        trajno = pt
        H = pl.zeros_like(F[rls])
        
        flx = (1/9.81)*cross[trajno]
        lamloc = NearestIndex(lonx,points[pt,0]) # nearest longitude grid point
        philoc = NearestIndex(eralat,points[pt,1]) # nearest latitude grid point
        dlam = eralam[lamloc] - eralam[lamloc-1] # longitude spacing
        dphi = eraphi[philoc-1] - eraphi[philoc] # latitude spacing
        
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
    
    Ks = pl.sum(Ki,axis=0) # sum over all points for this trajectory release
    # make it for the ERA-Interim grid:
    Kx = pl.zeros([eralat.size,eralon.size])
    Kx[:,:28] = Ks[:,540:]; Kx[:,484:] = Ks[:,:28]; Kx = Kx + Ks[:,28:540]
    Kz[rls] = Kx # stick into Kz

Fx = pl.zeros([NRLS,256,512])
Fx[:,:,:28] = F[:,:,540:]; Fx[:,:,484:] = F[:,:,:28]; Fx = Fx + F[:,:,28:540]
#Ktot = pl.sum(Kz,axis=0) # sum over all trajectory releases
eralam = pl.radians(eralon)
lat_half = HalfGrid(eraphi)
nlon = eralon.shape[0] # number of longitude points
delta_lambda = (2*pl.pi)/nlon
D = pl.zeros([NRLS,eraphi.size,eralam.size])
T = pl.zeros([NRLS,eraphi.size,eralam.size])

for rls in range(NRLS):
    for i in range(lat_half.size-1): # loops over 256
        dmu = pl.sin(lat_half[i]) - pl.sin(lat_half[i+1])
        for j in range(eralon.size): # loops over 512
            D[rls,i,j] = Kz[rls,i,j]*dmu*delta_lambda
            T[rls,i,j] = Fx[rls,i,j]*dmu*delta_lambda