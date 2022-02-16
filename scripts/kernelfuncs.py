# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 09:12:16 2017

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
#from trajfuncs import Haversine,NearestIndex
# previous line commented out as functions are now in this file

def NearestIndex(array_in,point_in):
    """Function to the the nearest index to a specified geographic co-ordinate from an
    array of latitude or longitude co-ordinates
    
    Args:
        array_in (array): longitude or latitude array
        point_in (float): longitude or latitude co-ordinate
    
    Returns:
        index (int): index of array which has value closest to point_in
    """
    index = pl.absolute(array_in-point_in).argmin()
    
    return index

def Haversine(pt1,pt2):
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
    phi1 = pt1[1]; phi2 = pt2[1]; dphi = phi2 - phi1
    lam1 = pt1[0]; lam2 = pt2[0]; dlam = lam2 - lam1
    
    # calculate dsigma:
    dsig = 2*pl.arcsin(pl.sqrt((pl.sin(dphi/2)**2 + pl.cos(phi1)*pl.cos(phi2)*pl.sin(dlam/2)**2)))
    
    R = 6.37e6 # radius of Earth (m)
    
    d = R*dsig # calculate distance
    
    return d

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
    K = pl.zeros_like(theta)
    
    for i in range(newlam.size): # loop over longitude subset
        for j in range(newphi.size): # loop over latitude subset
            x_grid = (newlam[i],newphi[j],a) # point from subset of ERA-I grid
            # calculate theta using dot product in spherical co-ordinates:
            theta[i,j] = pl.arccos(SphereDot(point,x_grid))
            r = a*theta[i,j] # value of r
            if r > R: # check if r in kernel
                r = 0.; K[i,j] = 0. # set r & K as zero if outside kernel
            else: # otherwise calculate value of kernel
                N = 1/(2*pl.pi*(1-(2*a/R)*pl.sin(R/a) + (2*a**2/R**2)*(1-pl.cos(R/a))))
                K[i,j] = N*(1-(r**2)/(R**2)) # kernel function
    
    return K

def HalfGrid(lat):
    """Function to create latitude half-grid
    
    Args:
        lat (array): latitude array in radians
    
    Returns:
        lath (array): half-grid latitude array in radians
    """
    # set up empty array, size one more than lat
    lath = pl.zeros([lat.shape[0]+1])
    # set first & last elements of lat_half seperately as pi/2 & -pi/2:
    lath[0] = (pl.pi)/2; lath[-1] = -(pl.pi)/2
    # loop over lat_half from index 1 to index -2:
    for i in range(1,lath.shape[0]-1): # loops over 256
        lath[i] = 0.5*(lat[i]+lat[i-1])
    
    return lath

def MainKernelFunc(catchcount,eralon,eralat,crossec,rlslabs,CAT,shed,bas,year,month):
    """Calculate the density of trajectory origins
    
    Args:
        catchcount (list): list of length equal to number of trajectory releases,
                            contains trajectory no., origin category, origin
                            timestep & lat/lon of origin for a particular
                            catchment
        eralon (array): longitude from ERA-Interim
        eralat (array): latitude from ERA-Interim
        crossec (array): vertical cross-section of moisture flux, not weighted
                         by dl
        rlslabs (array): lon,lat co-ordinates of trajectory release points
        CAT (string): origin category of trajectories (all, CAT I, CAT II)
        shed (string): release line of trajectories
        bas (string): origin catchment of trajectories
        year (string): year
        month (string): month
    
    Returns:
        D (array): sum of all kernels integrated over sphere for each trajectory
                   release (not normalized!)
        T (array): moisture flux-weighted density estimate integrated over sphere
                   for each trajectory release
    """
    R = 200000 # kernel radius
    a = 6.371*10**6 # radius of earth
    
    NRLS = len(catchcount) # number of trajectory releases
    NREL = rlslabs.shape[0] # number of release points
    dl = pl.zeros([NREL]) # empty arrays for distance weights
    for pt in range(1,NREL-1): # exclude first & last release points
        # use Haversine formula to calculate dl
        dl[pt] = 0.5*(Haversine(rlslabs[pt-1,:2],rlslabs[pt,:2]) + 
                        Haversine(rlslabs[pt,:2],rlslabs[pt+1,:2]))
    # use first 2 points & last 2 points for first & last dl
    dl[0] = Haversine(rlslabs[0,:2],rlslabs[1,:2])
    dl[-1] = Haversine(rlslabs[-2,:2],rlslabs[-1,:2])
    
    sfx = crossec*dl # vertical cross-section of moisture flux, weighted by dl
    
    # To deal with points near zero longitude, extend the longitude array
    # extended longitude array ranges from -320 to 380
    lonx = pl.zeros([eralon.size+56])
    lonx[:28] = -1*pl.flipud(eralon[1:29]); lonx[28:540] = eralon;
    lonx[540:] = eralon[:28] + 360
    
    eraphi = pl.radians(eralat); eralam = pl.radians(lonx) # lat,lon in radians
    
    # empty Kernel array for all releases:
    Kz = pl.zeros([NRLS,eralat.size,eralon.size])
    F = pl.zeros([NRLS,eraphi.size,eralam.size])
    NTRAJ = pl.zeros([NRLS])
    for rls in range(NRLS): # loop over trajectory releases
        if catchcount[rls].size == 0.: # if no origins in release...
            continue # ... skip to next iteration in loop
        # extract origin point from catchcount array, flip so lon-lat:
        points = catchcount[rls][:,3:]; points = pl.fliplr(points)
    
        # empty array for Kernel for this release for all points
        Ki = pl.zeros([points.shape[0],eralat.size,lonx.size])
        cross = sfx[rls].flatten() # flatten cross-section so 1D
    
        for pt in range(points.shape[0]): # loop over points
            if CAT == 'all':
                pass
            elif catchcount[rls][pt,1] != int(CAT):
                continue
    #        if catchcount[rls][pt,2] == 0.:
    #            continue
        
            NTRAJ[rls] = NTRAJ[rls] + 1
            trajno = catchcount[rls][pt,0]#; tj.append(trajno)
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
    
    # write NTRAJ to a text file:
    clusdir = '/net/glusterfs/scenario/users/np838619/traj/' + year +'/' + shed
    f = open(clusdir+'/NTRAJ/'+shed+'_'+bas+year+month+'_'+CAT+'.csv','w')
    pl.savetxt(f,NTRAJ.T)
    f.close()
    
    return D, T

def MainKernelFunc2(newlabs,eralon,eralat,crossec,rlslabs,CAT,shed,bas,year,month):
    """Calculate the density of trajectory origins
    
    Args:
        newlabs (array): array of dimensions NRLS x NPART x 4, contains info on
                        origin type, time and lat/lon co-ordinates
        eralon (array): longitude from ERA-Interim
        eralat (array): latitude from ERA-Interim
        crossec (array): vertical cross-section of moisture flux, not weighted
                         by dl
        rlslabs (array): lon,lat co-ordinates of trajectory release points
        CAT (string): origin category of trajectories (all, CAT I, CAT II)
        shed (string): release line of trajectories
        bas (string): origin catchment of trajectories
        year (string): year
        month (string): month
    
    Returns:
        D (array): sum of all kernels integrated over sphere for each trajectory
                   release (not normalized!)
        T (array): moisture flux-weighted density estimate integrated over sphere
                   for each trajectory release
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
    
    sfx = crossec*dl # vertical cross-section of moisture flux, weighted by dl
    
    # To deal with points near zero longitude, extend the longitude array
    # extended longitude array ranges from -320 to 380
    lonx = pl.zeros([eralon.size+56])
    lonx[:28] = -1*pl.flipud(eralon[1:29]); lonx[28:540] = eralon;
    lonx[540:] = eralon[:28] + 360
    
    eraphi = pl.radians(eralat); eralam = pl.radians(lonx) # lat,lon in radians
    
    # empty Kernel array for all releases:
    Kz = pl.zeros([NRLS,eralat.size,eralon.size])
    F = pl.zeros([NRLS,eraphi.size,eralam.size])
    NTRAJ = pl.zeros([NRLS])
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
        
            NTRAJ[rls] = NTRAJ[rls] + 1
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
    
    # write NTRAJ to a text file:
    clusdir = '/net/glusterfs/scenario/users/np838619/traj/' + year +'/' + shed
    f = open(clusdir+'/NTRAJ/'+shed+'_'+bas+year+month+'_'+CAT+'.csv','w')
    pl.savetxt(f,NTRAJ.T)
    f.close()
    
    return D, T

def MakeNCfile(trajfiles,eralon,eralat,shed):
    """Function to open up a netcdf file and create the latitude, longitude and
    time dimensions and variables in preparation for adding the density fields
    later.
    
    Args:
        trajfiles (list): namse of trajectory release files with date & time stamp
        eralon (array): longitude array from ERA-Interim
        eralat (array): latitude array from ERA-Interim
        shed (string): three letter abbreviation indicating catchment boundary
                       from which trajectories were released
    """
    year = trajfiles[0][12:16]
    month = trajfiles[0][16:18]
    clusdir = '/net/glusterfs/scenario/users/np838619/traj/' + year + '/'\
                                                            + shed + '/density/'
    ncfile = Dataset(clusdir+ 'ordens_' + shed + year + month + '.nc','w')
    
    lat_dim = ncfile.createDimension('lat',eralat.size)
    lat_in = ncfile.createVariable('lat',pl.float64,('lat',))
    lat_in.units = 'degrees_north'
    lat_in.long_name = 'latitude'
    lat_in[:] = eralat
    
    lon_dim = ncfile.createDimension('lon',eralon.size)
    lon_in = ncfile.createVariable('lon',pl.float64,('lon',))
    lon_in.units = 'degrees_east'
    lon_in.long_name = 'longitude'
    lon_in[:] = eralon
    
    time_dim = ncfile.createDimension('time',len(trajfiles))
    time = ncfile.createVariable('time',pl.float64,('time',))
    time.units = 'hours'
    time.long_name = 'release date/time'
    
    ncfile.close()

def SaveFields(trajfiles,densfield,fluxfield,shed,catch,cat):
    """Function to append origin density & flux density fields to a netcdf file
    which has already been created earlier.
    
    Args:
        trajfiles (list): name of trajectory release files with date & time stamp
        densfield (array): trajectory origin density field (unnormalized),
                            dimensions time x lat x lon, units trajectories/sr
        fluxfield (array): moisture flux origin density field, same dimensions
                            as above
        shed (string): catchment boundary from which trajectories have been released
        catch (string): catchment area of origin
        cat (string): catgeory of trajectory origins
    """
    # which category of trajectories? all, CAT 1 or CAT 2?
    if cat == 'all':
        field_ex = ' all'
    elif cat == '1':
        field_ex = ' CAT I'
    elif cat == '2':
        field_ex = ' CAT II'
    
    year = trajfiles[0][12:16]
    month = trajfiles[0][16:18]
    clusdir = '/net/glusterfs/scenario/users/np838619/traj/' + year + '/'\
                              + shed + '/density/'
    ncfile = Dataset(clusdir+'ordens_' + shed + year + month + '.nc','a')
    
    # write origin density field:
    D = ncfile.createVariable('origin density ' + shed + catch + field_ex,pl.float64,
                                                          ('time','lat','lon'))
    D.units = 'trajectories per steradian'
    D.standard_name = 'trajectory origin density ' + shed + catch + field_ex
    D[:,:,:] = densfield
    
    F = ncfile.createVariable('flux density ' + shed + catch + field_ex,pl.float64,
                                                          ('time','lat','lon'))
    F.units = 'Sv per steradian'
    F.standard_name = 'moisture flux origin density ' + shed + catch + field_ex
    F[:,:,:] = fluxfield
    
    ncfile.close()

def Kernels_netcdf(newlabs,atlcount,indcount,paccount,arccount,soucount,
                                   trajfiles,sf3,eralon,eralat,rlslabs,shed):
    """Function to call the kernel functions and write all the monthly data to
    netcdf files.
    
    Args:
        newlabs (array): origin category, timestep & co-ordinates of every
                         trajectory for every release
        atlcount (list): trajectories with Atlantic origin, contains trajectory
                         number, origin type/timestep/co-ordinates
        indcount (list): as above for the Indian Ocean
        paccount (list): as above for the Pacific Ocean
        indcount (list): as above for the Arctic Ocean
        soucount (list): as above for the Southern Ocean
        trajfiles (list): names of trajectory data files
        sf3 (array): cross-section of moisture flux on catchment boundary
        eralon (array): ERA-Interim longitude array
        eralat (array): ERA-Interim latitude array
        rlslabs (array): longitude, latitude points along the catchment boundary
        shed (string): catchment boundary from which trajectories are released
    """
    # Set up empty lists for each month of traj data files:
    janfiles = []; febfiles = []; marfiles = []; aprfiles = []; mayfiles = []
    junfiles = []; julfiles = []; augfiles = []; sepfiles = []; octfiles = []
    novfiles = []; decfiles = []
    
    # for 2014 1st Jan 00z 2015 is included and should go into December list
    for i in range(len(trajfiles)): # loop over trajfiles
        if trajfiles[i][16:18] == '01': # if month is January...
            janfiles.append(trajfiles[i]) # ... add to January list
        elif trajfiles[i][16:18] == '02': # if month is February...
            febfiles.append(trajfiles[i]) # etc.
        elif trajfiles[i][16:18] == '03':
            marfiles.append(trajfiles[i])
        elif trajfiles[i][16:18] == '04':
            aprfiles.append(trajfiles[i])
        elif trajfiles[i][16:18] == '05':
            mayfiles.append(trajfiles[i])
        elif trajfiles[i][16:18] == '06':
            junfiles.append(trajfiles[i])
        elif trajfiles[i][16:18] == '07':
            julfiles.append(trajfiles[i])
        elif trajfiles[i][16:18] == '08':
            augfiles.append(trajfiles[i])
        elif trajfiles[i][16:18] == '09':
            sepfiles.append(trajfiles[i])
        elif trajfiles[i][16:18] == '10':
            octfiles.append(trajfiles[i])
        elif trajfiles[i][16:18] == '11':
            novfiles.append(trajfiles[i])
        elif trajfiles[i][16:18] == '12':
            decfiles.append(trajfiles[i])

    #if trajfiles[-2][12:16] == '2014' and trajfiles[-1][12:16] == '2015':
    #    decfiles.append(trajfiles[-1])

    # need to use the lengths of each list to split up catchcount lists
    janlen = len(janfiles); feblen = len(febfiles); marlen = len(marfiles)
    aprlen = len(aprfiles); maylen = len(mayfiles); junlen = len(junfiles)
    jullen = len(julfiles); auglen = len(augfiles); seplen = len(sepfiles)
    octlen = len(octfiles); novlen = len(novfiles); declen = len(decfiles)
    
    # split the newlabs & catchcount arrays/lists into months:
    janlabs = newlabs[:janlen,:,:]; sf_jan = sf3[:janlen,:,:] # January
    jan_atl = atlcount[:janlen]; jan_ind = indcount[:janlen]
    jan_pac = paccount[:janlen]; jan_arc = arccount[:janlen]
    jan_sou = soucount[:janlen]
    oldlen = janlen; newlen = oldlen + feblen
    feblabs = newlabs[oldlen:newlen,:,:]; sf_feb = sf3[oldlen:newlen,:,:] # February
    feb_atl = atlcount[oldlen:newlen]; feb_ind = indcount[oldlen:newlen]
    feb_pac = paccount[oldlen:newlen]; feb_arc = arccount[oldlen:newlen]
    feb_sou = soucount[oldlen:newlen]
    oldlen = newlen; newlen = oldlen + marlen
    marlabs = newlabs[oldlen:newlen,:,:]; sf_mar = sf3[oldlen:newlen,:,:] # March
    mar_atl = atlcount[oldlen:newlen]; mar_ind = indcount[oldlen:newlen]
    mar_pac = paccount[oldlen:newlen]; mar_arc = arccount[oldlen:newlen]
    mar_sou = soucount[oldlen:newlen]
    oldlen = newlen; newlen = oldlen + aprlen
    aprlabs = newlabs[oldlen:newlen,:,:]; sf_apr = sf3[oldlen:newlen,:,:] # April
    apr_atl = atlcount[oldlen:newlen]; apr_ind = indcount[oldlen:newlen]
    apr_pac = paccount[oldlen:newlen]; apr_arc = arccount[oldlen:newlen]
    apr_sou = soucount[oldlen:newlen]
    oldlen = newlen; newlen = oldlen + maylen
    maylabs = newlabs[oldlen:newlen,:,:]; sf_may = sf3[oldlen:newlen,:,:] # May
    may_atl = atlcount[oldlen:newlen]; may_ind = indcount[oldlen:newlen]
    may_pac = paccount[oldlen:newlen]; may_arc = arccount[oldlen:newlen]
    may_sou = soucount[oldlen:newlen]
    oldlen = newlen; newlen = oldlen + junlen
    junlabs = newlabs[oldlen:newlen,:,:]; sf_jun = sf3[oldlen:newlen,:,:] # June
    jun_atl = atlcount[oldlen:newlen]; jun_ind = indcount[oldlen:newlen]
    jun_pac = paccount[oldlen:newlen]; jun_arc = arccount[oldlen:newlen]
    jun_sou = soucount[oldlen:newlen]
    oldlen = newlen; newlen = oldlen + jullen
    jullabs = newlabs[oldlen:newlen,:,:]; sf_jul = sf3[oldlen:newlen,:,:] # July
    jul_atl = atlcount[oldlen:newlen]; jul_ind = indcount[oldlen:newlen]
    jul_pac = paccount[oldlen:newlen]; jul_arc = arccount[oldlen:newlen]
    jul_sou = soucount[oldlen:newlen]
    oldlen = newlen; newlen = oldlen + auglen
    auglabs = newlabs[oldlen:newlen,:,:]; sf_aug = sf3[oldlen:newlen,:,:] # August
    aug_atl = atlcount[oldlen:newlen]; aug_ind = indcount[oldlen:newlen]
    aug_pac = paccount[oldlen:newlen]; aug_arc = arccount[oldlen:newlen]
    aug_sou = soucount[oldlen:newlen]
    oldlen = newlen; newlen = oldlen + seplen
    seplabs = newlabs[oldlen:newlen,:,:]; sf_sep = sf3[oldlen:newlen,:,:] # September
    sep_atl = atlcount[oldlen:newlen]; sep_ind = indcount[oldlen:newlen]
    sep_pac = paccount[oldlen:newlen]; sep_arc = arccount[oldlen:newlen]
    sep_sou = soucount[oldlen:newlen]
    oldlen = newlen; newlen = oldlen + octlen
    octlabs = newlabs[oldlen:newlen,:,:]; sf_oct = sf3[oldlen:newlen,:,:] # October
    oct_atl = atlcount[oldlen:newlen]; oct_ind = indcount[oldlen:newlen]
    oct_pac = paccount[oldlen:newlen]; oct_arc = arccount[oldlen:newlen]
    oct_sou = soucount[oldlen:newlen]
    oldlen = newlen; newlen = oldlen + novlen
    novlabs = newlabs[oldlen:newlen,:,:]; sf_nov = sf3[oldlen:newlen,:,:] # November
    nov_atl = atlcount[oldlen:newlen]; nov_ind = indcount[oldlen:newlen]
    nov_pac = paccount[oldlen:newlen]; nov_arc = arccount[oldlen:newlen]
    nov_sou = soucount[oldlen:newlen]
    oldlen = newlen; newlen = oldlen + declen
    declabs = newlabs[oldlen:newlen,:,:]; sf_dec = sf3[oldlen:newlen,:,:] # December
    dec_atl = atlcount[oldlen:newlen]; dec_ind = indcount[oldlen:newlen]
    dec_pac = paccount[oldlen:newlen]; dec_arc = arccount[oldlen:newlen]
    dec_sou = soucount[oldlen:newlen]
    
    yr = trajfiles[0][12:16]
    # January:
    mnt = janfiles[0][16:18]
    # all catchments:
    all_DA01, all_FA01 = MainKernelFunc2(janlabs,eralon,eralat,sf_jan,rlslabs,'all',shed,'all',yr,mnt)
    all_D101, all_F101 = MainKernelFunc2(janlabs,eralon,eralat,sf_jan,rlslabs,'1',shed,'all',yr,mnt)
    all_D201, all_F201 = MainKernelFunc2(janlabs,eralon,eralat,sf_jan,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA01, atl_FA01 = MainKernelFunc(jan_atl,eralon,eralat,sf_jan,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D101, atl_F101 = MainKernelFunc(jan_atl,eralon,eralat,sf_jan,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D201, atl_F201 = MainKernelFunc(jan_atl,eralon,eralat,sf_jan,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA01, ind_FA01 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D101, ind_F101 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D201, ind_F201 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA01, pac_FA01 = MainKernelFunc(jan_pac,eralon,eralat,sf_jan,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D101, pac_F101 = MainKernelFunc(jan_pac,eralon,eralat,sf_jan,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D201, pac_F201 = MainKernelFunc(jan_pac,eralon,eralat,sf_jan,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA01, arc_FA01 = MainKernelFunc(jan_arc,eralon,eralat,sf_jan,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D101, arc_F101 = MainKernelFunc(jan_arc,eralon,eralat,sf_jan,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D201, arc_F201 = MainKernelFunc(jan_arc,eralon,eralat,sf_jan,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA01, sou_FA01 = MainKernelFunc(jan_sou,eralon,eralat,sf_jan,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D101, sou_F101 = MainKernelFunc(jan_sou,eralon,eralat,sf_jan,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D201, sou_F201 = MainKernelFunc(jan_sou,eralon,eralat,sf_jan,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********January done**********'

    # February:
    mnt = febfiles[0][16:18]
    # all catchments:
    all_DA02, all_FA02 = MainKernelFunc2(feblabs,eralon,eralat,sf_feb,rlslabs,'all',shed,'all',yr,mnt)
    all_D102, all_F102 = MainKernelFunc2(feblabs,eralon,eralat,sf_feb,rlslabs,'1',shed,'all',yr,mnt)
    all_D202, all_F202 = MainKernelFunc2(feblabs,eralon,eralat,sf_feb,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA02, atl_FA02 = MainKernelFunc(feb_atl,eralon,eralat,sf_feb,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D102, atl_F102 = MainKernelFunc(feb_atl,eralon,eralat,sf_feb,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D202, atl_F202 = MainKernelFunc(feb_atl,eralon,eralat,sf_feb,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA02, ind_FA02 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D102, ind_F102 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D202, ind_F202 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA02, pac_FA02 = MainKernelFunc(feb_pac,eralon,eralat,sf_feb,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D102, pac_F102 = MainKernelFunc(feb_pac,eralon,eralat,sf_feb,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D202, pac_F202 = MainKernelFunc(feb_pac,eralon,eralat,sf_feb,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA02, arc_FA02 = MainKernelFunc(feb_arc,eralon,eralat,sf_feb,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D102, arc_F102 = MainKernelFunc(feb_arc,eralon,eralat,sf_feb,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D202, arc_F202 = MainKernelFunc(feb_arc,eralon,eralat,sf_feb,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA02, sou_FA02 = MainKernelFunc(feb_sou,eralon,eralat,sf_feb,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D102, sou_F102 = MainKernelFunc(feb_sou,eralon,eralat,sf_feb,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D202, sou_F202 = MainKernelFunc(feb_sou,eralon,eralat,sf_feb,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********February done**********'

    # March:
    mnt = marfiles[0][16:18]
    # all catchments:
    all_DA03, all_FA03 = MainKernelFunc2(marlabs,eralon,eralat,sf_mar,rlslabs,'all',shed,'all',yr,mnt)
    all_D103, all_F103 = MainKernelFunc2(marlabs,eralon,eralat,sf_mar,rlslabs,'1',shed,'all',yr,mnt)
    all_D203, all_F203 = MainKernelFunc2(marlabs,eralon,eralat,sf_mar,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA03, atl_FA03 = MainKernelFunc(mar_atl,eralon,eralat,sf_mar,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D103, atl_F103 = MainKernelFunc(mar_atl,eralon,eralat,sf_mar,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D203, atl_F203 = MainKernelFunc(mar_atl,eralon,eralat,sf_mar,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA03, ind_FA03 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D103, ind_F103 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D203, ind_F203 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA03, pac_FA03 = MainKernelFunc(mar_pac,eralon,eralat,sf_mar,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D103, pac_F103 = MainKernelFunc(mar_pac,eralon,eralat,sf_mar,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D203, pac_F203 = MainKernelFunc(mar_pac,eralon,eralat,sf_mar,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA03, arc_FA03 = MainKernelFunc(mar_arc,eralon,eralat,sf_mar,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D103, arc_F103 = MainKernelFunc(mar_arc,eralon,eralat,sf_mar,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D203, arc_F203 = MainKernelFunc(mar_arc,eralon,eralat,sf_mar,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA03, sou_FA03 = MainKernelFunc(mar_sou,eralon,eralat,sf_mar,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D103, sou_F103 = MainKernelFunc(mar_sou,eralon,eralat,sf_mar,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D203, sou_F203 = MainKernelFunc(mar_sou,eralon,eralat,sf_mar,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********March done**********'

    # April:
    mnt = aprfiles[0][16:18]
    # all catchments:
    all_DA04, all_FA04 = MainKernelFunc2(aprlabs,eralon,eralat,sf_apr,rlslabs,'all',shed,'all',yr,mnt)
    all_D104, all_F104 = MainKernelFunc2(aprlabs,eralon,eralat,sf_apr,rlslabs,'1',shed,'all',yr,mnt)
    all_D204, all_F204 = MainKernelFunc2(aprlabs,eralon,eralat,sf_apr,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA04, atl_FA04 = MainKernelFunc(apr_atl,eralon,eralat,sf_apr,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D104, atl_F104 = MainKernelFunc(apr_atl,eralon,eralat,sf_apr,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D204, atl_F204 = MainKernelFunc(apr_atl,eralon,eralat,sf_apr,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA04, ind_FA04 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D104, ind_F104 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D204, ind_F204 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA04, pac_FA04 = MainKernelFunc(apr_pac,eralon,eralat,sf_apr,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D104, pac_F104 = MainKernelFunc(apr_pac,eralon,eralat,sf_apr,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D204, pac_F204 = MainKernelFunc(apr_pac,eralon,eralat,sf_apr,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA04, arc_FA04 = MainKernelFunc(apr_arc,eralon,eralat,sf_apr,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D104, arc_F104 = MainKernelFunc(apr_arc,eralon,eralat,sf_apr,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D204, arc_F204 = MainKernelFunc(apr_arc,eralon,eralat,sf_apr,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA04, sou_FA04 = MainKernelFunc(apr_sou,eralon,eralat,sf_apr,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D104, sou_F104 = MainKernelFunc(apr_sou,eralon,eralat,sf_apr,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D204, sou_F204 = MainKernelFunc(apr_sou,eralon,eralat,sf_apr,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********April done**********'

    # May:
    mnt = mayfiles[0][16:18]
    # all catchments:
    all_DA05, all_FA05 = MainKernelFunc2(maylabs,eralon,eralat,sf_may,rlslabs,'all',shed,'all',yr,mnt)
    all_D105, all_F105 = MainKernelFunc2(maylabs,eralon,eralat,sf_may,rlslabs,'1',shed,'all',yr,mnt)
    all_D205, all_F205 = MainKernelFunc2(maylabs,eralon,eralat,sf_may,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA05, atl_FA05 = MainKernelFunc(may_atl,eralon,eralat,sf_may,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D105, atl_F105 = MainKernelFunc(may_atl,eralon,eralat,sf_may,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D205, atl_F205 = MainKernelFunc(may_atl,eralon,eralat,sf_may,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA05, ind_FA05 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D105, ind_F105 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D205, ind_F205 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA05, pac_FA05 = MainKernelFunc(may_pac,eralon,eralat,sf_may,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D105, pac_F105 = MainKernelFunc(may_pac,eralon,eralat,sf_may,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D205, pac_F205 = MainKernelFunc(may_pac,eralon,eralat,sf_may,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic;
    arc_DA05, arc_FA05 = MainKernelFunc(may_arc,eralon,eralat,sf_may,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D105, arc_F105 = MainKernelFunc(may_arc,eralon,eralat,sf_may,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D205, arc_F205 = MainKernelFunc(may_arc,eralon,eralat,sf_may,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA05, sou_FA05 = MainKernelFunc(may_sou,eralon,eralat,sf_may,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D105, sou_F105 = MainKernelFunc(may_sou,eralon,eralat,sf_may,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D205, sou_F205 = MainKernelFunc(may_sou,eralon,eralat,sf_may,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********May done**********'

    # June:
    mnt = junfiles[0][16:18]
    # all catchments:
    all_DA06, all_FA06 = MainKernelFunc2(junlabs,eralon,eralat,sf_jun,rlslabs,'all',shed,'all',yr,mnt)
    all_D106, all_F106 = MainKernelFunc2(junlabs,eralon,eralat,sf_jun,rlslabs,'1',shed,'all',yr,mnt)
    all_D206, all_F206 = MainKernelFunc2(junlabs,eralon,eralat,sf_jun,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA06, atl_FA06 = MainKernelFunc(jun_atl,eralon,eralat,sf_jun,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D106, atl_F106 = MainKernelFunc(jun_atl,eralon,eralat,sf_jun,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D206, atl_F206 = MainKernelFunc(jun_atl,eralon,eralat,sf_jun,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA06, ind_FA06 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D106, ind_F106 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D206, ind_F206 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA06, pac_FA06 = MainKernelFunc(jun_pac,eralon,eralat,sf_jun,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D106, pac_F106 = MainKernelFunc(jun_pac,eralon,eralat,sf_jun,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D206, pac_F206 = MainKernelFunc(jun_pac,eralon,eralat,sf_jun,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA06, arc_FA06 = MainKernelFunc(jun_arc,eralon,eralat,sf_jun,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D106, arc_F106 = MainKernelFunc(jun_arc,eralon,eralat,sf_jun,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D206, arc_F206 = MainKernelFunc(jun_arc,eralon,eralat,sf_jun,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA06, sou_FA06 = MainKernelFunc(jun_sou,eralon,eralat,sf_jun,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D106, sou_F106 = MainKernelFunc(jun_sou,eralon,eralat,sf_jun,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D206, sou_F206 = MainKernelFunc(jun_sou,eralon,eralat,sf_jun,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********June done**********'

    # July:
    mnt = julfiles[0][16:18]
    # all catchments:
    all_DA07, all_FA07 = MainKernelFunc2(jullabs,eralon,eralat,sf_jul,rlslabs,'all',shed,'all',yr,mnt)
    all_D107, all_F107 = MainKernelFunc2(jullabs,eralon,eralat,sf_jul,rlslabs,'1',shed,'all',yr,mnt)
    all_D207, all_F207 = MainKernelFunc2(jullabs,eralon,eralat,sf_jul,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA07, atl_FA07 = MainKernelFunc(jul_atl,eralon,eralat,sf_jul,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D107, atl_F107 = MainKernelFunc(jul_atl,eralon,eralat,sf_jul,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D207, atl_F207 = MainKernelFunc(jul_atl,eralon,eralat,sf_jul,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA07, ind_FA07 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D107, ind_F107 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D207, ind_F207 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA07, pac_FA07 = MainKernelFunc(jul_pac,eralon,eralat,sf_jul,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D107, pac_F107 = MainKernelFunc(jul_pac,eralon,eralat,sf_jul,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D207, pac_F207 = MainKernelFunc(jul_pac,eralon,eralat,sf_jul,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA07, arc_FA07 = MainKernelFunc(jul_arc,eralon,eralat,sf_jul,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D107, arc_F107 = MainKernelFunc(jul_arc,eralon,eralat,sf_jul,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D207, arc_F207 = MainKernelFunc(jul_arc,eralon,eralat,sf_jul,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA07, sou_FA07 = MainKernelFunc(jul_sou,eralon,eralat,sf_jul,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D107, sou_F107 = MainKernelFunc(jul_sou,eralon,eralat,sf_jul,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D207, sou_F207 = MainKernelFunc(jul_sou,eralon,eralat,sf_jul,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********July done**********'

    # August:
    mnt = augfiles[0][16:18]
    # all catchments:
    all_DA08, all_FA08 = MainKernelFunc2(auglabs,eralon,eralat,sf_aug,rlslabs,'all',shed,'all',yr,mnt)
    all_D108, all_F108 = MainKernelFunc2(auglabs,eralon,eralat,sf_aug,rlslabs,'1',shed,'all',yr,mnt)
    all_D208, all_F208 = MainKernelFunc2(auglabs,eralon,eralat,sf_aug,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA08, atl_FA08 = MainKernelFunc(aug_atl,eralon,eralat,sf_aug,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D108, atl_F108 = MainKernelFunc(aug_atl,eralon,eralat,sf_aug,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D208, atl_F208 = MainKernelFunc(aug_atl,eralon,eralat,sf_aug,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA08, ind_FA08 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D108, ind_F108 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D208, ind_F208 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA08, pac_FA08 = MainKernelFunc(aug_pac,eralon,eralat,sf_aug,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D108, pac_F108 = MainKernelFunc(aug_pac,eralon,eralat,sf_aug,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D208, pac_F208 = MainKernelFunc(aug_pac,eralon,eralat,sf_aug,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA08, arc_FA08 = MainKernelFunc(aug_arc,eralon,eralat,sf_aug,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D108, arc_F108 = MainKernelFunc(aug_arc,eralon,eralat,sf_aug,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D208, arc_F208 = MainKernelFunc(aug_arc,eralon,eralat,sf_aug,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA08, sou_FA08 = MainKernelFunc(aug_sou,eralon,eralat,sf_aug,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D108, sou_F108 = MainKernelFunc(aug_sou,eralon,eralat,sf_aug,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D208, sou_F208 = MainKernelFunc(aug_sou,eralon,eralat,sf_aug,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********August done**********'

    # September:
    mnt = sepfiles[0][16:18]
    # all catchments:
    all_DA09, all_FA09 = MainKernelFunc2(seplabs,eralon,eralat,sf_sep,rlslabs,'all',shed,'all',yr,mnt)
    all_D109, all_F109 = MainKernelFunc2(seplabs,eralon,eralat,sf_sep,rlslabs,'1',shed,'all',yr,mnt)
    all_D209, all_F209 = MainKernelFunc2(seplabs,eralon,eralat,sf_sep,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA09, atl_FA09 = MainKernelFunc(sep_atl,eralon,eralat,sf_sep,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D109, atl_F109 = MainKernelFunc(sep_atl,eralon,eralat,sf_sep,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D209, atl_F209 = MainKernelFunc(sep_atl,eralon,eralat,sf_sep,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA09, ind_FA09 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D109, ind_F109 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D209, ind_F209 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA09, pac_FA09 = MainKernelFunc(sep_pac,eralon,eralat,sf_sep,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D109, pac_F109 = MainKernelFunc(sep_pac,eralon,eralat,sf_sep,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D209, pac_F209 = MainKernelFunc(sep_pac,eralon,eralat,sf_sep,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA09, arc_FA09 = MainKernelFunc(sep_arc,eralon,eralat,sf_sep,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D109, arc_F109 = MainKernelFunc(sep_arc,eralon,eralat,sf_sep,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D209, arc_F209 = MainKernelFunc(sep_arc,eralon,eralat,sf_sep,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA09, sou_FA09 = MainKernelFunc(sep_sou,eralon,eralat,sf_sep,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D109, sou_F109 = MainKernelFunc(sep_sou,eralon,eralat,sf_sep,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D209, sou_F209 = MainKernelFunc(sep_sou,eralon,eralat,sf_sep,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********September done**********'

    # October:
    mnt = octfiles[0][16:18]
    # all catchments:
    all_DA10, all_FA10 = MainKernelFunc2(octlabs,eralon,eralat,sf_oct,rlslabs,'all',shed,'all',yr,mnt)
    all_D110, all_F110 = MainKernelFunc2(octlabs,eralon,eralat,sf_oct,rlslabs,'1',shed,'all',yr,mnt)
    all_D210, all_F210 = MainKernelFunc2(octlabs,eralon,eralat,sf_oct,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA10, atl_FA10 = MainKernelFunc(oct_atl,eralon,eralat,sf_oct,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D110, atl_F110 = MainKernelFunc(oct_atl,eralon,eralat,sf_oct,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D210, atl_F210 = MainKernelFunc(oct_atl,eralon,eralat,sf_oct,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA10, ind_FA10 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D110, ind_F110 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D210, ind_F210 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA10, pac_FA10 = MainKernelFunc(oct_pac,eralon,eralat,sf_oct,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D110, pac_F110 = MainKernelFunc(oct_pac,eralon,eralat,sf_oct,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D210, pac_F210 = MainKernelFunc(oct_pac,eralon,eralat,sf_oct,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA10, arc_FA10 = MainKernelFunc(oct_arc,eralon,eralat,sf_oct,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D110, arc_F110 = MainKernelFunc(oct_arc,eralon,eralat,sf_oct,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D210, arc_F210 = MainKernelFunc(oct_arc,eralon,eralat,sf_oct,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA10, sou_FA10 = MainKernelFunc(oct_sou,eralon,eralat,sf_oct,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D110, sou_F110 = MainKernelFunc(oct_sou,eralon,eralat,sf_oct,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D210, sou_F210 = MainKernelFunc(oct_sou,eralon,eralat,sf_oct,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********October done**********'

    # November:
    mnt = novfiles[0][16:18]
    # all catchments:
    all_DA11, all_FA11 = MainKernelFunc2(novlabs,eralon,eralat,sf_nov,rlslabs,'all',shed,'all',yr,mnt)
    all_D111, all_F111 = MainKernelFunc2(novlabs,eralon,eralat,sf_nov,rlslabs,'1',shed,'all',yr,mnt)
    all_D211, all_F211 = MainKernelFunc2(novlabs,eralon,eralat,sf_nov,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA11, atl_FA11 = MainKernelFunc(nov_atl,eralon,eralat,sf_nov,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D111, atl_F111 = MainKernelFunc(nov_atl,eralon,eralat,sf_nov,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D211, atl_F211 = MainKernelFunc(nov_atl,eralon,eralat,sf_nov,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA11, ind_FA11 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D111, ind_F111 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D211, ind_F211 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA11, pac_FA11 = MainKernelFunc(nov_pac,eralon,eralat,sf_nov,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D111, pac_F111 = MainKernelFunc(nov_pac,eralon,eralat,sf_nov,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D211, pac_F211 = MainKernelFunc(nov_pac,eralon,eralat,sf_nov,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA11, arc_FA11 = MainKernelFunc(nov_arc,eralon,eralat,sf_nov,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D111, arc_F111 = MainKernelFunc(nov_arc,eralon,eralat,sf_nov,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D211, arc_F211 = MainKernelFunc(nov_arc,eralon,eralat,sf_nov,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA11, sou_FA11 = MainKernelFunc(nov_sou,eralon,eralat,sf_nov,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D111, sou_F111 = MainKernelFunc(nov_sou,eralon,eralat,sf_nov,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D211, sou_F211 = MainKernelFunc(nov_sou,eralon,eralat,sf_nov,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********November done**********'

    # December:
    mnt = decfiles[0][16:18]
    # all catchments:
    all_DA12, all_FA12 = MainKernelFunc2(declabs,eralon,eralat,sf_dec,rlslabs,'all',shed,'all',yr,mnt)
    all_D112, all_F112 = MainKernelFunc2(declabs,eralon,eralat,sf_dec,rlslabs,'1',shed,'all',yr,mnt)
    all_D212, all_F212 = MainKernelFunc2(declabs,eralon,eralat,sf_dec,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA12, atl_FA12 = MainKernelFunc(dec_atl,eralon,eralat,sf_dec,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D112, atl_F112 = MainKernelFunc(dec_atl,eralon,eralat,sf_dec,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D212, atl_F212 = MainKernelFunc(dec_atl,eralon,eralat,sf_dec,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA12, ind_FA12 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D112, ind_F112 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D212, ind_F212 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA12, pac_FA12 = MainKernelFunc(dec_pac,eralon,eralat,sf_dec,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D112, pac_F112 = MainKernelFunc(dec_pac,eralon,eralat,sf_dec,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D212, pac_F212 = MainKernelFunc(dec_pac,eralon,eralat,sf_dec,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA12, arc_FA12 = MainKernelFunc(dec_arc,eralon,eralat,sf_dec,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D112, arc_F112 = MainKernelFunc(dec_arc,eralon,eralat,sf_dec,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D212, arc_F212 = MainKernelFunc(dec_arc,eralon,eralat,sf_dec,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA12, sou_FA12 = MainKernelFunc(dec_sou,eralon,eralat,sf_dec,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D112, sou_F112 = MainKernelFunc(dec_sou,eralon,eralat,sf_dec,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D212, sou_F212 = MainKernelFunc(dec_sou,eralon,eralat,sf_dec,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********December done**********'
    
    # Save these fields to netcdf files:
    # January:
    # all catchments:
    SaveFields(janfiles,all_DA01,all_FA01,shed,'all','all')
    SaveFields(janfiles,all_D101,all_F101,shed,'all','1')
    SaveFields(janfiles,all_D201,all_F201,shed,'all','2')
    # Atlantic
    SaveFields(janfiles,atl_DA01,atl_FA01,shed,'atl','all')
    SaveFields(janfiles,atl_D101,atl_F101,shed,'atl','1')
    SaveFields(janfiles,atl_D201,atl_F201,shed,'atl','2')
    # Indian
    SaveFields(janfiles,ind_DA01,ind_FA01,shed,'ind','all')
    SaveFields(janfiles,ind_D101,ind_F101,shed,'ind','1')
    SaveFields(janfiles,ind_D201,ind_F201,shed,'ind','2')
    # Pacific
    SaveFields(janfiles,pac_DA01,pac_FA01,shed,'pac','all')
    SaveFields(janfiles,pac_D101,pac_F101,shed,'pac','1')
    SaveFields(janfiles,pac_D201,pac_F201,shed,'pac','2')
    # Arctic
    SaveFields(janfiles,arc_DA01,arc_FA01,shed,'arc','all')
    SaveFields(janfiles,arc_D101,arc_F101,shed,'arc','1')
    SaveFields(janfiles,arc_D201,arc_F201,shed,'arc','2')
    # Southern
    SaveFields(janfiles,sou_DA01,sou_FA01,shed,'sou','all')
    SaveFields(janfiles,sou_D101,sou_F101,shed,'sou','1')
    SaveFields(janfiles,sou_D201,sou_F201,shed,'sou','2')
    
    # February:
    # all catchments:
    SaveFields(febfiles,all_DA02,all_FA02,shed,'all','all')
    SaveFields(febfiles,all_D102,all_F102,shed,'all','1')
    SaveFields(febfiles,all_D202,all_F202,shed,'all','2')
    # Atlantic
    SaveFields(febfiles,atl_DA02,atl_FA02,shed,'atl','all')
    SaveFields(febfiles,atl_D102,atl_F102,shed,'atl','1')
    SaveFields(febfiles,atl_D202,atl_F202,shed,'atl','2')
    # Indian
    SaveFields(febfiles,ind_DA02,ind_FA02,shed,'ind','all')
    SaveFields(febfiles,ind_D102,ind_F102,shed,'ind','1')
    SaveFields(febfiles,ind_D202,ind_F202,shed,'ind','2')
    # Pacific
    SaveFields(febfiles,pac_DA02,pac_FA02,shed,'pac','all')
    SaveFields(febfiles,pac_D102,pac_F102,shed,'pac','1')
    SaveFields(febfiles,pac_D202,pac_F202,shed,'pac','2')
    # Arctic
    SaveFields(febfiles,arc_DA02,arc_FA02,shed,'arc','all')
    SaveFields(febfiles,arc_D102,arc_F102,shed,'arc','1')
    SaveFields(febfiles,arc_D202,arc_F202,shed,'arc','2')
    # Southern
    SaveFields(febfiles,sou_DA02,sou_FA02,shed,'sou','all')
    SaveFields(febfiles,sou_D102,sou_F102,shed,'sou','1')
    SaveFields(febfiles,sou_D202,sou_F202,shed,'sou','2')
    
    # March:
    # all catchments:
    SaveFields(marfiles,all_DA03,all_FA03,shed,'all','all')
    SaveFields(marfiles,all_D103,all_F103,shed,'all','1')
    SaveFields(marfiles,all_D203,all_F203,shed,'all','2')
    # Atlantic
    SaveFields(marfiles,atl_DA03,atl_FA03,shed,'atl','all')
    SaveFields(marfiles,atl_D103,atl_F103,shed,'atl','1')
    SaveFields(marfiles,atl_D203,atl_F203,shed,'atl','2')
    # Indian
    SaveFields(marfiles,ind_DA03,ind_FA03,shed,'ind','all')
    SaveFields(marfiles,ind_D103,ind_F103,shed,'ind','1')
    SaveFields(marfiles,ind_D203,ind_F203,shed,'ind','2')
    # Pacific
    SaveFields(marfiles,pac_DA03,pac_FA03,shed,'pac','all')
    SaveFields(marfiles,pac_D103,pac_F103,shed,'pac','1')
    SaveFields(marfiles,pac_D203,pac_F203,shed,'pac','2')
    # Arctic
    SaveFields(marfiles,arc_DA03,arc_FA03,shed,'arc','all')
    SaveFields(marfiles,arc_D103,arc_F103,shed,'arc','1')
    SaveFields(marfiles,arc_D203,arc_F203,shed,'arc','2')
    # Southern
    SaveFields(marfiles,sou_DA03,sou_FA03,shed,'sou','all')
    SaveFields(marfiles,sou_D103,sou_F103,shed,'sou','1')
    SaveFields(marfiles,sou_D203,sou_F203,shed,'sou','2')
    
    # April:
    # all catchments:
    SaveFields(aprfiles,all_DA04,all_FA04,shed,'all','all')
    SaveFields(aprfiles,all_D104,all_F104,shed,'all','1')
    SaveFields(aprfiles,all_D204,all_F204,shed,'all','2')
    # Atlantic
    SaveFields(aprfiles,atl_DA04,atl_FA04,shed,'atl','all')
    SaveFields(aprfiles,atl_D104,atl_F104,shed,'atl','1')
    SaveFields(aprfiles,atl_D204,atl_F204,shed,'atl','2')
    # Indian
    SaveFields(aprfiles,ind_DA04,ind_FA04,shed,'ind','all')
    SaveFields(aprfiles,ind_D104,ind_F104,shed,'ind','1')
    SaveFields(aprfiles,ind_D204,ind_F204,shed,'ind','2')
    # Pacific
    SaveFields(aprfiles,pac_DA04,pac_FA04,shed,'pac','all')
    SaveFields(aprfiles,pac_D104,pac_F104,shed,'pac','1')
    SaveFields(aprfiles,pac_D204,pac_F204,shed,'pac','2')
    # Arctic
    SaveFields(aprfiles,arc_DA04,arc_FA04,shed,'arc','all')
    SaveFields(aprfiles,arc_D104,arc_F104,shed,'arc','1')
    SaveFields(aprfiles,arc_D204,arc_F204,shed,'arc','2')
    # Southern
    SaveFields(aprfiles,sou_DA04,sou_FA04,shed,'sou','all')
    SaveFields(aprfiles,sou_D104,sou_F104,shed,'sou','1')
    SaveFields(aprfiles,sou_D204,sou_F204,shed,'sou','2')
    
    # May:
    # all catchments:
    SaveFields(mayfiles,all_DA05,all_FA05,shed,'all','all')
    SaveFields(mayfiles,all_D105,all_F105,shed,'all','1')
    SaveFields(mayfiles,all_D205,all_F205,shed,'all','2')
    # Atlantic
    SaveFields(mayfiles,atl_DA05,atl_FA05,shed,'atl','all')
    SaveFields(mayfiles,atl_D105,atl_F105,shed,'atl','1')
    SaveFields(mayfiles,atl_D205,atl_F205,shed,'atl','2')
    # Indian
    SaveFields(mayfiles,ind_DA05,ind_FA05,shed,'ind','all')
    SaveFields(mayfiles,ind_D105,ind_F105,shed,'ind','1')
    SaveFields(mayfiles,ind_D205,ind_F205,shed,'ind','2')
    # Pacific
    SaveFields(mayfiles,pac_DA05,pac_FA05,shed,'pac','all')
    SaveFields(mayfiles,pac_D105,pac_F105,shed,'pac','1')
    SaveFields(mayfiles,pac_D205,pac_F205,shed,'pac','2')
    # Arctic
    SaveFields(mayfiles,arc_DA05,arc_FA05,shed,'arc','all')
    SaveFields(mayfiles,arc_D105,arc_F105,shed,'arc','1')
    SaveFields(mayfiles,arc_D205,arc_F205,shed,'arc','2')
    # Southern
    SaveFields(mayfiles,sou_DA05,sou_FA05,shed,'sou','all')
    SaveFields(mayfiles,sou_D105,sou_F105,shed,'sou','1')
    SaveFields(mayfiles,sou_D205,sou_F205,shed,'sou','2')
    
    # June:
    # all catchments:
    SaveFields(junfiles,all_DA06,all_FA06,shed,'all','all')
    SaveFields(junfiles,all_D106,all_F106,shed,'all','1')
    SaveFields(junfiles,all_D206,all_F206,shed,'all','2')
    # Atlantic
    SaveFields(junfiles,atl_DA06,atl_FA06,shed,'atl','all')
    SaveFields(junfiles,atl_D106,atl_F106,shed,'atl','1')
    SaveFields(junfiles,atl_D206,atl_F206,shed,'atl','2')
    # Indian
    SaveFields(junfiles,ind_DA06,ind_FA06,shed,'ind','all')
    SaveFields(junfiles,ind_D106,ind_F106,shed,'ind','1')
    SaveFields(junfiles,ind_D206,ind_F206,shed,'ind','2')
    # Pacific
    SaveFields(junfiles,pac_DA06,pac_FA06,shed,'pac','all')
    SaveFields(junfiles,pac_D106,pac_F106,shed,'pac','1')
    SaveFields(junfiles,pac_D206,pac_F206,shed,'pac','2')
    # Arctic
    SaveFields(junfiles,arc_DA06,arc_FA06,shed,'arc','all')
    SaveFields(junfiles,arc_D106,arc_F106,shed,'arc','1')
    SaveFields(junfiles,arc_D206,arc_F206,shed,'arc','2')
    # Southern
    SaveFields(junfiles,sou_DA06,sou_FA06,shed,'sou','all')
    SaveFields(junfiles,sou_D106,sou_F106,shed,'sou','1')
    SaveFields(junfiles,sou_D206,sou_F206,shed,'sou','2')
    
    # July:
    # all catchments:
    SaveFields(julfiles,all_DA07,all_FA07,shed,'all','all')
    SaveFields(julfiles,all_D107,all_F107,shed,'all','1')
    SaveFields(julfiles,all_D207,all_F207,shed,'all','2')
    # Atlantic
    SaveFields(julfiles,atl_DA07,atl_FA07,shed,'atl','all')
    SaveFields(julfiles,atl_D107,atl_F107,shed,'atl','1')
    SaveFields(julfiles,atl_D207,atl_F207,shed,'atl','2')
    # Indian
    SaveFields(julfiles,ind_DA07,ind_FA07,shed,'ind','all')
    SaveFields(julfiles,ind_D107,ind_F107,shed,'ind','1')
    SaveFields(julfiles,ind_D207,ind_F207,shed,'ind','2')
    # Pacific
    SaveFields(julfiles,pac_DA07,pac_FA07,shed,'pac','all')
    SaveFields(julfiles,pac_D107,pac_F107,shed,'pac','1')
    SaveFields(julfiles,pac_D207,pac_F207,shed,'pac','2')
    # Arctic
    SaveFields(julfiles,arc_DA07,arc_FA07,shed,'arc','all')
    SaveFields(julfiles,arc_D107,arc_F107,shed,'arc','1')
    SaveFields(julfiles,arc_D207,arc_F207,shed,'arc','2')
    # Southern
    SaveFields(julfiles,sou_DA07,sou_FA07,shed,'sou','all')
    SaveFields(julfiles,sou_D107,sou_F107,shed,'sou','1')
    SaveFields(julfiles,sou_D207,sou_F207,shed,'sou','2')
    
    # August:
    # all catchments:
    SaveFields(augfiles,all_DA08,all_FA08,shed,'all','all')
    SaveFields(augfiles,all_D108,all_F108,shed,'all','1')
    SaveFields(augfiles,all_D208,all_F208,shed,'all','2')
    # Atlantic
    SaveFields(augfiles,atl_DA08,atl_FA08,shed,'atl','all')
    SaveFields(augfiles,atl_D108,atl_F108,shed,'atl','1')
    SaveFields(augfiles,atl_D208,atl_F208,shed,'atl','2')
    # Indian
    SaveFields(augfiles,ind_DA08,ind_FA08,shed,'ind','all')
    SaveFields(augfiles,ind_D108,ind_F108,shed,'ind','1')
    SaveFields(augfiles,ind_D208,ind_F208,shed,'ind','2')
    # Pacific
    SaveFields(augfiles,pac_DA08,pac_FA08,shed,'pac','all')
    SaveFields(augfiles,pac_D108,pac_F108,shed,'pac','1')
    SaveFields(augfiles,pac_D208,pac_F208,shed,'pac','2')
    # Arctic
    SaveFields(augfiles,arc_DA08,arc_FA08,shed,'arc','all')
    SaveFields(augfiles,arc_D108,arc_F108,shed,'arc','1')
    SaveFields(augfiles,arc_D208,arc_F208,shed,'arc','2')
    # Southern
    SaveFields(augfiles,sou_DA08,sou_FA08,shed,'sou','all')
    SaveFields(augfiles,sou_D108,sou_F108,shed,'sou','1')
    SaveFields(augfiles,sou_D208,sou_F208,shed,'sou','2')
    
    # September:
    # all catchments:
    SaveFields(sepfiles,all_DA09,all_FA09,shed,'all','all')
    SaveFields(sepfiles,all_D109,all_F109,shed,'all','1')
    SaveFields(sepfiles,all_D209,all_F209,shed,'all','2')
    # Atlantic
    SaveFields(sepfiles,atl_DA09,atl_FA09,shed,'atl','all')
    SaveFields(sepfiles,atl_D109,atl_F109,shed,'atl','1')
    SaveFields(sepfiles,atl_D209,atl_F209,shed,'atl','2')
    # Indian
    SaveFields(sepfiles,ind_DA09,ind_FA09,shed,'ind','all')
    SaveFields(sepfiles,ind_D109,ind_F109,shed,'ind','1')
    SaveFields(sepfiles,ind_D209,ind_F209,shed,'ind','2')
    # Pacific
    SaveFields(sepfiles,pac_DA09,pac_FA09,shed,'pac','all')
    SaveFields(sepfiles,pac_D109,pac_F109,shed,'pac','1')
    SaveFields(sepfiles,pac_D209,pac_F209,shed,'pac','2')
    # Arctic
    SaveFields(sepfiles,arc_DA09,arc_FA09,shed,'arc','all')
    SaveFields(sepfiles,arc_D109,arc_F109,shed,'arc','1')
    SaveFields(sepfiles,arc_D209,arc_F209,shed,'arc','2')
    # Southern
    SaveFields(sepfiles,sou_DA09,sou_FA09,shed,'sou','all')
    SaveFields(sepfiles,sou_D109,sou_F109,shed,'sou','1')
    SaveFields(sepfiles,sou_D209,sou_F209,shed,'sou','2')
    
    # October:
    # all catchments:
    SaveFields(octfiles,all_DA10,all_FA10,shed,'all','all')
    SaveFields(octfiles,all_D110,all_F110,shed,'all','1')
    SaveFields(octfiles,all_D210,all_F210,shed,'all','2')
    # Atlantic
    SaveFields(octfiles,atl_DA10,atl_FA10,shed,'atl','all')
    SaveFields(octfiles,atl_D110,atl_F110,shed,'atl','1')
    SaveFields(octfiles,atl_D210,atl_F210,shed,'atl','2')
    # Indian
    SaveFields(octfiles,ind_DA10,ind_FA10,shed,'ind','all')
    SaveFields(octfiles,ind_D110,ind_F110,shed,'ind','1')
    SaveFields(octfiles,ind_D210,ind_F210,shed,'ind','2')
    # Pacific
    SaveFields(octfiles,pac_DA10,pac_FA10,shed,'pac','all')
    SaveFields(octfiles,pac_D110,pac_F110,shed,'pac','1')
    SaveFields(octfiles,pac_D210,pac_F210,shed,'pac','2')
    # Arctic
    SaveFields(octfiles,arc_DA10,arc_FA10,shed,'arc','all')
    SaveFields(octfiles,arc_D110,arc_F110,shed,'arc','1')
    SaveFields(octfiles,arc_D210,arc_F210,shed,'arc','2')
    # Southern
    SaveFields(octfiles,sou_DA10,sou_FA10,shed,'sou','all')
    SaveFields(octfiles,sou_D110,sou_F110,shed,'sou','1')
    SaveFields(octfiles,sou_D210,sou_F210,shed,'sou','2')
    
    # November:
    # all catchments:
    SaveFields(novfiles,all_DA11,all_FA11,shed,'all','all')
    SaveFields(novfiles,all_D111,all_F111,shed,'all','1')
    SaveFields(novfiles,all_D211,all_F211,shed,'all','2')
    # Atlantic
    SaveFields(novfiles,atl_DA11,atl_FA11,shed,'atl','all')
    SaveFields(novfiles,atl_D111,atl_F111,shed,'atl','1')
    SaveFields(novfiles,atl_D211,atl_F211,shed,'atl','2')
    # Indian
    SaveFields(novfiles,ind_DA11,ind_FA11,shed,'ind','all')
    SaveFields(novfiles,ind_D111,ind_F111,shed,'ind','1')
    SaveFields(novfiles,ind_D211,ind_F211,shed,'ind','2')
    # Pacific
    SaveFields(novfiles,pac_DA11,pac_FA11,shed,'pac','all')
    SaveFields(novfiles,pac_D111,pac_F111,shed,'pac','1')
    SaveFields(novfiles,pac_D211,pac_F211,shed,'pac','2')
    # Arctic
    SaveFields(novfiles,arc_DA11,arc_FA11,shed,'arc','all')
    SaveFields(novfiles,arc_D111,arc_F111,shed,'arc','1')
    SaveFields(novfiles,arc_D211,arc_F211,shed,'arc','2')
    # Southern
    SaveFields(novfiles,sou_DA11,sou_FA11,shed,'sou','all')
    SaveFields(novfiles,sou_D111,sou_F111,shed,'sou','1')
    SaveFields(novfiles,sou_D211,sou_F211,shed,'sou','2')
    
    # December:
    # all catchments:
    SaveFields(decfiles,all_DA12,all_FA12,shed,'all','all')
    SaveFields(decfiles,all_D112,all_F112,shed,'all','1')
    SaveFields(decfiles,all_D212,all_F212,shed,'all','2')
    # Atlantic
    SaveFields(decfiles,atl_DA12,atl_FA12,shed,'atl','all')
    SaveFields(decfiles,atl_D112,atl_F112,shed,'atl','1')
    SaveFields(decfiles,atl_D212,atl_F212,shed,'atl','2')
    # Indian
    SaveFields(decfiles,ind_DA12,ind_FA12,shed,'ind','all')
    SaveFields(decfiles,ind_D112,ind_F112,shed,'ind','1')
    SaveFields(decfiles,ind_D212,ind_F212,shed,'ind','2')
    # Pacific
    SaveFields(decfiles,pac_DA12,pac_FA12,shed,'pac','all')
    SaveFields(decfiles,pac_D112,pac_F112,shed,'pac','1')
    SaveFields(decfiles,pac_D212,pac_F212,shed,'pac','2')
    # Arctic
    SaveFields(decfiles,arc_DA12,arc_FA12,shed,'arc','all')
    SaveFields(decfiles,arc_D112,arc_F112,shed,'arc','1')
    SaveFields(decfiles,arc_D212,arc_F212,shed,'arc','2')
    # Southern
    SaveFields(decfiles,sou_DA12,sou_FA12,shed,'sou','all')
    SaveFields(decfiles,sou_D112,sou_F112,shed,'sou','1')
    SaveFields(decfiles,sou_D212,sou_F212,shed,'sou','2')
    
    print '*******netcdf files for ' + shed + ' written************'


def Kernels_netcdf_JFM(newlabs,atlcount,indcount,paccount,arccount,soucount,
                                   trajfiles,sf3,eralon,eralat,rlslabs,shed):
    """Function to call the kernel functions and write all the monthly data to
    netcdf files for January, February and March.
    
    Args:
        newlabs (array): origin category, timestep & co-ordinates of every
                         trajectory for every release
        atlcount (list): trajectories with Atlantic origin, contains trajectory
                         number, origin type/timestep/co-ordinates
        indcount (list): as above for the Indian Ocean
        paccount (list): as above for the Pacific Ocean
        indcount (list): as above for the Arctic Ocean
        soucount (list): as above for the Southern Ocean
        trajfiles (list): names of trajectory data files
        sf3 (array): cross-section of moisture flux on catchment boundary
        eralon (array): ERA-Interim longitude array
        eralat (array): ERA-Interim latitude array
        rlslabs (array): longitude, latitude points along the catchment boundary
        shed (string): catchment boundary from which trajectories are released
    """
    # Set up empty lists for each month of traj data files:
    janfiles = []; febfiles = []; marfiles = []
    
    # for 2014 1st Jan 00z 2015 is included and should go into December list
    for i in range(len(trajfiles)-1): # loop over trajfiles
        if trajfiles[i][16:18] == '01': # if month is January...
            janfiles.append(trajfiles[i]) # ... add to January list
        elif trajfiles[i][16:18] == '02': # if month is February...
            febfiles.append(trajfiles[i]) # etc.
        elif trajfiles[i][16:18] == '03':
            marfiles.append(trajfiles[i])
    
    # need to use the lengths of each list to split up catchcount lists
    janlen = len(janfiles); feblen = len(febfiles); marlen = len(marfiles)
    
    # split the newlabs & catchcount arrays/lists into months:
    janlabs = newlabs[:janlen,:,:]; sf_jan = sf3[:janlen,:,:] # January
    jan_atl = atlcount[:janlen]; jan_ind = indcount[:janlen]
    jan_pac = paccount[:janlen]; jan_arc = arccount[:janlen]
    jan_sou = soucount[:janlen]
    oldlen = janlen; newlen = oldlen + feblen
    feblabs = newlabs[oldlen:newlen,:,:]; sf_feb = sf3[oldlen:newlen,:,:] # February
    feb_atl = atlcount[oldlen:newlen]; feb_ind = indcount[oldlen:newlen]
    feb_pac = paccount[oldlen:newlen]; feb_arc = arccount[oldlen:newlen]
    feb_sou = soucount[oldlen:newlen]
    oldlen = newlen; newlen = oldlen + marlen
    marlabs = newlabs[oldlen:newlen,:,:]; sf_mar = sf3[oldlen:newlen,:,:] # March
    mar_atl = atlcount[oldlen:newlen]; mar_ind = indcount[oldlen:newlen]
    mar_pac = paccount[oldlen:newlen]; mar_arc = arccount[oldlen:newlen]
    mar_sou = soucount[oldlen:newlen]
    
    yr = trajfiles[0][12:16]
    # January:
    mnt = janfiles[0][16:18]
    # all catchments:
    all_DA01, all_FA01 = MainKernelFunc2(janlabs,eralon,eralat,sf_jan,rlslabs,'all',shed,'all',yr,mnt)
    all_D101, all_F101 = MainKernelFunc2(janlabs,eralon,eralat,sf_jan,rlslabs,'1',shed,'all',yr,mnt)
    all_D201, all_F201 = MainKernelFunc2(janlabs,eralon,eralat,sf_jan,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA01, atl_FA01 = MainKernelFunc(jan_atl,eralon,eralat,sf_jan,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D101, atl_F101 = MainKernelFunc(jan_atl,eralon,eralat,sf_jan,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D201, atl_F201 = MainKernelFunc(jan_atl,eralon,eralat,sf_jan,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA01, ind_FA01 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D101, ind_F101 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D201, ind_F201 = MainKernelFunc(jan_ind,eralon,eralat,sf_jan,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA01, pac_FA01 = MainKernelFunc(jan_pac,eralon,eralat,sf_jan,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D101, pac_F101 = MainKernelFunc(jan_pac,eralon,eralat,sf_jan,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D201, pac_F201 = MainKernelFunc(jan_pac,eralon,eralat,sf_jan,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA01, arc_FA01 = MainKernelFunc(jan_arc,eralon,eralat,sf_jan,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D101, arc_F101 = MainKernelFunc(jan_arc,eralon,eralat,sf_jan,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D201, arc_F201 = MainKernelFunc(jan_arc,eralon,eralat,sf_jan,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA01, sou_FA01 = MainKernelFunc(jan_sou,eralon,eralat,sf_jan,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D101, sou_F101 = MainKernelFunc(jan_sou,eralon,eralat,sf_jan,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D201, sou_F201 = MainKernelFunc(jan_sou,eralon,eralat,sf_jan,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********January done**********'

    # February:
    mnt = febfiles[0][16:18]
    # all catchments:
    all_DA02, all_FA02 = MainKernelFunc2(feblabs,eralon,eralat,sf_feb,rlslabs,'all',shed,'all',yr,mnt)
    all_D102, all_F102 = MainKernelFunc2(feblabs,eralon,eralat,sf_feb,rlslabs,'1',shed,'all',yr,mnt)
    all_D202, all_F202 = MainKernelFunc2(feblabs,eralon,eralat,sf_feb,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA02, atl_FA02 = MainKernelFunc(feb_atl,eralon,eralat,sf_feb,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D102, atl_F102 = MainKernelFunc(feb_atl,eralon,eralat,sf_feb,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D202, atl_F202 = MainKernelFunc(feb_atl,eralon,eralat,sf_feb,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA02, ind_FA02 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D102, ind_F102 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D202, ind_F202 = MainKernelFunc(feb_ind,eralon,eralat,sf_feb,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA02, pac_FA02 = MainKernelFunc(feb_pac,eralon,eralat,sf_feb,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D102, pac_F102 = MainKernelFunc(feb_pac,eralon,eralat,sf_feb,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D202, pac_F202 = MainKernelFunc(feb_pac,eralon,eralat,sf_feb,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA02, arc_FA02 = MainKernelFunc(feb_arc,eralon,eralat,sf_feb,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D102, arc_F102 = MainKernelFunc(feb_arc,eralon,eralat,sf_feb,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D202, arc_F202 = MainKernelFunc(feb_arc,eralon,eralat,sf_feb,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA02, sou_FA02 = MainKernelFunc(feb_sou,eralon,eralat,sf_feb,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D102, sou_F102 = MainKernelFunc(feb_sou,eralon,eralat,sf_feb,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D202, sou_F202 = MainKernelFunc(feb_sou,eralon,eralat,sf_feb,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********February done**********'

    # March:
    mnt = marfiles[0][16:18]
    # all catchments:
    all_DA03, all_FA03 = MainKernelFunc2(marlabs,eralon,eralat,sf_mar,rlslabs,'all',shed,'all',yr,mnt)
    all_D103, all_F103 = MainKernelFunc2(marlabs,eralon,eralat,sf_mar,rlslabs,'1',shed,'all',yr,mnt)
    all_D203, all_F203 = MainKernelFunc2(marlabs,eralon,eralat,sf_mar,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA03, atl_FA03 = MainKernelFunc(mar_atl,eralon,eralat,sf_mar,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D103, atl_F103 = MainKernelFunc(mar_atl,eralon,eralat,sf_mar,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D203, atl_F203 = MainKernelFunc(mar_atl,eralon,eralat,sf_mar,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA03, ind_FA03 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D103, ind_F103 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D203, ind_F203 = MainKernelFunc(mar_ind,eralon,eralat,sf_mar,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA03, pac_FA03 = MainKernelFunc(mar_pac,eralon,eralat,sf_mar,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D103, pac_F103 = MainKernelFunc(mar_pac,eralon,eralat,sf_mar,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D203, pac_F203 = MainKernelFunc(mar_pac,eralon,eralat,sf_mar,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA03, arc_FA03 = MainKernelFunc(mar_arc,eralon,eralat,sf_mar,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D103, arc_F103 = MainKernelFunc(mar_arc,eralon,eralat,sf_mar,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D203, arc_F203 = MainKernelFunc(mar_arc,eralon,eralat,sf_mar,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA03, sou_FA03 = MainKernelFunc(mar_sou,eralon,eralat,sf_mar,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D103, sou_F103 = MainKernelFunc(mar_sou,eralon,eralat,sf_mar,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D203, sou_F203 = MainKernelFunc(mar_sou,eralon,eralat,sf_mar,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********March done**********'
    
    # Save these fields to netcdf files:
    # January:
    # all catchments:
    SaveFields(janfiles,all_DA01,all_FA01,shed,'all','all')
    SaveFields(janfiles,all_D101,all_F101,shed,'all','1')
    SaveFields(janfiles,all_D201,all_F201,shed,'all','2')
    # Atlantic
    SaveFields(janfiles,atl_DA01,atl_FA01,shed,'atl','all')
    SaveFields(janfiles,atl_D101,atl_F101,shed,'atl','1')
    SaveFields(janfiles,atl_D201,atl_F201,shed,'atl','2')
    # Indian
    SaveFields(janfiles,ind_DA01,ind_FA01,shed,'ind','all')
    SaveFields(janfiles,ind_D101,ind_F101,shed,'ind','1')
    SaveFields(janfiles,ind_D201,ind_F201,shed,'ind','2')
    # Pacific
    SaveFields(janfiles,pac_DA01,pac_FA01,shed,'pac','all')
    SaveFields(janfiles,pac_D101,pac_F101,shed,'pac','1')
    SaveFields(janfiles,pac_D201,pac_F201,shed,'pac','2')
    # Arctic
    SaveFields(janfiles,arc_DA01,arc_FA01,shed,'arc','all')
    SaveFields(janfiles,arc_D101,arc_F101,shed,'arc','1')
    SaveFields(janfiles,arc_D201,arc_F201,shed,'arc','2')
    # Southern
    SaveFields(janfiles,sou_DA01,sou_FA01,shed,'sou','all')
    SaveFields(janfiles,sou_D101,sou_F101,shed,'sou','1')
    SaveFields(janfiles,sou_D201,sou_F201,shed,'sou','2')
    
    # February:
    # all catchments:
    SaveFields(febfiles,all_DA02,all_FA02,shed,'all','all')
    SaveFields(febfiles,all_D102,all_F102,shed,'all','1')
    SaveFields(febfiles,all_D202,all_F202,shed,'all','2')
    # Atlantic
    SaveFields(febfiles,atl_DA02,atl_FA02,shed,'atl','all')
    SaveFields(febfiles,atl_D102,atl_F102,shed,'atl','1')
    SaveFields(febfiles,atl_D202,atl_F202,shed,'atl','2')
    # Indian
    SaveFields(febfiles,ind_DA02,ind_FA02,shed,'ind','all')
    SaveFields(febfiles,ind_D102,ind_F102,shed,'ind','1')
    SaveFields(febfiles,ind_D202,ind_F202,shed,'ind','2')
    # Pacific
    SaveFields(febfiles,pac_DA02,pac_FA02,shed,'pac','all')
    SaveFields(febfiles,pac_D102,pac_F102,shed,'pac','1')
    SaveFields(febfiles,pac_D202,pac_F202,shed,'pac','2')
    # Arctic
    SaveFields(febfiles,arc_DA02,arc_FA02,shed,'arc','all')
    SaveFields(febfiles,arc_D102,arc_F102,shed,'arc','1')
    SaveFields(febfiles,arc_D202,arc_F202,shed,'arc','2')
    # Southern
    SaveFields(febfiles,sou_DA02,sou_FA02,shed,'sou','all')
    SaveFields(febfiles,sou_D102,sou_F102,shed,'sou','1')
    SaveFields(febfiles,sou_D202,sou_F202,shed,'sou','2')
    
    # March:
    # all catchments:
    SaveFields(marfiles,all_DA03,all_FA03,shed,'all','all')
    SaveFields(marfiles,all_D103,all_F103,shed,'all','1')
    SaveFields(marfiles,all_D203,all_F203,shed,'all','2')
    # Atlantic
    SaveFields(marfiles,atl_DA03,atl_FA03,shed,'atl','all')
    SaveFields(marfiles,atl_D103,atl_F103,shed,'atl','1')
    SaveFields(marfiles,atl_D203,atl_F203,shed,'atl','2')
    # Indian
    SaveFields(marfiles,ind_DA03,ind_FA03,shed,'ind','all')
    SaveFields(marfiles,ind_D103,ind_F103,shed,'ind','1')
    SaveFields(marfiles,ind_D203,ind_F203,shed,'ind','2')
    # Pacific
    SaveFields(marfiles,pac_DA03,pac_FA03,shed,'pac','all')
    SaveFields(marfiles,pac_D103,pac_F103,shed,'pac','1')
    SaveFields(marfiles,pac_D203,pac_F203,shed,'pac','2')
    # Arctic
    SaveFields(marfiles,arc_DA03,arc_FA03,shed,'arc','all')
    SaveFields(marfiles,arc_D103,arc_F103,shed,'arc','1')
    SaveFields(marfiles,arc_D203,arc_F203,shed,'arc','2')
    # Southern
    SaveFields(marfiles,sou_DA03,sou_FA03,shed,'sou','all')
    SaveFields(marfiles,sou_D103,sou_F103,shed,'sou','1')
    SaveFields(marfiles,sou_D203,sou_F203,shed,'sou','2')
    
    print '*******netcdf files for ' + shed + ' JFM written************'


def Kernels_netcdf_AMJ(newlabs,atlcount,indcount,paccount,arccount,soucount,
                                   trajfiles,sf3,eralon,eralat,rlslabs,shed):
    """Function to call the kernel functions and write all the monthly data to
    netcdf files for April, May and June.
    
    Args:
        newlabs (array): origin category, timestep & co-ordinates of every
                         trajectory for every release
        atlcount (list): trajectories with Atlantic origin, contains trajectory
                         number, origin type/timestep/co-ordinates
        indcount (list): as above for the Indian Ocean
        paccount (list): as above for the Pacific Ocean
        indcount (list): as above for the Arctic Ocean
        soucount (list): as above for the Southern Ocean
        trajfiles (list): names of trajectory data files
        sf3 (array): cross-section of moisture flux on catchment boundary
        eralon (array): ERA-Interim longitude array
        eralat (array): ERA-Interim latitude array
        rlslabs (array): longitude, latitude points along the catchment boundary
        shed (string): catchment boundary from which trajectories are released
    """
    # Set up empty lists for each month of traj data files:
    aprfiles = []; mayfiles = []; junfiles = []
    
    for i in range(len(trajfiles)): # loop over trajfiles
        if trajfiles[i][16:18] == '04':
            aprfiles.append(trajfiles[i])
        elif trajfiles[i][16:18] == '05':
            mayfiles.append(trajfiles[i])
        elif trajfiles[i][16:18] == '06':
            junfiles.append(trajfiles[i])
    
    aprlen = len(aprfiles); maylen = len(mayfiles); junlen = len(junfiles)
    x = pl.where(trajfiles==aprfiles[0])
    
    #oldlen = newlen; newlen = oldlen + aprlen
    aprlabs = newlabs[x[0][0]:x[0][0]+aprlen,:,:]; sf_apr = sf3[x[0][0]:x[0][0]+aprlen,:,:] # April
    apr_atl = atlcount[x[0][0]:x[0][0]+aprlen]; apr_ind = indcount[x[0][0]:x[0][0]+aprlen]
    apr_pac = paccount[x[0][0]:x[0][0]+aprlen]; apr_arc = arccount[x[0][0]:x[0][0]+aprlen]
    apr_sou = soucount[x[0][0]:x[0][0]+aprlen]
    oldlen = x[0][0]+aprlen; newlen = oldlen + maylen
    maylabs = newlabs[oldlen:newlen,:,:]; sf_may = sf3[oldlen:newlen,:,:] # May
    may_atl = atlcount[oldlen:newlen]; may_ind = indcount[oldlen:newlen]
    may_pac = paccount[oldlen:newlen]; may_arc = arccount[oldlen:newlen]
    may_sou = soucount[oldlen:newlen]
    oldlen = newlen; newlen = oldlen + junlen
    junlabs = newlabs[oldlen:newlen,:,:]; sf_jun = sf3[oldlen:newlen,:,:] # June
    jun_atl = atlcount[oldlen:newlen]; jun_ind = indcount[oldlen:newlen]
    jun_pac = paccount[oldlen:newlen]; jun_arc = arccount[oldlen:newlen]
    jun_sou = soucount[oldlen:newlen]
    
    yr = trajfiles[0][12:16]
    # April:
    mnt = aprfiles[0][16:18]
    # all catchments:
    all_DA04, all_FA04 = MainKernelFunc2(aprlabs,eralon,eralat,sf_apr,rlslabs,'all',shed,'all',yr,mnt)
    all_D104, all_F104 = MainKernelFunc2(aprlabs,eralon,eralat,sf_apr,rlslabs,'1',shed,'all',yr,mnt)
    all_D204, all_F204 = MainKernelFunc2(aprlabs,eralon,eralat,sf_apr,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA04, atl_FA04 = MainKernelFunc(apr_atl,eralon,eralat,sf_apr,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D104, atl_F104 = MainKernelFunc(apr_atl,eralon,eralat,sf_apr,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D204, atl_F204 = MainKernelFunc(apr_atl,eralon,eralat,sf_apr,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA04, ind_FA04 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D104, ind_F104 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D204, ind_F204 = MainKernelFunc(apr_ind,eralon,eralat,sf_apr,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA04, pac_FA04 = MainKernelFunc(apr_pac,eralon,eralat,sf_apr,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D104, pac_F104 = MainKernelFunc(apr_pac,eralon,eralat,sf_apr,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D204, pac_F204 = MainKernelFunc(apr_pac,eralon,eralat,sf_apr,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA04, arc_FA04 = MainKernelFunc(apr_arc,eralon,eralat,sf_apr,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D104, arc_F104 = MainKernelFunc(apr_arc,eralon,eralat,sf_apr,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D204, arc_F204 = MainKernelFunc(apr_arc,eralon,eralat,sf_apr,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA04, sou_FA04 = MainKernelFunc(apr_sou,eralon,eralat,sf_apr,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D104, sou_F104 = MainKernelFunc(apr_sou,eralon,eralat,sf_apr,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D204, sou_F204 = MainKernelFunc(apr_sou,eralon,eralat,sf_apr,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********April done**********'

    # May:
    mnt = mayfiles[0][16:18]
    # all catchments:
    all_DA05, all_FA05 = MainKernelFunc2(maylabs,eralon,eralat,sf_may,rlslabs,'all',shed,'all',yr,mnt)
    all_D105, all_F105 = MainKernelFunc2(maylabs,eralon,eralat,sf_may,rlslabs,'1',shed,'all',yr,mnt)
    all_D205, all_F205 = MainKernelFunc2(maylabs,eralon,eralat,sf_may,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA05, atl_FA05 = MainKernelFunc(may_atl,eralon,eralat,sf_may,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D105, atl_F105 = MainKernelFunc(may_atl,eralon,eralat,sf_may,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D205, atl_F205 = MainKernelFunc(may_atl,eralon,eralat,sf_may,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA05, ind_FA05 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D105, ind_F105 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D205, ind_F205 = MainKernelFunc(may_ind,eralon,eralat,sf_may,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA05, pac_FA05 = MainKernelFunc(may_pac,eralon,eralat,sf_may,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D105, pac_F105 = MainKernelFunc(may_pac,eralon,eralat,sf_may,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D205, pac_F205 = MainKernelFunc(may_pac,eralon,eralat,sf_may,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic;
    arc_DA05, arc_FA05 = MainKernelFunc(may_arc,eralon,eralat,sf_may,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D105, arc_F105 = MainKernelFunc(may_arc,eralon,eralat,sf_may,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D205, arc_F205 = MainKernelFunc(may_arc,eralon,eralat,sf_may,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA05, sou_FA05 = MainKernelFunc(may_sou,eralon,eralat,sf_may,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D105, sou_F105 = MainKernelFunc(may_sou,eralon,eralat,sf_may,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D205, sou_F205 = MainKernelFunc(may_sou,eralon,eralat,sf_may,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********May done**********'

    # June:
    mnt = junfiles[0][16:18]
    # all catchments:
    all_DA06, all_FA06 = MainKernelFunc2(junlabs,eralon,eralat,sf_jun,rlslabs,'all',shed,'all',yr,mnt)
    all_D106, all_F106 = MainKernelFunc2(junlabs,eralon,eralat,sf_jun,rlslabs,'1',shed,'all',yr,mnt)
    all_D206, all_F206 = MainKernelFunc2(junlabs,eralon,eralat,sf_jun,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA06, atl_FA06 = MainKernelFunc(jun_atl,eralon,eralat,sf_jun,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D106, atl_F106 = MainKernelFunc(jun_atl,eralon,eralat,sf_jun,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D206, atl_F206 = MainKernelFunc(jun_atl,eralon,eralat,sf_jun,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA06, ind_FA06 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D106, ind_F106 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D206, ind_F206 = MainKernelFunc(jun_ind,eralon,eralat,sf_jun,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA06, pac_FA06 = MainKernelFunc(jun_pac,eralon,eralat,sf_jun,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D106, pac_F106 = MainKernelFunc(jun_pac,eralon,eralat,sf_jun,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D206, pac_F206 = MainKernelFunc(jun_pac,eralon,eralat,sf_jun,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA06, arc_FA06 = MainKernelFunc(jun_arc,eralon,eralat,sf_jun,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D106, arc_F106 = MainKernelFunc(jun_arc,eralon,eralat,sf_jun,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D206, arc_F206 = MainKernelFunc(jun_arc,eralon,eralat,sf_jun,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA06, sou_FA06 = MainKernelFunc(jun_sou,eralon,eralat,sf_jun,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D106, sou_F106 = MainKernelFunc(jun_sou,eralon,eralat,sf_jun,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D206, sou_F206 = MainKernelFunc(jun_sou,eralon,eralat,sf_jun,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********June done**********'
    
    # April:
    # all catchments:
    SaveFields(aprfiles,all_DA04,all_FA04,shed,'all','all')
    SaveFields(aprfiles,all_D104,all_F104,shed,'all','1')
    SaveFields(aprfiles,all_D204,all_F204,shed,'all','2')
    # Atlantic
    SaveFields(aprfiles,atl_DA04,atl_FA04,shed,'atl','all')
    SaveFields(aprfiles,atl_D104,atl_F104,shed,'atl','1')
    SaveFields(aprfiles,atl_D204,atl_F204,shed,'atl','2')
    # Indian
    SaveFields(aprfiles,ind_DA04,ind_FA04,shed,'ind','all')
    SaveFields(aprfiles,ind_D104,ind_F104,shed,'ind','1')
    SaveFields(aprfiles,ind_D204,ind_F204,shed,'ind','2')
    # Pacific
    SaveFields(aprfiles,pac_DA04,pac_FA04,shed,'pac','all')
    SaveFields(aprfiles,pac_D104,pac_F104,shed,'pac','1')
    SaveFields(aprfiles,pac_D204,pac_F204,shed,'pac','2')
    # Arctic
    SaveFields(aprfiles,arc_DA04,arc_FA04,shed,'arc','all')
    SaveFields(aprfiles,arc_D104,arc_F104,shed,'arc','1')
    SaveFields(aprfiles,arc_D204,arc_F204,shed,'arc','2')
    # Southern
    SaveFields(aprfiles,sou_DA04,sou_FA04,shed,'sou','all')
    SaveFields(aprfiles,sou_D104,sou_F104,shed,'sou','1')
    SaveFields(aprfiles,sou_D204,sou_F204,shed,'sou','2')
    
    # May:
    # all catchments:
    SaveFields(mayfiles,all_DA05,all_FA05,shed,'all','all')
    SaveFields(mayfiles,all_D105,all_F105,shed,'all','1')
    SaveFields(mayfiles,all_D205,all_F205,shed,'all','2')
    # Atlantic
    SaveFields(mayfiles,atl_DA05,atl_FA05,shed,'atl','all')
    SaveFields(mayfiles,atl_D105,atl_F105,shed,'atl','1')
    SaveFields(mayfiles,atl_D205,atl_F205,shed,'atl','2')
    # Indian
    SaveFields(mayfiles,ind_DA05,ind_FA05,shed,'ind','all')
    SaveFields(mayfiles,ind_D105,ind_F105,shed,'ind','1')
    SaveFields(mayfiles,ind_D205,ind_F205,shed,'ind','2')
    # Pacific
    SaveFields(mayfiles,pac_DA05,pac_FA05,shed,'pac','all')
    SaveFields(mayfiles,pac_D105,pac_F105,shed,'pac','1')
    SaveFields(mayfiles,pac_D205,pac_F205,shed,'pac','2')
    # Arctic
    SaveFields(mayfiles,arc_DA05,arc_FA05,shed,'arc','all')
    SaveFields(mayfiles,arc_D105,arc_F105,shed,'arc','1')
    SaveFields(mayfiles,arc_D205,arc_F205,shed,'arc','2')
    # Southern
    SaveFields(mayfiles,sou_DA05,sou_FA05,shed,'sou','all')
    SaveFields(mayfiles,sou_D105,sou_F105,shed,'sou','1')
    SaveFields(mayfiles,sou_D205,sou_F205,shed,'sou','2')
    
    # June:
    # all catchments:
    SaveFields(junfiles,all_DA06,all_FA06,shed,'all','all')
    SaveFields(junfiles,all_D106,all_F106,shed,'all','1')
    SaveFields(junfiles,all_D206,all_F206,shed,'all','2')
    # Atlantic
    SaveFields(junfiles,atl_DA06,atl_FA06,shed,'atl','all')
    SaveFields(junfiles,atl_D106,atl_F106,shed,'atl','1')
    SaveFields(junfiles,atl_D206,atl_F206,shed,'atl','2')
    # Indian
    SaveFields(junfiles,ind_DA06,ind_FA06,shed,'ind','all')
    SaveFields(junfiles,ind_D106,ind_F106,shed,'ind','1')
    SaveFields(junfiles,ind_D206,ind_F206,shed,'ind','2')
    # Pacific
    SaveFields(junfiles,pac_DA06,pac_FA06,shed,'pac','all')
    SaveFields(junfiles,pac_D106,pac_F106,shed,'pac','1')
    SaveFields(junfiles,pac_D206,pac_F206,shed,'pac','2')
    # Arctic
    SaveFields(junfiles,arc_DA06,arc_FA06,shed,'arc','all')
    SaveFields(junfiles,arc_D106,arc_F106,shed,'arc','1')
    SaveFields(junfiles,arc_D206,arc_F206,shed,'arc','2')
    # Southern
    SaveFields(junfiles,sou_DA06,sou_FA06,shed,'sou','all')
    SaveFields(junfiles,sou_D106,sou_F106,shed,'sou','1')
    SaveFields(junfiles,sou_D206,sou_F206,shed,'sou','2')
    
    print '*******netcdf files for ' + shed + ' AMJ written************'

def Kernels_netcdf_JAS(newlabs,atlcount,indcount,paccount,arccount,soucount,
                                   trajfiles,sf3,eralon,eralat,rlslabs,shed):
    """Function to call the kernel functions and write all the monthly data to
    netcdf files for July, August and September.
    
    Args:
        newlabs (array): origin category, timestep & co-ordinates of every
                         trajectory for every release
        atlcount (list): trajectories with Atlantic origin, contains trajectory
                         number, origin type/timestep/co-ordinates
        indcount (list): as above for the Indian Ocean
        paccount (list): as above for the Pacific Ocean
        indcount (list): as above for the Arctic Ocean
        soucount (list): as above for the Southern Ocean
        trajfiles (list): names of trajectory data files
        sf3 (array): cross-section of moisture flux on catchment boundary
        eralon (array): ERA-Interim longitude array
        eralat (array): ERA-Interim latitude array
        rlslabs (array): longitude, latitude points along the catchment boundary
        shed (string): catchment boundary from which trajectories are released
    """
    # Set up empty lists for each month of traj data files:
    julfiles = []; augfiles = []; sepfiles = []
    
    for i in range(len(trajfiles)): # loop over trajfiles
        if trajfiles[i][16:18] == '07':
            julfiles.append(trajfiles[i])
        elif trajfiles[i][16:18] == '08':
            augfiles.append(trajfiles[i])
        elif trajfiles[i][16:18] == '09':
            sepfiles.append(trajfiles[i])

    jullen = len(julfiles); auglen = len(augfiles); seplen = len(sepfiles)
    x = pl.where(trajfiles==julfiles[0])
    
    #oldlen = newlen; newlen = oldlen + aprlen
    jullabs = newlabs[x[0][0]:x[0][0]+jullen,:,:]; sf_jul = sf3[x[0][0]:x[0][0]+jullen,:,:] # April
    jul_atl = atlcount[x[0][0]:x[0][0]+jullen]; jul_ind = indcount[x[0][0]:x[0][0]+jullen]
    jul_pac = paccount[x[0][0]:x[0][0]+jullen]; jul_arc = arccount[x[0][0]:x[0][0]+jullen]
    jul_sou = soucount[x[0][0]:x[0][0]+jullen]
    oldlen = jullen+x[0][0]; newlen = oldlen + auglen
    auglabs = newlabs[oldlen:newlen,:,:]; sf_aug = sf3[oldlen:newlen,:,:] # May
    aug_atl = atlcount[oldlen:newlen]; aug_ind = indcount[oldlen:newlen]
    aug_pac = paccount[oldlen:newlen]; aug_arc = arccount[oldlen:newlen]
    aug_sou = soucount[oldlen:newlen]
    oldlen = newlen; newlen = oldlen + seplen
    seplabs = newlabs[oldlen:newlen,:,:]; sf_sep = sf3[oldlen:newlen,:,:] # June
    sep_atl = atlcount[oldlen:newlen]; sep_ind = indcount[oldlen:newlen]
    sep_pac = paccount[oldlen:newlen]; sep_arc = arccount[oldlen:newlen]
    sep_sou = soucount[oldlen:newlen]
    
    yr = trajfiles[0][12:16]
    # July:
    mnt = julfiles[0][16:18]
    # all catchments:
    all_DA07, all_FA07 = MainKernelFunc2(jullabs,eralon,eralat,sf_jul,rlslabs,'all',shed,'all',yr,mnt)
    all_D107, all_F107 = MainKernelFunc2(jullabs,eralon,eralat,sf_jul,rlslabs,'1',shed,'all',yr,mnt)
    all_D207, all_F207 = MainKernelFunc2(jullabs,eralon,eralat,sf_jul,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA07, atl_FA07 = MainKernelFunc(jul_atl,eralon,eralat,sf_jul,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D107, atl_F107 = MainKernelFunc(jul_atl,eralon,eralat,sf_jul,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D207, atl_F207 = MainKernelFunc(jul_atl,eralon,eralat,sf_jul,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA07, ind_FA07 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D107, ind_F107 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D207, ind_F207 = MainKernelFunc(jul_ind,eralon,eralat,sf_jul,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA07, pac_FA07 = MainKernelFunc(jul_pac,eralon,eralat,sf_jul,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D107, pac_F107 = MainKernelFunc(jul_pac,eralon,eralat,sf_jul,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D207, pac_F207 = MainKernelFunc(jul_pac,eralon,eralat,sf_jul,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA07, arc_FA07 = MainKernelFunc(jul_arc,eralon,eralat,sf_jul,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D107, arc_F107 = MainKernelFunc(jul_arc,eralon,eralat,sf_jul,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D207, arc_F207 = MainKernelFunc(jul_arc,eralon,eralat,sf_jul,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA07, sou_FA07 = MainKernelFunc(jul_sou,eralon,eralat,sf_jul,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D107, sou_F107 = MainKernelFunc(jul_sou,eralon,eralat,sf_jul,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D207, sou_F207 = MainKernelFunc(jul_sou,eralon,eralat,sf_jul,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********July done**********'

    # August:
    mnt = augfiles[0][16:18]
    # all catchments:
    all_DA08, all_FA08 = MainKernelFunc2(auglabs,eralon,eralat,sf_aug,rlslabs,'all',shed,'all',yr,mnt)
    all_D108, all_F108 = MainKernelFunc2(auglabs,eralon,eralat,sf_aug,rlslabs,'1',shed,'all',yr,mnt)
    all_D208, all_F208 = MainKernelFunc2(auglabs,eralon,eralat,sf_aug,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA08, atl_FA08 = MainKernelFunc(aug_atl,eralon,eralat,sf_aug,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D108, atl_F108 = MainKernelFunc(aug_atl,eralon,eralat,sf_aug,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D208, atl_F208 = MainKernelFunc(aug_atl,eralon,eralat,sf_aug,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA08, ind_FA08 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D108, ind_F108 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D208, ind_F208 = MainKernelFunc(aug_ind,eralon,eralat,sf_aug,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA08, pac_FA08 = MainKernelFunc(aug_pac,eralon,eralat,sf_aug,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D108, pac_F108 = MainKernelFunc(aug_pac,eralon,eralat,sf_aug,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D208, pac_F208 = MainKernelFunc(aug_pac,eralon,eralat,sf_aug,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA08, arc_FA08 = MainKernelFunc(aug_arc,eralon,eralat,sf_aug,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D108, arc_F108 = MainKernelFunc(aug_arc,eralon,eralat,sf_aug,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D208, arc_F208 = MainKernelFunc(aug_arc,eralon,eralat,sf_aug,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA08, sou_FA08 = MainKernelFunc(aug_sou,eralon,eralat,sf_aug,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D108, sou_F108 = MainKernelFunc(aug_sou,eralon,eralat,sf_aug,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D208, sou_F208 = MainKernelFunc(aug_sou,eralon,eralat,sf_aug,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********August done**********'

    # September:
    mnt = sepfiles[0][16:18]
    # all catchments:
    all_DA09, all_FA09 = MainKernelFunc2(seplabs,eralon,eralat,sf_sep,rlslabs,'all',shed,'all',yr,mnt)
    all_D109, all_F109 = MainKernelFunc2(seplabs,eralon,eralat,sf_sep,rlslabs,'1',shed,'all',yr,mnt)
    all_D209, all_F209 = MainKernelFunc2(seplabs,eralon,eralat,sf_sep,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA09, atl_FA09 = MainKernelFunc(sep_atl,eralon,eralat,sf_sep,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D109, atl_F109 = MainKernelFunc(sep_atl,eralon,eralat,sf_sep,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D209, atl_F209 = MainKernelFunc(sep_atl,eralon,eralat,sf_sep,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA09, ind_FA09 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D109, ind_F109 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D209, ind_F209 = MainKernelFunc(sep_ind,eralon,eralat,sf_sep,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA09, pac_FA09 = MainKernelFunc(sep_pac,eralon,eralat,sf_sep,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D109, pac_F109 = MainKernelFunc(sep_pac,eralon,eralat,sf_sep,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D209, pac_F209 = MainKernelFunc(sep_pac,eralon,eralat,sf_sep,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA09, arc_FA09 = MainKernelFunc(sep_arc,eralon,eralat,sf_sep,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D109, arc_F109 = MainKernelFunc(sep_arc,eralon,eralat,sf_sep,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D209, arc_F209 = MainKernelFunc(sep_arc,eralon,eralat,sf_sep,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA09, sou_FA09 = MainKernelFunc(sep_sou,eralon,eralat,sf_sep,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D109, sou_F109 = MainKernelFunc(sep_sou,eralon,eralat,sf_sep,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D209, sou_F209 = MainKernelFunc(sep_sou,eralon,eralat,sf_sep,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********September done**********'
    
    # July:
    # all catchments:
    SaveFields(julfiles,all_DA07,all_FA07,shed,'all','all')
    SaveFields(julfiles,all_D107,all_F107,shed,'all','1')
    SaveFields(julfiles,all_D207,all_F207,shed,'all','2')
    # Atlantic
    SaveFields(julfiles,atl_DA07,atl_FA07,shed,'atl','all')
    SaveFields(julfiles,atl_D107,atl_F107,shed,'atl','1')
    SaveFields(julfiles,atl_D207,atl_F207,shed,'atl','2')
    # Indian
    SaveFields(julfiles,ind_DA07,ind_FA07,shed,'ind','all')
    SaveFields(julfiles,ind_D107,ind_F107,shed,'ind','1')
    SaveFields(julfiles,ind_D207,ind_F207,shed,'ind','2')
    # Pacific
    SaveFields(julfiles,pac_DA07,pac_FA07,shed,'pac','all')
    SaveFields(julfiles,pac_D107,pac_F107,shed,'pac','1')
    SaveFields(julfiles,pac_D207,pac_F207,shed,'pac','2')
    # Arctic
    SaveFields(julfiles,arc_DA07,arc_FA07,shed,'arc','all')
    SaveFields(julfiles,arc_D107,arc_F107,shed,'arc','1')
    SaveFields(julfiles,arc_D207,arc_F207,shed,'arc','2')
    # Southern
    SaveFields(julfiles,sou_DA07,sou_FA07,shed,'sou','all')
    SaveFields(julfiles,sou_D107,sou_F107,shed,'sou','1')
    SaveFields(julfiles,sou_D207,sou_F207,shed,'sou','2')
    
    # August:
    # all catchments:
    SaveFields(augfiles,all_DA08,all_FA08,shed,'all','all')
    SaveFields(augfiles,all_D108,all_F108,shed,'all','1')
    SaveFields(augfiles,all_D208,all_F208,shed,'all','2')
    # Atlantic
    SaveFields(augfiles,atl_DA08,atl_FA08,shed,'atl','all')
    SaveFields(augfiles,atl_D108,atl_F108,shed,'atl','1')
    SaveFields(augfiles,atl_D208,atl_F208,shed,'atl','2')
    # Indian
    SaveFields(augfiles,ind_DA08,ind_FA08,shed,'ind','all')
    SaveFields(augfiles,ind_D108,ind_F108,shed,'ind','1')
    SaveFields(augfiles,ind_D208,ind_F208,shed,'ind','2')
    # Pacific
    SaveFields(augfiles,pac_DA08,pac_FA08,shed,'pac','all')
    SaveFields(augfiles,pac_D108,pac_F108,shed,'pac','1')
    SaveFields(augfiles,pac_D208,pac_F208,shed,'pac','2')
    # Arctic
    SaveFields(augfiles,arc_DA08,arc_FA08,shed,'arc','all')
    SaveFields(augfiles,arc_D108,arc_F108,shed,'arc','1')
    SaveFields(augfiles,arc_D208,arc_F208,shed,'arc','2')
    # Southern
    SaveFields(augfiles,sou_DA08,sou_FA08,shed,'sou','all')
    SaveFields(augfiles,sou_D108,sou_F108,shed,'sou','1')
    SaveFields(augfiles,sou_D208,sou_F208,shed,'sou','2')
    
    # September:
    # all catchments:
    SaveFields(sepfiles,all_DA09,all_FA09,shed,'all','all')
    SaveFields(sepfiles,all_D109,all_F109,shed,'all','1')
    SaveFields(sepfiles,all_D209,all_F209,shed,'all','2')
    # Atlantic
    SaveFields(sepfiles,atl_DA09,atl_FA09,shed,'atl','all')
    SaveFields(sepfiles,atl_D109,atl_F109,shed,'atl','1')
    SaveFields(sepfiles,atl_D209,atl_F209,shed,'atl','2')
    # Indian
    SaveFields(sepfiles,ind_DA09,ind_FA09,shed,'ind','all')
    SaveFields(sepfiles,ind_D109,ind_F109,shed,'ind','1')
    SaveFields(sepfiles,ind_D209,ind_F209,shed,'ind','2')
    # Pacific
    SaveFields(sepfiles,pac_DA09,pac_FA09,shed,'pac','all')
    SaveFields(sepfiles,pac_D109,pac_F109,shed,'pac','1')
    SaveFields(sepfiles,pac_D209,pac_F209,shed,'pac','2')
    # Arctic
    SaveFields(sepfiles,arc_DA09,arc_FA09,shed,'arc','all')
    SaveFields(sepfiles,arc_D109,arc_F109,shed,'arc','1')
    SaveFields(sepfiles,arc_D209,arc_F209,shed,'arc','2')
    # Southern
    SaveFields(sepfiles,sou_DA09,sou_FA09,shed,'sou','all')
    SaveFields(sepfiles,sou_D109,sou_F109,shed,'sou','1')
    SaveFields(sepfiles,sou_D209,sou_F209,shed,'sou','2')
    
    print '*******netcdf files for ' + shed + ' JAS written************'

def Kernels_netcdf_OND(newlabs,atlcount,indcount,paccount,arccount,soucount,
                                   trajfiles,sf3,eralon,eralat,rlslabs,shed):
    """Function to call the kernel functions and write all the monthly data to
    netcdf files for October, November and December.
    
    Args:
        newlabs (array): origin category, timestep & co-ordinates of every
                         trajectory for every release
        atlcount (list): trajectories with Atlantic origin, contains trajectory
                         number, origin type/timestep/co-ordinates
        indcount (list): as above for the Indian Ocean
        paccount (list): as above for the Pacific Ocean
        indcount (list): as above for the Arctic Ocean
        soucount (list): as above for the Southern Ocean
        trajfiles (list): names of trajectory data files
        sf3 (array): cross-section of moisture flux on catchment boundary
        eralon (array): ERA-Interim longitude array
        eralat (array): ERA-Interim latitude array
        rlslabs (array): longitude, latitude points along the catchment boundary
        shed (string): catchment boundary from which trajectories are released
    """
    # Set up empty lists for each month of traj data files:
    octfiles = []; novfiles = []; decfiles = []
    
    for i in range(len(trajfiles)-1): # loop over trajfiles
        if trajfiles[i][16:18] == '10':
            octfiles.append(trajfiles[i])
        elif trajfiles[i][16:18] == '11':
            novfiles.append(trajfiles[i])
        elif trajfiles[i][16:18] == '12':
            decfiles.append(trajfiles[i])
    
    decfiles.append(trajfiles[-1])
    
    octlen = len(octfiles); novlen = len(novfiles); declen = len(decfiles)
    x = pl.where(trajfiles==octfiles[0])
    
    octlabs = newlabs[x[0][0]:x[0][0]+octlen,:,:]; sf_oct = sf3[x[0][0]:x[0][0]+octlen,:,:] # April
    oct_atl = atlcount[x[0][0]:x[0][0]+octlen]; oct_ind = indcount[x[0][0]:x[0][0]+octlen]
    oct_pac = paccount[x[0][0]:x[0][0]+octlen]; oct_arc = arccount[x[0][0]:x[0][0]+octlen]
    oct_sou = soucount[x[0][0]:x[0][0]+octlen]
    oldlen = octlen+x[0][0]; newlen = oldlen + novlen
    novlabs = newlabs[oldlen:newlen,:,:]; sf_nov = sf3[oldlen:newlen,:,:] # May
    nov_atl = atlcount[oldlen:newlen]; nov_ind = indcount[oldlen:newlen]
    nov_pac = paccount[oldlen:newlen]; nov_arc = arccount[oldlen:newlen]
    nov_sou = soucount[oldlen:newlen]
    oldlen = newlen; newlen = oldlen + declen
    declabs = newlabs[oldlen:newlen,:,:]; sf_dec = sf3[oldlen:newlen,:,:] # June
    dec_atl = atlcount[oldlen:newlen]; dec_ind = indcount[oldlen:newlen]
    dec_pac = paccount[oldlen:newlen]; dec_arc = arccount[oldlen:newlen]
    dec_sou = soucount[oldlen:newlen]
    
    yr = trajfiles[0][12:16]
    # October:
    mnt = octfiles[0][16:18]
    # all catchments:
    all_DA10, all_FA10 = MainKernelFunc2(octlabs,eralon,eralat,sf_oct,rlslabs,'all',shed,'all',yr,mnt)
    all_D110, all_F110 = MainKernelFunc2(octlabs,eralon,eralat,sf_oct,rlslabs,'1',shed,'all',yr,mnt)
    all_D210, all_F210 = MainKernelFunc2(octlabs,eralon,eralat,sf_oct,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA10, atl_FA10 = MainKernelFunc(oct_atl,eralon,eralat,sf_oct,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D110, atl_F110 = MainKernelFunc(oct_atl,eralon,eralat,sf_oct,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D210, atl_F210 = MainKernelFunc(oct_atl,eralon,eralat,sf_oct,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA10, ind_FA10 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D110, ind_F110 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D210, ind_F210 = MainKernelFunc(oct_ind,eralon,eralat,sf_oct,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA10, pac_FA10 = MainKernelFunc(oct_pac,eralon,eralat,sf_oct,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D110, pac_F110 = MainKernelFunc(oct_pac,eralon,eralat,sf_oct,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D210, pac_F210 = MainKernelFunc(oct_pac,eralon,eralat,sf_oct,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA10, arc_FA10 = MainKernelFunc(oct_arc,eralon,eralat,sf_oct,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D110, arc_F110 = MainKernelFunc(oct_arc,eralon,eralat,sf_oct,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D210, arc_F210 = MainKernelFunc(oct_arc,eralon,eralat,sf_oct,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA10, sou_FA10 = MainKernelFunc(oct_sou,eralon,eralat,sf_oct,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D110, sou_F110 = MainKernelFunc(oct_sou,eralon,eralat,sf_oct,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D210, sou_F210 = MainKernelFunc(oct_sou,eralon,eralat,sf_oct,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********October done**********'

    # November:
    mnt = novfiles[0][16:18]
    # all catchments:
    all_DA11, all_FA11 = MainKernelFunc2(novlabs,eralon,eralat,sf_nov,rlslabs,'all',shed,'all',yr,mnt)
    all_D111, all_F111 = MainKernelFunc2(novlabs,eralon,eralat,sf_nov,rlslabs,'1',shed,'all',yr,mnt)
    all_D211, all_F211 = MainKernelFunc2(novlabs,eralon,eralat,sf_nov,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA11, atl_FA11 = MainKernelFunc(nov_atl,eralon,eralat,sf_nov,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D111, atl_F111 = MainKernelFunc(nov_atl,eralon,eralat,sf_nov,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D211, atl_F211 = MainKernelFunc(nov_atl,eralon,eralat,sf_nov,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA11, ind_FA11 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D111, ind_F111 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D211, ind_F211 = MainKernelFunc(nov_ind,eralon,eralat,sf_nov,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA11, pac_FA11 = MainKernelFunc(nov_pac,eralon,eralat,sf_nov,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D111, pac_F111 = MainKernelFunc(nov_pac,eralon,eralat,sf_nov,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D211, pac_F211 = MainKernelFunc(nov_pac,eralon,eralat,sf_nov,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA11, arc_FA11 = MainKernelFunc(nov_arc,eralon,eralat,sf_nov,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D111, arc_F111 = MainKernelFunc(nov_arc,eralon,eralat,sf_nov,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D211, arc_F211 = MainKernelFunc(nov_arc,eralon,eralat,sf_nov,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA11, sou_FA11 = MainKernelFunc(nov_sou,eralon,eralat,sf_nov,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D111, sou_F111 = MainKernelFunc(nov_sou,eralon,eralat,sf_nov,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D211, sou_F211 = MainKernelFunc(nov_sou,eralon,eralat,sf_nov,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********November done**********'

    # December:
    mnt = decfiles[0][16:18]
    # all catchments:
    all_DA12, all_FA12 = MainKernelFunc2(declabs,eralon,eralat,sf_dec,rlslabs,'all',shed,'all',yr,mnt)
    all_D112, all_F112 = MainKernelFunc2(declabs,eralon,eralat,sf_dec,rlslabs,'1',shed,'all',yr,mnt)
    all_D212, all_F212 = MainKernelFunc2(declabs,eralon,eralat,sf_dec,rlslabs,'2',shed,'all',yr,mnt)
    # Atlantic:
    atl_DA12, atl_FA12 = MainKernelFunc(dec_atl,eralon,eralat,sf_dec,rlslabs,'all',shed,'atl',yr,mnt)
    atl_D112, atl_F112 = MainKernelFunc(dec_atl,eralon,eralat,sf_dec,rlslabs,'1',shed,'atl',yr,mnt)
    atl_D212, atl_F212 = MainKernelFunc(dec_atl,eralon,eralat,sf_dec,rlslabs,'2',shed,'atl',yr,mnt)
    # Indian:
    ind_DA12, ind_FA12 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'all',shed,'ind',yr,mnt)
    ind_D112, ind_F112 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'1',shed,'ind',yr,mnt)
    ind_D212, ind_F212 = MainKernelFunc(dec_ind,eralon,eralat,sf_dec,rlslabs,'2',shed,'ind',yr,mnt)
    # Pacific:
    pac_DA12, pac_FA12 = MainKernelFunc(dec_pac,eralon,eralat,sf_dec,rlslabs,'all',shed,'pac',yr,mnt)
    pac_D112, pac_F112 = MainKernelFunc(dec_pac,eralon,eralat,sf_dec,rlslabs,'1',shed,'pac',yr,mnt)
    pac_D212, pac_F212 = MainKernelFunc(dec_pac,eralon,eralat,sf_dec,rlslabs,'2',shed,'pac',yr,mnt)
    # Arctic:
    arc_DA12, arc_FA12 = MainKernelFunc(dec_arc,eralon,eralat,sf_dec,rlslabs,'all',shed,'arc',yr,mnt)
    arc_D112, arc_F112 = MainKernelFunc(dec_arc,eralon,eralat,sf_dec,rlslabs,'1',shed,'arc',yr,mnt)
    arc_D212, arc_F212 = MainKernelFunc(dec_arc,eralon,eralat,sf_dec,rlslabs,'2',shed,'arc',yr,mnt)
    # Southern:
    sou_DA12, sou_FA12 = MainKernelFunc(dec_sou,eralon,eralat,sf_dec,rlslabs,'all',shed,'sou',yr,mnt)
    sou_D112, sou_F112 = MainKernelFunc(dec_sou,eralon,eralat,sf_dec,rlslabs,'1',shed,'sou',yr,mnt)
    sou_D212, sou_F212 = MainKernelFunc(dec_sou,eralon,eralat,sf_dec,rlslabs,'2',shed,'sou',yr,mnt)
    #print '**********December done**********'
    
    # October:
    # all catchments:
    SaveFields(octfiles,all_DA10,all_FA10,shed,'all','all')
    SaveFields(octfiles,all_D110,all_F110,shed,'all','1')
    SaveFields(octfiles,all_D210,all_F210,shed,'all','2')
    # Atlantic
    SaveFields(octfiles,atl_DA10,atl_FA10,shed,'atl','all')
    SaveFields(octfiles,atl_D110,atl_F110,shed,'atl','1')
    SaveFields(octfiles,atl_D210,atl_F210,shed,'atl','2')
    # Indian
    SaveFields(octfiles,ind_DA10,ind_FA10,shed,'ind','all')
    SaveFields(octfiles,ind_D110,ind_F110,shed,'ind','1')
    SaveFields(octfiles,ind_D210,ind_F210,shed,'ind','2')
    # Pacific
    SaveFields(octfiles,pac_DA10,pac_FA10,shed,'pac','all')
    SaveFields(octfiles,pac_D110,pac_F110,shed,'pac','1')
    SaveFields(octfiles,pac_D210,pac_F210,shed,'pac','2')
    # Arctic
    SaveFields(octfiles,arc_DA10,arc_FA10,shed,'arc','all')
    SaveFields(octfiles,arc_D110,arc_F110,shed,'arc','1')
    SaveFields(octfiles,arc_D210,arc_F210,shed,'arc','2')
    # Southern
    SaveFields(octfiles,sou_DA10,sou_FA10,shed,'sou','all')
    SaveFields(octfiles,sou_D110,sou_F110,shed,'sou','1')
    SaveFields(octfiles,sou_D210,sou_F210,shed,'sou','2')
    
    # November:
    # all catchments:
    SaveFields(novfiles,all_DA11,all_FA11,shed,'all','all')
    SaveFields(novfiles,all_D111,all_F111,shed,'all','1')
    SaveFields(novfiles,all_D211,all_F211,shed,'all','2')
    # Atlantic
    SaveFields(novfiles,atl_DA11,atl_FA11,shed,'atl','all')
    SaveFields(novfiles,atl_D111,atl_F111,shed,'atl','1')
    SaveFields(novfiles,atl_D211,atl_F211,shed,'atl','2')
    # Indian
    SaveFields(novfiles,ind_DA11,ind_FA11,shed,'ind','all')
    SaveFields(novfiles,ind_D111,ind_F111,shed,'ind','1')
    SaveFields(novfiles,ind_D211,ind_F211,shed,'ind','2')
    # Pacific
    SaveFields(novfiles,pac_DA11,pac_FA11,shed,'pac','all')
    SaveFields(novfiles,pac_D111,pac_F111,shed,'pac','1')
    SaveFields(novfiles,pac_D211,pac_F211,shed,'pac','2')
    # Arctic
    SaveFields(novfiles,arc_DA11,arc_FA11,shed,'arc','all')
    SaveFields(novfiles,arc_D111,arc_F111,shed,'arc','1')
    SaveFields(novfiles,arc_D211,arc_F211,shed,'arc','2')
    # Southern
    SaveFields(novfiles,sou_DA11,sou_FA11,shed,'sou','all')
    SaveFields(novfiles,sou_D111,sou_F111,shed,'sou','1')
    SaveFields(novfiles,sou_D211,sou_F211,shed,'sou','2')
    
    # December:
    # all catchments:
    SaveFields(decfiles,all_DA12,all_FA12,shed,'all','all')
    SaveFields(decfiles,all_D112,all_F112,shed,'all','1')
    SaveFields(decfiles,all_D212,all_F212,shed,'all','2')
    # Atlantic
    SaveFields(decfiles,atl_DA12,atl_FA12,shed,'atl','all')
    SaveFields(decfiles,atl_D112,atl_F112,shed,'atl','1')
    SaveFields(decfiles,atl_D212,atl_F212,shed,'atl','2')
    # Indian
    SaveFields(decfiles,ind_DA12,ind_FA12,shed,'ind','all')
    SaveFields(decfiles,ind_D112,ind_F112,shed,'ind','1')
    SaveFields(decfiles,ind_D212,ind_F212,shed,'ind','2')
    # Pacific
    SaveFields(decfiles,pac_DA12,pac_FA12,shed,'pac','all')
    SaveFields(decfiles,pac_D112,pac_F112,shed,'pac','1')
    SaveFields(decfiles,pac_D212,pac_F212,shed,'pac','2')
    # Arctic
    SaveFields(decfiles,arc_DA12,arc_FA12,shed,'arc','all')
    SaveFields(decfiles,arc_D112,arc_F112,shed,'arc','1')
    SaveFields(decfiles,arc_D212,arc_F212,shed,'arc','2')
    # Southern
    SaveFields(decfiles,sou_DA12,sou_FA12,shed,'sou','all')
    SaveFields(decfiles,sou_D112,sou_F112,shed,'sou','1')
    SaveFields(decfiles,sou_D212,sou_F212,shed,'sou','2')
    
    print '*******netcdf files for ' + shed + ' OND written************'
