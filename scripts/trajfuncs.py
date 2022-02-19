# -*- coding: utf-8 -*-
"""
Created on Sat Jun 24 15:30:00 2017

@author: np838619
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap
from os import listdir
from os.path import isfile, join
import itertools
from netCDF4 import Dataset
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.path as mplPath
from scipy.stats import pearsonr
from scipy.interpolate import interp1d

exec(open('/home/users/np838619/PminusE_data/ERA_Int/functions.py').read())


def ListFiles(directory):
    """Function to list all the files in a directory and sort them.
	
	Args:
		directory (string): path of directory containing files

	Returns:
		list_of_files (array): sorted array of all filenames in the directory
    """
    list_of_files = [f for f in listdir(directory) if isfile(join(directory, f))]
    list_of_files = pl.sort(list_of_files)
    
    return list_of_files
    	
def ReadTxtFile(filename):
	"""Function to read a .txt file and make every line of the file an element of
	a list.

	Args:
		filename (string): path to .txt file
	
	Returns:
		flist (list): contents of file
	"""
	f = open(filename,'r')
	flist = []
	for line in f.readlines():
		flist.append(line.split())
	f.close()
	
	return flist
	

def TrajExtract(NSTEPS,data,flist,trajno,intervals):
    """Function to extract trajectory data.

	Args:
		NSTEPS (int): number of time steps
		data (int): number of fields for each trajectory
		flist (list): entire contents of trajectory file, each line is a list
				   within a list
		trajno (float): trajectory number

	Returns:
		trajdat (array): all trajectory data
    """
    trajdat = pl.zeros([NSTEPS,data])
    for step in range(NSTEPS):
        trajdat[step] = flist[18+step+trajno*(intervals+5)]
    
    return trajdat

def InitData(initfiles,pingdir):
    """Function to read in trajectory data from initialization files with the 
    trajectory data at t=0.
    
    Args:
        initfiles (array): filenames of files containing trajectory data at t=0
        pingdir (string): directory containing trajectory data files
        
    Returns:
        U (array): zonal velocity at t=0
        V (array): meridional velocity at t=0
    """
    global ReadTxtFile; global itertools
    U = []; V = []; SP = []
    for initfile in range(len(initfiles)): # loop over initfiles array
        # read in data file as list, each line is also a list
        flist = ReadTxtFile(pingdir + initfiles[initfile])
        
        if initfile == 0.: # info that doesn't change from file to file
            NPART = int(flist[2][-1]) # no. of trajectories
            NATTR = int(flist[3][-1]) + 2 # no. of attributes
            NCLUST = int(flist[6][-1]) # number of levels
            pointers = list(itertools.chain(*flist[8:10])) # traj no where each level starts
            pointers = pl.asarray(pointers,dtype='int')
            data = len(flist[16])-2 + NATTR # how much data for each trajectory
        
        initdata = pl.zeros([NPART,data+1]) # empty array for traj data
        
        # extract traj data from file:
        for traj in range(NPART): # loop over all trajectories
            initdata[traj,:-1] = flist[17+traj*6] # start from line 18 (zero indexing)
                                            # then increment by traj*5
            initdata[traj,-1] = flist[18+traj*6][0]
        
        # split data into sub-arrays
        lat = initdata[:,2]#; lon = initdata[:,3]
#        pres = initdata[:,4]; temp = initdata[:,5]
#        PV = initdata[:,6]; q = initdata[:,7]
#        height = initdata[:,8]; blz = initdata[:,9]
        u = initdata[:,10]; v = initdata[:,11]; sp = initdata[:,-1]
        
        NREL = NPART/NCLUST # number of release points
        # change shape of arrays to be no.of levels X no. of release points
        u = pl.reshape(u,(NCLUST,NREL)); v = pl.reshape(v,(NCLUST,NREL))
        lat = pl.reshape(lat,(NCLUST,NREL)); sp = pl.reshape(sp,(NCLUST,NREL))
        
        a = 6.37e6 # radius of Earth
        # u & v output in units 1/day, convert to m/s:
            # multiply by -1 as it is back trajectory:
        u = -1*u*a*pl.cos(pl.radians(lat))/86400; v = -1*v*a/86400
        
        U.append(u); V.append(v); SP.append(sp) # append to lists
    
    U = pl.asarray(U); V = pl.asarray(V); SP = pl.asarray(SP) # transform list to array

    return U, V, SP

def TrajPlot(m,lon,lat,originlabs):
    """Function to plot trajectories

	Args:
		lon (array): longitudes of trajectories at all time steps
		lat (array): latitudes of trajectories at all time steps
            originlabs (array): origin type; timestep of origin; lon/lat of origin
            
    """
    ortype = originlabs[0]
    if ortype == 1.: clr = 'b'
    elif ortype == 2.: clr = 'r'
    elif ortype == 0: return
    # Make map:
    if ortype == 1. or ortype == 2.:
        #m = Basemap(projection='nplaea',boundinglat=0.01,lon_0=270,resolution='l')
        #m.drawcoastlines(color='grey',linewidth=0.5,zorder=10)
        #x,y = m(lon,lat); m.plot(x,y,color=clr)
        xor,yor = m(originlabs[3],originlabs[2])
        m.plot(xor,yor,color=clr,marker='^',zorder=10)

def BoundaryLayer(blz,height,step,lat,lon):
    """Function to determine if and when trajectories have been in contact with the 
	boundary layer, and to calculate how many trajectories from one release have been
	in contact with the boundary layer.

	Args:
		blz (array): boundary layer heights of all trajectories at all time steps
		height (array): height of particle for all trajectories at all time steps
		step (array): time steps of all trajectories
		lat (array): latitudes along trajectory at all time steps
		lon (array): longitudes along trajectory at all time steps

	Returns:
		blcount (float): number of trajectories from release that have made
						contact with the boundary layer
		blprop (float): proportion of trajectories which have made contact with
						the boundary layer, expressed as a percentage
		blorigin (array): array containing trajectory label, timestep of BL
                                         origin, latitude/longitude of origin
    """
    NPART = blz.shape[0]
    blcount = 0.; blorigin = []
    for traj in range(NPART):
        # find where BL height exceeds particle height
        bl_inds = pl.where(blz[traj,:]>height[traj,:])
        if bl_inds[0].size > 0.: # if bl_inds size >0 then particle was in BL
            blcount = blcount + 1.
            # store traj no. & origin index, subtract 1 from index last in BL
            # as particle has still been in BL after this point
            if bl_inds[0][0] == 0.:
                blorigin.append((traj,step[traj,bl_inds[0][0]],#-1
                         lat[traj,bl_inds[0][0]],lon[traj,bl_inds[0][0]]))
            else:
                blorigin.append((traj,step[traj,bl_inds[0][0]-1],
                         lat[traj,bl_inds[0][0]-1],lon[traj,bl_inds[0][0]-1]))
        elif bl_inds[0].size == 0.:
            blcount = blcount
    
    blprop = (blcount/NPART)*100
    blorigin = pl.asarray(blorigin)    
    
    return blcount, blprop, blorigin

def BoundaryLayer2(blz,height,step,lat,lon,lsm,eralon,eralat):
    """Function to determine if and when trajectories have been in contact with the 
	boundary layer, and to calculate how many trajectories from one release have been
	in contact with the boundary layer.

	Args:
		blz (array): boundary layer heights of all trajectories at all time steps
		height (array): height of particle for all trajectories at all time steps
		step (array): time steps of all trajectories
		lat (array): latitudes along trajectory at all time steps
		lon (array): longitudes along trajectory at all time steps

	Returns:
		blcount (float): number of trajectories from release that have made
						contact with the boundary layer
		blprop (float): proportion of trajectories which have made contact with
						the boundary layer, expressed as a percentage
		blorigin (array): array containing trajectory label, timestep of BL
                                         origin, latitude/longitude of origin
    """
    NPART = blz.shape[0]
    blcount = 0.; blorigin = []
    for traj in range(NPART):
        # find where BL height exceeds particle height
        bl_inds = pl.where(blz[traj,:]>height[traj,:])
        if bl_inds[0].size > 0.: # if bl_inds size >0 then particle was in BL
            S = pl.zeros([bl_inds[0].size])
            for i in range(bl_inds[0].size):
                S[i] = BilinInterp((lon[traj,bl_inds[0][i]],lat[traj,bl_inds[0][i]]),
                                                            eralon,eralat,lsm)
            sea_inds = pl.where(S<0.5)
            if sea_inds[0].size > 0.:
                blcount = blcount + 1.
                #if bl_inds[0][sea_inds[0][0]] == 0.:
                blorigin.append((traj,step[traj,bl_inds[0][sea_inds[0][0]]],
                     lat[traj,bl_inds[0][sea_inds[0][0]]],
                                   lon[traj,bl_inds[0][sea_inds[0][0]]]))
                #else:
                #    blorigin.append((traj,step[traj,bl_inds[0][sea_inds[0][0]]-1],
                 #        lat[traj,bl_inds[0][sea_inds[0][0]]-1],
                  #                  lon[traj,bl_inds[0][sea_inds[0][0]]-1]))
            elif sea_inds[0].size == 0.:
                blcount = blcount
            # store traj no. & origin index, subtract 1 from index last in BL
            # as particle has still been in BL after this point
            #elif
            #    blorigin.append((traj,step[traj,bl_inds[0][0]],#-1
            #                 lat[traj,bl_inds[0][0]],lon[traj,bl_inds[0][0]]))
      #      elif sea_inds[0].size == 0.:
      #          if bl_inds[0][0] == 0.:
      #              blorigin.append((traj,step[traj,bl_inds[0][0]],#-1
      #                       lat[traj,bl_inds[0][0]],lon[traj,bl_inds[0][0]]))
      #          else:
      #              blorigin.append((traj,step[traj,bl_inds[0][0]-1],
      #                       lat[traj,bl_inds[0][0]-1],lon[traj,bl_inds[0][0]-1]))
        elif bl_inds[0].size == 0.:
            blcount = blcount
    
    blprop = (blcount/NPART)*100
    blorigin = pl.asarray(blorigin)    
    
    return blcount, blprop, blorigin

def OriginScatter(blorigin):
    """Function to plot the boundary layer origin on a map
    
    Args:
        blorigin (array): trajectory number, timestep of BL origin & latitude/
        longitude of origin location
        
    Returns:
        m (Basemap object): map of trajectory region
    """
    m = Basemap(projection='mill',resolution='l',llcrnrlat=-80.,urcrnrlat=80.,
                llcrnrlon=0.,urcrnrlon=360.,#lon_0=lon.mean(),lat_0=lat.mean(),
                lat_ts=20)
    #m.drawcoastlines()
    a,b = m(blorigin[:,3],blorigin[:,2]) # longitude, latitude
    m.scatter(a,b,c=blorigin[:,1]/24,cmap='viridis')#,linewidths=0.)
    cbar = m.colorbar(location='bottom')
    cbar.set_label('Time along back trajectories (days)',fontsize=16)
    # send the Basemap object to main code, draw coastlines in main code
    return m

def AttrScatter(originlabs,attr,ortype):
    """Function to to make a scatter plot on a map of trajectory origin locations
    coloured by an attribute. To be used as part of a subplot so object returned
    to main code to help make colour bar which will be the same for both maps.
    
    Args:
        originlabs (array): label denoting if trajectory has bl origin ('1'),
                            partial mixing origin ('2') or no origin ('0');
                            timestep of origin (int or None); lat/lon of origin
        attr (array): attribute along all trajectories at all timesteps
        ortype (int): origin label used to plot specific trajectories
    
    Returns:
        pmc (matplotlib.collections.PathCollection object): used to create
                        colour bar shared by subplots outside of function
    """
    if ortype == 1. or ortype == 2.: # check if valid origin type requested
        x = pl.where(originlabs[:,0]==ortype) # find relevant trajectories
        B = pl.zeros([x[0].size,3]) # make new array for attribute & lat/lon
        
        for trajno in range(x[0].size): # loop over relevant trajectories
            tn = x[0][trajno] # get trajectory number
            step = originlabs[tn,1] # get timestep of origin
            B[trajno,0] = attr[tn,step] # extract attribute value at origin
            B[trajno,1:] = originlabs[tn,2:] # extract lon/lat of origin
        
        # make world map:
        m = Basemap(projection='mill',resolution='l',llcrnrlat=-80.,urcrnrlat=80.,
                 llcrnrlon=0.,urcrnrlon=360.,lat_ts=20)
        m.drawcoastlines()
        a,b = m(B[:,2],B[:,1]) # convert lon.lat to projection co-ordinates
        # plot origin locations & colour by attribute
        pcm = m.scatter(a,b,c=B[:,0],cmap='viridis',norm=pl.Normalize(0,5))
    else: # if ortype entered wrong, raise exception
        raise Exception('Wrong origin type. Must be 1 or 2!')
    
    return pcm#cbar

def TrajDayCount(origin,NPART,NSTEPS):
    """Function to calculate the number & percentage of trajectories from one 
    release that have BL or PM origin at each time step.
    
    Args:
        origins (array): trajectory number, timestep (hours) of origin, lat/lon
                         of origin
        NPART (int): number of trajectories
        NSTEPS (int): number of time steps
    
    Returns:
        orprop (array): percentage of trajectories which have origin at each
                        timestep
    """
    # set up empty array for timestep of trajectory origin:
    orloc = pl.zeros([NPART,NSTEPS])
    for trajno in range(origin.shape[0]):
        # if a traj has origin, set orloc to 1 at the relevant timestep:
        orloc[origin[trajno,0],origin[trajno,1]/6] = 1.
    
    # number of new trajectories with origin at each timestep:
    newor = pl.sum(orloc,axis=0)
    # cumulative number of trajectories with origin at each timestep: 
    ordaycount = pl.zeros_like(newor)
    for step in range(len(newor)):
        ordaycount[step] = pl.sum(newor[:step+1])
    
    # calculate proportions of trajectories with origin at each timestep:
    orprop = (ordaycount/NPART)*100
    
    return orprop

def BLoriginFile(pingdir,BASETIME,filename,origin_data):
    """Function to write a file containing boundary layer origin data i.e.
    the trajectory number, time step of origin & latitude/longitude of origin.
    
    Args:
        pingdir (str): path to directory of trajectory data files
        BASETIME: start time of back trajectories
        filename (str): name of trajectory data file
        origin_data (array): boundary layer origin data
    """
    newfile = open(pingdir+'blorigin'+str(BASETIME)+'.txt','w')
    newfile.write('Trajectory number, time step of BL origin & latitude/longitude of'
    + ' origin of 10 day back trajectories which have boundary layer origin in'
    + ' file ' + filename + '\nTRAJECTORIES = ' + str(origin_data.shape[0]) + 
    '       BASETIME = ' + str(BASETIME) + '\n\n')
    newfile.write('TRAJNO     STEP     LAT     LON\n')
    #origin_data.tofile(newfile,sep="   ",format='%3f')
    pl.savetxt(newfile,origin_data,fmt='%3f',delimiter='    ')
    newfile.close()

def VMagnitude(phi,lambd,dt):
    """Function to calculate the magnitude of the 2D velocity vector.
	
	Args:
		phi (array): latitudes of trajectory locations
		lambd (array): longitudes of trajectory locations
		dt (int): timestep between trajectory locations
	
	Returns:
		V (float): magnitude of 2D velocity vector
    """
    dphi = pl.radians(phi[0]-phi[1])
    dlam = pl.radians(lambd[0]+360)-pl.radians(lambd[1])
    phi0 = pl.radians(phi[0])
    R = 6.37*10**6 # radius of the Earth
    V = pl.sqrt((R**2)*(dphi**2) + (R**2)*(pl.cos(phi0)**2)*(dlam**2))/(dt*60*60)
    
    return V

def TrajSegLabel(loc):
    """Function to assign each trajectory release point a label referring to which
	segment of the watershed it is from.

	Args:
		loc (string): stem of continental watershed e.g. NCA = North/Central America
	
	Returns:
             seglab (array): segment labels (integers) of each release point
             rlspts (array): trajectory release points
    """
    sheddir = '/home/users/np838619/Watershed/shed_defs/'
    
    #endlist = ReadTxtFile(sheddir + loc + '_clicks.txt')
    #endpts = pl.asarray(endlist[5:],dtype='float')
    
    rlslist = ReadTxtFile(sheddir + loc + '_traj_release.txt')
    rlspts = pl.asarray(rlslist[5:],dtype='float')
    
    seglab = [1] # first release point has to be on the first segment
    #segnos = pl.linspace(1,40,40)
    count = 1
    
    for rls in range(1,rlspts.shape[0]):
        if rlspts[rls,0] == rlspts[rls-1,0] and rlspts[rls,1] == rlspts[rls-1,1]:
            count = count + 1
            seglab.append(count)
        else:
            count = count
            seglab.append(count)
    seglab = pl.asarray(seglab)
    
    return seglab, rlspts

def RemoveRepeated(interp_pts,labels):
    """Function to remove repeated points at ends of segments from catchment 
    boundary - neccessary because net flux will be estimated incorrectly as it 
    would otherwise be calculated along too many points.
    
    Args:
        interp_pts (array): points along a catchment boundary, interpolated 
                            between segment end points
        labels (array): segment labels of each point along catchment boundary
    """
    repeat = []
    for r in range(1,interp_pts.shape[0]):
        if interp_pts[r,0] == interp_pts[r-1,0] and interp_pts[r,1] == interp_pts[r-1,1]:
            repeat.append(r)
    
    RP2 = []; L2 = []
    for r in range(interp_pts.shape[0]):
        if r not in repeat:
            RP2.append(interp_pts[r])
            L2.append(labels[r])
    rlspts = pl.asarray(RP2); labs = pl.asarray(L2)
    
    return rlspts, labs

def LocalCartesian(coords):
    """Function to convert latitude-longitude co-ordinates into local
	Cartesian co-ordinates.

	Args:
		co-ords (array): Nx2 array of co-ordinates in degrees

	Returns:
		loccar (array): Nx2 array of local Cartesian co-ordinates
    """
    coords = pl.radians(coords) # convert from degrees to radians
    loccar = pl.zeros_like(coords) # creat empty array same size as coords
    R = 6.37*10**6 # radius of the Earth
	# x = R*cos(lat)*lon
    loccar[:,0] = R*pl.sin(coords[:,1])*pl.cos(coords[:,0])
    loccar[:,1] = R*pl.sin(coords[:,0])*pl.sin(coords[:,1]) # y = R*lat
    
    return loccar

def SegLength(sheddir,loc):
    """Function to calculate segment lengths.
	
	Args:
		sheddir (string): path to directory with end points file
		loc (string): stem indicating which watershed is being used

	Returns:
		seglengths (array): length of segments in metres
    """
	# read in end points file:
    endlist = ReadTxtFile(sheddir + loc + '_clicks.txt')
    endpts = pl.asarray(endlist[5:],dtype='float') # convert list to array
    
    #endpts[:,0] = endpts[:,0] + 360.
    endpts = LocalCartesian(endpts) # convert from lon-lat to local Cartesian co-ords
    end2 = pl.zeros_like(endpts)
    #endpts = m(endpts[:,0],endpts[:,1])
    end2[:,0] = endpts[:,0]; end2[:,1] = endpts[:,1]
    
    seglenths = pl.zeros(len(end2)-1) # empty array for seg lengths
    for seg in range(seglenths.shape[0]):
		# calculate lengths with Pythagoras:
        seglenths[seg] = pl.sqrt((end2[seg+1,0]-end2[seg,0])**2 + 
                                    (end2[seg+1,1]-end2[seg,1])**2)
    
    return seglenths
            
def NormalVector(sheddir,loc):
    """Function to find the normal vectors to the line segments along a watershed.

	Args:
		sheddir (string): path to directory with end points file
		loc (string): stem indicating which watershed is being used

	Returns:
		n (array): unit normal vectors to line segments
    """
    endlist = ReadTxtFile(sheddir + loc + '_clicks.txt')
    endpts = pl.asarray(endlist[5:],dtype='float')
    endpts = Add360(endpts)
    
    # convert endpts array to local cartesian co-ordinates:
    #endpts[:,0] = endpts[:,0] + 360.
    endpts = pl.radians(endpts)#loccar = LocalCartesian(endpts)
    R = 6.37*10**6
    
    #loccar = pl.zeros_like(endpts)
    #loccar[:,0] = endpts[:,0]; loccar[:,1] = endpts[:,1]
    
    nhat = pl.zeros([endpts.shape[0]-1,2])
    for point in range(endpts.shape[0]-1):
        if endpts[point+1,1] == endpts[point,1]: # same lat. co-ordinate ...
            nhat[point] = pl.array([0,1]) # ... northward pointing normal
        elif endpts[point+1,0] == endpts[point,0]: # same lon. co-ordinate ...
            nhat[point] = pl.array([1,0]) # ... eastward pointing normal
        else:
            # dx = R*cos(lat)*d(lon)
            dx = R*pl.cos(endpts[point,1])*(endpts[point+1,0]-endpts[point,0])
            # dy = R*d(lat)
            dy = R*(endpts[point+1,1]-endpts[point,1])
            n = pl.array([-dy,dx])
            nhat[point] = n/pl.norm(n) # normalize

    return nhat

def TrajLab(seglabs,NCLUST,NPART):
    """Function to label each trajectory with a segment number

	Args:
		seglabs (array): segment labels
		NCLUST (int): number of vertical levels
		NPART (int) number of trajectories
    """
    NSEG = seglabs.max() # highest values is number of segments
    NREL = NPART/NCLUST # number of release points
    
    trajlab = pl.zeros([NCLUST,NREL]) # emptry array for traj labs
    
    for seg in range(1,NSEG+1): # loop over segments
        # where does current seg no. match seg label?
        labs = pl.where(seglabs==seg)
        for lev in range(trajlab.shape[0]): # loop over levels
			# assign seg label to trajectories
            trajlab[lev,labs[0][0]:labs[0][-1]+1] = seg
    
    return trajlab

def LevelSplit(var,NCLUST,NPART,NSTEPS,pointers):
    """Function to split variables up into vertical levels

	Args:
		var (array): variable requiring splitting
		NCLUST (int): number of vertical levels
		NPART (int): number of trajectories
		NSTEPS (int): number of time steps
		pointers (array): trajnos where level changes

	Returns:
		varsplt (array): variable split into vertical levels
    """
	# empty array for split variable:
    varsplt = pl.zeros([NCLUST,NPART/NCLUST,NSTEPS])
    pointers = pointers - 1. # change to 0 based indexing
    
    for lev in range(NCLUST-1): # loop over levels
        varsplt[lev] = var[pointers[lev]:pointers[lev+1]]
    varsplt[-1] = var[pointers[-1]:]
    
    return varsplt


def FluxCalc(rlslabs,u0,v0,surfp_tj,etalevs,normals,q_traj,pres_traj):
    """Function to calculate the vertically intehrated moisture flux at each 
    trajectory release point.
            (1/g)*|(uq,vq).n*dl*deta
    
    Args:
        rlslabs (array): longitude/latitude of trajectory release point and label
                        indicating which segment of watershed the point is on
        u0 (array): zonal velocity (m/s) at release point
        v0 (array): meridional velocity (m/s) at release point
        etalevs (array): model levels which trajectories are released from
        normals (array): normal vectors for each segment of watershed
        q_traj (array): specific humidity (kg/kg) along trajectories
        pres_traj (array): pressure (hPa) along trajectories
    
    Returns:
        shedflx (array): vertically integrated moisture flux at each trajectory
                        release point (kg/s)
    """
    global Haversine; global VIeta
    NPART = q_traj.shape[0] # get no. of trajectories from q_traj shape
    NREL = rlslabs.shape[0] # get no. of release points from rlslabs shape
    NCLUST = int(NPART/NREL) #  get no. of levels
    NSTEPS = q_traj.shape[1] # get no. of timesteps from q_traj shape
    
    # change shape of q & pressure arrays: no. levels X no. points X no. timesteps
    ps = pl.reshape(pres_traj,(NCLUST,NREL,NSTEPS))
    qs = pl.reshape(q_traj,(NCLUST,NREL,NSTEPS))
    
    
    # get zonal and meridional moisture fluxes on all levels at release points:
    uq = u0*qs[:,:,0]; vq = v0*qs[:,:,0]
    
    dl = pl.zeros([NREL]) # empty arrays for distance weights
    for pt in range(1,NREL-1): # exclude first & last release points
        # use Haversine formula to calculate dl
        dl[pt] = 0.5*(Haversine(rlslabs[pt-1,:2],rlslabs[pt,:2]) + 
                        Haversine(rlslabs[pt,:2],rlslabs[pt+1,:2]))
    # use first 2 points & last 2 points for first & last dl
    dl[0] = Haversine(rlslabs[0,:2],rlslabs[1,:2])
    dl[-1] = Haversine(rlslabs[-2,:2],rlslabs[-1,:2])
    
    fluxes_uv = pl.zeros([NREL,2]) # empty array for vertically integrated fluxes
    shedflx = pl.zeros([NREL]) # empty array for total fluxes
    
    for pt in range(NREL): # loop over release points
        # Integrate on eta levels
        fluxes_uv[pt,0] = VIeta(uq[:,pt],ps[:,pt,0],surfp_tj[pt],etalevs)
        fluxes_uv[pt,1] = VIeta(vq[:,pt],ps[:,pt,0],surfp_tj[pt],etalevs)
        seglab = rlslabs[pt,-1] - 1 # get label of segment
        # use dot product of moisture fluxes & normal vector, multiply by dl
        shedflx[pt] = pl.dot(fluxes_uv[pt],normals[seglab])*dl[pt]
    
    return shedflx

def FluxCalc2(rlslabs,u0,v0,surfp_tj,etalevs,normals,q_traj,pres_traj):
    """
    """
    global Haversine; global VIeta
    NPART = q_traj.shape[0] # get no. of trajectories from q_traj shape
    NREL = rlslabs.shape[0] # get no. of release points from rlslabs shape
    NCLUST = int(NPART/NREL) #  get no. of levels
    NSTEPS = q_traj.shape[1] # get no. of timesteps from q_traj shape
    
    # change shape of q & pressure arrays: no. levels X no. points X no. timesteps
    ps = pl.reshape(pres_traj,(NCLUST,NREL,NSTEPS))
    qs = pl.reshape(q_traj,(NCLUST,NREL,NSTEPS))
    
    # get zonal and meridional moisture fluxes on all levels at release points:
    uq = u0*qs[:,:,0]; vq = v0*qs[:,:,0]
    
    fluxes_uv = pl.zeros([NREL,2]) # empty array for vertically integrated fluxes
    shedflx = pl.zeros([NREL]) # empty array for total fluxes
    
    for pt in range(NREL): # loop over release points
        # Integrate on eta levels
        fluxes_uv[pt,0] = VIeta(uq[:,pt],ps[:,pt,0],surfp_tj[pt],etalevs)
        fluxes_uv[pt,1] = VIeta(vq[:,pt],ps[:,pt,0],surfp_tj[pt],etalevs)
        seglab = rlslabs[pt,-1] - 1 # get label of segment
        # use dot product of moisture fluxes & normal vector
        shedflx[pt] = pl.dot(fluxes_uv[pt],normals[seglab])

    return shedflx    

def FluxCalc3(rlslabs,u0,v0,surfp_tj,etalevs,normals,q_traj,pres_traj):
    """
    """
    global Haversine; global VIeta
    NPART = q_traj.shape[0] # get no. of trajectories from q_traj shape
    NREL = rlslabs.shape[0] # get no. of release points from rlslabs shape
    NCLUST = int(NPART/NREL) #  get no. of levels
    NSTEPS = q_traj.shape[1] # get no. of timesteps from q_traj shape
    
    # change shape of q & pressure arrays: no. levels X no. points X no. timesteps
    ps = pl.reshape(pres_traj,(NCLUST,NREL,NSTEPS))
    qs = pl.reshape(q_traj,(NCLUST,NREL,NSTEPS))
    
    # get zonal and meridional moisture fluxes on all levels at release points:
    uq = u0*qs[:,:,0]; vq = v0*qs[:,:,0]
    
    fluxes_uv = pl.zeros([NCLUST,NREL,2]) # empty array for vertically integrated fluxes
    shedflx = pl.zeros([NCLUST,NREL]) # empty array for total fluxes
    
    for pt in range(NREL): # loop over release points
        # Integrate on eta levels
        fluxes_uv[:,pt,0] = CrossSec(uq[:,pt],ps[:,pt,0],surfp_tj[0,pt],etalevs)
        fluxes_uv[:,pt,1] = CrossSec(vq[:,pt],ps[:,pt,0],surfp_tj[0,pt],etalevs)
        seglab = rlslabs[pt,-1] - 1 # get label of segment
        # use dot product of moisture fluxes & normal vector
        shedflx[:,pt] = pl.dot(fluxes_uv[:,pt],normals[seglab])

    return shedflx 

def VIpres(variable,surfp,pressure):
    """Function to integrate vertically in pressure co-ordinates.
                (1/g)|A*dp
    
    Args:
        variable (array): quantity to be integrated
        surfp (array): surface pressure (Pa)
        pressure (array): pressure levels (hPa)
    
    Returns:
        vint (float): vertical integral of a quantity
    """
    g = 9.81 # acceleration due to gravity (m/s^2)
    #erapres = erapres*100
    pressure = pressure*100 # change units of pressure to Pascals
    
    dp = pl.zeros([pressure.shape[0]]) # empty array for dp terms
    
    if surfp == pressure[0]: # if surface pressure is the same as lowest pressure level
        # lowest dp is half of difference between 2 lowest pressure levels:
        dp[0] = 0.5*(pressure[0] - pressure[1])
        # highest dp is half of difference between 2 highest pressure levels:
        dp[-1] = 0.5*(pressure[-2] - pressure[-1])
        for p in range(1,len(dp)-1):
            # rest of dp are half of difference between pressure levels on each side:
            dp[p] = 0.5*(pressure[p-1]-pressure[p+1])
    elif surfp > pressure[0]: # if surface pressure below lowest pressure level
        # lowest dp is surfp minus average of two lowest pressure levels:
        dp[0] = surfp - 0.5*(pressure[0]+pressure[1])
        # highest dp in half of difference between 2 highest pressure levels:
        dp[-1] = 0.5*(pressure[-2] - pressure[-1])
        for p in range(1,len(dp)-1):
            # rest of dp are half of difference between pressure levels on each side:
            dp[p] = 0.5*(pressure[p-1]-pressure[p+1])
    elif surfp < pressure[0]:# if surface pressure above lowest pressure level
        # find lowest pressure level above surface pressure
        b = pl.where(surfp > pressure); b = b[0][0]
        dp[:b] = 0. # all dp below this level are 0
        # dp at this level is surface pressure minus average of next two levels
        dp[b] = surfp - 0.5*(pressure[b]+pressure[b+1])
        # highest dp is half of difference between 2 highest levels
        dp[-1] = 0.5*(pressure[-2] - pressure[-1])
        for p in range(b+1,len(dp)-1):
            # rest of dp are half of difference between pressure levels on each side
            dp[p] = 0.5*(pressure[p-1]-pressure[p+1])
    
    X = pl.zeros([pressure.size]) # empty array, same size as pressure array
                
    for p in range(len(dp)): # loop over dp
        #a = pl.where(pressure[p]==erapres)
        #lev = a[0][0]
        X[p] = variable[p]*dp[p] # multiply variable by dp weight
    
    vint = (1/g)*pl.nansum(X,axis=0) # integrate over column
    
    return vint

def VIeta(variable,pressure,surfp,etalevs):
    """Function to integrate vertically in eta co-ordinates:
                (1/g)*|A*(dp/deta)*deta
    
    Args:
        variable (array): quantity to be integrated
        pressure (array): pressure levels (hPa)
        surfp (array): surface pressure (Pa)
        etalevs (array): model levels
    
    Returns:
        vint (float): vertical integral of a quantity
    """
    g = 9.81 # acceleration due to gravity (m/s^2)
    pressure = pressure*100 # change units of pressure to Pascals
    surfp = surfp*100
    deta = etalevs[0] - etalevs[1] # should be constant
    
    dp = pl.zeros([pressure.shape[0]]) # empty array for dp terms
    # deta is actually irrelevant as it is cancelled out (see function description)
    
    for p in range(1,len(dp)-1):
        # all but 1st and last dp are half of difference between pressure levels on either side
        dp[p] = 0.5*(pressure[p-1]-pressure[p+1])
    # lowest dp is half of difference between surface pressure and lowest pressure level
    dp[0] = 0.5*(surfp-pressure[0])
    # highest dp is half of difference between 2 highest pressure levels:
    dp[-1] =  0.5*(pressure[-2] - pressure[-1])
    
    X = pl.zeros_like(pressure) # empty array, same size as pressure array
    
    for lev in range(len(etalevs)): # loop over model levels
        X[lev] = variable[lev]*dp[lev] # mutliply variable by dp weight
    
    vint = (1/g)*pl.sum(X,axis=0) # integrate over column
    
    return vint

def CrossSec(variable,pressure,surfp,etalevs):
    """
    """
    pressure = pressure*100 # change units of pressure to Pascals
    surfp = surfp*100
    deta = etalevs[0] - etalevs[1] # should be constant
    
    dp = pl.zeros([pressure.shape[0]]) # empty array for dp terms
    # deta is actually irrelevant as it is cancelled out (see function description)
    
    for p in range(1,len(dp)-1):
        # all but 1st and last dp are half of difference between pressure levels on either side
        dp[p] = 0.5*(pressure[p-1]-pressure[p+1])
    # lowest dp is half of difference between surface pressure and lowest pressure level
    dp[0] = 0.5*(surfp-pressure[0])
    # highest dp is half of difference between 2 highest pressure levels:
    dp[-1] =  0.5*(pressure[-2] - pressure[-1])
    
    X = pl.zeros_like(pressure) # empty array, same size as pressure array
    
    for lev in range(len(etalevs)): # loop over model levels
        X[lev] = variable[lev]*dp[lev] # mutliply variable by dp weight
    
    return X

def SaveCrossSec(csec,year,shed):
    """
    """
    clusdir = '/glusterfs/scenario/users/np838619/traj/'
    ncfile = Dataset(clusdir+year+'/'+shed+'/CSEC/pres'+year+'_'+shed+'.nc','w')
    
    time_dim = ncfile.createDimension('releases',csec.shape[0])
    time = ncfile.createVariable('releases',pl.float64,('releases',))
    #time.units = 'hours'
    time.long_name = 'trajectory releases'
    
    height = ncfile.createDimension('levels',csec.shape[1])
    height = ncfile.createVariable('levels',pl.float64,('levels'))
    height.long_name = 'model levels'
    
    length = ncfile.createDimension('points',csec.shape[2])
    length = ncfile.createVariable('points',pl.float64,('points'))
    length.long_name = 'release points'
    
    C = ncfile.createVariable('pressure cross-section',pl.float64,
                                      ('releases','levels','points'))
    C.units = 'hPa'
    C.standard_name = 'pressure'
    C[:,:,:] = csec[:,:,:]
    
    ncfile.close()
    
    print '***********cross-sections written to netcdf file**************'
    
    return None

def BLorigin_flux(blorigin,var,NPART):
    """Function to eliminate the specific humidities of trajectories with no BL
    origin.
    
    Args:
        blorigin (array): trajectory numbers, time steps of BL origin and
                          latitude/longitude co-ordinates of BL origin location
        q (array): specific humidities at each time step for all trajectories
        NPART (int): number of trajectories
    
    Returns:
        q_bl (array): specific humidities of all trajectories with BL origin,
                      trajectories with no BL origin have zero specific humidity
    """
    var_bl = pl.zeros_like(var) # emptry array for trajectories with BL origin
    if blorigin.size == 0.:
        var_bl = var_bl
    else:
        for traj in range(NPART):
            if traj in blorigin[:,0]: #is traj number in blorigin array?
                var_bl[traj,:] = var[traj,:] # extract trajectories with BL origin
    
    return var_bl

def OriginFlux(originlabs,var):
    """Function to eliminate the values of an attribute along trajectories with
    no boundary layer or partial mixing origin. 
    
    Args:
        originlabs (array): label denoting if trajectory has bl origin ('1'),
                            partial mixing origin ('2') or no origin ('0');
                            timestep of origin (int or None); lat/lon of origin
        var (array): attribute along trajectory e.g. specific humidity or PV
    
    Returns
        var_or (array): same shape as var, but trajectories which have no origin
                        have been eliminated
    """
    # NPART = no. of trajectories, not neccessarily true NPART:
    NPART = originlabs.shape[0]
    var_or = pl.zeros_like(var) # set up empty array for var data to be kept
    
    for traj in range(NPART): # loop over trajectories
        if originlabs[traj,0] != 0.: # if label is nonzero...
            var_or[traj,:] = var[traj,:] # ... keep data along this trajectory
    
    return var_or

def LandSeaOrigin(lsm,eralat,eralon,blorigin):
    """
    """
    orland = 0.; orocean = 0.

    for i in range(int(blcount)):
        a = NearestIndex(eralat,blorigin[i,-2])
        b = NearestIndex(eralon,blorigin[i,-1])
        orig = lsm[a,b]
        if orig > 0.5:
            orland = orland + 1.
        if orig <= 0.5:
            orocean = orocean +1.
    print orland, ' trajectories have land origin'
    print orocean, ' trajectories have ocean origin'
    
def theta_e(q,temp,pres):
    """Function to calculate equivalent potential temperature using Bolton's
    formula.
    
    Args:
        q (float): specific humidity in kg/kg
        theta (float): potential temperature in K
        pres (float): pressure in hPa
    
    Returns:
        ept (float): equivalent potential temperature in K
    """
    # convert specific humidity to mass mixing ratio:
    w = q[:]/(1-q[:])
    
    # values less than 10e-10 should be set to 10e-10
    for step in range(len(w)):
        if w[step] < 10e-10:
            w[step] = 10e-10
    
    # calculate T from theta:
    theta = pl.zeros_like(temp)
    theta[:] = temp[:]*(1000./pres[:])**(0.2854*(1.0-0.28*w[:]))   
    
    frac1 = pl.zeros_like(theta); frac2 = pl.zeros_like(frac1)
    frac1[:] = 2840./(3.5*pl.log(temp[:]) - 
                        pl.log(100*pres[:]*w[:]/(0.622+0.378*w[:])) - 0.1998)
    frac2[:] = 3.376/(frac1[:] + 55.)
    ept = pl.zeros_like(frac2)
    ept[:] = theta[:]*pl.exp((frac2[:]-0.00254)*(10**3)*w[:]*(1+0.81*w[:]))
    
    return ept

def PartialMixing(equiv_traj,q_traj):
    """Function to determine if and where a trajectory experiences partial mixing.
    
    Args:
        equiv_traj (array): equivalent potential temperature (K) of air parcel
                            at all points along a trajectory
        q_traj (array): specific humidity (kg/kg) of air parcel at all points
                        along a trajectory
    
    Returns:
        maxloc (int or None): location of partial mixing
        PM (bool): variable stating whether or nor partial mixing occurs
    """
    PM = None # start with no partial mixing
    # set up empty arrays for change in theta e & q:
    dte = pl.zeros([len(equiv_traj)-1]); dq = pl.zeros([len(q_traj)-1])
    # calculate changes in theta_e & q backwards along trajectory:
    for step in range(len(dte)):
        dte[step] = equiv_traj[step] - equiv_traj[step+1]
        dq[step] = q_traj[step] - q_traj[step+1]
    
    # calculate all 3-timestep consecutive increases in theta_e:
    incloc = pl.where(dte>0)
    for i in range(len(incloc[0])-3):
        if incloc[0][i+1] == incloc[0][i]+1 and incloc[0][i+2] == incloc[0][i]+2:
            # total increase in theta_e must exceed 2.5 K
            if pl.sum(dte[incloc[0][i]:incloc[0][i]+3]) > 2.5:
                # dte_tot = pl.sum(dte[H[0][i]:H[0][i]+3])
                # find indices of consecutive increases in theta_e
                ind1 = incloc[0][i]; ind2 = incloc[0][i+1]; ind3 = incloc[0][i+2]
                break # no need to continue as most recent PM is required
        else: # if no consecutive increases in theta_e exist then no PM occurs
            PM = False
    
    if 'ind1' in locals(): # has ind1 been defined?
        #max_dte = pl.sum(dte[ind1:ind3+1])
        F = pl.where(dte==dte[ind1:ind3+1].max()) # where is the max increase?
        maxloc = F[0][0]
        if dq[maxloc] > 0:
            PM = True # if dq>0 with max increase in theta_e occurs then PM occurs
        else: # if dq<0 with max increase in theta_e then no PM occurs
            maxloc = None
            PM = False
    else: # if ind1 not defined then no PM occurs
        maxloc = None
        PM = False
    
    return maxloc, PM

def PMcalcs(equiv,q,step,lat,lon):
    """Function to calculate the number and percentage of trajectories from one
    release which have experienced partial mixing.
    
    Args:
        equiv (array): equivalent potential temperature (K) of air parcel at all 
                        points along all trajectories
        q (array): specific humidity (kg/kg) of air parcel at all points along 
                    all trajectories
        step (array): time steps (hours) of all trajectories
        lat (array): latitudes (degrees) along trajectory at all time steps
        lon (array): longitudes (degrees) along trajectory at all time steps
    
    Returns:
        pmcount (float): number of trajectories from release which have experienced 
                         partial mixing
        pmprop (float): proportion of trajectories which have experienced partial
                         mixing, expressed as a percentage
        pmorigin (array): array containing trajectory label, timestep of PM
                           origin, latitude/longitude of origin
    """
    global PartialMixing
    NPART = equiv.shape[0] # number of trajectories
    pmcount = 0.; pmorigin = [] # initialize pmcount as zero & make empty list
    
    for traj in range(NPART):
        # find location of PM and if it occurs:
        a,b = PartialMixing(equiv[traj],q[traj])
        if b == True: 
            pmcount = pmcount + 1. # if PM occurs, add 1 to pmcount
            # then add traj no., timestep & lat/lon to pmorigin list
            if a < 2: # 3 consecutive increases needed, so change a < 2 to 2:
                pmorigin.append((traj,step[traj,2],lat[traj,a],lon[traj,a]))
            else:
                pmorigin.append((traj,step[traj,a],lat[traj,a],lon[traj,a]))
        else: # if no PM occurs then pmcount doesn't change
            pmcount = pmcount
    
    pmprop = (pmcount/NPART)*100 # calculate % of trajectories which undergo PM
    pmorigin = pl.asarray(pmorigin) # convert pmorigin list to array
    
    return pmcount, pmprop, pmorigin

def LastOrigin(originlabs,stratloc):
    """
    """
    if originlabs.shape[0] != stratloc.shape[0]:
        raise Exception('lengths do not match in LastOrigin!!!!')
    
    NPART = originlabs.shape[0]    
    a = pl.isnan(stratloc)
    ol2 = pl.zeros_like(originlabs)
    
    for traj in range(NPART):
        if a[traj] == True:
            ol2[traj,0] = originlabs[traj,0]
            ol2[traj,1:] = originlabs[traj,1:]
        elif originlabs[traj,1] < stratloc[traj]:
            ol2[traj,0] = originlabs[traj,0]
            ol2[traj,1:] = originlabs[traj,1:]
        elif originlabs[traj,1] > stratloc[traj]:
            ol2[traj,0] = 3.
            ol2[traj,1] = stratloc[traj]
            ol2[traj,2:] = pl.float64('nan')
    
    return ol2

def OriginCheck(blorigin,pmorigin,NPART):
    """Function to determine whether a trajectory has boundary layer origin,
    partial mixing origin or unknown origin. Boundary layer origin takes
    precedence over partial mixing origin. If a trajectory makes contact with 
    the boundary layer at any point along trajectory then it assigned label '1'.
    If the trajectory has not made contact with the boundary layer then check
    for partial mixing and assign label '2' if this occurs. If neither origin
    type occurs then origin is unknown so assign label '0'.
    
    Args:
        blorigin (array): trajectory number, timestep of boundary layer origin,
                          lat/lon of boundary layer origin
        pmorigin (array): trajectory number, timestep of partial mixing origin,
                          lat/lon of partial mixing origin
        NPART (float): number of trajectories
    
    Returns:
        originlabs (array): origin labels of trajectories; timesteps of origins;
                            lat/lon of origins
    """
    originlabs = pl.zeros([NPART,4]) # empty array for origin labels, timestep
                                     # of origin & lat/lon of origin
    
    for traj in range(NPART): # check all trajectories
        if blorigin.size == 0. and pmorigin.size == 0.:
            pass
        elif traj in blorigin[:,0]: # if trajectory has bl origin,
            originlabs[traj,0] = 1.   # assign label '1'
            x = pl.where(blorigin[:,0]==traj)
            originlabs[traj,1] = blorigin[x[0][0],1]/6.
            originlabs[traj,2:] = blorigin[x[0][0],2:]
        elif pmorigin.size == 0: # check if there are any pm origins
            originlabs[traj,0] = 0.   # assign label '0'
            originlabs[traj,1:] = None
        # if trajectory has no bl origin, but pm origin,
        elif traj in pmorigin[:,0] and pmorigin.size > 0.:
            originlabs[traj,0] = 2.       # assign label '2'
            x = pl.where(pmorigin[:,0]==traj)
            originlabs[traj,1] = pmorigin[x[0][0],1]/6.
            originlabs[traj,2:] = pmorigin[x[0][0],2:]
        else: # or if no origin can be found, origin is unknown.
            originlabs[traj,0] = 0.   # assign label '0'
            originlabs[traj,1:] = None
    
    return originlabs

def OriginsDayCount(originlabs,NPART,NSTEPS):
    """Function to calculate the total number of trajectories with an origin 
    from either boundary layer or partial mixing at each timestep where
    boundary layer takes precedence over partial mixing.
    
    Args:
        originlabs (array): label denoting if trajectory has bl origin ('1'),
                            partial mixing origin ('2') or no origin ('0');
                            timestep of origin (int or None); lat/lon of origin
        NPART (int): number of trajectories
        NSTEPS (int): number of timesteps
    
    Returns:
        orprop (float): proportion of trajectories with either boundary layer
                        or partial mixing origin, expressed as a percentage
    """
    # set up empty array for timestep of trajectory origin:
    orloc = pl.zeros([NPART,NSTEPS])
    
    for trajno in range(NPART):
        # does trajectory have either bl or pm origin? 1 = bl, 2 = pm
        if originlabs[trajno,0] == 1. or originlabs[trajno,0] == 2.:
            # set orloc to 1 at timestep of origin
            orloc[trajno,originlabs[trajno,1]] = 1.
    
    # number of new trajectories with origin at each timestep:
    newor = pl.sum(orloc,axis=0)
    # cumulative number of trajectories with origin at each timestep: 
    ordaycount = pl.zeros_like(newor)
    for step in range(len(newor)):
        ordaycount[step] = pl.sum(newor[:step+1])
    
    # calculate proportions of trajectories with origin at each timestep:
    orprop = (ordaycount/NPART)*100
    
    return orprop

def StratStep(theta,PV,height):
    """Function to find the timestep at which a trajectory was last in the
    stratosphere. Checks for potential temperature at least 380K, potential
    vorticity at least 3PVU and height greater than 6000m to avoid high
    boundary layer PV.
    
    Args:
        theta (array): potential temperature (K) of air mass at all time
                       steps along all trajectories
        PV (array): potential vorticity (PVU) of air mass at all time steps
                    along all trajectories
        height (array): elevation (metres above surface) of air parcel at all
                        time steps along all trajectories
    
    Returns:
        stratloc (array): timesteps at which each trajectory was most recently 
                          in the stratosphere
    """
    NPART = theta.shape[0] # number of trajectories
    NSTEPS = theta.shape[1] # number of time steps
    stratloc = pl.zeros([NPART]) # empty array for most recent step in stratosphere
    for traj in range(NPART): # loop over trajectories
        for step in range(NSTEPS): # loop over steps
            # check theta and PV conditions
            if theta[traj,step] >= 380. or PV[traj,step] >= 3.:
                # check height condition
                if height[traj,step] > 6000.:
                    # at first step where condition is met, add step to stratloc
                    # use break to move to next trajectory once condition is met
                    stratloc[traj] = step; break
                else:
                    stratloc[traj] = None
            else:
                stratloc[traj] = None
    
    return stratloc

def StratOrigin(stratloc,lat,lon):
    """Function to make array of trajectories which have made contact with the
    stratosphere containing trajectory number, timestep of last contact with
    stratosphere and latidude/longitude of last contact with stratosphere.
    
    Args:
       stratloc (array): timesteps at which each trajectory was most recently 
                         in the stratosphere; 'nan' shows that trajectory did
                         not make contact with stratosphere
       lat (array): latitudes (degrees) of all trajectories at all timesteps
       lon (array): longitudes (degrees) of all trajectories at all timesteps
    
    Returns:
        stratlist (array): trajectory numbers of trajectories which have
                           stratospheric origin, with timesteps and lat/lon
    """
    check = pl.isnan(stratloc) #  find the nan entries
    strators = pl.where(check==False) # the non-nan entries are important here
    # make empty array, no. of trajectories in strat X 4:
    stratlist = pl.zeros([strators[0].size,4])
    
    for trajno in range(strators[0].size): # loop over trajectories in strat.
        stratlist[trajno,0] = strators[0][trajno] # trajectory number
        stratlist[trajno,1] = stratloc[strators[0][trajno]] # timestep
        stratlist[trajno,2] = lat[strators[0][trajno],stratloc[strators[0][trajno]]]
        stratlist[trajno,3] = lon[strators[0][trajno],stratloc[strators[0][trajno]]]
    
    return stratlist

def StratVar(strator,var):
    """Function to eliminate the values of an attribute along trajectories which
    do not have stratospheric origin i.e. no boundary layer or partial mixing 
    origin, trajectory in stratosphere at final timestep.
    
    Args:
        strator (array): trajectory numbers of trajectories which have
                         stratospheric origin
        var (array): attribute along trajectory e.g. specific humidity or PV
    
    Returns
        var_strat (array): same shape as var, but trajectories which have no
                            stratospheric origin are eliminated
    """
    var_strat = pl.zeros_like(var) # empty array for variable
    
    for trajno in range(len(strator)): # loop over traj numbers with strat origin
        # extract data for trajectories with stratospheric origin
        var_strat[strator[trajno,0]] = var[strator[trajno,0]]
    
    return var_strat

def AttrPlot(equiv,q,height,originlab):
    """Function to plot equivalent potential temperature, specific humidity and
    pressure along a trajectory with the timestep of origin also shown. Trajectories
    are coloured by their origin type and the matplotlib.lines.Line2D object is
    returned along with the origin label so the figure legend can be created
    in the main code.
    
    Args:
        equiv (array): equivalent potential temperature (K) at all timesteps
                       along a trajectory
        q (array): specific humidity (kg/kg) at all timesteps along a trajectory
        pres (array): pressure (hPa) at all timesteps along a trajectory
        originlab (array): label denoting if trajectory has bl origin ('1'),
                           partial mixing origin ('2') or no origin ('0');
                           timestep of origin (int or None); lat/lon of origin
   
   Returns:
       a[0] (object): matplotlib.lines.Line2D object to be used to make figure
                      legend
       originlab[0] (float): label signifying origin type of trajectory; '1'
                             refers to boundary layer origin, '2' refers to
                             partial mixing origin and '0' refers to unknown
                             origin
    """
    NSTEPS = equiv.size; step = pl.linspace(0,NSTEPS-1,NSTEPS)
    # change values of q less than 10e-10 to 10e-10:
    r = pl.where(q<10e-10)
    for i in range(r[0].size):
        q[r[0][i],r[1][i]] = 10e-10
    
    # if origin label is 1 or 2 (bl or pm) then use timestep as loc
    if originlab[1] > 0.: loc = originlab[1]
    else: loc = None # if origin label is 0 (unknown) there is no loc
    
    if originlab[0] == 1.: clr = 'b' # BL origin trajectories are blue
    elif originlab[0] == 2.: clr = 'r' # PM origin trajectories are red
    elif originlab[0] == 0: clr = 'w' # unknown origin trajectories are white
    
    # plot equivalent potential temperature:
    ax1 = pl.subplot(3,1,1)
    a = pl.plot(step/4,equiv,color=clr)#; pl.legend(loc=0,ncol=3)
    # if trajectory has an origin then plot as black triangle
    if loc != None:
        pl.plot(step[loc]/4,equiv[loc],marker='^',color='k',ms=7.5,zorder=5.)
    pl.ylabel('$\\theta_e$ (K)',fontsize=18)
    ax1.xaxis.set_major_formatter(pl.NullFormatter())
    
    # plot specific humidity:
    ax2 = pl.subplot(3,1,2)
    pl.plot(step/4,q*(10**3),color=clr)
    pl.ylabel('specific humidty (g/kg)',fontsize=18)
    if loc != None:
        pl.plot(step[loc]/4,q[loc]*(10**3),marker='^',color='k',ms=7.5,zorder=5.)
    ax2.xaxis.set_major_formatter(pl.NullFormatter())
    
    # plot pressure:
    ax3 = pl.subplot(3,1,3)
    pl.plot(step/4,height/1000.,color=clr)
    pl.ylabel('height (km)',fontsize=18); #pl.ylim(1050,100); 
    if loc != None:
        pl.plot(step[loc]/4,height[loc]/1000.,marker='^',color='k',ms=7.5,zorder=5.)
    ax3.set_xlabel('Time along trajectory (days)',fontsize=18)
    
    pl.subplots_adjust(top=0.95,bottom=0.06,wspace=0.20,hspace=0.08)
    
    return a[0], originlab[0]

def ZPDF(originlabs,height,ortype):
    """Function to produce histograms of trajectory heights at CAT I & CAT II
    origins.
    
    Args:
        originlabs (array): label denoting if trajectory has bl origin ('1'),
                           partial mixing origin ('2') or no origin ('0');
                           timestep of origin (int or None); lat/lon of origin
       height (array): heights above surface of trajectories at all timesteps
       ortype (int): label denoting if trajectory has bl origin ('1'),
                     partial mixing origin ('2') or no origin ('0')
    """
    # colours used in histogram and index of required colour:
    clrs = ['b','r']; ind = ortype-1
    if ortype == 1. or ortype == 2.: # chech if BL or PM origin
        x = pl.where(originlabs[:,0]==ortype) # which trajectories have origin?
        H = pl.zeros(x[0].size) # empty array, size no of trajectories with ortype
    
        for trajno in range(x[0].size): # loop over trajectories with this ortype
            tn = x[0][trajno] # trajectory number, trajno is array index
            step = originlabs[x[0][trajno],1] # timestep of origin
            H[trajno] = height[tn,step] # extract height of trajectory at origin
        
        pl.hist(H,color=clrs[ind])# make histogram of heights at origin
        pl.xlabel('heights (m)',fontsize=16); pl.ylabel('counts',fontsize=16)
        #pl.annotate(labs[ind],xy=(0.9,0.9),xycoords="figure points")
        pl.subplots_adjust(hspace=0.34,top=0.94)
    else: # if ortype entered wrong, raise exception
        raise Exception('Wrong origin type. Must be 1 or 2!')

def PrecPDF(originlabs,precip,ortype,bins):
    """Function to produce histograms of precipitation at origin location
    
    Args:
        originlabs (array): origin type and location of each trajectory
        precip (array): convective or stratiform precip along each trajectory
        ortype (float): 1 or 2, denoting CAT I or CAT II origin
        bins (int): number of bins for histogram
    """
    precip = precip#*1000
    clrs = ['b','r']; ind = ortype-1
    if ortype == 1. or ortype == 2.: # chech if BL or PM origin
        x = pl.where(originlabs[:,0]==ortype) # which trajectories have origin?
        H = pl.zeros(x[0].size) # empty array, size no of trajectories with ortype
    
        for trajno in range(x[0].size): # loop over trajectories with this ortype
            tn = x[0][trajno] # trajectory number, trajno is array index
            step = originlabs[x[0][trajno],1] # timestep of origin
            H[trajno] = precip[tn,step] # extract height of trajectory at origin
        
        pl.hist(H*1000,bins=bins,color=clrs[ind],normed=True,log=False)# make histogram of heights at origin
        pl.xlim(xmin=0); pl.ylim(0,1)
        #pl.xlabel('precipitation (mm)',fontsize=16)#; pl.ylabel('counts',fontsize=16)
        #pl.annotate(labs[ind],xy=(0.9,0.9),xycoords="figure points")
        pl.subplots_adjust(hspace=0.34,top=0.94)
    else: # if ortype entered wrong, raise exception
        raise Exception('Wrong origin type. Must be 1 or 2!')
    

def OrOrder(blz,height,stepno,lat,lon,equiv,q,stratloc,lsm,eralon,eralat):
    """Function to determine if a trajectory has boundary layer, partial mixing
    or stratospheric origin where stratospheric origin is prioritized over
    boundary layer and partial mixing origins i.e. a trajectory has stratospheric
    origin if it meets the criteria for stratospheric location most recently along
    the back trajectory. Boundary layer origin currently supersedes partial
    mixing but the function takes into account that a trajectory can fulfill
    both criteria. If a trajectory experiences partial mixing more recently than
    making contact with the boundary layer then it is initially labelled as having
    origin from partial mixing, then resassigned to boundary layer origin later.
    
    Args:
        blz (array): boundary layer heights (metres above surface) of all 
                     trajectories at all time steps
        height (array): heights (metres above surface) of trajectories at all
                        timesteps
        stepno (array): number of hours elapsed at each timestep for all
                        trajectories
        lat (array): latitudes (degrees) along all trajectories at all timesteps
        lon (array): longitudes (degrees) along all trajectories at all timesteps
        equiv (array): equivalent potential temperature (K) along all trajectories
                       at all timesteps
        q (array): specific humidity (kg/kg) along all trajectories at all timesteps
        stratloc (array): timesteps of each trajectory's last contact with the
                         stratosphere; trajectories which have not made contact
                         are labelled 'nan'
       
    Returns:
        blprop (array): proportions of trajectories with boundary layer origin
                        at each timestep expressed as a percentage
        pmprop (array): proportions of trajectories with partial mixing origin
                        at each timestep expressed as a percentage
        strprop (array): proportions of trajectories with stratospheric origin
                         at each timestep
    """
    global BoundaryLayer2; global PMcalcs
    # Take NPART & NSTEPS from dimensions of an attribute array:
    NPART = blz.shape[0]; NSTEPS = blz.shape[1]
    blors = pl.zeros([NPART,NSTEPS]) # empty array for origin labels
    pmors = pl.zeros_like(blors); pmdis = pl.zeros_like(blors) # discounted pmorigins
    strators = pl.zeros_like(blors); pmout = pl.zeros_like(blors) # reassigned pmorigins
    
    # Use BoundaryLayer & PMcalcs functions to find all trajectories that have
    # boundary layer & partial mixing origins & their origin timesteps:
    blc, blp, blorigin = BoundaryLayer2(blz,height,stepno,lat,lon,lsm,eralon,eralat)
    pmc, pmp, pmorigin = PMcalcs(equiv,q,stepno,lat,lon)
    
    for trajno in range(blorigin.shape[0]): # loop over BL origin trajectories
        # for all trajectories with BL origin, blors[traj,step] = 1
        blors[blorigin[trajno,0],blorigin[trajno,1]/6.] = 1.
    
    for trajno in range(pmorigin.shape[0]): # loop over PM origin trajectories
        # for all trajectories with PM origin, pmors[traj,step] = 1
        pmors[pmorigin[trajno,0],pmorigin[trajno,1]/6.] = 1.
    
    for traj in range(NPART): # loop over all trajectories
        a = pl.isnan(stratloc[traj]) # check if traj has strat origin
        if a == False: # if traj has strar origin ....
            # ... set strators[traj,step] = 1
            strators[traj,stratloc[traj]] = 1.

    # prioritize origins by stratosphere, BL, PM along back trajectory
    # make array for origin type labels; '1' for BL, '1' for PM & '-1' for strat
    # start with this array equal to blors array
    ORS = blors
    for traj in range(NPART): # loop over all trajectories
        # find where BL & PM origins occur along trajectory:
        x = pl.where(ORS[traj]==1.); y = pl.where(pmors[traj]==1.)
        if x[0].size == 1. and y[0].size == 1.: # if both BL & PM occur...
            ORS[traj] = ORS[traj] # ... prioritize BL
            if y[0][0] < x[0][0]: # ... but if PM is more recent than BL...
                pmdis[traj,y[0][0]] = 1. # .. label discounted PM origin
                pmout[traj,x[0][0]] = 1. # ... label where PM is reassigned 
        elif x[0].size == 1. and y[0].size == 0.: # if only BL origin occurs ...
             ORS[traj] = ORS[traj] # ... do nothing
        elif x[0].size == 0. and y[0].size == 1.: # if only PM occurs ...
            ORS[traj,y[0][0]] = 2. # ... label '2' for PM origin
    
    for traj in range(NPART): # loop over all trajectories
        # find where BL & strat origins occur along trajectory
        x = pl.where(ORS[traj]==1.); y = pl.where(strators[traj]==1.)
        if x[0].size == 1. and y[0].size == 1.: # if both BL & strat occur ...
            if y[0][0] < x[0][0]: # ... if strat is more recent...
                # ... remove BL label & add strat label at strat timestep
                ORS[traj,x[0][0]] = 0.; ORS[traj,y[0][0]] = -1.
            elif x[0][0] < y[0][0]: # ... if BL is more recent ...
                ORS[traj] = ORS[traj] # ... do nothing
        elif x[0].size == 0. and y[0].size == 1.: # if only strat occurs ...
            ORS[traj,y[0][0]] = -1. # .. label '-1' for strat origin
        elif x[0].size == 1. and y[0].size == 0.: # if only BL occurs ...
            ORS[traj] = ORS[traj] # ... do nothing
    
    for traj in range(NPART):
        # find where PM origins and strat origins occur along trajectory:
        x = pl.where(ORS[traj]==2.); y = pl.where(strators[traj]==1.)
        if x[0].size == 1. and y[0].size == 1.: # if both PM & strat occur ...
            if y[0][0] < x[0][0]: # ... if strat is more recent ...
                # ... remove PM label & add strat label at strat timestep
                ORS[traj,x[0][0]] = 0.; ORS[traj,y[0][0]] = -1.
            elif x[0][0] < y[0][0]: # ... if PM is more recent ...
                ORS[traj] = ORS[traj] # ... do nothing
        elif x[0].size == 0. and y[0].size == 1.: # if only strat occurs ...
            ORS[traj,y[0][0]] = -1. # .. label '-1' for strat origin
        elif x[0].size == 1. and y[0].size == 0.: # if only PM occurs ...
            ORS[traj] = ORS[traj] # ... do nothing
    
    for traj in range(NPART):
        nonzero = pl.where(ORS[traj]!=0.)
        if nonzero[0].size > 1.:
            ORS[traj,nonzero[0][1]] = 0.

    newlabs = pl.zeros([NPART,4])
    for traj in range(NPART):
        nonzero = pl.where(ORS[traj]!=0.)
        if nonzero[0].size > 0:
            newlabs[traj,0] = ORS[traj,nonzero[0][0]]
            newlabs[traj,1] = nonzero[0][0]
            newlabs[traj,2] = lat[traj,nonzero[0][0]]
            newlabs[traj,3] = lon[traj,nonzero[0][0]]
        else:
            newlabs[traj,1:] = None

    #  empty arrays for each origin type:
    bllocs = pl.zeros_like(ORS); pmlocs = pl.zeros_like(ORS)
    strlocs = pl.zeros_like(ORS)
    
    for traj in range(NPART): # loop over all trajectories
        for step in range(NSTEPS): # loop over all timesteps
            if ORS[traj,step] == 1.: # if BL origin ...
                bllocs[traj,step] = 1.      # .. assign to bllocs
            elif ORS[traj,step] == 2.: # if PM origin ...
                pmlocs[traj,step] = 1.      # ... assign to pmlocs
            elif ORS[traj,step] == -1: # if strat origin ...
                strlocs[traj,step] = 1.     # ... assign to strlocs

    # take sum of location arrays along all trajectories:
    blcount = pl.sum(bllocs,axis=0); pmcount = pl.sum(pmlocs,axis=0)
    stratcount = pl.sum(strlocs,axis=0)
    # do the same for discounted PM origins and reassignments arrays:
    discount = pl.sum(pmdis,axis=0); outcount = pl.sum(pmout,axis=0)
    # set up empty arrays for cumulative sums:
    blsums = pl.zeros([NSTEPS]); pmsums = pl.zeros([NSTEPS])
    strsums = pl.zeros([NSTEPS])
    
    # since PM origins are reassigned, the total number of trajectories with
    # PM origin changes when trajectory is reassigned to BL origin, so
    # add discounted origins & subtract reassigned origins
    pmcount = pmcount + discount - outcount
    
    for step in range(NSTEPS): # loop over all timesteps
        # find totals at each timestep for each origin type:
        blsums[step] = pl.sum(blcount[:step+1])
        strsums[step] = pl.sum(stratcount[:step+1])
        pmsums[step] = pl.sum(pmcount[:step+1])
        
    # convert totals to percentages:
    blprop = (blsums/NPART)*100; pmprop = (pmsums/NPART)*100
    strprop = (strsums/NPART)*100

    return blprop, pmprop, strprop, newlabs

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
    
    if relpt[-1] > lat[0]:
        latx = pl.zeros([lat.size+1]); latx[1:] = lat; latx[0] = 90.0
        lat = latx.copy()
        b = NearestIndex(lat,relpt[1])
        flux2 = pl.zeros([flux.shape[0]+1,flux.shape[1]])
        flux2[1:,:] = flux; flux2[0] = flux[0]
        flux = flux2
    elif relpt[-1] < lat[-1]:
        latx = pl.zeros([lat.size+1]); latx[:-1] = lat; latx[-1] = -90.0
        lat = latx.copy()
        b = NearestIndex(lat,relpt[1])
        flux2 = pl.zeros([flux.shape[0]+1,flux.shape[1]])
        flux2[:-1,:] = flux; flux2[-1] = flux[-1]
        flux = flux2
    else:
        b = NearestIndex(lat,relpt[1]) # nearest latitude index
    
    if relpt[0] == lon[a] and relpt[1] == lat[b]:
        return flux[b,a]
    elif relpt[0] == lon[a]:
        F = interp1D(relpt,lon,lat,flux)
        return F
    elif relpt[1] == lat[b]:
        F = interp1D(relpt,lon,lat,flux)
        return F
    
    if lon[a] < relpt[0]: # nearest lon west of relpt
        p1[0] = lon[a]; p3[0] = lon[a];  p2[0] = lon[a+1]; p4[0] = lon[a+1]
    elif lon[a] > relpt[0]: # nearest lon east of relpt
        p2[0] = lon[a]; p4[0] = lon[a]; p1[0] = lon[a-1]; p3[0] = lon[a-1]
        
    # does not take 0 meridian into account yet

    
    if lat[b] < relpt[1]: # nearest lat south of relpt
        p1[1] = lat[b]; p2[1] = lat[b]; p3[1] = lat[b-1]; p4[1] = lat[b-1]
    elif lat[b] > relpt[1]: # nearest lat north of relpt
        p3[1] = lat[b]; p4[1] = lat[b]; p1[1] = lat[b+1]; p2[1] = lat[b+1]
    #elif lat[b] == relpt[1]: # lat equal to relpt
    #    p3[1] = lat[b]; p4[1] = lat[b]; p1[1] = lat[b]; p2[1] = lat[b]
    
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

def interp1D(relpt,lon,lat,flux):
    """
    """
    # if release point longitude co-ordinate greater than max. ERA-Int longitude
    if relpt[0] > lon[-1]:
        a = -1 # take max ERA-Interim longitude as nearest longitude index
    else:
        a = NearestIndex(lon,relpt[0]) # nearest longitude index
    b = NearestIndex(lat,relpt[1]) # nearest latitude index
    
    # two IF things for is the lat or lon co-ordinate the problem
    if relpt[0] == lon[a]:
        if relpt[1] > lat[b]:
            X = (lat[b],lat[b-1]); Y = (flux[b,a],flux[b-1,a])
        elif relpt[1] < lat[b]:
            X = (lat[b+1],lat[b]); Y = (flux[b+1,a],flux[b,a])
        f = interp1d(X,Y); F = f(relpt[1])
    elif relpt[1] == lat[b]:
        if relpt[0] > lon[a]:
            if lon[a] == lon[-1]:
                X = (lon[a],lon[a+1]+360); Y = (flux[b,a],flux[b,0])
            else:
                X = (lon[a],lon[a+1]); Y = (flux[b,a],flux[b,a+1])
        elif relpt[0] < lon[a]:
            X = (lon[a-1],lon[a]); Y = (flux[b,a-1],flux[b,a])
        f = interp1d(X,Y); F = f(relpt[0])
    
    return F
    
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

def VarExtract(orlabs,var,label):
    """Function to extract a variable along a trajectory which has some origin
    type.
    
    Args:
        orlabs (array): origin labels of trajectories
        var (array): variable to be extracted
        label (float): 1 or 2, referring to to CAT I or CAT II origin
    
    Returns:
        var_or (array): variable along only trajectories with origin
    """
    var_or = pl.zeros_like(var) # empty array for extracted variable
    NPART = orlabs.shape[0] # number of trajectories
    for traj in range(NPART): # loop over number of trajectories
        if orlabs[traj,0] == label: # check if trajectory has correct origin type
            var_or[traj] = var[traj] # extract data for trajectory
    
    return var_or

def FluxEvol(blz,height,stepno,lat,lon,equiv,q,theta,PV,rlslabs,u0,v0,normals,
            pres,etalevs,mag,sp_tj,lsm,eralon,eralat):
    """Function to calculate how the magnitude of the cross-watershed moisture
    flux evolves back in time as more and more trajectories are assigned an
    origin type.
    
    Args:
        blz (array): boundary layer height (m) along all trajectories
        height (array): height above surface (m) of all trajectories at all timesteps
        stepno (array): timesteps in hours along back trajectory
        lat (array): latitude of all trajectories at all timesteps
        lon (array): longitude of all trajectories at all timesteps
        equiv (array): equivalent potential temperature (K) along all trajectories
        q (array): specific humidity (kg/kg) along all trajectories
        theta (array): dry potential temperature
        PV (array): potential vorticity (1 PVU = 1e-6*m^2Ks^1kg^-1) along all
                    trajectories
        rlslabs (array): longitude, latitude co-ordinates of trajectory release
                          points and segment label
        u0 (array): zonal velocity (m/s) of trajectory at t=0
        v0 (array): meridional velocity (m/s) of trajectory at t=0
        normals (array): normal vectors for each segment along the watershed
        pres (array): pressure (hPa) along all trajectories
        etalevs (array): model levels from trajectories were released
        mag (float): magnitude of the total cross-watershed moisture flux
        sp_tj (array): surface pressure (hPa) a t=0 for all trajectories
    
    Returns:
        flux_prop (array): proportion of the cross-watershed moisture flux
                            explained by all trajectories with origin at each
                            timestep
        bl_prop (array): proportion of the cross-watershed moisture flux explained
                         by trajectories with CAT I (boundary layer) origin at 
                         each timestep
        pm_prop (array): proportion of the cross-watershed moisture flux explained
                         by trajectories with CAT II (partial mixing) origin at
                         each timestep
        str_prop (array): proportion of the cross-watershed moisture flux 
                          explained by trajectories with stratospheric origin
                          at each timestep             
    """
    # get no. of trajectories & timesteps from shape of blz array:
    NPART = blz.shape[0]; NSTEPS = blz.shape[1]
    NREL = rlslabs.shape[0] # get no. of release points from shape of rlslabs array
    NCLUST = int(NPART/NREL) # no. of levels is no. trajectories / no. points
    # make empty arrays for flux proportions, size no. of steps:
    flux_prop = pl.zeros([NSTEPS]); bl_prop = pl.zeros([NSTEPS])
    pm_prop = pl.zeros([NSTEPS]); str_prop = pl.zeros([NSTEPS])
    
    # take absolute values of u, v and normals:
    u0 = pl.absolute(u0); v0 = pl.absolute(v0)
    normals = pl.absolute(normals)
    
    for step in range(NSTEPS): # loop over timesteps
        # find where all trajectories have last been in the stratosphere:
        sl = StratStep(theta[:,:step+1],PV[:,:step+1],height[:,:step+1])
        # get the flux proportion of each origin type at each timestep
        # and the array with the origin label, timestep and co-ordinates
        blprop, pmprop, strprop, originsteps = OrOrder(blz[:,:step+1],height[:,:step+1],
                                                   stepno[:,:step+1],lat[:,:step+1],
                                                    lon[:,:step+1],equiv[:,:step+1],
                                                    q[:,:step+1],sl,lsm,eralon,eralat)
    
        q_or = OriginFlux(originsteps,q) # extract q for all traj with origin
        # extract q for all traj with CAT I, CAT II or strat origin
        q_bl = VarExtract(originsteps,q,1); q_pm = VarExtract(originsteps,q,2)
        q_st = VarExtract(originsteps,q,-1)
        
        # Calculate vertically integrated moisture fluxes at each release point
        # for all trajectories with origin, CAT I, CAT II & strat origins
        stepflx = FluxCalc(rlslabs,u0,v0,sp_tj,etalevs,normals,q_or[:,:],
                                                 pres) # ALL
        blflx = FluxCalc(rlslabs,u0,v0,sp_tj,etalevs,normals,q_bl[:,:],
                                                 pres) # CAT I
        pmflx = FluxCalc(rlslabs,u0,v0,sp_tj,etalevs,normals,q_pm[:,:],
                                                 pres) # CAT II
        strflx = FluxCalc(rlslabs,u0,v0,sp_tj,etalevs,normals,q_st[:,:],
                                                  pres) # strat
       
        # Calculate the proportion of the magnitude of the flux explained by
        # all trajectories with origin, CAT I, CAT II and strat origins
        flux_prop[step] = ((pl.sum(pl.absolute(stepflx)))/(10**9))/mag # ALL
        bl_prop[step] = ((pl.sum(pl.absolute(blflx)))/(10**9))/mag # CAT I
        pm_prop[step] = ((pl.sum(pl.absolute(pmflx)))/(10**9))/mag # CAT II
        str_prop[step] = ((pl.sum(pl.absolute(strflx)))/(10**9))/mag # strat
    
    # turn these into percentages
    flux_prop = flux_prop*100; bl_prop = bl_prop*100
    pm_prop = pm_prop*100; str_prop = str_prop*100
   
    return flux_prop, bl_prop, pm_prop, str_prop

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

def CatchLine(NS_w,EW_n,NS_e,EW_s):
    """Function to create the catchment boundary of an ocean drainage basin
    using four lines dividing catchments.
                ONLY WORKS FOR ATLANTIC, INDIAN AND PACIFIC
                ARCTIC & SOUTHERN CATCHMENTS ARE BOUNDED BY
                A SINGLE LINE, THIS FUNCTION REQUIRES FOUR
                LINES TO CREATE A BOUNDARY!!!!
    
    Args:
        NS_w (array): line oriented north-south forming western boundary of catchment
        EW_n (array): line oriented east-west forming northern boundary of catchment
        NS_e (array): line oriented north-south forming eastern boundary of catchment
        EW_s (array): line oriented east-west forming southern boundary of catchment

    Returns:
        boundary (array): catchment boundary of ocean drainage basin, size sum
                          of the four component lines minus 4
    """
    # get the length of the zero-axis of each line:
    a = NS_w.shape[0]; b = EW_n.shape[0]; c = NS_e.shape[0]; d = EW_s.shape[0]
    # length of boundary is sum of the 4 lengths minus 4 as end points are repeated
    boundary = pl.zeros([a+b+c+d-4,2])
    
    # MAKE THE BOUNDARY:
    boundary[:b] = EW_n[:] # northern boundary
    boundary[b:c+b-1] = NS_e[1:] # eastern boundary, remove first point
    # southern boundary, reverse order and remove first point
    boundary[c+b-1:c+b+d-2] = pl.flipud(EW_s[:-1])
    # western boundary, reverse order, remove first and last points.
    boundary[c+b+d-2:] = pl.flipud(NS_w[1:-1])
    
    return boundary

def CatchPoly(boundary,m):
    """Function to make a polygon using a catchment boundary and map.
    
    Args:
        boundary (array): longitude, latitude co-ordinates (degrees) of the 
                          catchment boundary of an ocean drainage basin
        m (object): Basemap object defining the extent of a map
    
    Returns:
        poly: polygon defining the drainage basin of an ocean
    """
    global Polygon
    # empty array for map co-ordinates of boundary
    bnd_map = pl.zeros_like(boundary)
    # convert boundary to map co-ordinates
    bnd_map[:,0], bnd_map[:,1] = m(boundary[:,0],boundary[:,1])
    
    l = list(bnd_map) # Polygon class seems to like lists
    poly = Polygon(l,closed=True) # make the polygon of the drainage basin
                                  # closed=True means 1st point is the last one
    
    return poly

def OriginPartition(newlabs,lon_tj,lat_tj):
    """Function to split trajectory origins into ocean catchment areas.
    
    Args:
        newlabs (array): origin type (CAT I, CAT II, strat or none), timestep
                         & location for all trajectories
        lon_tj (array): longitude co-ordinates of all trajectories at all timesteps
        lat_tj (array): latitude co-ordinates of all trajectories at all timesteps
    
    Returns:
        atlcount (array): origin type, timestep & location for trajectories
                          with Atlantic origin
        indcount (array): as above but for the Indian Ocean
        paccount (array): as above but for the Pacific
        arccount (array): as above but for the Arctic
        soucount (array): as above but for the Southern Ocean
    """
    global Basemap; global CatchLine; global Add360; global CatchPoly
    global mplPath
    labels = pl.zeros([newlabs.shape[0]],dtype='S1') # why is this here?
    
    ##### NEED ALL LINES DEFINING SEGMENTS FOR TRAJECTORY RELEASE LINES #######
    sheddir = '/home/np838619/Watershed/shed_defs/' # directory with data
    
    endpts1 = pl.genfromtxt(sheddir+'Am_clicks.txt',skip_header=5) # Americas
    endpts2 = pl.genfromtxt(sheddir+'AfMe_clicks.txt',skip_header=5) # Africa
    endpts3 = pl.genfromtxt(sheddir+'EAA_clicks.txt',skip_header=5) # East Asia
    
    # Arctic lines:
    endpts4 = pl.genfromtxt(sheddir+'ArA_clicks.txt',skip_header=5) # Atlantic
    endpts5 = pl.genfromtxt(sheddir+'ArI_clicks.txt',skip_header=5) # Indian
    endpts6 = pl.genfromtxt(sheddir+'ArP_clicks.txt',skip_header=5) # Pacific
    
    # Southern Ocean lines
    endpts7 = pl.genfromtxt(sheddir+'SOA_clicks.txt',skip_header=5) # Atlantic
    endpts8 = pl.genfromtxt(sheddir+'SOI_clicks.txt',skip_header=5) # Indian
    endpts9 = pl.genfromtxt(sheddir+'SOP_clicks.txt',skip_header=5) # Pacific
    
    # Need traj release points for Arctic & Southern catchments because
    # end points lines means that latitude lines aren't correctly followed on
    # polar projections
    endpts10 = pl.genfromtxt(sheddir+'Ar_traj_release_new.txt',skip_header=5)
    endpts11 = pl.genfromtxt(sheddir+'SO_traj_release_new.txt',skip_header=5)
    
    ###### NEED 1 MAP FOR ATLANTIC/INDIAN CATCHMENTS, ANOTHER FOR PACIFIC #####
    ####### ALSO NEED SEPERATE MAPS FOR ARCTIC AND SOUTHERN CATCHMENTS ########
    m1 = Basemap(projection='cyl',llcrnrlon=-180,llcrnrlat=-80,urcrnrlon=180,urcrnrlat=80)
    m2 = Basemap(projection='cyl',llcrnrlon=0,llcrnrlat=-80,urcrnrlon=360,urcrnrlat=80)
    m3 = Basemap(projection='npaeqd',boundinglat=25,lon_0=270,round=True)
    m4 = Basemap(projection='spstere',boundinglat=-13,lon_0=270,round=True)

    ####### DEFINE POLYGONS FOR ATLANTIC, INDIAN AND PACIFIC CATCHMENT ########
    atl_bnd = CatchLine(endpts1,endpts4,endpts2,endpts7) # Atlantic boundary
    ind_bnd = CatchLine(endpts2,endpts5,endpts3,endpts8) # Indian boundary
    #e6 = endpts6; e1 = endpts1#; e1[:,0] = e1[:,0] + 360.
    # Negative lon co-ords in endpts1 & endpts6, correct this:
    e1 = Add360(endpts1); e6 = Add360(endpts6)
    #for i in range(e6.shape[0]):
    #    if e6[i,0] < 0.:
    #        e6[i,0] = e6[i,0] + 360.
    pac_bnd = CatchLine(endpts3,e6,e1,endpts9) # Pacific boundary
    # empty arrays for Arctic & Southern map co-ords, convert to map co-ords:
    arc_bnd = pl.zeros_like(endpts10); arc_bnd[:,0], arc_bnd[:,1] = m3(endpts10[:,0],endpts10[:,1])
    sou_bnd = pl.zeros_like(endpts11); sou_bnd[:,0], sou_bnd[:,1] = m4(endpts11[:,0],endpts11[:,1])
    
    # Make Polygons for Ocean catchments:
    atl_ply = CatchPoly(atl_bnd,m1); ind_ply = CatchPoly(ind_bnd,m1)
    pac_ply = CatchPoly(pac_bnd,m2); arc_ply = CatchPoly(arc_bnd,m3)
    sou_ply = CatchPoly(sou_bnd,m4)
    
    # MAKE EMPTY LISTS FOR EACH CATCHMENT AREA:
    atlcount = []; indcount = []; paccount = []; arccount = []; soucount = []
    strcount = []; unacount = []
    
    # LOOP OVER NEWLABS ARRAY
    for traj in range(newlabs.shape[0]):
        # IF LABEL = 1 OR 2
        if newlabs[traj,0] == 1. or newlabs[traj,0] == 2.:
            # NNED TO CHECK IF ORIGIN TIMESTEP=0 AS THIS MEANS THE ORIGIN CO-ORDS
            # ARE THE RELEASE POINT AND NOT IN A PARTICULAR CATCHMENT
            if newlabs[traj,1] == 0.:
                if lon_tj[traj,1] > 180.:
                    lonpt = lon_tj[traj,1] - 360.
                else:
                    lonpt = lon_tj[traj,1]
                # DEAL WITH THIS BY CHECKING WHICH CATCHMENT TRAJ WAS IN AT TIMESTEP=1
                # CONVERT ORIGIN CO-ORDINATES TO MAP CO-ORDINATES
                a,b = m1(lonpt,lat_tj[traj,1],inverse=False) # Atlantic/Indian
                c,d = m2(lon_tj[traj,1],lat_tj[traj,1],inverse=False) # Pacific
                e,f = m3(lon_tj[traj,1],lat_tj[traj,1],inverse=False) # Arctic
                g,h = m4(lon_tj[traj,1],lat_tj[traj,1],inverse=False) # Southern
            else: # ORIGIN TIMESTEP NOT EQUAL TO ZERO
                if newlabs[traj,-1] > 180.: # check if lon co-ord > 180
                    lonpt = newlabs[traj,-1] - 360.
                else:
                    lonpt = newlabs[traj,-1]
                a,b = m1(lonpt,newlabs[traj,-2],inverse=False) # Atlantic/Indian
                c,d = m2(newlabs[traj,-1],newlabs[traj,-2],inverse=False) # Pacific
                e,f = m3(newlabs[traj,-1],newlabs[traj,2],inverse=False) # Arctic
                g,h = m4(newlabs[traj,-1],newlabs[traj,2],inverse=False) # Southern
                
            ############# CHECK WHICH CATCHMENT EACH ORIGIN IS IN #############
            # create paths from the catchment boundaries
            atlPath = mplPath.Path(list(atl_bnd)); indPath = mplPath.Path(list(ind_bnd))
            pacPath = mplPath.Path(list(pac_bnd)); arcPath = mplPath.Path(list(arc_bnd))
            souPath = mplPath.Path(list(sou_bnd))
            # WHICH CATCHMENT WAS ORIGIN IN?
            # contains_point returns True if co-ords are inside given path
            X = atlPath.contains_point((a,b)); Y = indPath.contains_point((a,b))
            Z = pacPath.contains_point((c,d))
            W = arcPath.contains_point((e,f))
            V = souPath.contains_point((g,h))
            info = pl.zeros([5]) # empty array to append to catchcount lists
            info[0] = traj; info[1:] = newlabs[traj]
            if X == True: # ATLANTIC
                # IF ORIGIN IN THIS CATCHMENT, APPEND INFO FROM NEWLABS TO LIST
                atlcount.append(info)
            elif Y == True: # INDIAN
                indcount.append(info)
            elif Z == True: # PACIFIC
                paccount.append(info)
            elif W == True: # ARCTIC
                arccount.append(info)
            elif V == True: # SOUTHERN
                soucount.append(info)
        elif newlabs[traj,0] == -1.:
            info = pl.zeros([5])
            info[0] = traj; info[1:] = newlabs[traj]
            strcount.append(info)
        elif newlabs[traj,0] == 0:
            info = pl.zeros([5])
            info[0] = traj; info[1:] = newlabs[traj]
            unacount.append(info)

    # turn lists into arrays to output
    atlcount = pl.asarray(atlcount); indcount = pl.asarray(indcount)
    paccount = pl.asarray(paccount); arccount = pl.asarray(arccount)
    soucount = pl.asarray(soucount); strcount = pl.asarray(strcount)
    unacount = pl.asarray(unacount)
    
    # WHAT % OF TRAJECTORIES WITH ORIGIN IS EXPLAINED BY EACH CATCHMENT?
    return atlcount, indcount, paccount, arccount, soucount, strcount, unacount

def CatchScat(atlcount,indcount,paccount,arccount,soucount):
    """Function to make scatter plots of trajectory origins.
    
    Args:
        atlcount (array): origin type, timestep & location for trajectories
                          with Atlantic origin
        indcount (array): as above but for the Indian Ocean
        paccount (array): as above but for the Pacific
        arccount (array): as above but for the Arctic
        soucount (array): as above but for the Southern Ocean
    
    Returns:
        None.
    """
    ##### NEED ALL LINES DEFINING SEGMENTS FOR TRAJECTORY RELEASE LINES #######
                # SHOULD PROBABLY MAKE THIS ITS OWN FUNCTION #
    sheddir = '/home/np838619/Watershed/shed_defs/' # directory with data
    
    endpts1 = pl.genfromtxt(sheddir+'Am_clicks.txt',skip_header=5) # Americas
    endpts2 = pl.genfromtxt(sheddir+'AfMe_clicks.txt',skip_header=5) # Africa
    endpts3 = pl.genfromtxt(sheddir+'EAA_clicks.txt',skip_header=5) # East Asia
    
    # Arctic sections:
    endpts4 = pl.genfromtxt(sheddir+'ArA_clicks.txt',skip_header=5) # Atlantic
    endpts5 = pl.genfromtxt(sheddir+'ArI_clicks.txt',skip_header=5) # Indian
    endpts6 = pl.genfromtxt(sheddir+'ArP_clicks.txt',skip_header=5) # Pacific
    
    # Southern lines:
    endpts7 = pl.genfromtxt(sheddir+'SOA_clicks.txt',skip_header=5) # Atlantic
    endpts8 = pl.genfromtxt(sheddir+'SOI_clicks.txt',skip_header=5) # Indian
    endpts9 = pl.genfromtxt(sheddir+'SOP_clicks.txt',skip_header=5) # Pacific
    
    # Need traj release points for Arctic & Southern catchments because
    # end points lines means that latitude lines aren't correctly followed on
    # polar projections
    endpts10 = pl.genfromtxt(sheddir+'Ar_traj_release_new.txt',skip_header=5)
    endpts11 = pl.genfromtxt(sheddir+'SO_traj_release_new.txt',skip_header=5)
    
    ###### NEED 1 MAP FOR ATLANTIC/INDIAN CATCHMENTS, ANOTHER FOR PACIFIC #####
    ####### ALSO NEED SEPERATE MAPS FOR ARCTIC AND SOUTHERN CATCHMENTS ########
    m1 = Basemap(projection='cyl',llcrnrlon=-180,llcrnrlat=-80,urcrnrlon=180,urcrnrlat=80)
    m2 = Basemap(projection='cyl',llcrnrlon=0,llcrnrlat=-80,urcrnrlon=360,urcrnrlat=80)
    m3 = Basemap(projection='npaeqd',boundinglat=25,lon_0=270,round=True)
    m4 = Basemap(projection='spstere',boundinglat=-13,lon_0=270,round=True)

    ####### DEFINE POLYGONS FOR ATLANTIC, INDIAN AND PACIFIC CATCHMENTS #######
    atl_bnd = CatchLine(endpts1,endpts4,endpts2,endpts7) # Atlantic boundary
    ind_bnd = CatchLine(endpts2,endpts5,endpts3,endpts8) # Indian boundary
    # Negative lon co-ords in endpts1 & endpts6, correct this:
    e1 = Add360(endpts1); e6 = Add360(endpts6)
    pac_bnd = CatchLine(endpts3,e6,e1,endpts9) # Pacific boundary
    # empty arrays for Arctic & Southern map co-ords, convert to map co-ords:
    arc_bnd = pl.zeros_like(endpts10); arc_bnd[:,0], arc_bnd[:,1] = m3(endpts10[:,0],endpts10[:,1])
    sou_bnd = pl.zeros_like(endpts11); sou_bnd[:,0], sou_bnd[:,1] = m4(endpts11[:,0],endpts11[:,1])
    
    ### SET UP 3X2 PANEL PLOT - NOT SURE WHAT TO DO WITH FINAL PANEL THOUGH ###
    m = Basemap(projection='robin',lon_0=0.); my = Basemap(projection='robin',lon_0=180.)
    fig, ax = pl.subplots(3,2,figsize=(13,10))
    pl.subplot(321)
    m.drawcoastlines(); m.plot(atl_bnd[:,0],atl_bnd[:,1],latlon=True,color='b')
    x,y = m(atlcount[:,-1],atlcount[:,-2]); m.scatter(x,y,color='r')
    pl.subplot(322)
    m.drawcoastlines(); m.plot(ind_bnd[:,0],ind_bnd[:,1],latlon=True,color='b')
    x,y = m(indcount[:,-1],indcount[:,-2]); m.scatter(x,y,color='r')
    pl.subplot(323)
    my.drawcoastlines(); my.plot(pac_bnd[:,0],pac_bnd[:,1],latlon=True,color='b')
    x,y = my(paccount[:,-1],paccount[:,-2]); my.scatter(x,y,color='r')
    pl.subplot(324)
    m3.drawcoastlines(); m3.plot(arc_bnd[:,0],arc_bnd[:,1],latlon=False,color='b')
    x,y = m3(arccount[:,-1],arccount[:,-2]) ;m3.scatter(x,y,color='r')
    pl.subplot(325)
    m4.drawcoastlines(); m4.plot(sou_bnd[:,0],sou_bnd[:,1],latlon=False,color='b')
    #x,y = m4(soucount[:,-1],soucount[:,-2]); m4.scatter(x,y,color='r')
    pl.subplot(326)
    m.drawcoastlines()
    pl.tight_layout(); pl.subplots_adjust(wspace=-0.5)

def FluxPartition(catchcount,rlslabs,normals,etalevs,q_tj,surfp,u,v,pres_tj):
    """Function to extract the moisture flux with origin in a particular
    ocean catchment area. This is the version where the flux is weighted by dl.
    
    Args:
        catchcount (array): origin type, timestep & location for trajectories
                          with origin in a particular ocean catchment area
        rlslabs (array): longitude, latitude co-ordinates of trajectory release
                          points and segment label
        normals (array): normal vectors for each segment along the watershed
        etalevs (array): model levels from trajectories were released
        q_tj (array): specific humidity (kg/kg) along all trajectories
        surfp (array): surface pressure (hPa) a t=0 for all trajectories
        u (array): zonal velocity (m/s) of trajectory at t=0
        v (array): meridional velocity (m/s) of trajectory at t=0
        pres_tj (array): pressure (hPa) along all trajectories
    
    Returns:
        shedflx (array): vertically-integrated moisture flux in kg/s.
    """
    global FluxCalc
    NPART = q_tj.shape[0] # number of trajectories
    q_catch = pl.zeros_like(q_tj) # empty array for q with origin in this catchment
    
    if catchcount.size == 0.: # if no trajectories with origin here...
            q_catch[:] = 0. # ... set entire q array to zero
    else: # otherwise ...
        for trajno in range(NPART):  # ... loop over all trajectories ...
            if trajno in catchcount[:,0]: # ... if this traj has origin here ..
                q_catch[trajno,0] = q_tj[trajno,0] # ... save its q
            #else: # ... otherwise ...
             #   q_catch[trajno,0] = 0. # ... set its q to zero
    
    # calculate vertically-integrated moisture flux weighted by dl:
    shedflx = FluxCalc(rlslabs,u,v,surfp,etalevs,normals,q_catch,pres_tj)
    
    return shedflx

def FluxPart2(catchcount,rlslabs,normals,etalevs,q_tj,surfp,u,v,pres_tj):
    """Function to extract the moisture flux with origin in a particular
    ocean catchment area. This is the version where the flux is unweighted.
    
    Args:
        catchcount (array): origin type, timestep & location for trajectories
                          with origin in a particular ocean catchment area
        rlslabs (array): longitude, latitude co-ordinates of trajectory release
                          points and segment label
        normals (array): normal vectors for each segment along the watershed
        etalevs (array): model levels from trajectories were released
        q_tj (array): specific humidity (kg/kg) along all trajectories
        surfp (array): surface pressure (hPa) a t=0 for all trajectories
        u (array): zonal velocity (m/s) of trajectory at t=0
        v (array): meridional velocity (m/s) of trajectory at t=0
        pres_tj (array): pressure (hPa) along all trajectories
    
    Returns:
        shedflx (array): vertically-integrated moisture flux in kg/m/s.
    """
    global FluxCalc2
    NPART = q_tj.shape[0] # number of trajectories
    q_catch = pl.zeros_like(q_tj) # empty array for q with origin in this catchment
    
    if catchcount.size == 0.: # if no trajectories with origin here...
            q_catch[:] = 0. # ... set entire q array to zero
    else: # otherwise ...
        for trajno in range(NPART):  # ... loop over all trajectories ...
            if trajno in catchcount[:,0]: # ... if this traj has origin here ..
                q_catch[trajno,0] = q_tj[trajno,0] # ... save its q
            else: # ... otherwise ...
                q_catch[trajno,0] = 0. # ... set its q to zero
    
    # calculate unweighted vertically-integrated moisture flux;
    shedflx = FluxCalc2(rlslabs,u,v,surfp,etalevs,normals,q_catch,pres_tj)
    
    return shedflx

def FluxesLoop(catchcount,rlslabs,normals,etalevs,q,sp,u,v,pres):
    """Function to calculate the partitioned moisture flux for each release for
    a particular catchment. In this version fluxes are weighted by dl.
    
    Args:
        catchcount (array): origin type, timestep & location for trajectories
                          with origin in a particular ocean catchment area
        rlslabs (array): longitude, latitude co-ordinates of trajectory release
                          points and segment label
        normals (array): normal vectors for each segment along the watershed
        etalevs (array): model levels from trajectories were released
        q_tj (array): specific humidity (kg/kg) along all trajectories
        surfp (array): surface pressure (hPa) a t=0 for all trajectories
        u (array): zonal velocity (m/s) of trajectory at t=0
        v (array): meridional velocity (m/s) of trajectory at t=0
        pres_tj (array): pressure (hPa) along all trajectories
    
    Returns:
        catchflx (array): vertically-integrated moisture flux in kg/s for each
                          release
    """
    global FluxPartition
    # empty array for fluxes:
    catchflx = pl.zeros([len(catchcount),rlslabs.shape[0]])
    catchcount = pl.asarray(catchcount) # change from list to array
    
    for rls in range(len(catchcount)): #  loop of no. of releases
        # calculate vertically integrated moisture fluxes, weighted by dl
        catchflx[rls] = FluxPartition(catchcount[rls],rlslabs,normals,etalevs,
                                    q[rls],sp[rls,0],u[rls],v[rls],pres[rls])
    
    return catchflx

def FluxesLoop2(catchcount,rlslabs,normals,etalevs,q,sp,u,v,pres):
    """Function to calculate the partitioned moisture flux for each release for
    a particular catchment. In this version fluxes are unweighted.
    
    Args:
        catchcount (array): origin type, timestep & location for trajectories
                          with origin in a particular ocean catchment area
        rlslabs (array): longitude, latitude co-ordinates of trajectory release
                          points and segment label
        normals (array): normal vectors for each segment along the watershed
        etalevs (array): model levels from trajectories were released
        q_tj (array): specific humidity (kg/kg) along all trajectories
        surfp (array): surface pressure (hPa) a t=0 for all trajectories
        u (array): zonal velocity (m/s) of trajectory at t=0
        v (array): meridional velocity (m/s) of trajectory at t=0
        pres_tj (array): pressure (hPa) along all trajectories
    
    Returns:
        catchflx (array): vertically-integrated moisture flux in kg/s for each
                          release
    """
    global FluxPart2
    # empty array for fluxes:
    catflx = pl.zeros([len(catchcount),rlslabs.shape[0]])
    catchcount = pl.asarray(catchcount) # change from list to array
    
    for rls in range(len(catchcount)): #  loop of no. of releases
        # calculate vertically-integrated unweighted moisture fluxes
        catflx[rls] = FluxPart2(catchcount[rls],rlslabs,normals,etalevs,
                                    q[rls],sp[rls,0],u[rls],v[rls],pres[rls])
    
    return catflx


def CatchProps(catchcount,TJTOT):
    """Function to calculate the proportion of trajectories with origin in a
    particular catchment.
    
    Args:
        catchcount (array): origin type, timestep & location for trajectories
                          with origin in a particular ocean catchment area
        TJTOT (int): total number of trajectories
    
    Returns:
        catchprop (float): percentage of trajectories with origin in a particular
                            ocean catchment area
    """
    empty = [] # list for indices of empty arrays
    for rls in range(len(catchcount)): # loop over no. of releases
        if catchcount[rls].size == 0.: # if no traj have origin at this release
            empty.append(rls) # ... add this release index to list
    
    if len(empty) == 0.: # if list has zero length ...
        pass # .. move on
    else: # otherwise ...
        catchcount = pl.delete(catchcount,empty) # . delete the empty arrays
    
    if len(catchcount) == 0.: # if catchcount is now completely empty ...
        catchprop = 0. # ... proportion of trajectories is zero
        CAT1 = 0; CAT2 = 0 # .. proportions of CAT1 & CATII origins is zero
    else: # otherwise ...
        catchcount = pl.vstack(catchcount) # .. make 1 array
        catchprop = (catchcount.shape[0]/TJTOT)*100 # calculate proportion
        a = pl.where(catchcount[:,1]==1); b = pl.where(catchcount[:,1]==2)
        CAT1 = (a[0].size/TJTOT)*100; CAT2 = (b[0].size/TJTOT)*100
    
    return catchprop, CAT1, CAT2

def FluxCATFull(newlabs,rlslabs,normals,etalevs,q_tj,surfp,u,v,pres_tj):
    """
    """
    NPART = q_tj.shape[1]
    NRLS = newlabs.shape[0]
    NREL = rlslabs.shape[0]
    shedflx1 = pl.zeros([NRLS,NREL]); shedflx2 = pl.zeros_like(shedflx1)
    
    for rls in range(NRLS):
        q1 = pl.zeros_like(q_tj[rls]); q2 = pl.zeros_like(q1)
        for traj in range(NPART):
            if newlabs[rls,traj,0] == 1.:
                q1[traj,:] = q_tj[rls,traj,:]
            elif newlabs[rls,traj,0] == 2.:
                q2[traj,:] = q_tj[rls,traj,:]
            else:
                q1[traj,:] = 0.
                q2[traj,:] = 0.
        shedflx1[rls] = FluxCalc(rlslabs,u[rls],v[rls],surfp[rls,0],etalevs,
                                normals,q1,pres_tj[rls])
        shedflx2[rls] = FluxCalc(rlslabs,u[rls],v[rls],surfp[rls,0],etalevs,
                                normals,q2,pres_tj[rls])
    
    return shedflx1, shedflx2

def FluxCATFull2(newlabs,rlslabs,normals,etalevs,q_tj,surfp,u,v,pres_tj):
    """
    """
    NPART = q_tj.shape[1]
    NRLS = newlabs.shape[0]
    NREL = rlslabs.shape[0]
    shedflx1 = pl.zeros([NRLS,NREL]); shedflx2 = pl.zeros_like(shedflx1)
    
    for rls in range(NRLS):
        q1 = pl.zeros_like(q_tj[rls]); q2 = pl.zeros_like(q1)
        for traj in range(NPART):
            if newlabs[rls,traj,0] == 1.:
                q1[traj,:] = q_tj[rls,traj,:]
            elif newlabs[rls,traj,0] == 2.:
                q2[traj,:] = q_tj[rls,traj,:]
            else:
                q1[traj,:] = 0.
                q2[traj,:] = 0.
        shedflx1[rls] = FluxCalc2(rlslabs,u[rls],v[rls],surfp[rls,0],etalevs,
                                normals,q1,pres_tj[rls])
        shedflx2[rls] = FluxCalc2(rlslabs,u[rls],v[rls],surfp[rls,0],etalevs,
                                normals,q2,pres_tj[rls])
    
    return shedflx1, shedflx2

def FluxCATPart(newlabs,catchcount,rlslabs,normals,etalevs,q_tj,surfp,u,v,pres_tj):
    """Function to split the moisture flux with origin in a particular catchment
    into CAT I (boundary layer) & CAT II (partial mixing) contributions. This
    is the version for calculating the fluxes weighted by dl.
    
    Args:
        newlabs (array): origin type (CAT I, CAT II, strat or none), timestep
                         & location for all trajectories
        catchcount (array): origin type, timestep & location for trajectories
                          with origin in a particular ocean catchment area
        rlslabs (array): longitude, latitude co-ordinates of trajectory release
                          points and segment label
        normals (array): normal vectors for each segment along the watershed
        etalevs (array): model levels from trajectories were released
        q_tj (array): specific humidity (kg/kg) along all trajectories
        surfp (array): surface pressure (hPa) a t=0 for all trajectories
        u (array): zonal velocity (m/s) of trajectory at t=0
        v (array): meridional velocity (m/s) of trajectory at t=0
        pres_tj (array): pressure (hPa) along all trajectories
    
    Returns:
        shedflx1 (array): vertically-integrated moisture flux from CAT I origins
        shedflx2 (array): vertically-integrated moisture flux from CAT II origins
    """
    global FluxPartition
    NPART = q_tj.shape[1] # number of trajectories
    NRLS = newlabs.shape[0] # number of releases
    NREL = rlslabs.shape[0]
    shedflx1 = pl.zeros([NRLS,NREL]); shedflx2 = pl.zeros_like(shedflx1)
    
    for rls in range(NRLS):
        # empty arrays for CAT I & CAT II q:
        q_cat1 = pl.zeros_like(q_tj[rls]); q_cat2 = pl.zeros_like(q_cat1)
        for traj in range(NPART): # loop over trajectories
            if newlabs[rls,traj,0] == 1.: # if traj has CAT I origin ...
                q_cat1[traj,:] = q_tj[rls,traj,:] # ... save its q
            elif newlabs[rls,traj,0] == 2.: # if traj has CAT II origin ...
                q_cat2[traj,:] = q_tj[rls,traj,:] # ... save its q
            else: # otherwise ...
                q_cat1[traj,:] = 0. # set to zero
                q_cat2[traj,:] = 0. # set to zero
        shedflx1[rls] = FluxPartition(catchcount[rls],rlslabs,normals,etalevs,
                                q_cat1,surfp[rls,0],u[rls],v[rls],pres_tj[rls])
        shedflx2[rls] = FluxPartition(catchcount[rls],rlslabs,normals,etalevs,
                                q_cat2,surfp[rls,0],u[rls],v[rls],pres_tj[rls])
    
    return shedflx1, shedflx2

def FluxCATPart2(newlabs,catchcount,rlslabs,normals,etalevs,q_tj,surfp,u,v,pres_tj):
    """Function to split the moisture flux with origin in a particular catchment
    into CAT I (boundary layer) & CAT II (partial mixing) contributions. This 
    is the version for calculating the unweighted moisture flux.
    
    Args:
        newlabs (array): origin type (CAT I, CAT II, strat or none), timestep
                         & location for all trajectories
        catchcount (array): origin type, timestep & location for trajectories
                          with origin in a particular ocean catchment area
        rlslabs (array): longitude, latitude co-ordinates of trajectory release
                          points and segment label
        normals (array): normal vectors for each segment along the watershed
        etalevs (array): model levels from trajectories were released
        q_tj (array): specific humidity (kg/kg) along all trajectories
        surfp (array): surface pressure (hPa) a t=0 for all trajectories
        u (array): zonal velocity (m/s) of trajectory at t=0
        v (array): meridional velocity (m/s) of trajectory at t=0
        pres_tj (array): pressure (hPa) along all trajectories
    
    Returns:
        shedflx1 (array): vertically-integrated moisture flux from CAT I origins
        shedflx2 (array): vertically-integrated moisture flux from CAT II origins
    """
    global FluxPart2
    NPART = q_tj.shape[1] # number of trajectories
    NRLS = newlabs.shape[0] # number of releases
    NREL = rlslabs.shape[0]
    shedflx1 = pl.zeros([NRLS,NREL]); shedflx2 = pl.zeros_like(shedflx1)
    
    for rls in range(NRLS):
        # empty arrays for CAT I & CAT II q:
        q_cat1 = pl.zeros_like(q_tj[rls]); q_cat2 = pl.zeros_like(q_cat1)
        for traj in range(NPART): # loop over trajectories
            if newlabs[rls,traj,0] == 1.: # if traj has CAT I origin ...
                q_cat1[traj,:] = q_tj[rls,traj,:] # ... save its q
            elif newlabs[rls,traj,0] == 2.: # if traj has CAT II origin ...
                q_cat2[traj,:] = q_tj[rls,traj,:] # ... save its q
            else: # otherwise ...
                q_cat1[traj,:] = 0. # set to zero
                q_cat2[traj,:] = 0. # set to zero
        shedflx1[rls] = FluxPart2(catchcount[rls],rlslabs,normals,etalevs,
                                q_cat1,surfp[rls,0],u[rls],v[rls],pres_tj[rls])
        shedflx2[rls] = FluxPart2(catchcount[rls],rlslabs,normals,etalevs,
                                q_cat2,surfp[rls,0],u[rls],v[rls],pres_tj[rls])
    
    return shedflx1, shedflx2

def FluxCATLoop(newlabs,catchcount,rlslabs,normals,etalevs,q,surfp,u,v,pres):
    """
    """
    global FluxCatPart
    NRLS = len(catchcount); NREL = rlslabs.shape[0]
    CAT1 = pl.zeros([NRLS,NREL]); CAT2 = pl.zeros_like(CAT1)
    
    for rls in range(NRLS):
        CAT1[rls], CAT2[rls] = FluxCATPart(newlabs[rls],catchcount[rls],rlslabs,
                                            normals,etalevs,q[rls],surfp[rls,0],
                                            u[rls],v[rls],pres[rls])
    
    return CAT1, CAT2

def FluxCATLoop2(newlabs,catchcount,rlslabs,normals,etalevs,q,surfp,u,v,pres):
    """
    """
    global FluxCATPart2
    NRLS = len(catchcount); NREL = rlslabs.shape[0]
    CAT1 = pl.zeros([NRLS,NREL]); CAT2 = pl.zeros_like(CAT1)
    
    for rls in range(NRLS):
        CAT1[rls], CAT2[rls] = FluxCATPart2(newlabs[rls],catchcount[rls],rlslabs,
                                            normals,etalevs,q[rls],surfp[rls,0],
                                            u[rls],v[rls],pres[rls])
    
    return CAT1, CAT2

def CAT1LandSeaFull(lsm,newlabs,rlslabs,normals,etalevs,q_tj,surfp,u,v,pres_tj,
                    lon_tj,lat_tj,eralon,eralat):
    """
    """
    global FluxCalc
    NPART = q_tj.shape[1] # no. of trajectories
    NRLS = q_tj.shape[0] # no. of releases
    NREL = rlslabs.shape[0] # no. of pts
    shedflx_sea = pl.zeros([NRLS,NREL]); shedflx_lnd = pl.zeros_like(shedflx_sea)
    
    for rls in range(NRLS):
        q_land = pl.zeros_like(q_tj[rls]); q_sea = pl.zeros_like(q_tj[rls])
        for traj in range(NPART):
            if newlabs[rls,traj,0] != 1.:
                continue
            pt = (newlabs[rls,traj,-1],newlabs[rls,traj,-2])
            ts = newlabs[rls,traj,1]#; print traj
            #if ts.dtype == 'float64':
            #    continue
                #p2 = (lon_tj[rls,traj,ts],lat_tj[rls,traj,ts])
            
            x = BilinInterp(pt,eralon,eralat,lsm)
            #x2 = BilinInterp(p2,eralon,eralat,lsm)
            if x <= 0.5:
                q_sea[traj,:] = q_tj[rls,traj,:]
            #elif x2 <= 0.5:
            #    q_sea[traj,:] = q_tj[rls,traj,:]
            else:
                q_land[traj,:] = q_tj[rls,traj,:]
        shedflx_sea[rls] = FluxCalc(rlslabs,u[rls],v[rls],surfp[rls,0],
                                    etalevs,normals,q_sea,pres_tj[rls])
        shedflx_lnd[rls] = FluxCalc(rlslabs,u[rls],v[rls],surfp[rls,0],
                                    etalevs,normals,q_land,pres_tj[rls])
    
    return shedflx_sea, shedflx_lnd

def CAT1LandSea(lsm,catchcount,rlslabs,normals,etalevs,q_tj,surfp,u,v,pres_tj,
                lon_tj,lat_tj,eralon,eralat):
    """
    """
    global FluxPartition
    NPART = q_tj.shape[1] # no. of trajectories
    NRLS = q_tj.shape[0] # no. of releases
    NREL = rlslabs.shape[0] # no. of pts
    shedflx_sea = pl.zeros([NRLS,NREL]); shedflx_lnd = pl.zeros_like(shedflx_sea)
    
    for rls in range(NRLS): # loop over releases
        q_land = pl.zeros_like(q_tj[rls]); q_sea = pl.zeros_like(q_tj[rls])
        for trajno in range(catchcount[rls].shape[0]):
            tj = catchcount[rls][trajno,0]; ts = catchcount[rls][trajno,2]
            if catchcount[rls][trajno,1] == 1.:
                pt = (catchcount[rls][trajno,-1],catchcount[rls][trajno,-2])
#                if ts < q_tj.shape[-1]-1:
#                    p2 = (lon_tj[rls,tj,ts+1],lat_tj[rls,tj,ts+1])
#                elif ts == q_tj.shape[-1]-1:
#                    p2 = (lon_tj[rls,tj,ts],lat_tj[rls,tj,ts])
                x = BilinInterp(pt,eralon,eralat,lsm)
                #x2 = BilinInterp(p2,eralon,eralat,lsm)
                if x <= 0.5:
                    q_sea[tj,:] = q_tj[rls,tj,:]
                #elif x2 <= 0.5:
                #    q_sea[tj,:] = q_tj[rls,tj,:]
                else:
                    q_land[tj,:] = q_tj[rls,tj,:]
        shedflx_sea[rls] = FluxPartition(catchcount[rls],rlslabs,normals,etalevs,
                                q_sea,surfp[rls,0],u[rls],v[rls],pres_tj[rls])
        shedflx_lnd[rls] = FluxPartition(catchcount[rls],rlslabs,normals,etalevs,
                                q_land,surfp[rls,0],u[rls],v[rls],pres_tj[rls])
    
    return shedflx_sea, shedflx_lnd

def CAT1BL_full(lsm,newlabs,rlslabs,normals,etalevs,q_tj,surfp,u,v,pres_tj,
                eralon,eralat):
    """
    """
    global FluxCalc
    NPART = q_tj.shape[1] # no. of trajectories
    NRLS = q_tj.shape[0] # no. of releases
    NREL = rlslabs.shape[0] # no. of pts
    shedflx_BL = pl.zeros([NRLS,NREL])
    
    for rls in range(NRLS):
        q_BL = pl.zeros_like(q_tj[rls])
        for traj in range(NPART):
            pt = (newlabs[rls,traj,-1],newlabs[rls,traj,-2])
            if newlabs[rls,traj,0] == 1.:
                if newlabs[rls,traj,1] == 0.:
                #x = BilinInterp(pt,eralon,eralat,lsm)
                #if x > 0.5:
                    q_BL[traj,:] = q_tj[rls,traj,:]
        shedflx_BL[rls] = FluxCalc(rlslabs,u[rls],v[rls],surfp[rls,0],
                                    etalevs,normals,q_BL,pres_tj[rls])
      
    return shedflx_BL

def CAT1BL_part(lsm,catchcount,rlslabs,normals,etalevs,q_tj,surfp,u,v,pres_tj,
                eralon,eralat):
    """
    """
    global FluxPartition
    NPART = q_tj.shape[1] # no. of trajectories
    NRLS = q_tj.shape[0] # no. of releases
    NREL = rlslabs.shape[0] # no. of pts
    shedflx_BL = pl.zeros([NRLS,NREL])
    
    for rls in range(NRLS): # loop over releases
        q_BL = pl.zeros_like(q_tj[rls])
        for trajno in range(catchcount[rls].shape[0]):
            tj = catchcount[rls][trajno,0]
            pt = (catchcount[rls][trajno,-1],catchcount[rls][trajno,-2])
            if catchcount[rls][trajno,1] == 1.:
                if catchcount[rls][trajno,2] == 0.:
                #x = BilinInterp(pt,eralon,eralat,lsm)
                #if x > 0.5:
                    q_BL[tj,:] = q_tj[rls,tj,:]
        shedflx_BL[rls] = FluxPartition(catchcount[rls],rlslabs,normals,etalevs,
                                q_BL,surfp[rls,0],u[rls],v[rls],pres_tj[rls])

    return shedflx_BL


def DensityPlot(catchcount,rlspts):
    """Function to make a density plot from trajectory origin data.
    
    Args:
        catchcount (array): origin type, timestep & location for trajectories
                          with origin in a particular ocean catchment area
    
    Returns:
        None.
    """
    m = Basemap(projection='robin',lon_0=180.) # set up a map
    xpt,ypt = m(rlspts[:,0],rlspts[:,1])
    m.plot(xpt[:161],ypt[:161],color='k'); m.plot(xpt[161:],ypt[161:],color='k')
    nx = 100; ny = 50 # number of bin
    lon_bins = pl.linspace(catchcount[:,-1].min(),catchcount[:,-1].max(),nx+1)
    lat_bins = pl.linspace(catchcount[:,-2].min(),catchcount[:,-2].max(),ny+1)
    
    density, _, _ = pl.histogram2d(catchcount[:,-2],catchcount[:,-1],
                                   [lat_bins,lon_bins],normed=True)
    lon_bins2d, lat_bins2d = pl.meshgrid(lon_bins,lat_bins)
    xs,ys = m(lon_bins2d,lat_bins2d)
    
    # need to to do this for countour plotting
    density = pl.hstack((density,pl.zeros((density.shape[0],1))))
    density = pl.vstack((density,pl.zeros((density.shape[1]))))
    
    b = density[pl.nonzero(density)]
    levels = [0.0001,0.0004,0.0007,0.001,0.002]
    m.contourf(xs,ys,density,cmap='Accent',extend='max',levels=levels)
    m.drawcoastlines(); m.colorbar()
    
    return density

def MeanFluxesOutput(year,shed,means,cors,mags,cors2):
    """Function to write a text file with the mean net partitioned fluxes &  
    their magnitudes along with the correlations to the total flux.
    
    Args:
        year (string): from trajyear.py, runs at start of main code
        shed (string): three-letter string referring to the catchment boundary,
                        defined at the start of the main code
        means (array): mean net partioned fluxes
        cors (array): pearson correlation coefficients of the mean net partioned
                      fluxes to the total net flux
       mags (array): mean magnitudes of partitioned fluxes
       cors2 (array): pearson correlation coefficients of the mean partitioned
                      fluxes to the total flux magnitude
    """
    clusdir = '/glusterfs/scenario/users/np838619/traj/'
    f = open(clusdir + year +'/traj_' + shed + year + '.txt','w')
    # net fluxes:
    f.write('Mean Fluxes\n')
    pl.savetxt(f,means); f.write('Correlations\n')
    pl.savetxt(f,cors); f.write('\n')
    # Flux magnitudes
    f.write('Mean flux magnitudes\n')
    pl.savetxt(f,mags); f.write('\n'); f.write('Correlations\n')
    pl.savetxt(f,cors2); f.close()
    print 'Mean moisture flux data written to file'
    
    return None

def TrajPropsOutput(year,shed,catchprops,propsplit):
    """
    """
    clusdir = '/glusterfs/scenario/users/np838619/traj/'
    f = open(clusdir + year +'/props_' + shed + year + '.txt','w')
    f.write('Catchment proportions\n')
    pl.savetxt(f,catchprops)
    f.write('\n')
    f.write('Category 1 & 2 proportions')
    pl.savetxt(f,propsplit)
    f.close()
    print 'Trajectory origin proportions written to file'
    
    f.close()

def CATsplitOutput(year,shed,catsplits):
    """
    """
    clusdir = '/glusterfs/scenario/users/np838619/traj/'
    f = open(clusdir + year +'/CATsplit_' + shed + year + '.txt','w')
    f.write('Mean values of CAT I & CAT II moisture fluxes\n')
    pl.savetxt(f,catsplits)
    f.close()
    print 'Mean CAT I & CAT II moisture fluxes written to file'

def SaveTZero(shed,year,stamp,q,PV,T,zbl,pres,u0,v0,ps):
    """
    """
    clusdir = '/glusterfs/scenario/users/np838619/traj/'
    
    q = pl.reshape(q[:,0],(u0.shape[0],u0.shape[1]))
    PV = pl.reshape(PV[:,0],(u0.shape[0],u0.shape[1]))
    T = pl.reshape(T[:,0],(u0.shape[0],u0.shape[1]))
    zbl = pl.reshape(zbl[:,0],(u0.shape[0],u0.shape[1]))
    pres = pl.reshape(pres[:,0],(u0.shape[0],u0.shape[1]))
    
    fields = [q,PV,T,zbl,pres,u0,v0,ps]
    names = ['specific humidty t=0','potential vorticity t=0',
             'temperature t=0','boundary layer height t=0','pressure t=0',
             'u wind t=0','v wind t=0','surface pressure t=0']
    stand = ['q','PV','T','BLZ','pressure','u','v','psurf']
    units_list = ['kg/kg','K m**-2 kg**-1 s**-1','K','m','hPa','m s**-1',
                      'm s**-1','hPa']
    
    clusdir = '/glusterfs/scenario/users/np838619/traj/'
    ncfile = Dataset(clusdir+year+'/'+shed+'/tj_vertprofs_'+stamp[12:]+'_.nc','w')
    
    height = ncfile.createDimension('levels',q.shape[0])
    height = ncfile.createVariable('levels',pl.float64,('levels'))
    height.long_name = 'model levels'
    
    length = ncfile.createDimension('points',q.shape[1])
    length = ncfile.createVariable('points',pl.float64,('points'))
    length.long_name = 'release points'
    
    for i in range(len(fields)):
        C = ncfile.createVariable(names[i],pl.float64,
                                      ('levels','points'))
        C.units = units_list[i]
        C.standard_name = stand[i]
        C[:,:] = fields[i][:,:]
    
    ncfile.close()
    
    return None