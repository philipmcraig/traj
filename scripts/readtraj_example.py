# -*- coding: utf-8 -*-
"""
This script accompanies the example trajectory output data for the following 
paper submitted to Journal of Geophysical Research: Atmospheres

Craig, P.M., Ferreira, D. & Methven,J. (2023) <title undecided>

and as part of the dataset

Craig, P.M. (2023) Back trajectories released from ocean catchment boundaries

submitted to the University of Reading Research Repository (include doi later).

The trajectory initialization files & the full output files are read by this 
code for a user-specified catchment boundary. The boundaries are as follows:
    amr: Americas
    afr: Africa
    eaa: East Asia & Australia (referred to as South-East Asia in the paper)
    ara: Arctic Atlantic
    ari: Arctic Indian
    arp: Arctic Pacific
    soa: Southern Atlantic
    soi: Southern Indian
    sop: Southern Pacific

Data are extracted from the files & stored as numpy arrays.
The following metadata are output by the code:
    initial release time
    number of trajectories per file
    total number of trajectories
    number of vertical levels
    number of timesteps along the trajectories
    length of trajectories in days
    number of attributes interpolated along each trajectory & their names

@author: Philip Craig
"""

from __future__ import division
import pylab as pl
import itertools
import glob as glob

def InitData(initfiles):
    """Function to read in trajectory data from initialization files with the 
    trajectory data at t=0.
    
    Args:
        initfiles (array): filenames of files containing trajectory data at t=0
        
    Returns:
        U (array): zonal velocity at t=0
        V (array): meridional velocity at t=0
    """
    # empty lists for u, v & sp data
    U = []; V = []; SP = []
    for initfile in range(len(initfiles)): # loop over initfiles array
        # read in data file as list, each line is also a list
        flist = ReadTxtFile(initfiles[initfile])
        
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
                                            # then increment by traj*6
            initdata[traj,-1] = flist[18+traj*6][0]
        
        # split data into sub-arrays
        lat = initdata[:,2]#; lon = initdata[:,3]
        u = initdata[:,10] # zonal velocity
        v = initdata[:,11] # meridional velocity
        sp = initdata[:,-1] # surface pressure
        
        NREL = NPART/NCLUST # number of release points
        # change shape of arrays to be no.of levels X no. of release points
        u = pl.reshape(u,(NCLUST,int(NREL))); v = pl.reshape(v,(NCLUST,int(NREL)))
        lat = pl.reshape(lat,(NCLUST,int(NREL))); sp = pl.reshape(sp,(NCLUST,int(NREL)))
        
        a = 6.37e6 # radius of Earth
        # u & v output in units 1/day, convert to m/s:
            # multiply by -1 as it is back trajectory:
        u = -1*u*a*pl.cos(pl.radians(lat))/86400; v = -1*v*a/86400
        
        U.append(u); V.append(v); SP.append(sp) # append to lists
    
    # transform lists to arrays
    U = pl.asarray(U)
    V = pl.asarray(V)
    SP = pl.asarray(SP)

    return U, V, SP

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
    # loop over steps along trajectory
    for step in range(NSTEPS):
        trajdat[step] = flist[18+step+trajno*(intervals+5)]
        # data starts at line 19 in file, start loop there
        # end loop at line 19 plus no. of intervals
        # each element from the lists extracted becomes a row in the trajdat array
    
    return trajdat

shed = 'amr' # specific catchment boundary (e.g. amr = Americas, afr = Africa)
datadir = '/path/to/files/' # specify directory containing files

# glob finds all filenames in a directory
trajfiles = glob.glob(datadir+'utraj-df_'+shed+'*')
trajfiles = pl.asarray(trajfiles)
initfiles = glob.glob(datadir+'utraj-tr_'+shed+'*')
initfiles = pl.asarray(initfiles)

# extract data from initialization files
    #zonal & meridional velocities + surface pressure at t=0
u0, v0, sp_tj = InitData(initfiles)#,datadir)

# set up empty lists for attributes interpolated along trajectories
lat = [] # latitude
lon = [] # longitude
pres = [] # pressure
temp = [] # temperature
PV = [] # potential vorticity
q = [] # specific humidity
height = [] # height above surface
blz = [] # boundary layer height

# loop over trajectory output files (only 1 per boundary in this example)
for trajfile in range(trajfiles.shape[0]):
    flist = ReadTxtFile(trajfiles[trajfile]) # open traj data file
    
    BASETIME = int(flist[0][-1]) # start time of back trajectory
    if trajfile == 0.:
        print 'INITIAL RELEASE TIME:', BASETIME
        NPART = int(flist[3][-1]) # no. of trajectories
        NATTR = int(flist[4][-1]) # no. of attributes
        NCLUST = int(flist[7][-1]) # no. of vertical levels
        pointers = list(itertools.chain(*flist[9:11])) # traj no where each level starts
        pointers = pl.asarray(pointers,dtype='int') # turn into array of integers
        DINT = int(flist[2][3]) # size of data interval in hours
        TIMESTEPS = int(flist[2][-2]) # no. of time steps in an interval
        
        intervals = int(flist[15][4]) # no. of intervals stated in line 16 of file, 5th element of list
        NSTEPS = intervals + 1 #  no. of steps = no. of intervals (plus one)
        data = len(flist[17])-2 + NATTR # no. of rows needed
    		    # some headers in line 18 of file, remove MB units & attributes headers
    		    # add no. of attributes, line 5 in file, final element in list
    
    # empty array for trajectory data, shape no. of steps by no. of data fields
    alltraj = pl.zeros([NPART,NSTEPS,data])
    
    # step, hours, lat, lon, P, dry PT, PV, q, height, BL height
    for traj in range(NPART):
        alltraj[traj] = TrajExtract(NSTEPS,data,flist,traj,intervals)


    	# extract all the variables from trajdat array & name them:
    if trajfile == 0.:
        stepno = alltraj[:,:,0] # step numbers
        hours = alltraj[:,:,1] # hours along trajectory
    lat.append(alltraj[:,:,2]) # latitudes
    lon.append(alltraj[:,:,3]) # longitudes
    pres.append(alltraj[:,:,4]) # pressures
    temp.append(alltraj[:,:,5]) # temperatures
    PV.append(alltraj[:,:,6]) # potential vorticities
    q.append(alltraj[:,:,7]) # specific humidities
    height.append(alltraj[:,:,8]) # heights
    blz.append(alltraj[:,:,-1]) # boundary layer heights

TJTOT = NPART*len(trajfiles)

print NPART, 'TRAJECTORIES PER RELEASE', int(len(trajfiles)), 'RELEASES'
print 'TOTAL NO. OF TRAJECTORIES = ', TJTOT
print NCLUST, 'LEVELS'
print NSTEPS, 'STEPS PER TRAJECTORY, ', int((NSTEPS-1)/4), 'DAY TRAJECTORIES' 
print NATTR, 'ATTRIBUTES: temperature, potential vorticity, specific humidity,',\
            'height, boundary layer height'

# convert lists of arrays to arrays, shape=(len(trajfiles),NPART,NSTEPS)
lat = pl.asarray(lat) # latitude
lon = pl.asarray(lon) # longitude
pres = pl.asarray(pres) # pressure
temp = pl.asarray(temp) # temperature
PV = pl.asarray(PV) # potential vorticity
q = pl.asarray(q) # specific humidity
height = pl.asarray(height) # height
blz = pl.asarray(blz) # boundary layer height