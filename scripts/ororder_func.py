# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 15:32:56 2017

@author: np838619
"""

def OrOrder(blz,height,stepno,lat,lon,equiv,q,stratloc):
    """
    """
    # Take NPART & NSTEPS from dimensions of an attribute array:
    NPART = blz.shape[0]; NSTEPS = blz.shape[1]
    blors = pl.zeros([NPART,NSTEPS]) # empty array for origin labels
    pmors = pl.zeros_like(blors); pmdis = pl.zeros_like(blors) # discounted pmorigins
    strators = pl.zeros_like(blors); pmout = pl.zeros_like(blors) # reassigned pmorigins
    
    # Use BoundaryLayer & PMcalcs functions to find all trajectories that have
    # boundary layer & partial mixing origins & their origin timesteps:
    blc, blp, blorigin = BoundaryLayer(blz,height,stepno,lat,lon)
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

    return blprop, pmprop, strprop, newlabs#, ORS