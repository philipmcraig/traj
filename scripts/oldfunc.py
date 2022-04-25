# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 15:48:35 2016

@author: np838619
"""

def TrajDayCount(NPART,NSTEPS,blz,height):
    """Function to calculate the number & percentage of trajectories from one 
    release that have BL origin at each time step.
    
    Args:
        NPART (int): number of trajectories
        NSTEPS (int): number of time steps
        blz (array): boundary layer heights of all trajectories at all time steps
        height (array): height of particle for all trajectories at all time steps
    
    Returns:
        blprop (array): percentage of trajectories which have BL origin at each
                        time step
    """
    blcount = pl.zeros([NPART,NSTEPS])
    blsums = []
    for trajno in range(NPART):
        #blcount = 0.
        for step in range(NSTEPS):
            bl_ind = pl.where(blz[trajno,0:step+1]>height[trajno,0:step+1])
            if bl_ind[0].size > 0:
                blcount[trajno,step] = blcount[trajno,step] + 1.
            elif bl_ind[0].size == 0.:
                blcount[trajno,step] = blcount[trajno,step]
    
    blsums = pl.sum(blcount,axis=0)
    blprop = (blsums/NPART)*100  
    
    return blprop

def PMDayCount(equiv,q):
    """Function to calculate the percentage of trajectories with partial mixing
    origin at each time step.
    
    Args:
        equiv (array): equivalent potential temperature (K) of air parcel at all 
                        points along all trajectories
        q (array): specific humidity (kg/kg) of air parcel at all points along 
                    all trajectories
    
    Returns:
        pmprop (array): percentages of trajectories which have partial mixing
                        origin at all time steps
    """
    # number of trajtories & timesteps:
    NPART = equiv.shape[0]; NSTEPS = equiv.shape[1]
    pmcount = pl.zeros([NPART,NSTEPS]) # start with empty matrix for pmcount
    
    for traj in range(NPART):
        for step in range(NSTEPS):
            # find location of PM and if it occurs:
            a,b = PartialMixing(equiv[traj,0:step+1],q[traj,0:step+1])
            if b == True: 
                pmcount[traj,step] = 1. # if PM occurs, pmcount element = 1
            else: # if no PM occurs then pmcount stays 0
                pmcount[traj,step] = 0.
    
    # total no. of trajectories that experienced PM at each timestep:
    pmsums = pl.sum(pmcount,axis=0)
    # calculate % of trajectories which experienced PM eat each timestep:
    pmprop = (pmsums/NPART)*100
    
    return pmprop