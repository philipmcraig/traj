# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 14:41:13 2017

@author: np838619
"""

def OriginCheck(blorigin,pmorigin,storigin,NPART):
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
    if blorigin.size == 0. and pmorigin.size == 0. and storigin.size == 0.:
            return originlabs
    else: 
        for traj in range(NPART): # check all trajectories
            if traj in blorigin[:,0]: # if trajectory has bl origin,
                x = pl.where(blorigin[:,0]==traj); x = x[0][0]
                if traj in storigin[:,0]:
                    y = pl.where(storigin[:,0]==traj); y = y[0][0]
                    if blorigin[x,1]/6. < storigin[y,1]:
                        originlabs[traj,0] = 1.   # assign label '1'
                        originlabs[traj,1] = blorigin[x,1]/6.
                        originlabs[traj,2:] = blorigin[x,2:]
                    elif storigin[y,1] < blorigin[x,1]/6.:
                        originlabs[traj,0] = -1.   # assign label '-1'
                        originlabs[traj,1] = storigin[y,1]
                        originlabs[traj,2:] = storigin[y,2:]
                else:
                    originlabs[traj,0] = 1.   # assign label '1'
                    originlabs[traj,1] = blorigin[x,1]/6.
                    originlabs[traj,2:] = blorigin[x,2:]
            elif pmorigin.size == 0: # check if there are any pm origins
                originlabs[traj,0] = 0.   # assign label '0'
                originlabs[traj,1:] = None
            elif traj in pmorigin[:,0]:
                x = pl.where(pmorigin[:,0]==traj); x = x[0][0]
                if traj in storigin[:,0]:
                    y = pl.where(storigin[:,0]==traj); y = y[0][0]
                    if pmorigin[x,1]/6. < storigin[y,1]:
                        originlabs[traj,0] = 2.   # assign label '2'
                        originlabs[traj,1] = pmorigin[x,1]/6.
                        originlabs[traj,2:] = pmorigin[x,2:]
                    elif storigin[y,1] < pmorigin[x,1]/6.:
                        originlabs[traj,0] = -1.   # assign label '-1'
                        originlabs[traj,1] = storigin[y,1]
                        originlabs[traj,2:] = storigin[y,2:]
                else:
                    originlabs[traj,0] = 2.   # assign label '1'
                    originlabs[traj,1] = pmorigin[x,1]/6.
                    originlabs[traj,2:] = pmorigin[x,2:]
            elif traj in storigin[:,0]:
                y = pl.where(storigin[:,0]==traj); y = y[0][0]
                originlabs[traj,0] = -1.   # assign label '-1'
                originlabs[traj,1] = storigin[y,1]
                originlabs[traj,2:] = storigin[y,2:]
            else: # or if no origin can be found, origin is unknown.
                originlabs[traj] = originlabs[traj]
    
    return originlabs