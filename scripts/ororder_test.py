# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 10:15:41 2017

@author: np838619
"""

import pylab as pl

NPART = blz[0].shape[0]; NSTEPS = blz[0].shape[1]
NREL = rlslabs.shape[0] # get no. of release points from shape of rlslabs array
NCLUST = int(NPART/NREL) # no. of levels is no. trajectories / no. points
# make empty arrays for flux proportions, size no. of steps:
flux_prop = pl.zeros([NSTEPS]); bl_prop = pl.zeros([NSTEPS])
pm_prop = pl.zeros([NSTEPS]); str_prop = pl.zeros([NSTEPS])

# take absolute values of u, v and normals:
u0 = pl.absolute(u0[0]); v0 = pl.absolute(v0[0])
normals = pl.absolute(normals)

for step in range(NSTEPS): # loop over timesteps
    # find where all trajectories have last been in the stratosphere:
    sl = StratStep(theta[0,:,:step+1],PV[0,:,:step+1],height[0,:,:step+1])
    # get the flux proportion of each origin type at each timestep
    # and the array with the origin label, timestep and co-ordinates
    blprop, pmprop, strprop, originsteps = OrOrder(blz[0,:,:step+1],height[0,:,:step+1],
                                               stepno[:,:step+1],lat[0,:,:step+1],
                                                lon[0,:,:step+1],equiv[0,:,:step+1],
                                                q[0,:,:step+1],sl,lsm,eralon,eralat)

    q_or = OriginFlux(originsteps,q[0]) # extract q for all traj with origin
    # extract q for all traj with CAT I, CAT II or strat origin
    q_bl = VarExtract(originsteps,q[0],1); q_pm = VarExtract(originsteps,q[0],2)
    q_st = VarExtract(originsteps,q[0],-1)
    
    # Calculate vertically integrated moisture fluxes at each release point
    # for all trajectories with origin, CAT I, CAT II & strat origins
    stepflx = FluxCalc(rlslabs,u0,v0,sp_tj[0,0],etalevs,normals,q_or[:,:],
                                             pres[0]) # ALL
    blflx = FluxCalc(rlslabs,u0,v0,sp_tj[0,0],etalevs,normals,q_bl[:,:],
                                             pres[0]) # CAT I
    pmflx = FluxCalc(rlslabs,u0,v0,sp_tj[0,0],etalevs,normals,q_pm[:,:],
                                             pres[0]) # CAT II
    strflx = FluxCalc(rlslabs,u0,v0,sp_tj[0,0],etalevs,normals,q_st[:,:],
                                              pres[0]) # strat
   
    # Calculate the proportion of the magnitude of the flux explained by
    # all trajectories with origin, CAT I, CAT II and strat origins
    flux_prop[step] = ((pl.sum(pl.absolute(stepflx)))/(10**9))/mag[0] # ALL
    bl_prop[step] = ((pl.sum(pl.absolute(blflx)))/(10**9))/mag[0] # CAT I
    pm_prop[step] = ((pl.sum(pl.absolute(pmflx)))/(10**9))/mag[0] # CAT II
    str_prop[step] = ((pl.sum(pl.absolute(strflx)))/(10**9))/mag[0] # strat

# turn these into percentages
flux_prop = flux_prop*100; bl_prop = bl_prop*100
pm_prop = pm_prop*100; str_prop = str_prop*100