# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 12:30:14 2017

@author: np838619
"""

import pylab as pl

NPART = blz[0].shape[0]
blcount = 0.; blorigin = []
for traj in range(NPART):
    # find where BL height exceeds particle height
    bl_inds = pl.where(blz[0,traj,:]>height[0,traj,:])
    if bl_inds[0].size > 0.: # if bl_inds size >0 then particle was in BL
        S = pl.zeros([bl_inds[0].size])
        for i in range(bl_inds[0].size):
            S[i] = BilinInterp((lon[0,traj,bl_inds[0][i]],lat[0,traj,bl_inds[0][i]]),eralon,eralat,lsm)
        sea_inds = pl.where(S<=0.5)
        if sea_inds[0].size > 0.:
            blcount = blcount + 1.
            #if bl_inds[0][sea_inds[0][0]] == 0.:
                #print traj
            blorigin.append((traj,stepno[traj,bl_inds[0][sea_inds[0][0]]],
                 lat[0,traj,bl_inds[0][sea_inds[0][0]]],lon[0,traj,bl_inds[0][sea_inds[0][0]]]))
            #else:
            #    blorigin.append((traj,stepno[traj,bl_inds[0][sea_inds[0][0]]-1],
            #         lat[0,traj,bl_inds[0][sea_inds[0][0]]-1],lon[0,traj,bl_inds[0][sea_inds[0][0]]-1]))
        elif sea_inds[0].size == 0.:
            blcount = blcount
        # store traj no. & origin index, subtract 1 from index last in BL
        # as particle has still been in BL after this point
        #elif
        #    blorigin.append((traj,step[traj,bl_inds[0][0]],#-1
        #                 lat[traj,bl_inds[0][0]],lon[traj,bl_inds[0][0]]))
        #elif sea_inds[0].size == 0.:    
        #    if bl_inds[0][0] == 0.:
        #        blorigin.append((traj,stepno[traj,bl_inds[0][0]],#-1
        #                 lat[0,traj,bl_inds[0][0]],lon[0,traj,bl_inds[0][0]]))
        #    else:
        #        blorigin.append((traj,stepno[traj,bl_inds[0][0]-1],
        #                 lat[0,traj,bl_inds[0][0]-1],lon[0,traj,bl_inds[0][0]-1]))
    elif bl_inds[0].size == 0.:
        blcount = blcount

blprop = (blcount/NPART)*100
blorigin = pl.asarray(blorigin)