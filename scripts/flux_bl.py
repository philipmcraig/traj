# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 15:04:28 2017

@author: np838619
"""

stepflx = pl.zeros([NSTEPS])

#blcount,blprop,blsteps = BoundaryLayer(blz[0,:,:1],height[0,:,:1],stepno[:,:1],lat[0,:,:1],lon[0,:,:1])
#stepzero = blsteps
#q_bl = BLorigin_flux(stepzero,q[0],NPART)
#flx = FluxCalc(rlslabs,u,v,surfp,etalevs,normals,q_bl[:,:],pres[0],eralon,eralat,erapres)
#stepflx[0] = pl.sum(pl.absolute(flx))/(10**9)

for step in range(NSTEPS):
    blcount,blprop,blsteps = BoundaryLayer(blz[0,:,:step+1],height[0,:,:step+1],
                                           stepno[:,:step+1],lat[0,:,:step+1],lon[0,:,:step+1])
#    X = []
#    if step == 1.:
#        for i in range(blsteps.shape[0]):
#            if blsteps[i,0] not in stepzero[:,0]:
#                X.append(i)
#    elif step > 1.:
#        for i in range(blsteps.shape[0]):
#            if blsteps[i,0] not in prevstep[:,0]:
#                X.append(i)
#    if len(X) > 0:
#        Y = pl.zeros([len(X),4])
#        for i in range(len(X)):
#            Y[i] = blsteps[X[i]]
#    else:
#        Y = pl.zeros([1,4])
    
    q_bl = BLorigin_flux(blsteps,q[0],NPART)
    
    flx = FluxCalc(rlslabs,u,v,surfp,etalevs,normals,q_bl[:,:],pres[0],eralon,eralat,erapres)
    stepflx[step] = pl.sum(pl.absolute(flx))/(10**9)
    
    #prevstep = blsteps