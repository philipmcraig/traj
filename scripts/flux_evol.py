# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 10:45:24 2017

@author: np838619
"""

def VarExtract(orlabs,var,label):
    """
    """
    var_or = pl.zeros_like(var)
    NPART = orlabs.shape[0]
    for traj in range(NPART):
        if orlabs[traj,0] == label:
            var_or[traj] = var[traj]
    
    return var_or

NREL = rlslabs.shape[0]
flux_prop = pl.zeros([NSTEPS]); bl_prop = pl.zeros([NSTEPS]); str_prop = pl.zeros([NSTEPS])
pm_prop = pl.zeros([NSTEPS])
B = pl.zeros([NSTEPS]); S = pl.zeros_like(B); O = pl.zeros_like(B); Q = pl.zeros([NSTEPS])
ps = pl.reshape(pres[0],(NCLUST,NREL,NSTEPS))
mag = pl.sum(pl.absolute(shedflx))/(10**9)
um = pl.absolute(u0[0]); vm = pl.absolute(v0[0]); nm = pl.absolute(normals)

magflx = FluxCalc(rlslabs,um,vm,sp_tj[0,0],etalevs,nm,q[0],pres[0],eralon,eralat)
mag = pl.sum(magflx)/(10**9)

for step in range(NSTEPS):
#    blcount,blprop,blsteps = BoundaryLayer(blz[0,:,:step+1],height[0,:,:step+1],
#                                           stepno[:,:step+1],lat[0,:,:step+1],lon[0,:,:step+1])
#    pmcount,pmprop,pmsteps = PMcalcs(equiv[0,:,:step+1],q[0,:,:step+1],stepno[:,:step+1],
#                                     lat[0,:,:step+1],lon[0,:,:step+1])
    sl = StratStep(theta[0,:,:step+1],PV[0,:,:step+1],height[0,:,:step+1])
    #ststeps = StratOrigin(sl,lat[0,:,:step+1],lon[0,:,:step+1])
#    originsteps = OriginCheck(blsteps,pmsteps,ststeps,NPART)
    blprop, pmprop, strprop, originsteps = OrOrder(blz[0,:,:step+1],height[0,:,:step+1],
                                                   stepno[:,:step+1],lat[0,:,:step+1],
                                                    lon[0,:,:step+1],equiv[0,:,:step+1],
                                                    q[0,:,:step+1],sl)
#    X = []
#    if pmsteps.shape[0] > 0.:
#        for i in range(pmsteps.shape[0]):
#            if pmsteps[i,0] not in blsteps[:,0]:
#                X.append(i)
#        if len(X) > 0:
#            Y = pl.zeros([len(X),4])
#            for i in range(len(X)):
#                Y[i] = pmsteps[X[i]]
#        else:
#            Y = pl.zeros([1,4]); Y[:] = pl.float64('nan')
#    else:
#        Y = pl.zeros([1,4]); Y[:] = pl.float64('nan')
    q_or = OriginFlux(originsteps,q[0])
#    q_bl = BLorigin_flux(blsteps,q[0],NPART); q_pm = BLorigin_flux(Y,q[0],NPART)
#    q_st = StratVar(ststeps,q[0])
    q_bl = VarExtract(originsteps,q[0],1); q_pm = VarExtract(originsteps,q[0],2)
    q_st = VarExtract(originsteps,q[0],-1)
    stepflx = FluxCalc(rlslabs,um,vm,sp_tj[0,0],etalevs,nm,q_or[:,:],pres[0],eralon,eralat)
    blflx = FluxCalc(rlslabs,um,vm,sp_tj[0,0],etalevs,nm,q_bl[:,:],pres[0],eralon,eralat)
    pmflx = FluxCalc(rlslabs,um,vm,sp_tj[0,0],etalevs,nm,q_pm[:,:],pres[0],eralon,eralat)
    strflx = FluxCalc(rlslabs,um,vm,sp_tj[0,0],etalevs,nm,q_st[:,:],pres[0],eralon,eralat)
    flux_prop[step] = ((pl.sum(pl.absolute(stepflx)))/(10**9))/mag
    bl_prop[step] = ((pl.sum(pl.absolute(blflx)))/(10**9))/mag
    pm_prop[step] = ((pl.sum(pl.absolute(pmflx)))/(10**9))/mag
    str_prop[step] = ((pl.sum(pl.absolute(strflx)))/(10**9))/mag
    
#    qx = q_or - q_bl
#    q_step = pl.zeros([NREL]); qs = pl.reshape(qx,(NCLUST,NREL,NSTEPS))
#    for pt in range(NREL):
#        sp_rel = BilinInterp(rlslabs[pt,:2],eralon,eralat,surfp)
#        q_step[pt] = VIeta(qs[:,pt,0],ps[:,pt,0],sp_rel,etalevs)
#    Q[step] = pl.sum(q_step)
    #O[step] = pl.sum(pl.sum(q_or[:,0]))
    #B[step] = pl.sum(pl.sum(q_bl[:,0]))
    #S[step] = pl.sum(pl.absolute(stepflx))/(10**9)#pl.sum(pl.sum(q_st[:,0]))

#pl.plot(Q)
pl.plot(flux_prop*100)
pl.plot(bl_prop*100)
pl.plot(str_prop*100)
pl.plot(pm_prop*100)
pl.ylim(0,100)