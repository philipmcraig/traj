# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 14:13:48 2017

@author: np838619
"""


NPART = q.shape[1] # no. of trajectories
NRLS = q.shape[0] # no. of releases
NREL = rlslabs.shape[0] # no. of pts
shedflx_BL = pl.zeros([NRLS,NREL])

catchcount = atlcount
for rls in range(NRLS):
    q_BL = pl.zeros_like(q[rls])
    for trajno in range(catchcount[rls].shape[0]):
        #pt = (newlabs[rls,traj,-1],newlabs[rls,traj,-2])
        tj = catchcount[rls][trajno,0]
        pt = (catchcount[rls][trajno,-1],catchcount[rls][trajno,-2])
        if catchcount[rls][trajno,1] == 1.:
            if catchcount[rls][trajno,2] == 0.:
                print trajno
           #x = BilinInterp(pt,eralon,eralat,lsm)
            #if x > 0.5:
                q_BL[tj,:] = q[rls,tj,:]
    shedflx_BL[rls] = FluxPartition(catchcount[rls],rlslabs,normals,etalevs,
                                q_BL,sp_tj[rls,0],u0[rls],v0[rls],pres[rls])
    #shedflx_BL[rls] = FluxCalc(rlslabs,u0[rls],v0[rls],sp_tj[rls,0],
    #                            etalevs,normals,q_BL,pres[rls])