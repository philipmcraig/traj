# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 17:23:02 2017

@author: np838619
"""

NPART = q.shape[1] # no. of trajectories
NRLS = q.shape[0] # no. of releases
NREL = rlslabs.shape[0] # no. of pts
shedflx_sea = pl.zeros([NRLS,NREL]); shedflx_lnd = pl.zeros_like(shedflx_sea)

catchcount = atlcount
for rls in range(NRLS): # loop over releases
    q_land = pl.zeros_like(q[rls]); q_sea = pl.zeros_like(q[rls])
    for traj in range(NPART):
#        tj = newlabs[rls,traj,0]; ts = catchcount[rls,traj,2]
        if newlabs[rls,traj,0] != 1.:
            continue
        pt = (newlabs[rls,traj,-1],newlabs[rls,traj,-2])
        ts = newlabs[rls,traj,1]#; print traj
        #if ts.dtype == 'float64':
        #    continue
        #p2 = (lon[rls,traj,ts+1],lat[rls,traj,ts+1])
#            pt = (catchcount[rls][trajno,-1],catchcount[rls][trajno,-2])
#            if ts < q[rls].shape[-1]-1:
#                    p2 = (lon[rls,tj,ts+1],lat[rls,tj,ts+1])
#            elif ts == q_tj[rls].shape[-1]-1:
#                    p2 = (lon[rls,tj,ts],lat[rls,tj,ts])
            #p2 = (lon[rls,traj,ts+1],lat[rls,traj,ts+1])
        x = BilinInterp(pt,eralon,eralat,lsm)
            #x2 = BilinInterp(p2,eralon,eralat,lsm)
        if x <= 0.5:
            q_sea[traj,:] = q[rls,traj,:]
           # elif x2 <= 0.5:
                #print tj
            #    q_sea[traj,:] = q[rls,traj,:]
        else:
            print traj, x#, x2
            q_land[traj,:] = q[rls,traj,:]
    shedflx_sea[rls] = FluxPartition(catchcount[rls],rlslabs,normals,etalevs,
                            q_sea,sp_tj[rls,0],u0[rls],v0[rls],pres[rls])
    shedflx_lnd[rls] = FluxPartition(catchcount[rls],rlslabs,normals,etalevs,
                            q_land,sp_tj[rls,0],u0[rls],v0[rls],pres[rls])