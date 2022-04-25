# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 10:39:37 2017

@author: np838619
"""

def FluxCalc(rlslabs,u_era,v_era,surfp_era,etalevs,normals,q_traj,pres_traj,eralon,eralat):
    """
    """
    NPART = q_traj.shape[0]; NREL = rlslabs.shape[0]; NCLUST = NPART/NREL
    NSTEPS = q_traj.shape[1]
    sp_rel = pl.zeros([NREL])
    u_rel = pl.zeros([NCLUST,NREL]); v_rel = pl.zeros_like(u_rel)
    
    ps = pl.reshape(pres_traj,(NCLUST,NREL,NSTEPS))
    qs = pl.reshape(q_traj,(NCLUST,NREL,NSTEPS))
    
    for pt in range(NREL):
        sp_rel[pt] = BilinInterp(rlslabs[pt,:2],eralon,eralat,surfp_era)
        for lev in range(NCLUST):
            L = NearestIndex(erapres,ps[lev,pt,0])
            u_rel[lev,pt] = BilinInterp(rlslabs[pt,:2],eralon,eralat,u_era[L])
            v_rel[lev,pt] = BilinInterp(rlslabs[pt,:2],eralon,eralat,v_era[L])
    
    uq = u_rel*qs[:,:,0]; vq = v_rel*qs[:,:,0]
    
    dl = pl.zeros([NREL])
    for pt in range(1,NREL-1):
        dl[pt] = 0.5*(Haversine(rlslabs[pt-1,:2],rlslabs[pt,:2]) + 
                        Haversine(rlslabs[pt,:2],rlslabs[pt+1,:2]))
    dl[0] = Haversine(rlslabs[0,:2],rlslabs[1,:2])
    dl[-1] = Haversine(rlslabs[-2,:2],rlslabs[-1,:2])
    
    fluxes_uv = pl.zeros([NREL,2])
    shedflx = pl.zeros([NREL])
    
    for pt in range(NREL):
        fluxes_uv[pt,0] = VIeta(uq[:,pt],ps[:,pt,0],sp_rel[pt],etalevs)
        fluxes_uv[pt,1] = VIeta(vq[:,pt],ps[:,pt,0],sp_rel[pt],etalevs)
        seglab = rlslabs[pt,-1] - 1
        shedflx[pt] = pl.dot(fluxes_uv[pt],normals[seglab])*dl[pt]
    
    return shedflx