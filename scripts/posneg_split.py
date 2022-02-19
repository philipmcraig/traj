# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 15:13:54 2017

@author: np838619
"""

def FluxCalc_split(rlslabs,u0,v0,surfp_tj,etalevs,normals,q_traj,pres_traj):
    """
    """
    NPART = q_traj.shape[0] # get no. of trajectories from q_traj shape
    NREL = rlslabs.shape[0] # get no. of release points from rlslabs shape
    NCLUST = int(NPART/NREL) #  get no. of levels
    NSTEPS = q_traj.shape[1] # get no. of timesteps from q_traj shape
    
    # change shape of q & pressure arrays: no. levels X no. points X no. timesteps
    ps = pl.reshape(pres_traj,(NCLUST,NREL,NSTEPS))
    qs = pl.reshape(q_traj,(NCLUST,NREL,NSTEPS))
    
    
    # get zonal and meridional moisture fluxes on all levels at release points:
    uq = u0*qs[:,:,0]; vq = v0*qs[:,:,0]
    
    dl = pl.zeros([NREL]) # empty arrays for distance weights
    for pt in range(1,NREL-1): # exclude first & last release points
        # use Haversine formula to calculate dl
        dl[pt] = 0.5*(Haversine(rlslabs[pt-1,:2],rlslabs[pt,:2]) + 
                        Haversine(rlslabs[pt,:2],rlslabs[pt+1,:2]))
    # use first 2 points & last 2 points for first & last dl
    dl[0] = Haversine(rlslabs[0,:2],rlslabs[1,:2])
    dl[-1] = Haversine(rlslabs[-2,:2],rlslabs[-1,:2])
    
    fluxes_uv = pl.zeros([NREL,2]) # empty array for vertically integrated fluxes
    shedflx1 = pl.zeros([NREL]) # empty array for positive fluxes
    shedflx2 = pl.zeros([NREL]) # empty array for negative fluxes
    
    for pt in range(NREL): # loop over release points
        # Integrate on eta levels
        fluxes_uv[pt,0] = VIeta(uq[:,pt],ps[:,pt,0],surfp_tj[pt],etalevs)
        fluxes_uv[pt,1] = VIeta(vq[:,pt],ps[:,pt,0],surfp_tj[pt],etalevs)
        seglab = rlslabs[pt,-1] - 1 # get label of segment
        # use dot product of moisture fluxes & normal vector, multiply by dl
        flux = pl.dot(fluxes_uv[pt],normals[seglab])*dl[pt]
        if flux > 0.:
            shedflx1[pt] = flux
        elif flux < 0.:
            shedflx2[pt] = flux
    
    return shedflx1, shedflx2

def FluxEvol_split(blz,height,stepno,lat,lon,equiv,q,theta,PV,rlslabs,u0,v0,normals,
            pres,etalevs,sp_tj,lsm,eralon,eralat):
    """
    """
     # get no. of trajectories & timesteps from shape of blz array:
    NPART = blz.shape[0]; NSTEPS = blz.shape[1]
    NREL = rlslabs.shape[0] # get no. of release points from shape of rlslabs array
    NCLUST = int(NPART/NREL) # no. of levels is no. trajectories / no. points
    # make empty arrays for flux proportions, size no. of steps:
    flux_prop = pl.zeros([NSTEPS]); bl_prop = pl.zeros([NSTEPS])
    pm_prop = pl.zeros([NSTEPS]); str_prop = pl.zeros([NSTEPS])
    stepflx1 = pl.zeros([NSTEPS]); stepflx2 = pl.zeros([NSTEPS])
    
    # take absolute values of u, v and normals:
    #u0 = pl.absolute(u0); v0 = pl.absolute(v0)
    #normals = pl.absolute(normals)
    
    for step in range(NSTEPS): # loop over timesteps
        # find where all trajectories have last been in the stratosphere:
        sl = StratStep(theta[:,:step+1],PV[:,:step+1],height[:,:step+1])
        # get the flux proportion of each origin type at each timestep
        # and the array with the origin label, timestep and co-ordinates
        blprop, pmprop, strprop, originsteps = OrOrder(blz[:,:step+1],height[:,:step+1],
                                                   stepno[:,:step+1],lat[:,:step+1],
                                                    lon[:,:step+1],equiv[:,:step+1],
                                                    q[:,:step+1],sl,lsm,eralon,eralat)
    
        q_or = OriginFlux(originsteps,q) # extract q for all traj with origin
        # extract q for all traj with CAT I, CAT II or strat origin
        #q_bl = VarExtract(originsteps,q,1); q_pm = VarExtract(originsteps,q,2)
        #q_st = VarExtract(originsteps,q,-1)
        
        # Calculate vertically integrated moisture fluxes at each release point
        # for all trajectories with origin, CAT I, CAT II & strat origins
        flxpos, flxneg = FluxCalc_split(rlslabs,u0,v0,sp_tj,etalevs,normals,q_or[:,:],
                                                 pres) # ALL
        stepflx1[step] = pl.sum(flxpos)
        stepflx2[step] = pl.sum(flxneg)
        #blflx = FluxCalc(rlslabs,u0,v0,sp_tj,etalevs,normals,q_bl[:,:],
        #                                         pres) # CAT I
        #pmflx = FluxCalc(rlslabs,u0,v0,sp_tj,etalevs,normals,q_pm[:,:],
        #                                         pres) # CAT II
        #strflx = FluxCalc(rlslabs,u0,v0,sp_tj,etalevs,normals,q_st[:,:],
         #                                         pres) # strat
       
        # Calculate the proportion of the magnitude of the flux explained by
        # all trajectories with origin, CAT I, CAT II and strat origins
        #flux_prop[step] = (pl.absolute(stepflx))/(10**9))#/mag # ALL
        #bl_prop[step] = ((pl.sum(pl.absolute(blflx)))/(10**9))/mag # CAT I
        #pm_prop[step] = ((pl.sum(pl.absolute(pmflx)))/(10**9))/mag # CAT II
        #str_prop[step] = ((pl.sum(pl.absolute(strflx)))/(10**9))/mag # strat
    
    # turn these into percentages
    #flux_prop = flux_prop*100#; bl_prop = bl_prop*100
    #pm_prop = pm_prop*100; str_prop = str_prop*100
   
    return stepflx1, stepflx2#, bl_prop, pm_prop, str_prop

#a,b = FluxCalc_split(rlslabs,u0[0],v0[0],sp_tj[0,0],etalevs,normals,q[0],pres[0])

c = pl.zeros([len(trajfiles),NSTEPS]); d = pl.zeros([len(trajfiles),NSTEPS])

for rls in range(len(trajfiles)):
    c[rls],d[rls] = FluxEvol_split(blz[rls],height[rls],stepno,lat[rls],lon[rls],
                equiv[rls],q[rls],theta[rls],PV[rls],rlslabs,u0[rls],v0[rls],
                normals,pres[rls],etalevs,sp_tj[rls,0],lsm,eralon,eralat)

pos_mn = pl.mean(c,axis=0); neg_mn = pl.mean(d,axis=0)
H = pl.array([pos_mn,neg_mn])

#steps = pl.linspace(0,120,121)*6/24
#pl.plot(steps,pos_mn/(10**9),color='b',lw=2.,label='+ve')
#pl.plot(steps,neg_mn/(10**9),color='g',lw=2.,label='-ve')
#pl.plot(steps,(pos_mn+neg_mn)/(10**9),color='r',lw=2.,label='net')
#pl.legend(loc=0)
#pl.ylabel('Sv'); pl.xlabel('days')

f = open('/home/np838619/Trajectory/fluxevol/'+shed+'/posneg_Jul10.csv','w')
pl.savetxt(f,H.T,fmt='%9.5f',delimiter=' ')
f.close()