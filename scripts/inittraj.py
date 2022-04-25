# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 14:17:02 2017

@author: np838619
""" 

for i in range(len(initfiles)):
    flist = ReadTxtFile(pingdir + initfiles[0])
    
    BASETIME = int(flist[0][-1])
    NPART = int(flist[2][-1])
    NATTR = int(flist[3][-1]) + 2
    NCLUST = int(flist[6][-1])
    pointers = list(itertools.chain(*flist[8:10])) # traj no where each level starts
    pointers = pl.asarray(pointers,dtype='int')
    data = len(flist[16])-2 + NATTR

    initdata = pl.zeros([NPART,data])
    
    for traj in range(NPART):
        initdata[traj] = flist[17+traj*5]
    
    lat = initdata[:,2]; lon = initdata[:,3]
    pres = initdata[:,4]; temp = initdata[:,5]
    PV = initdata[:,6]; q = initdata[:,7]
    height = initdata[:,8]; blz = initdata[:,9]
    u = initdata[:,10]; v = initdata[:,11]

NREL = NPART/NCLUST
q = pl.reshape(q,(NCLUST,NREL)); u = pl.reshape(u,(NCLUST,NREL))
v = pl.reshape(v,(NCLUST,NREL)); lat = pl.reshape(lat,(NCLUST,NREL))
pres = pl.reshape(pres,(NCLUST,NREL)); temp = pl.reshape(temp,(NCLUST,NREL))
PV = pl.reshape(PV,(NCLUST,NREL))

a = 6.37e6
U = -1*u*(a*pl.cos(pl.radians(lat)))
V = -1*v*a

nc1 = Dataset('/home/np838619/Watershed/ggas201307201200.nc','r')
surfp = nc1.variables['SP'][:]
eralon = nc1.variables['longitude'][:]; eralat = nc1.variables['latitude'][:]
nc1.variables

nc2 = Dataset('/home/np838619/Watershed/ggap201307201200.nc','r')
uera = nc2.variables['U'][:]; vera = nc2.variables['V'][:]
erapres = nc2.variables['p'][:]; qera = nc2.variables['Q'][:]
PVera = nc2.variables['PV'][:]; Tera = nc2.variables['T'][:]
nc2.close()

surfp = pl.squeeze(surfp); uera = pl.squeeze(uera); vera = pl.squeeze(vera)
qera = pl.squeeze(qera); PVera = pl.squeeze(PVera); Tera = pl.squeeze(Tera)


u_shed = pl.zeros([erapres.size,NREL]); v_shed = pl.zeros_like(u_shed)
q_shed = pl.zeros_like(u_shed); sp_shed = pl.zeros([NREL])
PV_shed = pl.zeros_like(u_shed); T_shed = pl.zeros_like(u_shed)

for pt in range(NREL):
    sp_shed[pt] = BilinInterp(rlslabs[pt,:2],eralon,eralat,surfp)
    for lev in range(erapres.size):
        u_shed[lev,pt] = BilinInterp(rlslabs[pt,:2],eralon,eralat,uera[lev])
        v_shed[lev,pt] = BilinInterp(rlslabs[pt,:2],eralon,eralat,vera[lev])
        q_shed[lev,pt] = BilinInterp(rlslabs[pt,:2],eralon,eralat,qera[lev])
        PV_shed[lev,pt] = BilinInterp(rlslabs[pt,:2],eralon,eralat,PVera[lev])
        T_shed[lev,pt] = BilinInterp(rlslabs[pt,:2],eralon,eralat,Tera[lev])
        