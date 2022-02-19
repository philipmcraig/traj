# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 14:02:37 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
from netCDF4 import Dataset
from matplotlib.lines import Line2D

def P2Eta(psurf):
    """
    """
    trajdir = '/home/np838619/Trajectory/'
    p0 = 10**5
    psurf = psurf
    coeffs = pl.genfromtxt(trajdir + 'coeffs.txt')
    a = coeffs[:,1]; b = coeffs[:,2]
    
    pL = a + b*(psurf/p0)
    
    return pL

def EtaToP(a,b,psurf):
    """
    """
    p0 = 10**5
    A = pl.zeros([len(a)-1]); B = pl.zeros_like(A)
    for i in range(len(A)):
        A[i] = 0.5*(a[i+1]+a[i])
        B[i] = 0.5*(b[i+1]+b[i])
    
    P = A*p0 + B*psurf
    
    return P

def LocalCartesian(coords):
    """Function to convert latitude-longitude co-ordinates into local
	Cartesian co-ordinates.

	Args:
		co-ords (array): Nx2 array of co-ordinates in degrees

	Returns:
		loccar (array): Nx2 array of local Cartesian co-ordinates
    """
    coords = pl.radians(coords) # convert from degrees to radians
    loccar = pl.zeros_like(coords) # creat empty array same size as coords
    R = 6.37*10**6 # radius of the Earth
	# x = R*cos(lat)*lon
    loccar[:,0] = R*pl.sin(coords[:,1])*pl.cos(coords[:,0])#coords[:,0]#
    loccar[:,1] = R*pl.sin(coords[:,1])*pl.sin(coords[:,0]) #coords[:,1] y = R*lat
    
    return loccar

def TrajSegLabel(loc):
    """Function to assign each trajectory release point a label referring to which
	segment of the watershed it is from.

	Args:
		loc (string): stem of continental watershed e.g. NCA = North/Central America
	
	Returns:
             seglab (array): segment labels (integers) of each release point
             rlspts (array): trajectory release points
    """
    sheddir = '/home/np838619/Watershed/shed_defs/'
    
    rlspts = pl.genfromtxt(sheddir+loc+'_traj_release.txt',skip_header=5)
    
    seglab = [1] # first release point has to be on the first segment
    count = 1
    
    for rls in range(1,rlspts.shape[0]):
        if rlspts[rls,0] == rlspts[rls-1,0] and rlspts[rls,1] == rlspts[rls-1,1]:
            count = count + 1
            seglab.append(count)
        else:
            count = count
            seglab.append(count)
    seglab = pl.asarray(seglab)
    
    return seglab#, rlspts

def MidPts(rlspts):
    """
    """
    #rlspts[:,0] = rlspts[:,0] + 360.
    #loccar = LocalCartesian(rlspts)
    #m = Basemap(llcrnrlon=-180,llcrnrlat=-80.,urcrnrlon=180.,urcrnrlat=80.,
                #lat_ts=20.,resolution='l',projection='mill',suppress_ticks=True)
    #pl.zeros_like(rlspts)
    #loccar[:,0], loccar[:,1] = m(rlspts[:,0],rlspts[:,1])
    
    dl = pl.zeros([rlspts.shape[0]-1])
    for pt in range(rlspts.shape[0]-1):
        dl[pt] = Haversine(rlspts[pt],rlspts[pt+1])
    
    midpt = pl.zeros([rlspts.shape[0]])
    
    for pt in range(1,rlspts.shape[0]-1):
        midpt[pt] = 0.5*(dl[pt]+dl[pt-1])
    midpt[0] = dl[0]; midpt[-1] = dl[-1]
    
    return midpt

def RR2(labels,loc):
    rlspts = pl.genfromtxt(sheddir+loc+'_traj_release.txt',skip_header=5)
    repeat = []
    for r in range(1,rlspts.shape[0]):
        if rlspts[r,0] == rlspts[r-1,0] and rlspts[r,1] == rlspts[r-1,1]:
            repeat.append(r)
    
    L2 = []#; RP2 = [] 
    for r in range(rlspts.shape[0]):
        if r not in repeat:
        #RP2.append(interp_pts[r])
            L2.append(labels[r])
    labs = pl.asarray(L2)# rlspts = pl.asarray(RP2)
    
    return labs#,rlspts, labs

def MidFlux(rlspts,lon,lat,zon,mer):
    """
    """
    flux_uv = pl.zeros_like(rlspts); 
    for i in range(len(flux_uv)):
        a = NearestIndex(lon,rlspts[i,0]); b = NearestIndex(lat,rlspts[i,1])
        if lon[a] == rlspts[i,0]:
            flux_uv[i,0] = zon[b,a]; flux_uv[i,1] = mer[b,a]
        else:
            flux_uv[i,0] = BilinInterp(rlspts[i],lon,lat,zon)
            flux_uv[i,1] = BilinInterp(rlspts[i],lon,lat,mer)
    
    midflux = pl.zeros_like(rlspts)
    for f in range(1,midflux.shape[0]-1):
        midflux[f,0] = (flux_uv[f,0]+flux_uv[f-1,0])/2
        midflux[f,1] = (flux_uv[f,1]+flux_uv[f-1,1])/2
    midflux[0,:] = flux_uv[0,:]
    midflux[-1,:] = flux_uv[-1,:]
    
    return midflux

def NormalFlux(midflux,labs,nhat):
    """
    """
    FdotN = pl.zeros([midflux.shape[0]])
    for i in range(FdotN.shape[0]):
        segno = labs[i,0] - 1
        FdotN[i] = pl.dot(midflux[i],nhat[segno])
        #if segno not in labs:
        #    pass
        #else:
            #print i, segno
            #FdotN[i] = pl.dot(midflux[i],nhat[segno])
    
    return FdotN

def ShedFluxes(sheddir,loc,lon,lat,zon,mer):
    """
    """
    endpts = pl.genfromtxt(sheddir+loc+'_clicks.txt',skip_header=5)
    labs = TrajSegLabel(loc)
    l2 = pl.zeros([labs[0].size,3]); l2[:,0] = labs[0]; l2[:,1:] = labs[1]
    labs = l2.copy()
    nhat = NormalVector(sheddir,loc)
    labs = RR2(labs,loc)
    rlspts = pl.genfromtxt(sheddir+loc+'_traj_release_new.txt',skip_header=5)
    rlspts = Add360(rlspts)
    midpt = MidPts(rlspts)
    midflux = MidFlux(rlspts,lon,lat,zon,mer)
    FdotN = NormalFlux(midflux,labs,nhat)
    fluxes = FdotN[:]#*midpt[:]/(10**9)
    
    return fluxes

def FC2(rlslabs,u0,v0,surfp_tj,etalevs,normals,q_traj,pres_traj):
    """Function to calculate the vertically intehrated moisture flux at each 
    trajectory release point.
            (1/g)*|(uq,vq).n*dl*deta
    
    Args:
        rlslabs (array): longitude/latitude of trajectory release point and label
                        indicating which segment of watershed the point is on
        u0 (array): zonal velocity (m/s) at release point
        v0 (array): meridional velocity (m/s) at release point
        etalevs (array): model levels which trajectories are released from
        normals (array): normal vectors for each segment of watershed
        q_traj (array): specific humidity (kg/kg) along trajectories
        pres_traj (array): pressure (hPa) along trajectories
    
    Returns:
        shedflx (array): vertically integrated moisture flux at each trajectory
                        release point (kg/s)
    """
    global Haversine; global VIeta
    NPART = q_traj.shape[0] # get no. of trajectories from q_traj shape
    NREL = rlslabs.shape[0] # get no. of release points from rlslabs shape
    NCLUST = int(NPART/NREL) #  get no. of levels
    NSTEPS = q_traj.shape[1] # get no. of timesteps from q_traj shape
    
    # change shape of q & pressure arrays: no. levels X no. points X no. timesteps
    ps = pres_traj.copy()#pl.reshape(pres_traj,(NCLUST,NREL,NSTEPS))
    qs = q_traj.copy()#pl.reshape(q_traj,(NCLUST,NREL,NSTEPS))
    
    
    # get zonal and meridional moisture fluxes on all levels at release points:
    uq = u0*qs[:,:]; vq = v0*qs[:,:]
    
    dl = pl.zeros([NREL]) # empty arrays for distance weights
    for pt in range(1,NREL-1): # exclude first & last release points
        # use Haversine formula to calculate dl
        dl[pt] = 0.5*(Haversine(rlslabs[pt-1,:2],rlslabs[pt,:2]) + 
                        Haversine(rlslabs[pt,:2],rlslabs[pt+1,:2]))
    # use first 2 points & last 2 points for first & last dl
    dl[0] = Haversine(rlslabs[0,:2],rlslabs[1,:2])
    dl[-1] = Haversine(rlslabs[-2,:2],rlslabs[-1,:2])
    
    fluxes_uv = pl.zeros([NREL,2]) # empty array for vertically integrated fluxes
    shedflx = pl.zeros([NREL]) # empty array for total fluxes
    
    for pt in range(NREL): # loop over release points
        # Integrate on eta levels
        fluxes_uv[pt,0] = VIeta(uq[:,pt],ps[:,pt],surfp_tj[pt],etalevs)
        fluxes_uv[pt,1] = VIeta(vq[:,pt],ps[:,pt],surfp_tj[pt],etalevs)
        seglab = rlslabs[pt,-1] - 1 # get label of segment
        # use dot product of moisture fluxes & normal vector, multiply by dl
        shedflx[pt] = pl.dot(fluxes_uv[pt],normals[seglab])#*dl[pt]
    
    return shedflx, dl

pl.close('all')

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())
exec(open('/home/np838619/Trajectory/trajfuncs.py').read())

clusdir = '/glusterfs/scenario/users/np838619/traj/'
trajdir = '/home/np838619/Trajectory/'
eradir = '/glusterfs/scenario/users/np838619/ERA/'
sheddir = '/home/np838619/Watershed/shed_defs/'

trajlevs = pl.linspace(0.95,0.15,17)

#coeffs = pl.genfromtxt(trajdir + 'coeffs.txt')
#a = coeffs[:,1]; b = coeffs[:,2]
filenames = PrintFiles(clusdir+'2014/amr/','tj_')
loc = 'Am'; shed = 'amr'
year = filenames[0][13:17]

rlspts = pl.genfromtxt(sheddir+loc+'_traj_release_new.txt',skip_header=5)
rlspts = Add360(rlspts)

ncfile = Dataset(clusdir+'2014/amr/'+filenames[0],'r')
ps_tj = ncfile.variables['surface pressure t=0'][:]
q_tj = ncfile.variables['specific humidty t=0'][:]
pres_tj = ncfile.variables['pressure t=0'][:]
pv_tj = ncfile.variables['potential vorticity t=0'][:]
tmp_tj = ncfile.variables['temperature t=0'][:]
u_tj = ncfile.variables['u wind t=0'][:]
v_tj = ncfile.variables['v wind t=0'][:]
blz_tj = ncfile.variables['boundary layer height t=0'][:]
ncfile.close()

qu_tj = u_tj*q_tj; qv_tj = v_tj*q_tj

erafile1 = Dataset(eradir+'ggap201407011200.nc','r')
eralat = erafile1.variables['latitude'][:]
eralon = erafile1.variables['longitude'][:]
erapres = erafile1.variables['p'][:]
eraq = erafile1.variables['Q'][:]
erapv = erafile1.variables['PV'][:]
eratmp = erafile1.variables['T'][:]
erau = erafile1.variables['U'][:]
erav = erafile1.variables['V'][:]
erafile1.close()
eraq = pl.squeeze(eraq); erapv = pl.squeeze(erapv*(10**6))
eratmp = pl.squeeze(eratmp); erau = pl.squeeze(erau)
erav = pl.squeeze(erav)

era_qu = eraq*erau; era_qv = eraq*erav

erafile2 = Dataset(eradir+'ggas201407011200.nc','r')
eraps = erafile2.variables['SP'][:]
erafile2.close()
eraps = pl.squeeze(eraps)

erafile3 = Dataset(eradir+'ggfs201407010012.nc','r')
erablz = erafile3.variables['BLH'][:]
erafile3.close()
erablz = pl.squeeze(erablz)

sp_era = pl.zeros([rlspts.shape[0]]); blz_era = pl.zeros_like(sp_era)
q_era = pl.zeros([erapres.size,rlspts.shape[0]])
pv_era = pl.zeros([erapres.size,rlspts.shape[0]])
tmp_era = pl.zeros_like(pv_era); u_era = pl.zeros_like(pv_era)
v_era = pl.zeros_like(u_era)
#qu_era = pl.zeros_like(u_era); qv_era = pl.zeros_like(v_era)

for pt in range(rlspts.shape[0]):
    sp_era[pt] =  BilinInterp(rlspts[pt],eralon,eralat,eraps)
    blz_era[pt] = BilinInterp(rlspts[pt],eralon,eralat,erablz)
    for lev in range(erapres.size):
        q_era[lev,pt] = BilinInterp(rlspts[pt],eralon,eralat,eraq[lev])
        pv_era[lev,pt] = BilinInterp(rlspts[pt],eralon,eralat,erapv[lev])
        tmp_era[lev,pt] = BilinInterp(rlspts[pt],eralon,eralat,eratmp[lev])
        u_era[lev,pt] = BilinInterp(rlspts[pt],eralon,eralat,erau[lev])
        v_era[lev,pt] = BilinInterp(rlspts[pt],eralon,eralat,erav[lev])
        #qu_era[lev,pt] = BilinInterp(rlspts[pt],eralon,eralat,era_qu[lev])
        #qv_era[lev,pt] = BilinInterp(rlspts[pt],eralon,eralat,era_qv[lev])

fig, ax = pl.subplots(2,3,figsize=(14,8))
P = 45
print rlspts[45,0], rlspts[45,1]

ax1 = pl.subplot(231)
ax1.plot(q_tj[:,P]*1000,pres_tj[:,P],color='b',lw=3,label='trajectory')
ax1.plot(q_era[:,P]*1000,erapres,color='r',lw=2,ls='--',label='ERA-Interim')
pl.ylim(1000,0); pl.title('(a) $q$',fontsize=18); ax1.grid(axis='both')
pl.xlabel('g kg$^{-1}$',fontsize=15); pl.ylabel('hPa',fontsize=15,labelpad=-5)
ax1.legend(loc=1)

ax2 = pl.subplot(232)
ax2.plot(tmp_tj[:,P],pres_tj[:,P],color='b',lw=3)
ax2.plot(tmp_era[:,P],erapres,color='r',lw=2,ls='--')
pl.ylim(1000,0); pl.title('(b) $T$',fontsize=18)
pl.xlabel('$K$',fontsize=15); ax2.grid(axis='both')
ax2.yaxis.set_major_formatter(pl.NullFormatter())

ax3 = pl.subplot(233)
ax3.plot(pv_tj[:,P],pres_tj[:,P],color='b',lw=3)
ax3.plot(pv_era[:,P],erapres,color='r',ls='--',lw=2)
pl.ylim(1000,0); pl.xlim(-1,1); pl.title('(c) PV',fontsize=18)
pl.xlabel('PVU',fontsize=15); ax3.grid(axis='both')
ax3.yaxis.set_major_formatter(pl.NullFormatter())

ax4 = pl.subplot(234)
ax4.plot(u_tj[:,P],pres_tj[:,P],color='b',lw=3)
ax4.plot(u_era[:,P],erapres,color='r',lw=2,ls='--')
ax4.plot(v_tj[:,P],pres_tj[:,P],color='g',lw=3)
ax4.plot(v_era[:,P],erapres,color='orange',lw=2,ls='--')
pl.ylim(1000,0); pl.title('(d) $u,\;v$',fontsize=18); pl.xlim(-40,40)
pl.xlabel('m s$^{-1}$',fontsize=15); ax4.grid(axis='both')
pl.ylabel('hPa',fontsize=15,labelpad=-5); pl.xlim(-15,15)
l1 = ax4.legend([Line2D([0],[0],color='b',lw=3),Line2D([0],[0],color='r',ls='--',lw=2)],
            [r'$u_{traj}$',r'$u_{ERA}$'],loc=3,handletextpad=0.05)
l2 = ax4.legend([Line2D([0],[0],color='g',lw=3),Line2D([0],[0],color='orange',ls='--',lw=2)],
            [r'$v_{traj}$',r'$v_{ERA}$'],loc=4,handletextpad=0.05)
pl.gca().add_artist(l1)

ax5 = pl.subplot(235)
ax5.plot(ps_tj[0,:],color='b',lw=3)
ax5.plot(sp_era/100,color='r',lw=2,ls='--')
pl.ylim(1000,0); pl.title('(e) $\\tilde{p}_{surf}$',fontsize=18)
ax5.grid(axis='both'); pl.xlim(0,rlspts.shape[0])
ax5.yaxis.set_major_formatter(pl.NullFormatter())
a = pl.arange(0,rlspts.shape[0],30); b = pl.around(rlspts[a,1],decimals=0)
ax5.set_xticks(a); ax5.set_xticklabels(b.astype('int'),fontsize=14)
pl.xlabel('latitude',fontsize=16)

ax6 = pl.subplot(236)
ax6.plot(blz_tj[0,:]/1000,color='b',lw=3)
ax6.plot(blz_era/1000,color='r',lw=2,ls='--')
pl.title('(f) $z_{BL}$',fontsize=18); pl.xlim(0,rlspts.shape[0])
pl.ylabel('km',fontsize=15); ax6.grid(axis='both')
a = pl.arange(0,rlspts.shape[0],30); b = pl.around(rlspts[a,1],decimals=0)
ax6.set_xticks(a); ax6.set_xticklabels(b.astype('int'),fontsize=14)
pl.xlabel('latitude',fontsize=16)

pl.tight_layout()
pl.savefig(clusdir+year+'/'+shed+'/datacomp_'+filenames[0][13:-4]+'.png')
#pl.show()

# vertically integrate qu_tj & qv_tj:
QU_tj = pl.zeros([rlspts.shape[0]]); QV_tj = pl.zeros_like(QU_tj)
for pt in range(rlspts.shape[0]):
    QU_tj[pt] = VIeta(qu_tj[:,pt],pres_tj[:,pt],ps_tj[0,pt],trajlevs)
    QV_tj[pt] = VIeta(qv_tj[:,pt],pres_tj[:,pt],ps_tj[0,pt],trajlevs)


labs, rlspts = TrajSegLabel(loc)
rlspts, labs = RemoveRepeated(rlspts,labs)
rlslabs = pl.zeros([len(labs),3])
rlslabs[:,:2] = rlspts[:,:]; rlslabs[:,2] = labs[:]
rlslabs[:,:2] = Add360(rlslabs[:,:2])
normals = NormalVector(sheddir,loc)
F1, dl = FC2(rlslabs,u_tj,v_tj,ps_tj[0],trajlevs,normals,q_tj,pres_tj)

# vertically integrate era_qu & era_qv:
QU_era = pl.zeros([eralat.size,eralon.size]); QV_era = pl.zeros_like(QU_era)
for i in range(eralat.size):
    for j in range(eralon.size):
        QU_era[i,j] = VIpres(era_qu[:,i,j],eraps[i,j],erapres)
        QV_era[i,j] = VIpres(era_qv[:,i,j],eraps[i,j],erapres)

F2 = ShedFluxes(sheddir,loc,eralon,eralat,QU_era,QV_era)

int1 = pl.sum(F1*dl)/(10**9); int2 = pl.sum(F2*dl)/(10**9)

fig = pl.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(F1,color='b',lw=3,label='trajectory')
ax.plot(F2,color='r',lw=2,ls='--',label='ERA-Interim')
ax.legend(loc=3); ax.grid(axis='y')
a = pl.arange(0,rlspts.shape[0],20); b = pl.around(rlspts[a,1],decimals=0)
ax.set_xticks(a); ax.set_xticklabels(b.astype('int'),fontsize=14)
pl.xlim(0,rlspts.shape[0])
pl.xlabel('latitude',fontsize=16)
pl.ylabel('kg m$^{-1}$ s$^{-1}$',fontsize=16,labelpad=-5)
ax.tick_params(axis='x',pad=7)
ax.tick_params(axis='y',labelsize=14)
ax.annotate('$\int\mathbf{Q}\cdot\mathbf{\hat{n}}_{traj}$= '+'{:.2f}'.format(int1)+' Sv',
            xy=(5,-250),fontsize=16,color='b')
ax.annotate('$\int\mathbf{Q}\cdot\mathbf{\hat{n}}_{ERA}$= '+'{:.2f}'.format(int2)+' Sv',
            xy=(5,-350),fontsize=16,color='r')
#pl.savefig(clusdir+year+'/'+shed+'/profcomp_'+filenames[0][13:-4]+'.png')
#pl.show()
