# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 10:01:45 2016

@author: np838619
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.basemap import Basemap
from os import listdir
from os.path import isfile, join
import itertools
from netCDF4 import Dataset
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.path as mplPath
from scipy.stats import pearsonr 
import timeit
from trajfuncs import *
from kernelfuncs import *

start_time = timeit.default_timer()
exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())
pl.close('all')
#datadir = '/export/pinggwo/data-01/np838619/traj/Am30_out07/'
shed = 'amr'; year='2010'
#exec(open('/home/np838619/Trajectory/trajyear.py').read())
datadir = '/glusterfs/scenario/users/np838619/traj/'+year+'/'+shed+'/TD/'
#datadir = '/net/glusterfs/scenario/users/np838619/tjexp/'+shed+'/wint/TD/'

trajfiles = PrintFiles(datadir,'utraj-df')
trajfiles = pl.sort(trajfiles)#; trajfiles = pl.flipud(trajfiles)
initfiles = PrintFiles(datadir,'utraj-tr')
initfiles = pl.sort(initfiles)#; initfiles = pl.flipud(initfiles)

#trajfiles = trajfiles[:10]; initfiles = initfiles[:10]

u0, v0, sp_tj = InitData(initfiles,datadir)

# set up empty lists for BL & PM origins:
# evolution of % of trajectories with either origin type appended for each release
#bl_days = []; pm_days = []; st_days = []; origin_days = []
lat = []; lon = []; pres = []; temp = []; PV =[]; q = []; height = []; stratp = [];
convp = []; blz = []

for trajfile in range(trajfiles.shape[0]):
    flist = ReadTxtFile(datadir + trajfiles[trajfile]) # open traj data file
    
    BASETIME = int(flist[0][-1]) # start time of back trajectory
    if trajfile == 0.:
        print 'INITIAL RELEASE TIME:', BASETIME
        NPART = int(flist[3][-1]) # no. of trajectories
        NATTR = int(flist[4][-1]) # no. of attributes
        NCLUST = int(flist[7][-1]) # no. of vertical levels
        pointers = list(itertools.chain(*flist[9:11])) # traj no where each level starts
        pointers = pl.asarray(pointers,dtype='int') # turn into array of integers
        DINT = int(flist[2][3]) # size of data interval in hours
        TIMESTEPS = int(flist[2][-2]) # no. of time steps in an interval
        
        intervals = int(flist[15][4]) # no. of intervals stated in line 16 of file, 5th element of list
        NSTEPS = intervals + 1 #  no. of steps = no. of intervals (plus one)
        data = len(flist[17])-2 + NATTR # no. of rows needed
    		    # some headers in line 18 of file, remove MB units & attributes headers
    		    # add no. of attributes, line 5 in file, final element in list
    
    	# empty array for trajectory data, shape no. of steps by no. of data fields:
    alltraj = pl.zeros([NPART,NSTEPS,data])
    
    	# step, hours, lat, lon, P, dry PT, PV, q, height, stratiform precip, convective precip, BL height
    for traj in range(NPART):
        alltraj[traj] = TrajExtract(NSTEPS,data,flist,traj,intervals)
    	    #for j in range(NSTEPS): # loop over no. of intervals
    	     #   trajdat[i,j] = flist[18+i*45] # data starts at line 19 in file, start loop there
    		                    # end loop at line 19 plus no. of intervals
    				            # each element from the lists extracted becomes a row in the trajdat array
    	#flist = []
    	# extract all the variables from trajdat array & name them:
    if trajfile == 0.:
        stepno = alltraj[:,:,0]; hours = alltraj[:,:,1]
    lat.append(alltraj[:,:,2]); lon.append(alltraj[:,:,3]);
    pres.append(alltraj[:,:,4]); temp.append(alltraj[:,:,5]);
    PV.append(alltraj[:,:,6]); q.append(alltraj[:,:,7]);
    height.append(alltraj[:,:,8])#; stratp.append(alltraj[:,:,9]);
    #convp.append(alltraj[:,:,10]);
    blz.append(alltraj[:,:,-1])
    #stepno = alltraj[:,:,0]; hours = alltraj[:,:,1]; lat = alltraj[:,:,2]; lon = alltraj[:,:,3];
    #pres = alltraj[:,:,4]; temp = alltraj[:,:,5]; PV = alltraj[:,:,6]; q = alltraj[:,:,7];
    #height = alltraj[:,:,8]; stratp = alltraj[:,:,9]; convp = alltraj[:,:,10]; blz = alltraj[:,:,11]

TJTOT = NPART*len(trajfiles)

print NPART, 'TRAJECTORIES PER RELEASE', int(len(trajfiles)), 'RELEASES'
#print 'TRAJECTORIES RELEASED EVERY', int(trajfiles[0][-4:-2]) - int(trajfiles[1][-4:-2]), 'DAYS AT', trajfiles[1][-2:],'Z'
print 'TOTAL NO. OF TRAJECTORIES = ', TJTOT
print NCLUST, 'LEVELS'
print NSTEPS, 'STEPS PER TRAJECTORY, ', int((NSTEPS-1)/4), 'DAY TRAJECTORIES' 
print NATTR, 'ATTRIBUTES: temperature, potential vorticity, specific humidity,',\
            'height, boundary layer height'

del alltraj
lat = pl.asarray(lat); lon = pl.asarray(lon); pres = pl.asarray(pres)
temp = pl.asarray(temp); PV = pl.asarray(PV); q = pl.asarray(q)
height = pl.asarray(height); stratp = pl.asarray(stratp)
convp = pl.asarray(convp); blz = pl.asarray(blz)

#lon[:,0] = lon[:,0]  + 360. # will need to change to if statement or something
theta = temp*(1000/pres)**(0.2854) # dry theta
w = q/(1-q) # mixing ratio

equiv = pl.zeros_like(q)
stratloc = pl.zeros([len(trajfiles),NPART])
for trajfile in range(len(trajfiles)):
    stratloc[trajfile] = StratStep(theta[trajfile],PV[trajfile],height[trajfile])
    for traj in range(NPART):
        equiv[trajfile,traj] = theta_e(q[trajfile,traj],temp[trajfile,traj],
                                                        pres[trajfile,traj])

	#TrajPlot(NPART,NSTEPS,lon,lat)
	#pl.show()
	#pl.savefig('/export/pinggwo/data-01/np838619/traj/trajplot.png')

#for trajfile in range(len(trajfiles)):
    

nc1 = Dataset('/home/np838619/Watershed/ggis198901010000.nc','r')
eralon = nc1.variables['longitude'][:]; eralat = nc1.variables['latitude'][:]
lsm = nc1.variables['LSM'][:]
nc1.close()
lsm = pl.squeeze(lsm)

bl_days = []; pm_days = []; st_days = []; origin_days = []; nl = []
strator = []
for trajfile in range(len(trajfiles)):
    blprops, pmprops, stprops, newlabs = OrOrder(blz[trajfile],height[trajfile],
                                                 stepno,lat[trajfile],lon[trajfile],
                                                 equiv[trajfile],q[trajfile],
                                                stratloc[trajfile],lsm,eralon,eralat)
    bl_days.append(blprops); pm_days.append(pmprops); st_days.append(stprops)
    totprops = blprops + pmprops + stprops; origin_days.append(totprops)
    nl.append(newlabs)
    strator.append(StratOrigin(stratloc[trajfile],lat[trajfile],lon[trajfile]))

bl_days = pl.asarray(bl_days); bl_days = pl.squeeze(bl_days)
pm_days = pl.asarray(pm_days); pm_days = pl.squeeze(pm_days)
st_days = pl.asarray(st_days); st_days = pl.squeeze(st_days)
origin_days = pl.asarray(origin_days); origin_days = pl.squeeze(origin_days)
newlabs = pl.asarray(nl)

#tj_mn = pl.mean(origin_days,axis=0); bl_mn = pl.mean(bl_days,axis=0)
#pm_mn = pl.mean(pm_days,axis=0); st_mn = pl.mean(st_days,axis=0)
#tjprops = pl.array([tj_mn,bl_mn,pm_mn,st_mn])

#f = open('/home/np838619/Trajectory/fluxevol/'+shed+'/jul2010_tj_'+shed+'.csv','w')
#f.write('Total  CATI  CATII  Stratosphere\n')
#pl.savetxt(f,tjprops.T,fmt='%9.5f',delimiter=' ')
#f.close()
#print '****Trajectory proportions written to file*****'

##fig, ax = pl.subplots()
#pl.figure(2)
"""pl.plot(stepno[0]/24,bl_days,label='CAT I origin',color='b',linewidth=2)
pl.plot(stepno[0]/24,pm_days,label='CAT II origin',color='g',linewidth=2)
pl.plot(stepno[0]/24,st_days,label='Strat. origin',color='purple',linewidth=2)
pl.plot(stepno[0]/24,origin_days,label='Total',color='r',linewidth=2)
pl.xlabel('Time along trajectory (days)',fontsize=22)
#ax.set_xticklabels()
pl.ylabel('% of trajectories',fontsize=22); pl.legend(loc=2)"""
#origin_mean = pl.mean(origin_days,axis=0)

sheddir = '/home/np838619/Watershed/shed_defs/'
loc = 'Am'
labs, rlspts = TrajSegLabel(loc)
 
rlspts, labs = RemoveRepeated(rlspts,labs)
rlslabs = pl.zeros([len(labs),3])
rlslabs[:,:2] = rlspts[:]; rlslabs[:,2] = labs[:]
rlslabs[:,:2] = Add360(rlslabs[:,:2])


normals = NormalVector(sheddir,loc)
NSEG = labs.max()


#sp_rel = pl.zeros([rlslabs.shape[0]]); mslp_rel = pl.zeros([rlslabs.shape[0]])
#for pt in range(rlslabs.shape[0]):
#    sp_rel[pt] = BilinInterp(rlslabs[pt],eralon,eralat,surfp)
#    mslp_rel[pt] = BilinInterp(rlslabs[pt],eralon,eralat,mslp)

#pl.plot(sp_rel,label='psurf ERA')
#pl.plot(sp_tj[0,0]*100,label='psurf traj')
#pl.plot(mslp_rel,label='mslp ERA')
#pl.legend()

#nc2 = Dataset('/home/np838619/Watershed/ggap201307201200.nc','r')
#u = nc2.variables['U'][:]; v = nc2.variables['V'][:]
#erapres = nc2.variables['p'][:]
#nc2.close()

#surfp = pl.squeeze(surfp)#; u = pl.squeeze(u); v = pl.squeeze(v)
etalevs = pl.linspace(0.95,0.15,17)

shedflx = pl.zeros([len(trajfiles),rlslabs.shape[0]]) # flux weighted by dl
sf2 = pl.zeros_like(shedflx) # flux unweighted by dl
sf3 = pl.zeros([len(trajfiles),NCLUST,rlslabs.shape[0]])
magflx = pl.zeros_like(shedflx)
atlcount = []; indcount = []; paccount = []; arccount = []; soucount = []
strcount = []; unacount = []
for rls in range(len(trajfiles)):
    shedflx[rls] = FluxCalc(rlslabs,u0[rls],v0[rls],sp_tj[rls,0],etalevs,normals,q[rls],pres[rls])
    sf2[rls] = FluxCalc2(rlslabs,u0[rls],v0[rls],sp_tj[rls,0],etalevs,normals,q[rls],pres[rls])
    sf3[rls] = FluxCalc3(rlslabs,u0[rls],v0[rls],sp_tj[rls],etalevs,normals,q[rls],pres[rls])
    magflx[rls] = FluxCalc(rlslabs,pl.absolute(u0[rls]),pl.absolute(v0[rls]),sp_tj[rls,0],
                      etalevs,pl.absolute(normals),q[rls],pres[rls])
    atc, inc, pac, aic, soc, stc, uac = OriginPartition(newlabs[rls],
                                                    lon[rls],lat[rls])
    atlcount.append(atc)
    indcount.append(inc)
    paccount.append(pac)
    arccount.append(aic)
    soucount.append(soc)
    strcount.append(stc)
    unacount.append(uac)

crosssec = pl.mean(sf3,axis=0)
pl.pcolormesh(crosssec); pl.colorbar()

fluxtot = pl.sum(shedflx,axis=1)/(10**9); print pl.mean(fluxtot), ' Sv'
flxprof = pl.mean(sf2,axis=0) # profile along boundary should be unweighted flux

#for rls in range(len(trajfiles)):
    
#mag = pl.sum(magflx,axis=1)/(10**9)
#print pl.mean(mag), ' Sv'

ordens, tjdens = MainKernelFunc2(newlabs,eralon,eralat,sf3,rlslabs,'all')

#fxmag = pl.zeros([len(trajfiles),NSTEPS]); blmag = pl.zeros_like(fxmag)
#pmmag = pl.zeros_like(fxmag); stmag = pl.zeros_like(fxmag)
#
#for rls in range(len(trajfiles)):
#    fxmag[rls],blmag[rls],pmmag[rls],stmag[rls] = FluxEvol(blz[rls],height[rls],
#                            stepno,lat[rls],lon[rls],equiv[rls],q[rls],theta[rls],
#                            PV[rls],rlslabs,u0[rls],v0[rls],normals,pres[rls],etalevs,
#                            mag[rls],sp_tj[rls,0,:],lsm,eralon,eralat)
#
#fx_mn = pl.mean(fxmag,axis=0); bl_mn = pl.mean(blmag,axis=0)
#pm_mn = pl.mean(pmmag,axis=0); st_mn = pl.mean(stmag,axis=0)
#allmags = pl.array([fx_mn,bl_mn,pm_mn,st_mn])
#f = open('/home/np838619/Trajectory/fluxevol/'+shed+'/jul2010_fx_'+shed+'.csv','w')
#f.write('Total  CATI  CATII  Stratosphere\n')
#pl.savetxt(f,allmags.T,fmt='%9.5f',delimiter=' ')
#f.close()
#print '****Flux proportions written to file****'

#flux_mag,fbl_mag,fpm_mag,fst_mag = FluxEvol(blz[0],height[0],stepno,lat[0],lon[0],equiv[0],
#                                    q[0],theta[0],PV[0],rlslabs,u0[0],v0[0],normals,
#                                    pres[0],etalevs,mag[0],sp_tj[0,0],lsm,eralon,eralat)

#pl.figure(2)
#pl.plot(stepno[0]/24,origin_days[0],color='r',label='% of trajectories',
#                                                                linewidth=2)
#pl.plot(stepno[0]/24,flux_mag,color='k',label='% of flux',linewidth=2)
#pl.plot(stepno[0]/24,fbl_mag,color='b',label='CAT I contribution',
#                                                        ls='--',linewidth=2)
#pl.plot(stepno[0]/24,fpm_mag,color='g',label='CAT II contribution',
#                                                        ls='--',linewidth=2)
#pl.plot(stepno[0]/24,fst_mag,color='purple',label='Strat. contribution',
#                                                        ls='--',linewidth=2)
#pl.xlabel('time along trajectories (days)',fontsize=22)
#pl.ylabel('%',fontsize=22); pl.legend(loc=(0.575,0.270),ncol=1); pl.ylim(0,100)

#atlcount = []; indcount = []; paccount = []; arccount = []; soucount = []
#strcount = []; unacount = []
#for trajfile in range(len(trajfiles)):
#    atc, inc, pac, aic, soc, stc, uac = OriginPartition(newlabs[trajfile],
#                                                    lon[trajfile],lat[trajfile])
#    atlcount.append(atc)
#    indcount.append(inc)
#    paccount.append(pac)
#    arccount.append(aic)
#    soucount.append(soc)
#    strcount.append(stc)
#    unacount.append(uac)

#normals = pl.absolute(normals)
#u0 = pl.absolute(u0); v0 = pl.absolute(v0)

#atlflx = FluxesLoop(atlcount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres)
#pacflx = FluxesLoop(paccount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres)
#indflx = FluxesLoop(indcount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres)
#arcflx = FluxesLoop(arccount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres)
#souflx = FluxesLoop(soucount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres)
#strflx = FluxesLoop(strcount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres)
#unaflx = FluxesLoop(unacount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres)
#
#atlflx2 = FluxesLoop2(atlcount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres)
#pacflx2 = FluxesLoop2(paccount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres)
#indflx2 = FluxesLoop2(indcount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres)
#arcflx2 = FluxesLoop2(arccount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres)
#souflx2 = FluxesLoop2(soucount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres)
#strflx2 = FluxesLoop2(strcount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres)
#unaflx2 = FluxesLoop2(unacount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres)
#
#atlsum = pl.sum(atlflx,axis=1)/(10**9); atllin = pl.mean(atlflx2,axis=0)
#pacsum = pl.sum(pacflx,axis=1)/(10**9); paclin = pl.mean(pacflx2,axis=0)
#indsum = pl.sum(indflx,axis=1)/(10**9); indlin = pl.mean(indflx2,axis=0)
#arcsum = pl.sum(arcflx,axis=1)/(10**9); arclin = pl.mean(arcflx2,axis=0)
#sousum = pl.sum(souflx,axis=1)/(10**9); soulin = pl.mean(souflx2,axis=0)
#strsum = pl.sum(strflx,axis=1)/(10**9); strlin = pl.mean(strflx2,axis=0)
#unasum = pl.sum(unaflx,axis=1)/(10**9); unalin = pl.mean(unaflx2,axis=0)

#catchints = zip(fluxtot,atlsum,indsum,pacsum,arcsum,sousum,strsum,unasum)
#f = open('/home/np838619/Trajectory/fluxpart/'+year+'/'+shed+'/catchints.csv','w')
#pl.savetxt(f,catchints); f.close()
#
#meanprofiles = zip(flxprof,atllin,indlin,paclin,arclin,soulin,strlin,unalin)
#f = open('/home/np838619/Trajectory/fluxpart/'+year+'/'+shed+'/meanprofs.csv','w')
#pl.savetxt(f,meanprofiles); f.close()

#print 'Atlantic = ', pl.mean(atlsum), ' Sv, correlation = ', pearsonr(fluxtot,atlsum)
#print 'Pacific = ', pl.mean(pacsum), ' Sv, correlation = ', pearsonr(fluxtot,pacsum)
#print 'Indian = ', pl.mean(indsum), ' Sv, correlation = ', pearsonr(fluxtot,indsum)
#print 'Arctic = ', pl.mean(arcsum), ' Sv, correlation = ', pearsonr(fluxtot,arcsum)
#print 'Southern = ', pl.mean(sousum), ' Sv, correlation = ', pearsonr(fluxtot,sousum)
#print 'Stratosphere = ', pl.mean(strsum), ' Sv, correlation = ', pearsonr(fluxtot,strsum)
#print 'No Origin = ', pl.mean(unasum), ' Sv, correlation = ', pearsonr(fluxtot,unasum)
#means = pl.array([pl.mean(atlsum),pl.mean(indsum),pl.mean(pacsum),
#                  pl.mean(arcsum),pl.mean(sousum),pl.mean(strsum),pl.mean(unasum)])
#pl.savetxt(f,means); f.write('\n'); f.write('Correlations\n')

#atlcor = pearsonr(fluxtot,atlsum); indcor = pearsonr(fluxtot,indsum)
#paccor = pearsonr(fluxtot,pacsum); arccor = pearsonr(fluxtot,arcsum)
#soucor = pearsonr(fluxtot,sousum); strcor = pearsonr(fluxtot,strsum)
#unacor = pearsonr(fluxtot,unasum)
#
#cors = pl.array([atlcor[0],indcor[0],paccor[0],arccor[0],soucor[0],strcor[0],unacor[0]])

# write mean flux data to file:
#MeanFluxesOutput(year,shed,means,cors,mags,cors2)

savedir = '/home/np838619/Trajectory/fluxpart/'+year+'/'+shed+'/'
#pl.figure(1,figsize=(15,15))
#exec(open('/home/np838619/Trajectory/plot_fluxpart.py').read())
#pl.subplots_adjust(top=0.85,bottom=0.11)
#pl.savefig(savedir+'flux_mags_'+shed+'_'+year+'.png')

#pl.figure(2)
#pl.plot(flxprof,color='k',linewidth=2,label='$F_{tot}$')
#pl.plot(paclin,color='g',ls='--',linewidth=2,label='$F_{Pac}$')
#pl.plot(atllin,color='r',ls='--',linewidth=2,label='$F_{Atl}$')
##pl.plot(arclin,color='deeppink',ls='--',linewidth=2,label='$F_{Arc}$')
#pl.ylabel('kg m$^{-1}$ s$^{-1}$',fontsize=25); pl.xlim(0,rlslabs.shape[0])
#pl.legend(loc=0,frameon=False,fontsize=18,handlelength=1.6)
#pl.subplots_adjust(right=0.92,left=0.14)
#pl.savefig(savedir+'flux_prfs_'+shed'_'+year+'.png')

#atlcount = pl.vstack(atlcount); indcount = pl.vstack(indcount); paccount = pl.vstack(paccount)
#arccount = pl.vstack(arccount); soucount = pl.vstack(soucount); strcount = pl.vstack(strcount)
#unacount = pl.vstack(unacount)
#
#atlprop = CatchProps(atlcount,TJTOT)#; print atlprop[0], ' % TRAJECTORIES WITH ATLANTIC ORIGIN'
#indprop = CatchProps(indcount,TJTOT)#; print indprop[0], ' % TRAJECTORIES WITH INDIAN ORIGIN'
#pacprop = CatchProps(paccount,TJTOT)#; print pacprop[0], ' % TRAJECTORIES WITH PACIFIC ORIGIN'
#arcprop = CatchProps(arccount,TJTOT)#; print arcprop[0], ' % TRAJECTORIES WITH ARCTIC ORIGIN'
#souprop = CatchProps(soucount,TJTOT)#; print souprop[0], ' % TRAJECTORIES WITH SOUTHERN ORIGIN'
#strprop = CatchProps(strcount,TJTOT)#; print strprop[0], ' % TRAJECTORIES WITH STRATOSPHERIC ORIGIN'
#unaprop = CatchProps(unacount,TJTOT)#; print unaprop[0], ' % TRAJECTORIES WITH NO ORIGIN'

#catchprops = pl.array([atlprop[0],indprop[0],pacprop[0],arcprop[0],souprop[0],
#                       strprop[0],unaprop[0]])
#propsplit = pl.zeros([5,2])
#propsplit[0,:] = atlprop[1:]; propsplit[1,:] = indprop[1:]
#propsplit[2,:] = pacprop[1:]; propsplit[3,:] = arcprop[1:]
#propsplit[4,:] = souprop[1:]
#TrajPropsOutput(year,shed,catchprops,propsplit)

#totCAT1, totCAT2 = FluxCATFull(newlabs,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres)
#atlCAT1, atlCAT2 = FluxCATPart(newlabs,atlcount,rlslabs,normals,etalevs,q,
#                                                           sp_tj,u0,v0,pres)
#pacCAT1, pacCAT2 = FluxCATPart(newlabs,paccount,rlslabs,normals,etalevs,q,
#                                                           sp_tj,u0,v0,pres)
#indCAT1, indCAT2 = FluxCATPart(newlabs,indcount,rlslabs,normals,etalevs,q,
#                                                           sp_tj,u0,v0,pres)
#arcCAT1, arcCAT2 = FluxCATPart(newlabs,arccount,rlslabs,normals,etalevs,q,
#                                                           sp_tj,u0,v0,pres)
#souCAT1, souCAT2 = FluxCATPart(newlabs,soucount,rlslabs,normals,etalevs,q,
#                                                           sp_tj,u0,v0,pres)
#
#totC1s = pl.sum(totCAT1,axis=1)/(10**9); totC2s = pl.sum(totCAT2,axis=1)/(10**9)
#atlC1s = pl.sum(atlCAT1,axis=1)/(10**9); atlC2s = pl.sum(atlCAT2,axis=1)/(10**9)
#indC1s = pl.sum(indCAT1,axis=1)/(10**9); indC2s = pl.sum(indCAT2,axis=1)/(10**9)
#pacC1s = pl.sum(pacCAT1,axis=1)/(10**9); pacC2s = pl.sum(pacCAT2,axis=1)/(10**9)
#arcC1s = pl.sum(arcCAT1,axis=1)/(10**9); arcC2s = pl.sum(arcCAT2,axis=1)/(10**9)
#souC1s = pl.sum(souCAT1,axis=1)/(10**9); souC2s = pl.sum(souCAT2,axis=1)/(10**9)
#
#CAT1ints = zip(totC1s,atlC1s,indC1s,pacC1s,arcC1s,souC1s)
#CAT2ints = zip(totC2s,atlC2s,indC2s,pacC2s,arcC2s,souC2s)
#f = open('/home/np838619/Trajectory/fluxpart/'+year+'/'+shed+'/CAT1ints.csv','w')
#pl.savetxt(f,CAT1ints)
#f.close()
#f = open('/home/np838619/Trajectory/fluxpart/'+year+'/'+shed+'/CAT2ints.csv','w')
#pl.savetxt(f,CAT2ints)
#f.close()

#totC1m = pl.mean(totC1s); totC2m = pl.mean(totC2s)
#atlC1m = pl.mean(atlC1s); atlC2m = pl.mean(atlC2s)
#indC1m = pl.mean(indC1s); indC2m = pl.mean(indC2s)
#pacC1m = pl.mean(pacC1s); pacC2m = pl.mean(pacC2s)
#arcC1m = pl.mean(arcC1s); arcC2m = pl.mean(arcC2s)
#souC1m = pl.mean(souC1s); souC2m = pl.mean(souC2s)
#
#totCT1u, totCT2u = FluxCATFull2(newlabs,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres)
#atlCT1u, atlCT2u = FluxCATPart2(newlabs,atlcount,rlslabs,normals,etalevs,q,
#                                                           sp_tj,u0,v0,pres)
#pacCT1u, pacCT2u = FluxCATPart2(newlabs,paccount,rlslabs,normals,etalevs,q,
#                                                           sp_tj,u0,v0,pres)
#indCT1u, indCT2u = FluxCATPart2(newlabs,indcount,rlslabs,normals,etalevs,q,
#                                                           sp_tj,u0,v0,pres)
#arcCT1u, arcCT2u = FluxCATPart2(newlabs,arccount,rlslabs,normals,etalevs,q,
#                                                           sp_tj,u0,v0,pres)
#souCT1u, souCT2u = FluxCATPart2(newlabs,soucount,rlslabs,normals,etalevs,q,
#                                                           sp_tj,u0,v0,pres)
#
#totC1p = pl.mean(totCT1u,axis=0); totC2p = pl.mean(totCT2u,axis=0)
#atlC1p = pl.mean(atlCT1u,axis=0); atlC2p = pl.mean(atlCT2u,axis=0)
#pacC1p = pl.mean(pacCT1u,axis=0); pacC2p = pl.mean(pacCT2u,axis=0)
#indC1p = pl.mean(indCT1u,axis=0); indC2p = pl.mean(indCT2u,axis=0)
#arcC1p = pl.mean(arcCT1u,axis=0); arcC2p = pl.mean(arcCT2u,axis=0)
#souC1p = pl.mean(souCT1u,axis=0); souC2p = pl.mean(souCT2u,axis=0)
#
#CAT1profiles = zip(totC1p,atlC1p,indC1p,pacC1p,arcC1p,souC1p)
#CAT2profiles = zip(totC2p,atlC2p,indC2p,pacC2p,arcC2p,souC2p)
#f = open('/home/np838619/Trajectory/fluxpart/'+year+'/'+shed+'/CAT1profs.csv','w')
#pl.savetxt(f,CAT1profiles); f.close()
#f = open('/home/np838619/Trajectory/fluxpart/'+year+'/'+shed+'/CAT2profs.csv','w')
#pl.savetxt(f,CAT2profiles); f.close()


#catsplits = pl.zeros([6,2])
#catsplits[0,0] = totC1m; catsplits[0,1] = totC2m
#catsplits[1,0] = atlC1m; catsplits[1,1] = atlC2m
#catsplits[2,0] = indC1m; catsplits[2,1] = indC2m
#catsplits[3,0] = pacC1m; catsplits[3,1] = pacC2m
#catsplits[4,0] = arcC1m; catsplits[4,1] = arcC2m
#catsplits[5,0] = souC1m; catsplits[5,1] = souC2m
##CATsplitOutput(year,shed,catsplits)
#
#totC1_sea, totC1_land = CAT1LandSeaFull(lsm,newlabs,rlslabs,normals,etalevs,q,
#                                        sp_tj,u0,v0,pres,lon,lat,eralon,eralat)
#atlC1_sea, atlC1_land = CAT1LandSea(lsm,atlcount,rlslabs,normals,etalevs,q,sp_tj,
#                                  u0,v0,pres,lon,lat,eralon,eralat)
#pacC1_sea, pacC1_land = CAT1LandSea(lsm,paccount,rlslabs,normals,etalevs,q,sp_tj,
#                                  u0,v0,pres,lon,lat,eralon,eralat)
#indC1_sea, indC1_land = CAT1LandSea(lsm,indcount,rlslabs,normals,etalevs,q,sp_tj,
#                                  u0,v0,pres,lon,lat,eralon,eralat)
#arcC1_sea, arcC1_land = CAT1LandSea(lsm,arccount,rlslabs,normals,etalevs,q,sp_tj,
#                                  u0,v0,pres,lon,lat,eralon,eralat)
#souC1_sea, souC1_land = CAT1LandSea(lsm,soucount,rlslabs,normals,etalevs,q,sp_tj,
#                                  u0,v0,pres,lon,lat,eralon,eralat)
#
#totC1_bl = CAT1BL_full(lsm,newlabs,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres,
#                                                              eralon,eralat)
#atlC1_bl = CAT1BL_part(lsm,atlcount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres,
#                                                              eralon,eralat)
#pacC1_bl = CAT1BL_part(lsm,paccount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres,
#                                                              eralon,eralat)
#indC1_bl = CAT1BL_part(lsm,indcount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres,
#                                                              eralon,eralat)
#arcC1_bl = CAT1BL_part(lsm,arccount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres,
#                                                              eralon,eralat)
#souC1_bl = CAT1BL_part(lsm,soucount,rlslabs,normals,etalevs,q,sp_tj,u0,v0,pres,
#                                                              eralon,eralat)
#
#totC1s_sea = pl.sum(totC1_sea,axis=1)/(10**9); totC1s_lnd = pl.sum(totC1_land,axis=1)/(10**9)
#atlC1s_sea = pl.sum(atlC1_sea,axis=1)/(10**9); atlC1s_lnd = pl.sum(atlC1_land,axis=1)/(10**9)
#pacC1s_sea = pl.sum(pacC1_sea,axis=1)/(10**9); pacC1s_lnd = pl.sum(pacC1_land,axis=1)/(10**9)
#indC1s_sea = pl.sum(indC1_sea,axis=1)/(10**9); indC1s_lnd = pl.sum(indC1_land,axis=1)/(10**9)
#arcC1s_sea = pl.sum(arcC1_sea,axis=1)/(10**9); arcC1s_lnd = pl.sum(arcC1_land,axis=1)/(10**9)
#souC1s_sea = pl.sum(souC1_sea,axis=1)/(10**9); souC1s_lnd = pl.sum(souC1_land,axis=1)/(10**9)
#
#CAT1ints_sea = zip(totC1s_sea,atlC1s_sea,pacC1s_sea,indC1s_sea,arcC1s_sea,souC1s_sea)
#CAT1ints_lnd = zip(totC1s_lnd,atlC1s_lnd,pacC1s_lnd,indC1s_lnd,arcC1s_lnd,souC1s_lnd)
#f = open('/home/np838619/Trajectory/fluxpart/'+year+'/'+shed+'/CAT1sea.csv','w')
#pl.savetxt(f,CAT1ints_sea); f.close()
#f = open('/home/np838619/Trajectory/fluxpart/'+year+'/'+shed+'/CAT1land.csv','w')
#pl.savetxt(f,CAT1ints_lnd); f.close()

#totC1s_bl = pl.sum(totC1_bl,axis=1); atlC1s_bl = pl.sum(atlC1_bl,axis=1)
#pacC1s_bl = pl.sum(pacC1_bl,axis=1); indC1s_bl = pl.sum(indC1_bl,axis=1)
#arcC1s_bl = pl.sum(arcC1_bl,axis=1); souC1s_bl = pl.sum(souC1_bl,axis=1)

#CAT1ints_bl = zip(totC1s_bl,atlC1s,pacC1s_bl,indC1s_bl,arcC1s_bl,souC1s_bl)
#f = open('/home/np838619/Trajectory/fluxpart/'+year+'/'+shed+'/CAT1_BL.csv','w')
#pl.savetxt(f,CAT1ints_bl); f.close()

elapsed = timeit.default_timer() - start_time
print elapsed

#fig,ax = pl.subplots(1,2,figsize=(12,6))
#ax1 = pl.subplot(121)
#pl.plot(atllin,color='r',linewidth=2.,label='$F_{Atl}$')
#pl.plot(atlC1p,color='r',linewidth=2.,ls='--',label='CAT I')
#pl.plot(atlC2p,color='r',linewidth=2.,ls=':',label='CAT II')
#pl.xlim(0,rlslabs.shape[0]); pl.ylabel('kg m$^{-1}$ s$^{-1}$',fontsize=25)
#pl.text(3,-225,'Atlantic',fontsize=22)
#ax2 = pl.subplot(122)
#pl.plot(paclin,color='g',linewidth=2.,label='$F_{Pac}$')
#pl.plot(pacC1p,color='g',linewidth=2.,ls='--',label='CAT I')
#pl.plot(pacC2p,color='g',linewidth=2.,ls=':',label='CAT II')
#ax2.yaxis.set_major_formatter(pl.NullFormatter())
#pl.xlim(0,rlslabs.shape[0]); pl.ylim(-250,50)
#pl.text(3,-225,'Pacific',fontsize=22)
#pl.legend(loc=4,fontsize=22)
#pl.subplots_adjust(wspace=0.07,left=0.1,right=0.97)
#
#fig,ax = pl.subplots(1,2,figsize=(12,6))
#ax1 = pl.subplot(121)
#pl.plot(totC1p,color='k',linewidth=2.,label='$F_{Tot}$')
#pl.plot(atlC1p,color='r',linewidth=2.,ls='--',label='$F_{Atl}$')
#pl.plot(pacC1p,color='g',linewidth=2.,ls='--',label='$F_{Pac}$')
#pl.xlim(0,rlslabs.shape[0]); pl.ylabel('kg m$^{-1}$ s$^{-1}$',fontsize=25)
#pl.text(3,-225,'CAT I',fontsize=22)
#ax2 = pl.subplot(122)
#pl.plot(totC2p,color='k',linewidth=2.,label='$F_{Tot}$')
#pl.plot(atlC2p,color='r',linewidth=2.,ls='--',label='$F_{Atl}$')
#pl.plot(pacC2p,color='g',linewidth=2.,ls='--',label='$F_{Pac}$')
#ax2.yaxis.set_major_formatter(pl.NullFormatter())
#pl.xlim(0,rlslabs.shape[0]); pl.ylim(-250,50)
#pl.text(3,-225,'CAT II',fontsize=22)
#pl.legend(loc=4,fontsize=22)
#pl.subplots_adjust(wspace=0.07,left=0.1,right=0.97)

#CatchScat(atlcount,indcount,paccount,arccount,soucount)
#pl.figure(1); density = DensityPlot(atlcount,rlslabs)
#pl.figure(2); density = DensityPlot(paccount,rlslabs)
#pl.figure(3); density = DensityPlot(indcount,rlslabs)
#pl.figure(4); density = DensityPlot(arccount,rlslabs)
#pl.figure(5); density = DensityPlot(soucount,rlslabs)
# Read ERA-Interim mask
"""maskfile = Dataset('/home/np838619/PminusE_data/ERA_Int/direct_EmP.nc','r')
lsm = maskfile.variables['LSM'][:]
eralat = maskfile.variables['lat'][:]
eralon = maskfile.variables['lon'][:]
maskfile.close()
lsm = pl.squeeze(lsm)""" # get rid of 1-D axis

#LandSeaOrigin(lsm,eralat,eralon,blorigin)

#q_or = pl.zeros([len(allfiles),NPART,NSTEPS])
#q_strat = pl.zeros_like(q_or)
#for trajfile in range(len(allfiles)):
## find specific humidities of trajectories with ...
#    q_or[trajfile] = OriginFlux(newlabs[trajfile],q[trajfile]) # ... BL or PM origin
#    q_strat[trajfile] = StratVar(strator[trajfile],q[trajfile]) # ... stratospheric origin


#start_time = timeit.default_timer()
#a,b,c = NewFunc(blz,height,stepno,lat,lon,equiv,q,stratloc)
#elapsed = timeit.default_timer() - start_time
#print elapsed

#g = [7, 84,200, 262,335,449,515, 631,1051] # trajectories time series to plot
#k = [] # empty list for matplotlib.lines.Line2D object and trajectory label
#ax, fig = pl.subplots(3,1,figsize=(12,10))#ax, fig = pl.subplots(2,1)
#for i in range(len(g)):
#    trajno = g[i]
#    z = AttrPlot(equiv[-1,trajno],q[-1,trajno],height[-1,trajno],newlabs[-1,trajno])
#    k.append(z) # add matplotlib.lines.Line2D object & trajectory label to list
#B = pl.zeros([2],dtype=object) # 1 entry will be blue line, the other red
#for i in range(len(k)):
#    if B[0] == 0.: # if 1st array entry is empty ...
#        if k[i][1] == 1.: # ... & if 2nd entry of list tuple is label '1' ...
#            B[0] = k[i][0] # ... add matplotlib.lines.Line2D object to array
#    elif B[1] == 0.: # ... if 2nd array entry is empty ... 
#        if k[i][1] == 2.: # ... & if 2nd entry of list tuple is label '2' ...
#            B[1] = k[i][0] # ... add matplotlib.lines.Line2D object to array
#    elif B[0] != 0. and B[1] != 0.:
#        break
#    
#ax.legend((B[0],B[1]),('CAT I','CAT II'),loc=(0.18,0.57))
#pl.annotate('Trajectories: ' + str(g),xy=(0.15,0.3),xycoords='figure fraction')

# histograms of traj height at origin:
#ax, fig = pl.subplots(2,1)
#pl.subplot(211); ZPDF(pl.reshape(newlabs,(len(allfiles)*NPART,4)),pl.reshape(height,(len(allfiles)*NPART,NSTEPS)),1)
#ax.text(0.8,0.85,'CAT I',fontsize=16)
#pl.subplot(212); ZPDF(pl.reshape(newlabs,(len(allfiles)*NPART,4)),pl.reshape(height,(len(allfiles)*NPART,NSTEPS)),2)
#ax.text(0.8,0.4,'CAT II',fontsize=16)
#ax.clear

#ax,fig = pl.subplots(2,1)
#pl.subplot(211); PrecPDF(pl.reshape(newlabs,(len(allfiles)*NPART,4)),pl.reshape(stratp,(len(allfiles)*NPART,NSTEPS)),1,20)
#ax.text(0.8,0.85,'CAT I',fontsize=16)
#pl.subplot(212); PrecPDF(pl.reshape(newlabs,(len(allfiles)*NPART,4)),pl.reshape(stratp,(len(allfiles)*NPART,NSTEPS)),2,20)
#ax.text(0.8,0.4,'CAT II',fontsize=16); pl.xlabel('precipitation (mm)',fontsize=16);
#pl.subplots_adjust(top=0.92,hspace=0.17)
#pl.suptitle('Stratiform precipitation',fontsize=23)
#
#ax,fig = pl.subplots(2,1)
#pl.subplot(211); PrecPDF(pl.reshape(newlabs,(len(allfiles)*NPART,4)),pl.reshape(convp,(len(allfiles)*NPART,NSTEPS)),1,20)
#ax.text(0.8,0.85,'CAT I',fontsize=16)
#pl.subplot(212); PrecPDF(pl.reshape(newlabs,(len(allfiles)*NPART,4)),pl.reshape(convp,(len(allfiles)*NPART,NSTEPS)),2,20)
#ax.text(0.8,0.4,'CAT II',fontsize=16); pl.xlabel('precipitation (mm)',fontsize=16)
#pl.subplots_adjust(top=0.92,hspace=0.17)
#pl.suptitle('Convective precipitation',fontsize=23)
#
##pl.figure(5)
#ax,fig = pl.subplots(2,1,figsize=(8,8))
#pl.subplot(211)
#pcm = AttrScatter(pl.reshape(newlabs,(len(allfiles)*NPART,4)),
#                  pl.reshape(stratp*1000,(len(allfiles)*NPART,NSTEPS)),1)
#pl.title('CAT I origin',fontsize=20)
##pl.figure(6)
#pl.subplot(212)
#pcm = AttrScatter(pl.reshape(newlabs,(len(allfiles)*NPART,4)),
#                  pl.reshape(stratp*1000,(len(allfiles)*NPART,NSTEPS)),2)
#pl.title('CAT II origin',fontsize=20)
#f=pl.gcf(); colax = f.add_axes([0.2,0.1,0.6,0.03])
#clb = pl.colorbar(pcm,cax=colax,orientation='horizontal',extend='max')
#clb.set_label('Stratiform precipitation at origin (mm)',fontsize=20)
#clb.ax.tick_params(labelsize=14)
#pl.subplots_adjust(top=0.95,bottom=0.15,hspace=0.12)