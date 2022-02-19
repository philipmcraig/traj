# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 19:17:35 2018

@author: np838619
"""

import pylab as pl
from scipy import stats
from netCDF4 import Dataset

def PlotTitle(shed):
    
    if shed == 'amr':
        title = 'Americas'; loc = 'Am'#; c = 'b'
    elif shed == 'afr':
        title = 'Africa'; loc = 'AfMe'#; c = 'k'
    elif shed == 'eaa':
        title = 'South-East Asia'; loc = 'EAA'#; c = 'g'
    elif shed == 'ara':
        title = 'Arctic Atlantic'; loc = 'ArA'#; c = 'r'
    elif shed == 'ari':
        title = 'Arctic Indian'; loc = 'ArI'#; c = 'r'
    elif shed == 'arp':
        title = 'Arctic Pacific'; loc = 'ArP'#; c = 'r'
    elif shed == 'soa':
        title = 'Southern Atlantic'; loc = 'SOA'#; c = 'darkgoldenrod'
    elif shed == 'soi':
        title = 'Southern Indian'; loc = 'SOI'#; c = 'darkgoldenrod'
    elif shed == 'sop':
        title = 'Southern Pacific'; loc = 'SOP'#; c = 'darkgoldenrod'
    
    return title, loc

pl.close('all')
clusdir = '/storage/shared/glusterfs/scenario/users/qx911590/traj/'
homedir = '/home/users/qx911590/np838619/Trajectory/fluxpart/'
sheddir = '/home/users/qx911590/np838619/Watershed/shed_defs/'

sheds = ['amr','afr','eaa','ara','ari','arp','soa','soi','sop']
rlspts = []

catchints = pl.zeros([len(sheds),8])
catchsem = pl.zeros_like(catchints)
profiles = []
c1int = pl.zeros([len(sheds),6]); c2int = pl.zeros_like(c1int)
c1sem = pl.zeros_like(c1int); c2sem = pl.zeros_like(c2int)
csecs = []; cs_sm = []; pres = []

c1prfs = []; c2prfs = []

# loop over sheds:
for sd in range(len(sheds)):

    shed = sheds[sd]; title, loc = PlotTitle(shed)
    if shed == 'afr':
        ec = 'r'
    else:
        ec = 'k'
    datadir = homedir# +'/'+ shed + '/'
    
    rlspts.append(pl.genfromtxt(sheddir+loc+'_traj_release_new.txt',skip_header=5))
    
    # read in some files here:
    # total & partioned moisture fluxes for each trajectory release in Sv:
    INTS = pl.genfromtxt(datadir+'annmeans_catchints_'+shed+'.csv',delimiter=',',skip_header=1)
    catchints[sd] = INTS[-1,1:]; catchsem[sd] = stats.sem(INTS[:-1,1:],axis=0)
    # mean profiles of total & partioned moisture fluxes in kg/m/s:
    PRFS = pl.genfromtxt(datadir+'profiles_annmean_'+shed+'.csv',delimiter=' ',skip_header=1)
    profiles.append(PRFS)
        # total & partitioned CAT1/2 moisture fluxes for each release in Sv
    CSPL = pl.genfromtxt(datadir+'annmeans_CATsplit_'+shed+'.csv',delimiter=',',skip_header=2)
    c1int[sd] = CSPL[-1,1::2]; c2int[sd] = CSPL[-1,2::2]
    c1sem[sd] = stats.sem(CSPL[:-1,1::2],axis=0)
    c2sem[sd]= stats.sem(CSPL[:-1,2::2],axis=0)
    
    #CAT1ints = pl.genfromtxt(datadir+'CAT1ints.csv',delimiter=' ',skip_header=1)
    #CAT2ints = pl.genfromtxt(datadir+'CAT2ints.csv',delimiter=' ',skip_header=1)
        # mean profiles of total & partioned moisture fluxes in kg/m/s
    C1P = pl.genfromtxt(datadir+'CAT1prf_annmean_'+shed+'.csv',delimiter=' ',skip_header=1)
    C2P = pl.genfromtxt(datadir+'CAT2prf_annmean_'+shed+'.csv',delimiter=' ',skip_header=1)
    c1prfs.append(C1P); c2prfs.append(C2P)
    
    CROSS = pl.zeros([5,17,rlspts[sd].shape[0]]); PRES = pl.zeros_like(CROSS)
    years = ['2010','2011','2012','2013','2014']
    for yr in range(len(years)):
        path = clusdir + years[yr] + '/' + sheds[sd] + '/CSEC/'
        nc1 = Dataset(path+'CSEC'+years[yr]+'_'+shed+'.nc','r')
        var = nc1.variables['flux cross-section'][:]
        nc1.close()
        CROSS[yr] = pl.mean(var,axis=0)/9.81
        nc2 = Dataset(path+'pres'+years[yr]+'_'+shed+'.nc','r')
        var = nc2.variables['pressure cross-section'][:]
        nc2.close()
        PRES[yr] = pl.mean(var,axis=0)
        
    csecs.append(pl.mean(CROSS,axis=0)); cs_sm.append(stats.sem(CROSS,axis=0))
    pres.append(pl.mean(PRES,axis=0))
   

title = ['(a) Americas','(b) Africa','(c) South-East Asia',
         '(d) Arctic Atlantic','(e) Arctic Indian','(f) Arctic Pacific',
         '(g) Southern Atlantic','(h) Southern Indian','(i) Southern Pacific']
fluxes = ['Total','Atlantic','Indian','Pacific','Arctic','Southern',
                          'Strat.','None']

###############################################################################
fig, ax = pl.subplots(3,3,figsize=(13,11))
y = pl.arange(len(fluxes))
xl = [(-0.4,0.2),(-0.4,0.2),(-0.4,0.6),(-0.2,0.2),(-0.2,0.2),(-0.2,0.2),
      (-0.4,0.2),(-0.4,0.2),(-0.8,0.4)]
error_kw=dict(elinewidth=2,ecolor='k',capsize=5,capthick=2)

for i in range(9):
    axx = pl.subplot(3,3,i+1)
    br = axx.barh(y,catchints[i,:],align='center',ec='k',xerr=catchsem[i],
        error_kw=error_kw)
    br[0].set_color('grey'); br[-1].set_color('w')
    br[1].set_color('r'); br[2].set_color('b'); br[3].set_color('g')
    br[4].set_color('deeppink'); br[5].set_color('darkgoldenrod')
    br[6].set_color('skyblue')
    br[1].set_edgecolor('k'); br[2].set_edgecolor('k'); br[3].set_edgecolor('k')
    br[4].set_edgecolor('k'); br[5].set_edgecolor('k'); br[6].set_edgecolor('k')
    br[-1].set_edgecolor('k'); br[0].set_edgecolor('k')
    pl.title(title[i],fontsize=18)
    pl.xlim(xl[i][0],xl[i][1])
    pl.xticks(fontsize=14); axx.set_yticks(y)
    axx.grid(axis='x')
    if i in (0,3,6):
        axx.set_yticklabels(fluxes,fontsize=15)
    else:
        axx.yaxis.set_major_formatter(pl.NullFormatter())
    if i in (3,4,5):
        axx.set_xticklabels(['$-$0.2','','$-$0.1','',0.0,'',0.1,'',0.2])
    if i > 5:
        pl.xlabel('Sv',fontsize=18)
    if i == 0:
        axx.legend(br[:4],fluxes[:4],ncol=1,fontsize=14,columnspacing=0.5,loc=2)
    if i == 2:
        axx.legend(br[4:],fluxes[4:],ncol=1,fontsize=14,columnspacing=0.5,loc=1)
    if i == 4:
        axx.arrow(-0.03,6.5,-0.1,0,head_width=0.4,head_length=0.03,width=0.1,color='k')
        axx.annotate('westward/\n southward',xy=(-0.15,6.8),fontsize=13)
        axx.arrow(0.03,6.5,0.1,0,head_width=0.4,head_length=0.03,width=0.1,color='k')
        axx.annotate('eastward/\n northward',xy=(0.05,6.8),fontsize=13)
pl.tight_layout()
#pl.savefig(homedir+'/annmean_plots/panels/fluxpart_panels.png')
###############################################################################

###############################################################################
eraprofs = pl.genfromtxt('/home/users/qx911590/np838619/Watershed/profiles_1014.csv',delimiter=',')
eraprofs = eraprofs.T
fig, ax = pl.subplots(3,3,figsize=(13,11))

P = []
for i in range(9):
    axx = pl.subplot(3,3,i+1)
    p0 = axx.plot(profiles[i][:,0],color='k',lw=2,label='total',zorder=2)
    EI = axx.plot(eraprofs[i,:len(rlspts[i])],color='grey',lw=1.5,zorder=1)
    if i == 0:
        p1 = axx.plot(profiles[i][:,1],color='r',lw=2,ls='--',label='Atlantic',zorder=2)
        axx.plot(profiles[i][:,3],color='g',lw=2,ls='--',label='Pacific',zorder=2)
        a = pl.arange(0,155,20); b = pl.around(rlspts[i][a,1],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        P.append(p0[0]); P.append(p1[0])
        pl.ylabel('kg/m/s',fontsize=16); pl.ylim(-250,150)
    elif i == 1:
        p2 = axx.plot(profiles[i][:,2],color='b',lw=2,ls='--',label='Indian',zorder=2)
        axx.plot(profiles[i][:,1],color='r',lw=2,ls='--',label='Atlantic',zorder=2)
        a = pl.arange(0,214,30); b = pl.around(rlspts[i][a,1],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        P.append(p2[0]); axx.legend(P,fluxes[:3],loc=4)
        pl.xlabel('latitude',fontsize=14); pl.ylim(-250,150)
    elif i == 2:
        p3 = axx.plot(profiles[i][:,3],color='g',lw=2,ls='--',label='Pacific',zorder=2)
        axx.plot(profiles[i][:,2],color='b',lw=2,ls='--',label='Indian',zorder=2)
        a = pl.arange(0,165,20); b = pl.around(rlspts[i][a,1],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        P.append(p3[0]); pl.ylim(-250,150)
        axx.legend([EI[0]],['ERA-Interim'],loc=3)
    elif i == 3:
        axx.plot(profiles[i][:,3],color='g',lw=2,ls='--',label='Pacific',zorder=2)
        axx.plot(profiles[i][:,1],color='r',lw=2,ls='--',label='Atlantic',zorder=2)
        p4 = axx.plot(profiles[i][:,4],color='deeppink',lw=2,ls='--',label='Arctic',zorder=2)
        a = pl.arange(0,197,30); b = pl.around(rlspts[i][a,0],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        P.append(p4[0])#; axx.annotate(title[i],(90,-55),fontsize=16)
        pl.ylabel('kg/m/s',fontsize=16); pl.ylim(-60,60)
    elif i == 4:
        axx.plot(profiles[i][:,2],color='b',lw=2,ls='--',label='Indian',zorder=2)
        axx.plot(profiles[i][:,1],color='r',lw=2,ls='--',label='Atlantic',zorder=2)
        axx.plot(profiles[i][:,4],color='deeppink',lw=2,ls='--',label='Arctic',zorder=2)
        a = pl.arange(0,103,20); b = pl.around(rlspts[i][a,0],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        pl.xlabel('longitude',fontsize=14); pl.ylim(-60,60)
        #axx.annotate(title[i],(50,-55),fontsize=16)
    elif i == 5:
        axx.plot(profiles[i][:,3],color='g',lw=2,ls='--',label='Pacific',zorder=2)
        axx.plot(profiles[i][:,1],color='r',lw=2,ls='--',label='Atlantic',zorder=2)
        #axx.plot(profiles[i][:,2],color='b',lw=2,ls='--',label='Indian')
        axx.plot(profiles[i][:,4],color='deeppink',lw=2,ls='--',label='Arctic',zorder=2)
        a = pl.arange(0,220,30); b = pl.around(rlspts[i][a,0],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        axx.legend(P[3:],fluxes[3:5],loc=4); pl.ylim(-60,60)
        #axx.annotate(title[i],(110,-55),fontsize=16)
    elif i == 6:
        axx.plot(profiles[i][:,3],color='g',lw=2,ls='--',label='Pacific',zorder=2)
        axx.plot(profiles[i][:,1],color='r',lw=2,ls='--',label='Atlantic',zorder=2)
        p5 = axx.plot(profiles[i][:,5],color='darkgoldenrod',lw=2,ls='--',label='Southern',zorder=2)
        a = pl.arange(0,185,30); b = pl.around(rlspts[i][a,0],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        P.append(p5[0])
        pl.ylabel('kg/m/s',fontsize=16); pl.ylim(-200,100)
    elif i == 7:
        axx.plot(profiles[i][:,2],color='b',lw=2,ls='--',label='Indian',zorder=2)
        axx.plot(profiles[i][:,1],color='r',lw=2,ls='--',label='Atlantic',zorder=2)
        axx.plot(profiles[i][:,5],color='darkgoldenrod',lw=2,ls='--',label='Southern',zorder=2)
        a = pl.arange(0,184,30); b = pl.around(rlspts[i][a,0],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        pl.xlabel('longitude',fontsize=14); pl.ylim(-200,100)
        axx.legend(P[5:],fluxes[5:6],loc=2)
    elif i == 8:
        axx.plot(profiles[i][:,2],color='b',lw=2,ls='--',label='Indian',zorder=2)
        axx.plot(profiles[i][:,3],color='g',lw=2,ls='--',label='Pacific',zorder=2)
        axx.plot(profiles[i][:,5],color='darkgoldenrod',lw=2,ls='--',label='Southern',zorder=2)
        R = rlspts[i][:,0].copy(); R[63:] = R[63:] - 360.
        a = pl.arange(0,212,30); b = pl.around(R[a],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        pl.ylim(-200,100)
    
    if i not in (0,3,6):
        axx.yaxis.set_major_formatter(pl.NullFormatter())
    if i == 3:
        axx.arrow(60,-10,0,-40,head_width=8,head_length=5,width=0.3,color='k')
        axx.annotate('westward/\n southward',xy=(62,-35),fontsize=13)
        axx.arrow(60,10,0,40,head_width=8,head_length=5,width=0.3,color='k')
        axx.annotate('eastward/\n northward',xy=(2,25),fontsize=13)
    
    pl.xlim(0,rlspts[i].shape[0])
    axx.grid(axis='y')
    axx.tick_params(axis='y',pad=5); axx.tick_params(axis='x',pad=7)
    pl.title(title[i],fontsize=18)

pl.tight_layout()
pl.subplots_adjust(wspace=0.08,hspace=0.30)
#pl.savefig(homedir+'/annmean_plots/panels/meanprfs_panels.png')
###############################################################################

###############################################################################
fig, ax = pl.subplots(3,3,figsize=(13,11))
y = pl.arange(len(fluxes[:6]))*5
h = pl.zeros_like(y);  h[:] = 2

for i in range(9):
    axx = pl.subplot(3,3,i+1)
    b1 = axx.barh(y-1,c1int[i,:],height=h,align='center',color='w',hatch='/',
        lw=2.,xerr=c1sem[i],error_kw=error_kw)
    b2 = axx.barh(y+1,c2int[i,:],height=h,align='center',color='w',hatch='x',
        lw=2.,xerr=c2sem[i],error_kw=error_kw)
    b1[0].set_edgecolor('grey'); b2[0].set_edgecolor('grey')
    b1[1].set_edgecolor('r'); b2[1].set_edgecolor('r')
    b1[2].set_edgecolor('b'); b2[2].set_edgecolor('b')
    b1[3].set_edgecolor('g'); b2[3].set_edgecolor('g')
    b1[4].set_edgecolor('deeppink'); b2[4].set_edgecolor('deeppink')
    b1[5].set_edgecolor('darkgoldenrod'); b2[5].set_edgecolor('darkgoldenrod')
    
    if i < 3:
        pl.xlim(-0.3,0.3)
    elif i in (3,4,5):
        pl.xlim(-0.1,0.1)
    elif i in (6,7):
        pl.xlim(-0.3,0.2)
    else:
        pl.xlim(-0.5,0.2)
    
    if i > 5:
        pl.xlabel('Sv',fontsize=16)
    
    if i in (0,3,6):
        axx.set_yticklabels(fluxes[:-2],fontsize=15)
    else:
        axx.yaxis.set_major_formatter(pl.NullFormatter())
    
    if i == 0:
        axx.legend(b1[:3],fluxes[:3],ncol=1,fontsize=14,columnspacing=0.5,loc=2)
    elif i == 1:
        axx.legend([b2[0],b1[0]],['CAT II','CAT I'],ncol=1,fontsize=14,loc=2)
    elif i == 7:
        axx.legend(b1[3:],fluxes[3:-2],ncol=1,fontsize=14,loc=2)
    
    if i == 4:
        axx.arrow(-0.01,25,-0.06,0,head_width=0.8,head_length=0.01,width=0.1,color='k')
        axx.annotate('westward/\n southward',xy=(-0.08,20),fontsize=13)
        axx.arrow(0.01,25,0.06,0,head_width=0.8,head_length=0.01,width=0.1,color='k')
        axx.annotate('eastward/\n northward',xy=(0.02,20),fontsize=13)
    
    pl.grid(axis='x'); axx.set_yticks(y); pl.xticks(fontsize=13)
    pl.title(title[i],fontsize=18)

pl.tight_layout()
#pl.savefig(homedir+'/annmean_plots/panels/catints_panels.png')
###############################################################################

###############################################################################
fig, ax = pl.subplots(3,3,figsize=(13,11))

P = []
for i in range(9):
    axx = pl.subplot(3,3,i+1)
    p0 = axx.plot(c2prfs[i][:,0],color='k',lw=2,label='total')
    if i == 0:
        p1 = axx.plot(c2prfs[i][:,1],color='r',lw=2,ls='--',label='Atlantic')
        axx.plot(c2prfs[i][:,3],color='g',lw=2,ls='--',label='Pacific')
        a = pl.arange(0,155,20); b = pl.around(rlspts[i][a,1],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        P.append(p0[0]); P.append(p1[0])
        pl.ylabel('kg/m/s',fontsize=16); pl.ylim(-250,150)
    elif i == 1:
        p2 = axx.plot(c2prfs[i][:,2],color='b',lw=2,ls='--',label='Indian')
        axx.plot(c2prfs[i][:,1],color='r',lw=2,ls='--',label='Atlantic')
        a = pl.arange(0,214,30); b = pl.around(rlspts[i][a,1],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        P.append(p2[0])#; axx.legend(P,fluxes[:3],loc=4)
        pl.xlabel('latitude',fontsize=14); pl.ylim(-250,150)
    elif i == 2:
        p3 = axx.plot(c2prfs[i][:,3],color='g',lw=2,ls='--',label='Pacific')
        axx.plot(c2prfs[i][:,2],color='b',lw=2,ls='--',label='Indian')
        a = pl.arange(0,165,20); b = pl.around(rlspts[i][a,1],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        P.append(p3[0]); pl.ylim(-250,150)
    elif i == 3:
        axx.plot(c2prfs[i][:,3],color='g',lw=2,ls='--',label='Pacific')
        axx.plot(c2prfs[i][:,1],color='r',lw=2,ls='--',label='Atlantic')
        p4 = axx.plot(c2prfs[i][:,4],color='deeppink',lw=2,ls='--',label='Arctic')
        a = pl.arange(0,197,30); b = pl.around(rlspts[i][a,0],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        P.append(p4[0])#; axx.annotate(title[i],(90,-55),fontsize=16)
        pl.ylabel('kg/m/s',fontsize=16); pl.ylim(-60,60)
    elif i == 4:
        axx.plot(c2prfs[i][:,2],color='b',lw=2,ls='--',label='Indian')
        axx.plot(c2prfs[i][:,1],color='r',lw=2,ls='--',label='Atlantic')
        axx.plot(c2prfs[i][:,4],color='deeppink',lw=2,ls='--',label='Arctic')
        a = pl.arange(0,103,20); b = pl.around(rlspts[i][a,0],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        pl.xlabel('longitude',fontsize=14); pl.ylim(-60,60)
    elif i == 5:
        axx.plot(c2prfs[i][:,3],color='g',lw=2,ls='--',label='Pacific')
        axx.plot(c2prfs[i][:,1],color='r',lw=2,ls='--',label='Atlantic')
        #axx.plot(profiles[i][:,2],color='b',lw=2,ls='--',label='Indian')
        axx.plot(c2prfs[i][:,4],color='deeppink',lw=2,ls='--',label='Arctic')
        a = pl.arange(0,220,30); b = pl.around(rlspts[i][a,0],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        #axx.legend(P[3:],fluxes[3:5],loc=2)
        pl.ylim(-60,60)
    elif i == 6:
        axx.plot(c2prfs[i][:,3],color='g',lw=2,ls='--',label='Pacific')
        axx.plot(c2prfs[i][:,1],color='r',lw=2,ls='--',label='Atlantic')
        p5 = axx.plot(c2prfs[i][:,5],color='darkgoldenrod',lw=2,ls='--',label='Southern')
        a = pl.arange(0,185,30); b = pl.around(rlspts[i][a,0],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        P.append(p5[0])
        pl.ylabel('kg/m/s',fontsize=16); pl.ylim(-200,100)
    elif i == 7:
        axx.plot(c2prfs[i][:,2],color='b',lw=2,ls='--',label='Indian')
        axx.plot(c2prfs[i][:,1],color='r',lw=2,ls='--',label='Atlantic')
        axx.plot(c2prfs[i][:,5],color='darkgoldenrod',lw=2,ls='--',label='Southern')
        a = pl.arange(0,184,30); b = pl.around(rlspts[i][a,0],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        pl.xlabel('longitude',fontsize=14); pl.ylim(-200,100)
        axx.legend(P,fluxes[:6],loc=3,ncol=2,columnspacing=0.5)
    elif i == 8:
        axx.plot(c2prfs[i][:,2],color='b',lw=2,ls='--',label='Indian')
        axx.plot(c2prfs[i][:,3],color='g',lw=2,ls='--',label='Pacific')
        axx.plot(c2prfs[i][:,5],color='darkgoldenrod',lw=2,ls='--',label='Southern')
        R = rlspts[i][:,0].copy(); R[63:] = R[63:] - 360.
        a = pl.arange(0,212,30); b = pl.around(R[a],decimals=0)
        axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
        pl.ylim(-200,100)

    if i not in (0,3,6):
        axx.yaxis.set_major_formatter(pl.NullFormatter())
    if i == 3:
        axx.arrow(60,-10,0,-40,head_width=8,head_length=5,width=0.3,color='k')
        axx.annotate('westward/\n southward',xy=(62,-35),fontsize=13)
        axx.arrow(60,10,0,40,head_width=8,head_length=5,width=0.3,color='k')
        axx.annotate('eastward/\n northward',xy=(2,25),fontsize=13)

    pl.xlim(0,rlspts[i].shape[0])
    pl.grid(axis='y')
    pl.title(title[i],fontsize=18)

pl.suptitle('CAT II $\overline{\mathbf{Q}}\cdot\mathbf{\hat{n}}$ profiles',
            fontsize=18,y=1.00)
pl.tight_layout()
pl.subplots_adjust(top=0.94,hspace=0.31)
#pl.savefig(homedir+'/annmean_plots/panels/cat2prfs_panels.png')
###############################################################################

###############################################################################
fig, ax = pl.subplots(3,3,figsize=(13,11))
levels = [-30,-20,-10,-5,-2,-1,#[-350,-300,-250,-200,-150,-100,-50,-25,-5,
          1,2,5,10,20,30]# 5,25,50,100,150,200,250,300,350]
         
eta = pl.linspace(0.95,0.15,17)
for i in range(9):
    axx = pl.subplot(3,3,i+1)
    cs = axx.contourf(csecs[i],cmap='seismic',norm=pl.Normalize(-30,30),
                 levels=levels,extend='both')
    axx.contour(csecs[i],colors='lightgrey',norm=pl.Normalize(-30,30),levels=[0],
                cmap=None,lw=0.7)
    cn = axx.contour(pres[i],colors='grey',norm=pl.Normalize(100,1000),cmap=None,
                levels=pl.linspace(100,1000,10),lw=0.7,linestyles='--')
    axx.clabel(cn,pl.linspace(100,1000,10).astype('int'),fontsize=10,fmt="%.0f")
    axx.set_yticklabels(eta[::2])
    pl.title(title[i],fontsize=16)
    
    if i in (0,2):
        a = pl.arange(0,rlspts[i].shape[0],20); b = pl.around(rlspts[i][a,1],decimals=0)
    elif i == 1:
        a = pl.arange(0,rlspts[i].shape[0],30); b = pl.around(rlspts[i][a,1],decimals=0)
    elif i in (3,5,6,7):
        a = pl.arange(0,rlspts[i].shape[0],30); b = pl.around(rlspts[i][a,0],decimals=0)
    elif i == 4:
        a = pl.arange(0,rlspts[i].shape[0],20); b = pl.around(rlspts[i][a,0],decimals=0)
    elif i == 8:
        R = rlspts[i][:,0].copy(); R[63:] = R[63:] - 360.
        a = pl.arange(0,rlspts[i].shape[0],30); b = pl.around(R[a],decimals=0)
    
    axx.set_xticks(a); axx.set_xticklabels(b.astype('int'),fontsize=13)
    
    if i in (0,3,6):
        pl.ylabel('$\eta$',fontsize=18)
    else:
        axx.yaxis.set_major_formatter(pl.NullFormatter())
    
    if i == 1:
        pl.xlabel('latitude',fontsize=16,labelpad=0.01)
    elif i in (4,7):
        pl.xlabel('longitude',fontsize=16,labelpad=0.01)

f = pl.gcf()
colax = f.add_axes([0.12,0.05,0.8,0.03])                   
cb = pl.colorbar(cs,orientation='horizontal',cax=colax)
cb.set_ticks(levels)
ticklabs = pl.asarray(levels)#; ticklabs = ticklabs/1000
cb.set_ticklabels(ticklabs)
cb.ax.tick_params(labelsize=13)
cb.update_ticks()#; cb.ax.set_aspect(0.09)
cb.set_label('kg/m/s',fontsize=20,labelpad=0.01)

pl.tight_layout(); pl.subplots_adjust(bottom=0.13,top=0.97,hspace=0.31)
#pl.savefig(homedir+'/annmean_plots/panels/csecs_panels.png')