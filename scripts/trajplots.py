# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 11:59:43 2017

@author: np838619
"""

import pylab as pl
from scipy import stats

def PlotTitle(shed):
    
    if shed == 'amr':
        title = 'Americas'; c = 'b'
    elif shed == 'afr':
        title = 'Africa'; c = 'k'
    elif shed == 'eaa':
        title = 'South-East Asia'; c = 'g'
    elif shed == 'ara':
        title = 'Arctic Atlantic'; c = 'r'
    elif shed == 'ari':
        title = 'Arctic Indian'; c = 'r'
    elif shed == 'arp':
        title = 'Arctic Pacific'; c = 'r'
    elif shed == 'soa':
        title = 'Southern Atlantic'; c = 'darkgoldenrod'
    elif shed == 'soi':
        title = 'Southern Indian'; c = 'darkgoldenrod'
    elif shed == 'sop':
        title = 'Southern Pacific'; c = 'darkgoldenrod'
    
    return title, c

pl.close('all')
clusdir = '/glusterfs/scenario/users/np838619/traj/'
homedir = '/home/np838619/Trajectory/fluxpart/'
year = '2014'
shed = 'sop'; title, c = PlotTitle(shed)
if shed == 'afr':
    ec = 'r'
else:
    ec = 'k'
datadir = homedir + year + '/'+ shed + '/'

rlspts = pl.genfromtxt('/home/np838619/Watershed/shed_defs/AfMe_traj_release_new.txt',skip_header=5)

# read in some files here:
    # total & partioned moisture fluxes for each trajectory release in Sv
catchints = pl.genfromtxt(datadir+'catchints.csv',delimiter=' ',skip_header=1)
    # mean profiles of total & partioned moisture fluxes in kg/m/s
meanprofs = pl.genfromtxt(datadir+'meanprofs.csv',delimiter=' ',skip_header=1)
    # total & partitioned CAT1/2 moisture fluxes for each release in Sv
CAT1ints = pl.genfromtxt(datadir+'CAT1ints.csv',delimiter=' ',skip_header=1)
CAT2ints = pl.genfromtxt(datadir+'CAT2ints.csv',delimiter=' ',skip_header=1)
    # mean profiles of total & partioned moisture fluxes in kg/m/s
CAT1profs = pl.genfromtxt(datadir+'CAT1profs.csv',delimiter=' ',skip_header=1)
CAT2profs = pl.genfromtxt(datadir+'CAT2profs.csv',delimiter=' ',skip_header=1)
    # total & partitioned CAT1 land/sea moisture fluxes for each release in Sv
CAT1lnd = pl.genfromtxt(datadir+'CAT1land.csv',delimiter=' ')
CAT1sea = pl.genfromtxt(datadir+'CAT1sea.csv',delimiter=' ')
    # total & partitioned BL at t=0 moisture fluxes for each release in Sv
CAT1bl = pl.genfromtxt(datadir+'CAT1_BL.csv',delimiter=' ')/(10**9)

#pl.figure(1) # total & partioned moisture fluxes for each release


#pl.figure(2) # mean profiles of total & partitioned moisture fluxes

#pl.figure(3) # total & partioned CAT1/2 moisture fluxes for each release

#pl.figure(4) # mean profiles of total & partitioned CAT1/2 moisture fluxes

fig, ax = pl.subplots()
fluxes = ['Total','Atlantic','Indian','Pacific','Arctic','Southern',
                          'Strat.','None']
catchsem = stats.sem(catchints,axis=0)
error_kw=dict(elinewidth=2,ecolor=ec,capsize=5,capthick=2)
y = pl.arange(len(fluxes))
ax.barh(y,pl.mean(catchints[:,:],axis=0),align='center',color=c,xerr=catchsem,
        error_kw=error_kw)
#pl.axvline(x=0,ls='--',color='k')
ax.set_yticks(y); pl.xticks(fontsize=18)
ax.set_yticklabels(fluxes,fontsize=22)
pl.xlabel('Sv',fontsize=22,labelpad=2)
pl.subplots_adjust(left=0.18,right=0.97)
ax.grid(axis='x'); pl.xlim(-0.8,0.4)
#ax.set_xticks(pl.linspace(-0.1,0.1,11))
ax.tick_params(axis='x',labelsize=16)
#ax.set_xticklabels([-0.6,-0.3,0,0.3])
pl.title(title+' '+year, fontsize=22)
#pl.savefig(homedir+year+'/'+shed+'/fluxpart_'+shed+year+'.png')

fig,ax = pl.subplots()
y = pl.arange(len(fluxes[:6]))*5
w1 = pl.mean(CAT1ints[:],axis=0)
w2 = pl.mean(CAT2ints[:],axis=0)
h = pl.zeros_like(y);  h[:] = 2
C1sem = stats.sem(CAT1ints,axis=0)
C2sem = stats.sem(CAT2ints,axis=0)
ax.barh(y-1,pl.mean(CAT1ints[:,:],axis=0),height=h,align='center',edgecolor=c,
        linewidth=2.,color='w',hatch='/',label='CAT I',xerr=C1sem,error_kw=error_kw)
ax.barh(y+1,pl.mean(CAT2ints[:,:],axis=0),height=h,align='center',edgecolor=c,
        linewidth=2.,color='w',hatch='x',label='CAT II',xerr=C2sem,error_kw=error_kw)
pl.legend(loc=2,fontsize=22)
ax.set_yticks(y[:6]); pl.xticks(fontsize=15)
ax.set_yticklabels(fluxes[:6],fontsize=22)
ax.grid(axis='x'); pl.xlim(-0.8,0.4)
pl.xlabel('Sv',fontsize=22,labelpad=2)
pl.subplots_adjust(left=0.18,right=0.97)
pl.title(title+' '+year, fontsize=22)
#pl.savefig(homedir+year+'/'+shed+'/catsplits_'+shed+year+'.png')

fig,ax = pl.subplots()
pl.plot(meanprofs[:,0],color='k',lw=2.,label='Total')
pl.plot(meanprofs[:,2],color='b',lw=2.,ls='--',label='Indian')
pl.plot(meanprofs[:,1],color='r',lw=2.,ls='--',label='Atlantic')
#pl.plot(meanprofs[:,5],color='darkgoldenrod',lw=2.,ls='--',label='Southern')
pl.grid(axis='y'); pl.yticks(fontsize=15)
pl.xlim(0,meanprofs.shape[0])
a = pl.arange(0,214,30); b = pl.around(rlspts[a,1],decimals=0)
ax.set_xticks(a); ax.set_xticklabels(b.astype('int'),fontsize=16)
ax.tick_params(axis='x', which='major', pad=6)
pl.legend(loc=3,fontsize=22); pl.xlabel('latitude',fontsize=18,labelpad=2)
pl.ylim(-250,150); pl.ylabel('kg/m/s',fontsize=22,labelpad=5)
pl.subplots_adjust(right=0.95,left=0.17)
pl.title(title+' '+year,fontsize=22)
#pl.savefig(homedir+year+'/'+shed+'/fluxprofs_'+shed+year+'.png')


fig, ax = pl.subplots(1,2,figsize=(12,6))
ax1 = pl.subplot(121)
pl.plot(CAT1profs[:,0],color='k',lw=2.,label='Total')
pl.plot(CAT1profs[:,2],color='b',lw=2.,ls='--',label='Indian')
pl.plot(CAT1profs[:,1],color='r',lw=2.,ls='--',label='Atlantic')
#pl.plot(CAT1profs[:,5],color='darkgoldenrod',lw=2.,ls='--',label='Southern')
pl.grid(axis='y')
pl.yticks(pl.linspace(-200,120,9),fontsize=14); pl.yticks(fontsize=15)
pl.ylabel('kg/m/s',fontsize=22)
pl.xlim(0,CAT1profs.shape[0]); pl.ylim(-200,120)
ax1.set_xticks(a); ax1.set_xticklabels(b.astype('int'),fontsize=16)
ax1.tick_params(axis='x', which='major', pad=6)
pl.text(3,-190,'CAT I',fontsize=22)
########################################################
ax2 = pl.subplot(122)
pl.plot(CAT2profs[:,0],color='k',lw=2.,label='Total')
pl.plot(CAT2profs[:,2],color='b',lw=2.,ls='--',label='Indian')
pl.plot(CAT2profs[:,1],color='r',lw=2.,ls='--',label='Atlantic')
#pl.plot(CAT2profs[:,5],color='darkgoldenrod',lw=2.,ls='--',label='Southern')
pl.grid(axis='y'); pl.yticks(pl.linspace(-200,120,9))
ax2.yaxis.set_major_formatter(pl.NullFormatter())
pl.xlim(0,CAT1profs.shape[0]); pl.ylim(-200,120)
ax2.set_xticks(a); ax2.set_xticklabels(b.astype('int'),fontsize=16)
ax2.tick_params(axis='x', which='major', pad=6)
pl.text(3,-190,'CAT II',fontsize=22)
pl.legend(loc=4,fontsize=22)
pl.subplots_adjust(wspace=0.08,left=0.10,right=0.92)
pl.suptitle(title+' '+year,fontsize=22,y=0.95)
fig.text(0.5, 0.01, 'latitude', ha='center',fontsize=18)
#pl.savefig(homedir+year+'/'+shed+'/catprofs_'+shed+year+'.png')

fig, ax = pl.subplots()
y = pl.arange(len(fluxes[:6]))*15
h = pl.zeros_like(y);  h[:] = 4
lnd_sem = stats.sem(CAT1lnd,axis=0); sea_sem = stats.sem(CAT1sea,axis=0)
bl_sem = stats.sem(CAT1bl,axis=0)
ax.barh(y-4,pl.mean(CAT1bl,axis=0),height=h,align='center',edgecolor=c,fill=False,
        lw=2.,label='BL at t=0',xerr=bl_sem,error_kw=error_kw)
ax.barh(y,pl.mean(CAT1lnd,axis=0),height=h,align='center',edgecolor=c,
        linewidth=2.,color='w',hatch='\\',label='land',xerr=lnd_sem,error_kw=error_kw)
ax.barh(y+4,pl.mean(CAT1sea,axis=0),height=h,align='center',edgecolor=c,
        linewidth=2.,color='w',hatch='|',label='sea',xerr=sea_sem,error_kw=error_kw)
pl.legend(loc=2,fontsize=22)
ax.set_yticks(y[:6]); pl.xticks(fontsize=16)
ax.set_yticklabels(fluxes[:6],fontsize=22)
ax.grid(axis='x'); pl.xlim(-0.8,0.4)
pl.xlabel('Sv',fontsize=22,labelpad=2)
pl.subplots_adjust(left=0.18,right=0.97)
pl.title(title+' '+year, fontsize=22)#; pl.xlim(-0.15,0.15)
#pl.savefig(homedir+year+'/'+shed+'/CAT1_landsea_'+shed+year+'.png')#