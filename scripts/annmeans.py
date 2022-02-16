# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 15:20:10 2017

@author: np838619
"""

import pylab as pl
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
from trajfuncs import Haversine, Add360
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
shed = 'ara'; title, c = PlotTitle(shed)
if shed == 'afr':
    ec = 'r'
else:
    ec = 'k'

rlspts = pl.genfromtxt('/home/np838619/Watershed/shed_defs/ArA_traj_release_new.txt',skip_header=5)
#rlspts = Add360(rlspts)
#dist = pl.zeros([rlspts.shape[0]-1]); cdis = pl.zeros_like(dist)
#for i in range(rlspts.shape[0]-1):
#    dist[i] = Haversine(rlspts[i],rlspts[i+1])
#    cdis[i] = pl.sum(dist[:i+1])
#cdis = cdis/1000.

# annual means of partioned fluxes:
fluxpart = pl.genfromtxt(homedir+'annmeans_catchints_'+shed+'.csv',
                                                 skip_header=1,delimiter=',')
means = fluxpart[-1,1:]
catchstd = pl.std(fluxpart[:-1,1:],axis=0)
catchsem = stats.sem(fluxpart[:-1,1:],axis=0)
error_kw=dict(elinewidth=2,ecolor=ec,capsize=5,capthick=2)

fig,ax = pl.subplots()
fluxes = ['Total','Atlantic','Indian','Pacific','Arctic','Southern',
                          'Strat.','None']
y = pl.arange(len(fluxes))
ax.barh(y,means,align='center',color=c,xerr=catchsem,
        error_kw=error_kw)
ax.set_yticks(y); pl.xticks(fontsize=16)
ax.set_yticklabels(fluxes,fontsize=22)
pl.xlabel('Sv',fontsize=22)
pl.subplots_adjust(left=0.18,right=0.97)
ax.grid(axis='x'); pl.xlim(-0.8,0.4)
pl.title(title, fontsize=22)
#pl.savefig(homedir+'annmean_plots/fluxpart_'+shed+'_annmean'+'.png')


CATsplit = pl.genfromtxt(homedir+'annmeans_CATsplit_'+shed+'.csv',skip_header=2,
                                                                 delimiter=',')
CAT1means = CATsplit[-1,1::2]
CAT2means = CATsplit[-1,2::2]
C1sem = stats.sem(CATsplit[:-1,1::2],axis=0)#C1std = pl.std(CATsplit[:-1,1::2],axis=0)
C2sem = stats.sem(CATsplit[:-1,2::2],axis=0)#C2std = pl.std(CATsplit[:-1,2::2],axis=0)

fig,ax = pl.subplots()
y = pl.arange(len(fluxes[:6]))*5
h = pl.zeros_like(y);  h[:] = 2
ax.barh(y-1,CAT1means,height=h,align='center',edgecolor=c,
        linewidth=2.,color='w',hatch='/',label='CAT I',xerr=C1sem,error_kw=error_kw)
ax.barh(y+1,CAT2means,height=h,align='center',edgecolor=c,
        linewidth=2.,color='w',hatch='x',label='CAT II',xerr=C2sem,error_kw=error_kw)
pl.legend(loc=2,fontsize=22)
ax.set_yticks(y[:6]); pl.xticks(fontsize=16)
ax.set_yticklabels(fluxes[:6],fontsize=22)
ax.grid(axis='x'); pl.xlim(-0.5,0.2)
pl.xlabel('Sv',fontsize=22)
pl.subplots_adjust(left=0.18,right=0.97)
pl.title(title, fontsize=22)
#pl.savefig(homedir+'annmean_plots/catsplits_'+shed+'_annmean'+'.png')


sscyc = pl.genfromtxt(homedir+'season_cycs_'+shed+'.csv',delimiter=' ',skip_header=1)

fig,ax = pl.subplots()
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
pl.plot(sscyc[:,0],color='k',lw=2.,label='$F_{Tot}$')
pl.plot(sscyc[:,1],color='r',lw=2.,ls='--',label='$F_{Atl}$')
pl.plot(sscyc[:,2],color='b',lw=2.,ls='--',label='$F_{Ind}$')
pl.plot(sscyc[:,3],color='g',lw=2.,ls='--',label='$F_{Pac}$')
pl.plot(sscyc[:,4],color='deeppink',lw=2.,ls='--',label='$F_{Arc}$')
pl.plot(sscyc[:,5],color='darkgoldenrod',lw=2.,ls='--',label='$F_{Sou}$')
pl.plot(sscyc[:,6],color='indigo',lw=2.,ls='--',label='$F_{Str}$')
pl.plot(sscyc[:,7],color='saddlebrown',lw=2.,ls='--',label='$F_{None}$')
pl.xlim(0,11); pl.ylim(-0.9,0.9); pl.ylabel('Sv',fontsize=22)
ax.set_xticks(pl.arange(0,12)) ;ax.set_xticklabels(months,fontsize=18)
pl.yticks([-0.9,-0.6,-0.3,0,0.3,0.6,0.9],fontsize=16)
ax.tick_params(axis='x', which='major', pad=12); pl.grid(axis='y')
pl.legend(loc=9,ncol=4,columnspacing=0.5,fontsize=20)
pl.title(title,fontsize=22)
#pl.savefig(homedir+'sscyc_'+shed+'.png')

profiles = pl.genfromtxt(homedir+'profiles_annmean_'+shed+'.csv',delimiter=' ',
                             skip_header=1)

fig,ax = pl.subplots()
ax.plot(profiles[:,0],color='k',lw=2.,label='Total')
ax.plot(profiles[:,3],color='g',lw=2.,ls='--',label='Pacific')
ax.plot(profiles[:,1],color='r',lw=2.,ls='--',label='Atlantic')
#ax.plot(profiles[:,5],color='darkgoldenrod',lw=2.,ls='--',label='Southern')
a = pl.arange(0,197,30); b = pl.around(rlspts[a,1],decimals=0)
ax.set_xticks(a); ax.set_xticklabels(b.astype('int'),fontsize=16)
ax.tick_params(axis='x', which='major', pad=6)
pl.grid(axis='y'); pl.xlabel('latitude',fontsize=18,labelpad=-1)
pl.xlim(0,profiles.shape[0]); pl.ylabel('kg/m/s',fontsize=22)
pl.legend(loc=3,fontsize=20,ncol=1,columnspacing=0.5); pl.yticks(fontsize=14)
pl.title(title,fontsize=22)
#pl.savefig(homedir+'annmean_plots/annmean_prof_'+shed+'.png')


CAT1prof = pl.genfromtxt(homedir+'CAT1prf_annmean_'+shed+'.csv',delimiter=' ',
                             skip_header=1)
CAT2prof = pl.genfromtxt(homedir+'CAT2prf_annmean_'+shed+'.csv',delimiter=' ',
                             skip_header=1)

fig,ax = pl.subplots(1,2,figsize=(12,6))
ax1 = pl.subplot(121)
ax1.plot(CAT1prof[:,0],color='k',lw=2.,label='Total')
ax1.plot(CAT1prof[:,3],color='g',lw=2.,ls='--',label='Pacific')
ax1.plot(CAT1prof[:,1],color='r',lw=2.,ls='--',label='Atlantic')
ax1.plot(CAT1prof[:,4],color='deeppink',lw=2.,ls='--',label='Arctic')
pl.ylim(-200,120); pl.yticks(fontsize=14); pl.grid(axis='y')
pl.xlim(0,profiles.shape[0]); pl.ylabel('kg/m/s',fontsize=22,labelpad=-5)
ax1.set_xticks(a); ax1.set_xticklabels(b.astype('int'),fontsize=16)
ax1.tick_params(axis='x', which='major', pad=6)
pl.text(3,-190,'CAT I',fontsize=22)
ax2 = pl.subplot(122)
ax2.plot(CAT2prof[:,0],color='k',lw=2.,label='Total')
ax2.plot(CAT2prof[:,3],color='g',lw=2.,ls='--',label='Pacific')
ax2.plot(CAT2prof[:,1],color='r',lw=2.,ls='--',label='Atlantic')
ax2.plot(CAT2prof[:,4],color='deeppink',lw=2.,ls='--',label='Arctic')
pl.grid(axis='y')#; pl.axhline(y=0,ls='--',color='k')
pl.xlim(0,profiles.shape[0]); pl.ylim(-200,120)
ax2.yaxis.set_major_formatter(pl.NullFormatter())
ax2.set_xticks(a); ax2.set_xticklabels(b.astype('int'),fontsize=16)
ax2.tick_params(axis='x', which='major', pad=6)
pl.text(3,-190,'CAT II',fontsize=22)
pl.legend(loc=4,fontsize=20)
pl.subplots_adjust(wspace=0.08,left=0.10,right=0.92)
pl.suptitle(title,fontsize=22,y=0.95)
fig.text(0.5, 0.01, 'longitude', ha='center',fontsize=18)
p#l.savefig(homedir+'annmean_plots/annmean_catprofs_'+shed+'.png')


CAT1_lnd = pl.genfromtxt(homedir+'annmeans_CAT1land_'+shed+'.csv',delimiter=',',
                         skip_header=1)
CAT1_sea = pl.genfromtxt(homedir+'annmeans_CAT1sea_'+shed+'.csv',delimiter=',',
                         skip_header=1)
CAT1_bl = pl.genfromtxt(homedir+'annmeans_CAT1bl_'+shed+'.csv',delimiter=',',
                        skip_header=1)#/(10**9)

lnd_mns = CAT1_lnd[-1,1:]; lnd_sem = stats.sem(CAT1_lnd[:-1,1:],axis=0)
sea_mns = CAT1_sea[-1,1:]; sea_sem = stats.sem(CAT1_sea[:-1,1:],axis=0)
bl_mns = CAT1_bl[-1,1:]; bl_sem = stats.sem(CAT1_bl[:-1,1:],axis=0)

fig,ax = pl.subplots()
y = pl.arange(len(fluxes[:6]))*15
h = pl.zeros_like(y);  h[:] = 4
ax.barh(y-4,CAT1_bl[-1,1:],height=h,align='center',edgecolor=c,fill=False,
        lw=2.,label='BL at t=0',xerr=bl_sem,error_kw=error_kw)
ax.barh(y,CAT1_lnd[-1,1:],height=h,align='center',edgecolor=c,
        lw=2.,color='w',hatch='\\',label='land',xerr=lnd_sem,error_kw=error_kw)
ax.barh(y+4,CAT1_sea[-1,1:],height=h,align='center',edgecolor=c,
        lw=2.,color='w',hatch='|',label='sea',xerr=sea_sem,error_kw=error_kw)
pl.legend(loc=2,fontsize=22)
ax.set_yticks(y[:6]); pl.xticks(fontsize=15)
ax.set_yticklabels(fluxes[:6],fontsize=22)
ax.grid(axis='x')#; pl.xlim(-0.25,0.15)
pl.xlabel('Sv',fontsize=22,labelpad=0)
pl.subplots_adjust(left=0.18,right=0.96)
pl.title(title, fontsize=22); pl.xlim(-0.5,0.2)
#pl.savefig(homedir+'annmean_plots/annmean_CAT1_landsea_'+shed+'.png')