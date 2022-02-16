# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 09:54:39 2017

@author: np838619
"""

import pylab as pl

pl.close('all')
shed = 'sop'
time = 'jul2010'; time2 = 'Jul10'
#datadir = '/net/glusterfs/scenario/users/np838619/tjexp/'+shed+'/summ/TD/'
trajdir = '/home/np838619/Trajectory/fluxevol/'

tj_data = pl.genfromtxt(trajdir+shed+'/'+time+'_tj_'+shed+'.csv',skip_header=1)

tj_tot = tj_data[:,0] # total magnitude of moisture flux
tj_CAT1 = tj_data[:,1] # magnitude of CATI moisture flux
tj_CAT2 = tj_data[:,2] # magnitude of CATII moisture flux
tj_strat = tj_data[:,3] # magntiude of stratospheric moisture flux

fx_data = pl.genfromtxt(trajdir+shed+'/'+time+'_fx_'+shed+'.csv',skip_header=1)

fx_tot = fx_data[:,0] # total magnitude of moisture flux
fx_CAT1 = fx_data[:,1] # magnitude of CATI moisture flux
fx_CAT2 = fx_data[:,2] # magnitude of CATII moisture flux
fx_strat = fx_data[:,3] # magntiude of stratospheric moisture flux

steps = pl.linspace(0,120,121)*6/24

fig,ax = pl.subplots(1,2,figsize=(14,6))

lw = 3
ax1 = pl.subplot(121)
ax1.plot(steps,tj_tot,color='r',lw=lw,label='total')
ax1.plot(steps,tj_CAT1,color='b',ls='-',lw=lw,label='CAT I')
ax1.plot(steps,tj_CAT2,color='g',ls='-',lw=lw,label='CAT II')
ax1.plot(steps,tj_strat,color='purple',ls='-',lw=lw,label='strat.')
pl.yticks(pl.linspace(0,100,11),fontsize=18)
pl.xticks(fontsize=18)
ax1.legend(loc=2,fontsize=19,ncol=2,columnspacing=0.6)
pl.grid(axis='y'); pl.ylim(0,100)
pl.ylabel('%',fontsize=22,labelpad=0.5)
pl.xlabel('days',fontsize=22,labelpad=0.5)
ax1.annotate('(a) trajectories',(0.33,0.84),xycoords='figure fraction',fontsize=18)

ax2 = pl.subplot(122)
#ax2.plot(steps,tj_tot,color='r',lw=2,label='total traj.')
ax2.plot(steps,fx_tot,color='r',lw=lw,label='total')
ax2.plot(steps,fx_CAT1,color='b',ls='-',lw=lw,label='CAT I')
ax2.plot(steps,fx_CAT2,color='g',ls='-',lw=lw,label='CAT II')
ax2.plot(steps,fx_strat,color='purple',ls='-',lw=lw,label='strat.')
pl.yticks(pl.linspace(0,100,11),fontsize=18)
pl.xticks(fontsize=18)
#ax2.legend(loc=8,ncol=2,fontsize=19,columnspacing=0.5)
pl.grid(axis='y')
ax2.yaxis.set_major_formatter(pl.NullFormatter())
pl.xlabel('days',fontsize=22,labelpad=0.5)
ax2.annotate('(b) flux magnitude',(0.54,0.84),xycoords='figure fraction',fontsize=18)

pl.suptitle('January 2011',fontsize=22)
pl.subplots_adjust(left=0.07,right=0.93,wspace=0.15)
#pl.savefig(trajdir+shed+'/traj_flux_evol_'+shed+'201101.png')

posneg = pl.genfromtxt(trajdir+shed+'/'+'posneg_'+time2+'.csv')
pl.figure(2)
pl.plot(steps,posneg[:,0]/(10**9),lw=lw,label='+ve')
pl.plot(steps,posneg[:,1]/(10**9),lw=lw,label='-ve')
pl.plot(steps,pl.sum(posneg,axis=1)/(10**9),lw=lw,label='net',color='k')
pl.legend(loc=0,fontsize=18,ncol=3); pl.grid(axis='y')
pl.xlabel('days',fontsize=22,labelpad=-0.5); pl.xticks(fontsize=18)
pl.ylabel('Sv',fontsize=22,labelpad=-5); pl.yticks(fontsize=18)
pl.title('July 2010',fontsize=22)
#pl.savefig(trajdir+shed+'/'+shed+'_posneg_201007.png')