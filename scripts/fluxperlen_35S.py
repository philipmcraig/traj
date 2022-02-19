# -*- coding: utf-8 -*-
"""
Created on Mon May  7 15:19:06 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
from trajfuncs import BilinInterp, Haversine, Add360

def GetLoc(shed):
    
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
    
    return loc

pl.close('all')
clusdir = '/glusterfs/scenario/users/np838619/traj/'
homedir = '/home/np838619/Trajectory/fluxpart/'
sheddir = '/home/np838619/Watershed/shed_defs/'

sheds = ['soa','soi','sop']
profiles = []
rlspts = []
dl = []

for sd in range(len(sheds)):

    shed = sheds[sd]; loc = GetLoc(shed)
    rp = pl.genfromtxt(sheddir+loc+'_traj_release_new.txt',skip_header=5)
    rlspts.append(rp)
    PRFS = pl.genfromtxt(homedir+'profiles_annmean_'+shed+'.csv',delimiter=' ',skip_header=1)
    profiles.append(PRFS)
    
    D = pl.zeros([rp.shape[0]])
    for pt in range(1,rp.shape[0]-1): # exclude first & last release points
        # use Haversine formula to calculate dl
        D[pt] = 0.5*(Haversine(rp[pt-1,:],rp[pt,:]) + 
                        Haversine(rp[pt,:],rp[pt+1,:]))
    # use first 2 points & last 2 points for first & last dl
    D[0] = Haversine(rp[0,:],rp[1,:])
    D[-1] = Haversine(rp[-2,:],rp[-1,:])
    
    dl.append(D)


soa_len = pl.sum(dl[0][97:])
soi_len = pl.sum(dl[1][:111])
sop_len = pl.sum(dl[2][31:181])

soa_flx = profiles[0]*dl[0][:,None]
soi_flx = profiles[1]*dl[1][:,None]
sop_flx = profiles[2]*dl[2][:,None]

soa_int = pl.sum(soa_flx[97:],axis=0)/(10**9)
soi_int = pl.sum(soi_flx[:111],axis=0)/(10**9)
sop_int = pl.sum(sop_flx[31:181],axis=0)/(10**9)

soa_ave = soa_int/soa_len
soi_ave = soi_int/soi_len
sop_ave = sop_int/sop_len