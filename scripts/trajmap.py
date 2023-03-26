# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 18:52:15 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
#from mpl_toolkits.basemap import Basemap
import cartopy
import cartopy.crs as ccrs
from netCDF4 import Dataset

def ShedPlot(sheddir):
    """
    """
    m = Basemap(projection='robin',lon_0=-180.,resolution='l')
    m.drawcoastlines(linewidth=0.4,color='lightgray',zorder=1)
    m.drawcountries(color='lightgray',zorder=1)
    
    #sheddir = '/home/np838619/Watershed/shed_defs/'
    filex = '_traj_release_new.txt'
    
    Am = pl.genfromtxt(sheddir + 'Am' + filex,skip_header=5)
    AfMe = pl.genfromtxt(sheddir + 'AfMe' + filex,skip_header=5)
    EAA = pl.genfromtxt(sheddir + 'EAA' + filex,skip_header=5)
    Ar = pl.genfromtxt(sheddir + 'Ar' + filex,skip_header=5)
    SO = pl.genfromtxt(sheddir + 'SO' + filex,skip_header=5)
    NAs = pl.genfromtxt(sheddir + 'NAs' + filex,skip_header=5)
    
    lw = 3
    m.plot(Am[:,0],Am[:,1],latlon=True,color='b',linewidth=lw)
    m.plot(AfMe[:,0],AfMe[:,1],latlon=True,color='k',linewidth=lw)
    m.plot(EAA[:,0],EAA[:,1],latlon=True,color='g',linewidth=lw)
    m.plot(NAs[:,0],NAs[:,1],latlon=True,color='maroon',ls='--',linewidth=2)
    
    F = pl.zeros_like(SO); F[:,0], F[:,1] = m(SO[:,0],SO[:,1])
    m.plot(F[:246,0], F[:246,1],color='darkgoldenrod',linewidth=lw)
    m.plot(F[246:555,0], F[246:555,1],color='darkgoldenrod',linewidth=lw)
    m.plot(F[555:,0],F[555:,1],color='darkgoldenrod',linewidth=lw)
    
    m.plot(Ar[:157,0],Ar[:157,1],latlon=True,linewidth=lw,color='r',zorder=2)
    m.plot(Ar[157:,0],Ar[157:,1],latlon=True,linewidth=lw,color='r',zorder=2)
    
    isect1 = pl.where(Am[0]==Ar); isect1=(Ar[isect1[0][0],0],Ar[isect1[0][0],1])
    isect2 = pl.where(AfMe[0]==Ar); isect2=(Ar[isect2[0][0],0],Ar[isect2[0][0],1])
    isect3 = pl.where(EAA[0]==Ar); isect3=(Ar[isect3[0][0],0],Ar[isect3[0][0],1])
    m.plot(isect1[0],isect1[1],'r.',markersize=13,mew=3,latlon=True)
    m.plot(isect2[0],isect2[1],'r.',markersize=13,mew=3,latlon=True)
    m.plot(isect3[0],isect3[1],'r.',markersize=13,mew=3,latlon=True)
    
    isect4 = pl.where(Am[-1]==SO); isect4=(SO[isect4[0][0],0],SO[isect4[0][0],1])
    isect5 = pl.where(AfMe[-1]==SO); isect5=(SO[isect5[0][0],0],SO[isect5[0][0],1])
    isect6 = pl.where(EAA[-1]==SO); isect6=(SO[isect6[0][0],0],SO[isect6[0][0],1])
    m.plot(isect4[0],isect4[1],color='darkgoldenrod',marker='.',markersize=13,mew=3,latlon=True)
    m.plot(isect5[0],isect5[1],color='darkgoldenrod',marker='.',markersize=13,mew=3,latlon=True)
    m.plot(isect6[0],isect6[1],color='darkgoldenrod',marker='.',markersize=13,mew=3,latlon=True)
    
    return m


def DrawArrow(m,catchints,shed_lon,shed_lat,shedx,shedy,dx,dy,hw,hl,w,size,cols):
    """
    """
    #cols = ['grey','r','b','g','deeppink','darkgoldenrod','skyblue','k']
    lw = 2
    x1,y1 = m(shed_lon,shed_lat)
    L = []
    
    for i in range(catchints.shape[0]):
        if catchints[i] == 0.0:
            pass
        else:
            if catchints[i] > 0:
                A = 1
            elif catchints[i] < 0:
                A = -1
            if i == catchints.shape[0]-1 and cols[i] == 'k':
                fc = 'white'; fill = False; lw = 2
            else:
                fc = cols[i]; fill = True
            ar = pl.arrow(x1[i],y1[i],A*dx,A*dy,fc=fc, ec=cols[i], lw = lw, 
                     head_width=hw,head_length=hl,width=w,zorder=3,
                                         fill=fill,head_starts_at_zero=False)
            pl.annotate('{:.2f}'.format(A*catchints[i]),xy=(shedx[i],shedy[i]),
                         xycoords='axes fraction',color=cols[i],size=size,zorder=3)
            L.append(ar)
    
    return L

def FluxArrows(m,catchints):
    """
    """
    cols = ['grey','r','b','g','deeppink','darkgoldenrod','skyblue','k']
    cols = pl.array(cols)
    lw = 2
    catchints = pl.around(catchints,2)
    
    amr_lon = [-102,-98,-115,-95,-95,-75,-85,-83]#[50,45,40,35,30,25,20,15]
    amr_lat = [30,22,27.4,14,14.9,8.6,2.3,6]#pl.linspace(40,-4,8)
    amr_x = [0.635,0.65,0.72,0.775,0.78,0.72,0.81,0.81]
    amr_y= [0.68,0.625,0.66,0.58,0.58,0.54,0.51,0.53]#pl.linspace(0.75,0.45,8)
    afr_lon = [40,30,43,30,43,30,28,10]
    afr_lat = [30,19,8,-3,4.9,-1.4,-7.7,-14]#pl.linspace(30,-14,8)
    afr_x = [0.055,0.13,0.05,0.125,0.125,0.05,0.125,0.065]
    afr_y = [0.68,0.61,0.54,0.475,0.52,0.485,0.445,0.385]#pl.linspace(0.7,0.4,8)
    eaa_lon = [93,93,92,109,95,96,135,137]
    eaa_lat = [28,17,6,-5,2.9,-3.4,-9.7,-14]#pl.linspace(28,-16,8)
    eaa_x = [0.31,0.305,0.30,0.23,0.305,0.31,0.42,0.425]
    eaa_y = [0.67,0.60,0.53,0.465,0.51,0.47,0.43,0.41]#pl.arange(0.67,0.36,-0.04)
    ara_lon1 = [-90,-58]; ara_lat1 = [70,80]#,70,70,70,70]
    ara_x1 = [0.675,0.705]; ara_y1 = [0.86,0.90]#,0.81,0.81,0.81,0.81]
    ara_lon2 = [26,29,41]; ara_lat2 = [67,56,47]#,70,70,70,70]
    ara_x2 = [0.21,0.18,0.185]; ara_y2 = [0.895,0.835,0.78]
    ari_lon = [49.6,59.55,69.5,68,79.45,77.5,83.8,89.4]#pl.linspace(49,89,8)
    ari_lat = [31,24,31,28,37,33,28,28]
    ari_x = [0.155,0.165,0.20,0.19,0.23,0.22,0.235,0.25]#pl.linspace(0.14,0.3,8)
    ari_y = [0.74,0.70,0.74,0.725,0.66,0.755,0.725,0.725]
    arp_lon = [140,160,155,180,200,200,215,220]#pl.linspace(125,230,8)
    arp_lat = [55,66,65,63,68,68,72,71]
    arp_x = [0.40,0.45,0.44,0.49,0.53,0.53,0.56,0.57]
    arp_y = [0.89,0.835,0.83,0.935,0.84,0.84,0.87,0.86]
    soa_lon = [-40,-30,-35,-20,-30,-10,-13,-10]#pl.linspace(-62,-5,8)
    soa_lat = [-32,-32,-32,-37,-32,-37,-37,-37]
    soa_x = [0.86,0.89,0.84,0.905,0.91,0.93,0.925,0.95]#pl.linspace(0.8,0.95,8)
    soa_y = [0.24,0.24,0.24,0.32,0.24,0.32,0.32,0.32]
    soi_lon = [40,57.5,75,69.3,80.7,92.5,103.6,110]#pl.linspace(35,115,8)
    soi_lat = [-32,-32,-32,-32,-32,-37,-32,-32]
    soi_x = [0.12,0.165,0.21,0.195,0.23,0.26,0.285,0.305]
    soi_y = [0.24,0.24,0.24,0.24,0.24,0.32,0.24,0.24]
    sop_lon = [180,184.3,200,220,227.1,240,255.7,260]#pl.linspace(170,270,8)
    sop_lat = [-32,-32,-37,-32,-37,-37,-32,-32]
    sop_x = [0.49,0.50,0.54,0.60,0.61,0.65,0.69,0.70]
    sop_y = [0.24,0.24,0.32,0.24,0.32,0.32,0.24,0.24]

    L1=DrawArrow(m,catchints[0],amr_lon,amr_lat,amr_x,amr_y,1000000,0,350000,
              350000,8000,20,cols)
    L2=DrawArrow(m,catchints[1],afr_lon,afr_lat,afr_x,afr_y,1000000,0,350000,
              350000,8000,20,cols)
    L3=DrawArrow(m,catchints[2],eaa_lon,eaa_lat,eaa_x,eaa_y,1000000,0,350000,
              350000,8000,20,cols)
    L4=DrawArrow(m,catchints[3,3:5],ara_lon1,ara_lat1,ara_x1,ara_y1,0,550000,
              250000,250000,7000,15,cols[3:5])
    L5=DrawArrow(m,catchints[3,[0,1,-1]],ara_lon2,ara_lat2,ara_x2,ara_y2,550000,0,
              250000,250000,7000,15,cols[[0,1,-1]])
    L6=DrawArrow(m,catchints[4],ari_lon,ari_lat,ari_x,ari_y,0,550000,
              250000,250000,7000,15,cols)
    L7=DrawArrow(m,catchints[5],arp_lon,arp_lat,arp_x,arp_y,0,550000,
              250000,250000,7000,15,cols)
    L8=DrawArrow(m,catchints[6],soa_lon,soa_lat,soa_x,soa_y,0,550000,
              250000,250000,7000,15,cols)
    L9=DrawArrow(m,catchints[7],soi_lon,soi_lat,soi_x,soi_y,0,550000,
              250000,250000,7000,15,cols)
    L10=DrawArrow(m,catchints[8],sop_lon,sop_lat,sop_x,sop_y,0,550000,
              250000,250000,7000,15,cols)
    
    L = [L1[0],L1[1],L2[2],L2[3],L4[1],L9[-2],L1[-1]]
    labs = ['total','Atlantic','Indian','Pacific','Arctic','Southern',
            'no origin']
    leg = pl.legend(L,labs,loc=3,ncol=4,columnspacing=0.7,fontsize=18,
                    title='arrows')
    pl.setp(leg.get_title(),fontsize='20')
    
    return None

def PmEboxes(m,catchints):
    """
    """
    atl = catchints[0,1:] - catchints[1,1:] - catchints[3,1:] + catchints[6,1:]
    ind = catchints[1,1:] - catchints[2,1:] - catchints[4,1:] + catchints[7,1:]
    pac = catchints[2,1:] - catchints[0,1:] - catchints[5,1:] + catchints[8,1:]
    arc = catchints[3,1:] + catchints[4,1:] + catchints[5,1:]
    sou = -catchints[6,1:] - catchints[7,1:] - catchints[8,1:]
    
    atl_sum = pl.sum(atl); ind_sum = pl.sum(ind); pac_sum = pl.sum(pac)
    arc_sum = pl.sum(arc); sou_sum = pl.sum(sou)
    
    bdict = {'facecolor':'white', 'alpha':0.5, 'pad':5, 'fill':'True',
             'color':'w', 'ec':'k'}
    
    a,b = m(180,15)
    pl.text(a,b,'{:.2f}'.format(pac_sum),bbox=bdict,fontsize=35)
    a,b = m(310,15)
    pl.text(a,b,'{:.2f}'.format(atl_sum),bbox=bdict,fontsize=35)
    a,b = m(55,-13)
    pl.text(a,b,'{:.2f}'.format(ind_sum),bbox=bdict,fontsize=35)
    a,b = m(100,70)
    pl.text(a,b,'{:.2f}'.format(arc_sum),bbox=bdict,fontsize=35)
    a,b = m(180,-55)
    pl.text(a,b,'{:.2f}'.format(sou_sum),bbox=bdict,fontsize=35)
    
    return None

def InOutPlot(m,IN_OUT):
    """
    """
    IN_OUT = pl.absolute(IN_OUT)
    
    a,b = m(307,3)
    pl.text(a,b,'  IN\n'+'{:.2f}'.format(IN_OUT[0,0]),fontsize=22)
    a,b = m(322,3)
    pl.text(a,b,'OUT\n'+'{:.2f}'.format(IN_OUT[0,1]),fontsize=22)
    
    a,b = m(50,-25)
    pl.text(a,b,'  IN\n'+'{:.2f}'.format(IN_OUT[1,0]),fontsize=22)
    a,b = m(65,-25)
    pl.text(a,b,'OUT\n'+'{:.2f}'.format(IN_OUT[0,1]),fontsize=22)
    
    a,b = m(177,3)
    pl.text(a,b,'  IN\n'+'{:.2f}'.format(IN_OUT[2,0]),fontsize=22)
    a,b = m(192,3)
    pl.text(a,b,'OUT\n'+'{:.2f}'.format(IN_OUT[2,1]),fontsize=22)
    
    a,b = m(65,67)
    pl.text(a,b,'  IN\n'+'{:.2f}'.format(IN_OUT[3,0]),fontsize=22)
    a,b = m(83,67)
    pl.text(a,b,'OUT\n'+'{:.2f}'.format(IN_OUT[3,1]),fontsize=22)
    
    a,b = m(174,-68)
    pl.text(a,b,'  IN\n'+'{:.2f}'.format(IN_OUT[4,0]),fontsize=22)
    a,b = m(194,-68)
    pl.text(a,b,'OUT\n'+'{:.2f}'.format(IN_OUT[4,1]),fontsize=22)
    
    return None

exec(open('C:/Users/phili/GitHub/traj/scripts/trajfuncs.py').read())

pl.close('all')
clusdir = '/glusterfs/scenario/users/np838619/traj/'
#homedir = '/home/np838619/Trajectory/fluxpart/'
sheddir = 'C:/Users/phili/GitHub/watershed/shed_defs/'#'/home/np838619/Watershed/shed_defs/'

pl.figure(figsize=(23.5,13))
m = ShedPlot(sheddir); pl.tight_layout()

sheds = ['amr','afr','eaa','ara','ari','arp','soa','soi','sop']
catchints = pl.zeros([len(sheds),8])

# loop over sheds:
for sd in range(len(sheds)):

    shed = sheds[sd]
    datadir = 'C:/Users/phili/GitHub/traj/fluxpart/'#homedir
    INTS = pl.genfromtxt(datadir+'annmeans_catchints_'+shed+'.csv',delimiter=',',skip_header=1)
    catchints[sd] = INTS[-1,1:]

atl = catchints[0,1:] - catchints[1,1:] - catchints[3,1:] + catchints[6,1:]
ind = catchints[1,1:] - catchints[2,1:] - catchints[4,1:] + catchints[7,1:]
pac = catchints[2,1:] - catchints[0,1:] - catchints[5,1:] + catchints[8,1:]
arc = catchints[3,1:] + catchints[4,1:] + catchints[5,1:]
sou = -catchints[6,1:] - catchints[7,1:] - catchints[8,1:]

basins = pl.array([atl,ind,pac,arc,sou])

IN_OUT = pl.zeros([5,2])
for i in range(5):
    for j in range(7):
        if basins[i,j] > 0:
            IN_OUT[i,0] = IN_OUT[i,0] + basins[i,j]
        elif basins[i,j] < 0:
            IN_OUT[i,1] = IN_OUT[i,1] + basins[i,j]

FluxArrows(m,catchints)
PmEboxes(m,catchints)
InOutPlot(m,IN_OUT)
pl.tight_layout()
#pl.savefig(homedir+'annmean_plots/panels/fluxpart_map_new.png')

atl_area = 0.746e14
ind_area = 0.450e14
pac_area = 1.405e14

# IN/OUT per unit area:
atl_pua = IN_OUT[0]/atl_area; print atl_pua
ind_pua = IN_OUT[1]/ind_area; print ind_pua
pac_pua = IN_OUT[2]/pac_area; print pac_pua


rlspts1 = pl.genfromtxt(sheddir+'Am_traj_release_new.txt',skip_header=5)
rlspts2 = pl.genfromtxt(sheddir+'AfMe_traj_release_new.txt',skip_header=5)
rlspts3 = pl.genfromtxt(sheddir+'EAA_traj_release_new.txt',skip_header=5)
rlspts4 = pl.genfromtxt(sheddir+'ArA_traj_release_new.txt',skip_header=5)
rlspts5 = pl.genfromtxt(sheddir+'ArI_traj_release_new.txt',skip_header=5)
rlspts6 = pl.genfromtxt(sheddir+'ArP_traj_release_new.txt',skip_header=5)
rlspts7 = pl.genfromtxt(sheddir+'SOA_traj_release_new.txt',skip_header=5)
rlspts8 = pl.genfromtxt(sheddir+'SOI_traj_release_new.txt',skip_header=5)
rlspts9 = pl.genfromtxt(sheddir+'SOP_traj_release_new.txt',skip_header=5)

atl_bnd = CatchLine(rlspts1,rlspts4,rlspts2,rlspts7)
ind_bnd = CatchLine(rlspts2,rlspts5,rlspts3,rlspts8)
pac_bnd = CatchLine(rlspts3,rlspts6,rlspts1,rlspts9)

atlbnd_len = pl.zeros([atl_bnd.shape[0]])
for i in range(atl_bnd.shape[0]-1):
    atlbnd_len[i] = Haversine(atl_bnd[i],atl_bnd[i+1])
atlbnd_len[-1] = Haversine(atl_bnd[-1],atl_bnd[0])

indbnd_len = pl.zeros([ind_bnd.shape[0]])
for i in range(ind_bnd.shape[0]-1):
    indbnd_len[i] = Haversine(ind_bnd[i],ind_bnd[i+1])
indbnd_len[-1] = Haversine(ind_bnd[-1],ind_bnd[0])

pacbnd_len = pl.zeros([pac_bnd.shape[0]])
for i in range(pac_bnd.shape[0]-1):
    pacbnd_len[i] = Haversine(pac_bnd[i],pac_bnd[i+1])
pacbnd_len[-1] = Haversine(pac_bnd[-1],pac_bnd[0])

