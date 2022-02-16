# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 10:36:30 2018

@author: np838619
"""

import cartopy
import cartopy.crs as ccrs
from matplotlib.lines import Line2D

pl.close('all')
fig,ax = pl.subplots(2,2,figsize=(20,12))
proj = ccrs.PlateCarree(central_longitude=-70)
#ax = pl.subplot(projection=ccrs.PlateCarree(central_longitude=110))

#amr_ext = [-180,180,80,-80]
afr_ext = [210,-160,80,-80]
#eaa_ext = [-210,150,80,-80]
titles = ['(a) Atlantic origin, winter', '(b) Indian origin, winter',
          '(c) Atlantic origin, summer', '(d) Indian origin, summer']

a = 0; b = 4
#catch1 = atlcount
#catch2 = indcount
catchcount = [atlcount[0],indcount[0],atlcount[1],indcount[1]]
#L1 = []; L2 = []; L3 = []; L4 = []; L5 = []
LINES = [Line2D([0], [0], color='grey'),Line2D([0], [0], color='darkgoldenrod'),
         Line2D([0], [0], color='g'),Line2D([0], [0], color='r'),
            Line2D([0], [0], color='b')]
LABS = ['$\\tilde{p}_{t=0}<600$ hPa','$600<\\tilde{p}_{t=0}<700$ hPa',
            '$700<\\tilde{p}_{t=0}<800$ hPa','$800<\\tilde{p}_{t=0}<900$ hPa',
                                                '$\\tilde{p}_{t=0}>900$ hPa']

for p in range(4):
    axx = pl.subplot(2,2,p+1,projection=proj)
    axx.coastlines()
    axx.set_extent([-180,180,80,-80],ccrs.PlateCarree())
    
    if p < 2:
        j = 0
    else:
        j = 1
    
    for i in range(catchcount[p].shape[0]):
        if len(rlspts)*a <= catchcount[p][i,0] < len(rlspts)*b:
            if pres[j,catchcount[p][i,0],0] < 600.:
                c = 'grey'; zo = 2
            elif 600 < pres[j,catchcount[p][i,0],0] < 700:
                c = 'darkgoldenrod'; zo = 2
            elif 700 < pres[j,catchcount[p][i,0],0] < 800:
                c = 'g'; zo = 1
            elif 800 < pres[j,catchcount[p][i,0],0] < 900:
                c = 'r'; zo = 1
            elif pres[j,catchcount[p][i,0],0] > 900:
                c= 'b'; zo = 1
        
            if catchcount[p][i,1] == 1:
                OR1 = axx.plot(lon[j,catchcount[p][i,0],catchcount[p][i,2]],
                        lat[j,catchcount[p][i,0],catchcount[p][i,2]],
                        color='k',marker='x',mew=0.8,transform=ccrs.Geodetic(),
                            zorder=3,ls='None')
            elif catchcount[p][i,1] == 2:
                OR2 = axx.plot(lon[j,catchcount[p][i,0],catchcount[p][i,2]],
                        lat[j,catchcount[p][i,0],catchcount[p][i,2]],
                        color='k',marker='o',mew=0.8,zorder=3,ms=5,
                        transform=ccrs.Geodetic(),ls='None')
            
            axx.plot(lon[j,catchcount[p][i,0],:],lat[j,catchcount[p][i,0],:],
                    transform=ccrs.Geodetic(),color=c,zorder=zo)
           
    gl = axx.gridlines(xlocs=[180,-120,-60,0,60,120,180],
                       ylocs=[-90,-60,-30,0,30,60,90],draw_labels=True)
    gl.xlabels_top = False
    
    g2 = axx.gridlines(xlocs=[-180,-120],ylocs=[-90,-60,-30,0,30,60,90],
                       draw_labels=False)
    
    if p in (1,3):
        gl.ylabels_left = False
    elif p in (0,2):
        gl.ylabels_right = False

    axx.plot(rlspts[:,0],rlspts[:,1],color='skyblue',transform=ccrs.PlateCarree(),
             zorder=3,lw=1.5)
    pl.title(titles[p],fontsize=18)

leg1 = fig.legend(LINES,LABS,loc=4,ncol=3,fontsize=20,columnspacing=0.5)
for legobj in leg1.legendHandles:
    legobj.set_linewidth(3.0)

leg2 = fig.legend([OR1[0],OR2[0]],['CAT I','CAT II'],loc=(0.01,0.025),fontsize=20)

pl.tight_layout()
pl.subplots_adjust(left=0.03,right=0.97,wspace=0.07,hspace=-0.3)