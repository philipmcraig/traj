# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 10:36:30 2018

@author: np838619
"""

import cartopy
import cartopy.crs as ccrs
from matplotlib.lines import Line2D

pl.close('all')
fig = pl.figure(figsize=(16, 10))
ax = pl.axes(projection=ccrs.PlateCarree(central_longitude=30))
#ax = pl.subplot(projection=ccrs.PlateCarree(central_longitude=110))
ax.coastlines()
#amr_ext = [-180,180,80,-80]
afr_ext = [210,-160,80,-80]
#eaa_ext = [-210,150,80,-80]
ax.set_extent(afr_ext,ccrs.PlateCarree())
pl.tight_layout()

a = 0; b = 4
catchcount = indcount
#L1 = []; L2 = []; L3 = []; L4 = []; L5 = []
LINES = [Line2D([0], [0], color='grey'),Line2D([0], [0], color='darkgoldenrod'),
         Line2D([0], [0], color='g'),Line2D([0], [0], color='r'),
            Line2D([0], [0], color='b')]
LABS = ['$\\tilde{p}_{t=0}<600$ hPa','$600<\\tilde{p}_{t=0}<700$ hPa',
            '$700<\\tilde{p}_{t=0}<800$ hPa','$800<\\tilde{p}_{t=0}<900$ hPa',
                                                '$\\tilde{p}_{t=0}>900$ hPa']

for i in range(catchcount[0].shape[0]):
    if 155*a < catchcount[0][i,0] < 155*b:
        if pres[0,catchcount[0][i,0],0] < 600.:
            c = 'grey'
        elif 600 < pres[0,catchcount[0][i,0],0] < 700:
            c = 'darkgoldenrod'
        elif 700 < pres[0,catchcount[0][i,0],0] < 800:
            c = 'g'
        elif 800 < pres[0,catchcount[0][i,0],0] < 900:
            c = 'r' 
        elif pres[0,catchcount[0][i,0],0] > 900:
            c= 'b' 
        if catchcount[0][i,1] == 1:
            ax.plot(lon[0,catchcount[0][i,0],catchcount[0][i,2]],
                    lat[0,catchcount[0][i,0],catchcount[0][i,2]],
                    color='k',marker='x',mew=0.8,transform=ccrs.Geodetic(),zorder=2)
        elif catchcount[0][i,1] == 2:
            ax.plot(lon[0,catchcount[0][i,0],catchcount[0][i,2]],
                    lat[0,catchcount[0][i,0],catchcount[0][i,2]],
                    color='k',marker='o',mew=0.8,zorder=2,ms=5,
                    transform=ccrs.Geodetic(),)
        ax.plot(lon[0,catchcount[0][i,0],:],lat[0,catchcount[0][i,0],:],
                color=c,transform=ccrs.Geodetic(),zorder=1)
                
#        if pres[0,catchcount[0][i,0],0] < 600.:
#            L1.append(T[0])
#        elif 600 < pres[0,catchcount[0][i,0],0] < 700:
#            L2.append(T[0])
#        elif 700 < pres[0,catchcount[0][i,0],0] < 800:
#            L3.append(T[0])
#        elif 800 < pres[0,catchcount[0][i,0],0] < 900:
#            L4.append(T[0])
#        elif pres[0,catchcount[0][i,0],0] > 900:
#            L5.append(T[0])
           
gl = ax.gridlines(ylocs=[0],draw_labels=True)

ax.plot(rlspts[:,0],rlspts[:,1],color='skyblue',transform=ccrs.PlateCarree(),
        zorder=3)
leg = ax.legend(LINES,LABS,loc=4)
for legobj in leg.legendHandles:
    legobj.set_linewidth(2.0)