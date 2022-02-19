# -*- coding: utf-8 -*-
"""
Created on Tue May 23 18:53:58 2017

@author: np838619
"""

from __future__ import division
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.basemap import Basemap
from matplotlib.collections import PolyCollection

pl.close('all')

m = Basemap(projection='cyl',llcrnrlat=-45,llcrnrlon=-155,urcrnrlat=60,urcrnrlon=-25)

fig = pl.figure()
ax = fig.gca(projection='3d')#ax = Axes3D(fig)
ax.add_collection3d(m.drawcoastlines(linewidth=0.5))
ax.add_collection3d(m.drawcountries(linewidth=0.5))

ax.set_axis_off()
ax.azim = 330
ax.elev=20
ax.dist = 7.5

amline = pl.genfromtxt('/home/np838619/Watershed/shed_defs/Am_traj_release_new.txt',skip_header=5)
line = pl.zeros([57,3])
line[:,0] = pl.array([ 266.79,  266.31,  265.81,  265.39,  265.6 ,  266.45,  267.65,
        269.22,  270.85,  272.57,  274.32,  276.05,  277.85,  279.57,
        281.1 ,  282.01,  281.77,  281.76,  282.08,  282.64,  283.46,
        284.64,  286.24,  288.09,  290.06,  291.68,  293.06,  294.33,
        295.38,  296.21,  297.04,  297.67,  297.9 ,  297.98,  298.06,
        297.72,  296.81,  295.27,  293.18,  290.15,  285.85,  280.78,
        275.7 ,  271.12,  266.19,  261.57,  255.73,  250.02,  244.61,
        238.09,  230.59,  220.93,  211.02,  203.61,  202.87,  206.93,
        210.69])
line[:,1] = pl.array([ 15.987 ,  14.696 ,  13.858 ,  13.006 ,  11.927 ,  10.982 ,
        10.059 ,   9.14  ,   8.5239,   8.1371,   7.8332,   7.4721,
         7.2029,   7.1174,   6.8787,   6.2905,   6.6885,   8.0356,
         8.934 ,   9.6416,  10.416 ,  11.111 ,  11.389 ,  11.716 ,
        12.684 ,  13.611 ,  14.411 ,  15.306 ,  16.004 ,  16.527 ,
        17.26  ,  18.257 ,  19.202 ,  19.961 ,  20.876 ,  22.022 ,
        23.266 ,  24.623 ,  25.928 ,  28.011 ,  31.107 ,  35.214 ,
        39.484 ,  42.626 ,  44.566 ,  45.776 ,  46.67  ,  46.836 ,
        45.803 ,  42.995 ,  39.616 ,  37.998 ,  38.43  ,  42.611 ,
        50.808 ,  57.869 ,  61.958 ])
line[:,2] = pl.array([ 526.36,  532.83,  523.07,  506.7 ,  498.1 ,  495.51,  491.67,
        472.92,  451.28,  431.07,  406.92,  397.57,  395.21,  386.9 ,
        442.8 ,  651.88,  802.11,  875.61,  900.4 ,  914.55,  920.03,
        921.13,  895.51,  932.5 ,  958.13,  964.34,  973.92,  977.51,
        974.85,  972.96,  971.9 ,  966.85,  962.4 ,  958.32,  952.07,
        931.81,  919.31,  903.5 ,  868.03,  829.86,  738.08,  643.64,
        570.69,  537.53,  526.37,  512.38,  492.54,  484.04,  478.41,
        515.15,  568.77,  545.83,  521.96,  455.64,  385.65,  361.9 ,
        341.18])
line[:,2] = line[:,2]/2000; line[0,2] = 0.03*9
#for i in range(1,100):
#    line[i,2] = random.uniform(0,0.34)

#m.plot(amline[:,0],amline[:,1],linewidth=3)
zs = 0.0
for i in range(17):
    zs = zs + 0.03
    ax.plot3D(amline[:,0],amline[:,1],zs=zs,zdir='z',color='b')

ax.scatter(line[0,0]-360,line[0,1],line[0,2],lw=1,color='r',marker='x')


ax.plot3D(line[:,0]-360,line[:,1],line[:,2],zdir='z',color='r',linewidth=2)
ax.plot3D([-70,-70],[-25,-25],[0.01,0.5],zdir='z',color='k',linewidth=2.5)
ax.plot3D([-70,-70],[-25,-24],[0.5,0.47],zdir='z',color='k',linewidth=2)
ax.plot3D([-70,-70],[-25,-26],[0.5,0.47],zdir='z',color='k',linewidth=2)
ax.plot3D([-70,-70],[-26,-24],[0.47,0.47],zdir='z',color='k',linewidth=2)
ax.text3D(-80,-29,0.25,'$z$',fontsize=30)

ax.plot3D([-39,-39],[0,45],[0,0],color='k',linewidth=2.5)
ax.plot3D([-39,-36],[45,42],[0,0],color='k',linewidth=2.5)
ax.plot3D([-39,-42],[45,42],[0,0],color='k',linewidth=2.5)
ax.plot3D([-36,-42],[42,42],[0,0],color='k',linewidth=2.5)
ax.text3D(-43,18,0,'N',fontsize=30)