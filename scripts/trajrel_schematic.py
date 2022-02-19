# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 13:45:13 2018

@author: np838619
"""

from __future__ import division
import pylab as pl
from scipy.stats import norm
import random

def TidyUp(axx):
    """
    """
    axx.spines['right'].set_color('none')
    axx.spines['top'].set_color('none')
    axx.spines['left'].set_color('none')
    axx.spines['bottom'].set_color('none')
    axx.xaxis.set_major_formatter(pl.NullFormatter())
    axx.yaxis.set_major_formatter(pl.NullFormatter())
    pl.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        left='off',
        top='off',         # ticks along the top edge are off
        right='off',
        labelbottom='off') # labels along the bottom edge are off

    return None

pl.close('all')
fig = pl.figure(figsize=(12,5))
ax = fig.add_subplot(1,1,1)
#ax = pl.figure(1)
TidyUp(ax)
pl.xlim(0,1); pl.ylim(0,1)

ax.axhline(y=0.05,xmax=0.5,color='b',lw=3)
ax.annotate('OCEAN',(0.2,-0.001),color='b',fontsize=20)

x = pl.linspace(-0.35,0.39,100)
ax.plot(x+0.85,-x**2+0.17,color='brown',lw=3)
ax.annotate('LAND',(0.8,0.05),color='brown',fontsize=20)

ax.axhline(y=0.15,xmax=0.5,color='r',ls='--',lw=2)
ax.axvline(x=0.5,ymin=0.15,ymax=0.22,color='r',ls='--',lw=2)
ax.plot(x+0.85,-x**2+0.35,color='r',lw=2,ls='--')
ax.annotate('ECMWF BL',(0.025,0.075),color='k',fontsize=18)

levs = pl.linspace(0.2,0.95,17)
xmax = 0.85131313131313124; ymax = 0.1699982756861545
ax.axvline(x=xmax,ymin=ymax,color='k',ls='--',lw=2)
for i in range(len(levs)):
    ax.axhline(y=levs[i],xmin=xmax-0.015,xmax=xmax+0.015,color='k',lw=2)

circ1 = pl.Circle((0.25,0.215),0.05,color='grey',ec='None')
circ2 = pl.Circle((0.22,0.245),0.03,color='grey',ec='None')
circ3 = pl.Circle((0.22,0.185),0.03,color='grey',ec='None')
circ4 = pl.Circle((0.20,0.215),0.03,color='grey',ec='None')
circ5 = pl.Circle((0.28,0.245),0.03,color='grey',ec='None')
circ6 = pl.Circle((0.25,0.265),0.03,color='grey',ec='None')
circ7 = pl.Circle((0.27,0.195),0.04,color='grey',ec='None')
circ8 = pl.Circle((0.25,0.185),0.03,color='grey',ec='None')

ax.add_artist(circ1); ax.add_artist(circ2); ax.add_artist(circ3)
ax.add_artist(circ4); ax.add_artist(circ5); ax.add_artist(circ6)
ax.add_artist(circ7); ax.add_artist(circ8)

xp = pl.linspace(xmax,0,500)
yp = pl.zeros_like(xp); yp[0] = levs[1]
yp[1:240] = pl.sin(50*xp[1:240])/100 + 0.25
yp[240:420] = -pl.cos(50*xp[240:420])/16 + 0.185
yp[380:] = pl.sin(50*xp[380:])/100 + 0.235
ax.annotate(' ',xy=(xp[45],yp[45]),xytext=(xp[55],yp[55]),size=40,
                arrowprops={'arrowstyle':'->','color':'b','lw':'2'})

ax.plot(xp,yp,color='b',lw=2,label='CAT I')
ax.axvline(xp[267],ymin=0.05,ymax=0.15,lw=2,ls='--',color='b')

yp = pl.zeros_like(xp); yp[0] = levs[5]
yp[1:280] = -pl.sin(50*xp[1:280])/100 + 0.425
yp[280:440] = -pl.cos(25*xp[280:440])/10 + 0.325
yp[440:] = pl.sin(50*xp[440:])/100 + 0.415
ax.annotate(' ',xy=(xp[45],yp[45]),xytext=(xp[55],yp[55]),size=40,
                arrowprops={'arrowstyle':'->','color':'g','lw':'2'})

ax.plot(xp,yp,color='g',lw=2,label='CAT II')
ax.axvline(xp[335],ymin=0.05,ymax=yp[335],color='g',lw=2,ls='--')

pl.axhline(y=0.96,ls='--',lw=2,color='purple')
yp = pl.zeros_like(xp); yp[0] = levs[15]
yp[1:150] = -pl.sin(50*xp[1:150])/100 + 0.89
yp[150:] = -pl.cos(10*xp[150:])/10 + 0.995
ax.annotate(' ',xy=(xp[45],yp[45]),xytext=(xp[55],yp[55]),size=40,
                arrowprops={'arrowstyle':'->','color':'purple','lw':'2'})

ax.plot(xp,yp,color='purple',lw=2,label='stratosphere')
ax.annotate('stratosphere',xy=(0.6,0.98),fontsize=18,color='purple')

yp = pl.zeros_like(xp); yp[0] = levs[9]
yp[1:] = -pl.sin(50*xp[1:])/100 + 0.61
ax.annotate(' ',xy=(xp[45],yp[45]),xytext=(xp[55],yp[55]),size=40,
                arrowprops={'arrowstyle':'->','color':'k','lw':'2'})

ax.plot(xp,yp,color='k',lw=2,label='no origin')

ax.legend(loc=(0.01,0.7),ncol=2,columnspacing=0.5)

ax.arrow(0.95,0.5,0.0,0.1,color='k')
ax.annotate('z',xy=(0.96,0.54),fontsize=18)

pl.tight_layout()
#pl.savefig('/home/np838619/Trajectory/trajrel_schematic.png')