# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 11:28:29 2017

@author: np838619

Code to plot the integrated net and partitioned moisture fluxes of a catchment
boundary over a year, month etc. Runs from readtraj_*.py
"""

pl.plot(fluxtot,color='k',linewidth=2,label='$F_{tot}$')
pl.plot(atlsum,color='r',ls='--',linewidth=2,label='$F_{Atl}$')
pl.plot(pacsum,color='g',ls='--',linewidth=2,label='$F_{Pac}$')
pl.plot(indsum,color='b',ls='--',linewidth=2,label='$F_{Ind}$')
pl.plot(arcsum,color='deeppink',ls='--',linewidth=2,label='$F_{Arc}$')
pl.plot(sousum,color='darkgoldenrod',ls='--',linewidth=2,label='$F_{Sou}$')
pl.plot(strsum,color='indigo',ls='--',linewidth=2,label='$F_{Str}$')
pl.plot(unasum,color='saddlebrown',ls='--',linewidth=2,label='$F_{None}$')
pl.ylabel('Sv',fontsize=25); pl.xlabel('days',fontsize=25)
pl.xlim(0,len(trajfiles))
pl.legend(ncol=8,loc=(0.01,1.01),labelspacing=0.5,columnspacing=0.5,frameon=True,
          fontsize=20,handlelength=1.6)
pl.subplots_adjust(top=0.85,bottom=0.11)