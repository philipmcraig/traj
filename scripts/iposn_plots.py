# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 16:06:10 2016

@author: np838619
"""

import pylab

q8 = pl.genfromtxt('/home/np838619/Trajectory/qout_I8.txt')
p8 = pl.genfromtxt('/home/np838619/Trajectory/pout_I8.txt')
q3 = pl.genfromtxt('/home/np838619/Trajectory/qout_I3.txt')
p3 = pl.genfromtxt('/home/np838619/Trajectory/pout_I3.txt')

pl.plot(q8[:,0],p8[:,0],label='I8(t=0)',color='b',ls='-')
pl.plot(q8[:,1],p8[:,1],label='I8(t=1)',color='b',ls='--')
pl.plot(q8[:,2],p8[:,2],label='I8(t=2)',color='b',ls=':')
pl.plot(q3[:,0],p3[:,0],label='I3(t=0)',color='r',ls='-')
pl.plot(q3[:,1],p3[:,1],label='I3(t=1)',color='r',ls='--')
pl.plot(q3[:,2],p3[:,2],label='I3(t=2)',color='r',ls=':')
pl.legend(loc=0)
pl.ylim(1100,0)