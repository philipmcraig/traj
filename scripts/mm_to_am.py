# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 18:01:33 2017

@author: np838619
"""

import pylab as pl

years = ['2010','2011','2012','2013','2014']
shed = 'sop'
homedir = '/home/np838619/Trajectory/fluxpart/'

data1 = pl.genfromtxt(homedir+years[0]+'/'+shed+'/catchmeans.csv',delimiter=' ')#/(10**9)
#                     skip_header=1)
data2 = pl.genfromtxt(homedir+years[1]+'/'+shed+'/catchmeans.csv',delimiter=' ')#/(10**9)
#                    skip_header=1)
data3 = pl.genfromtxt(homedir+years[2]+'/'+shed+'/catchmeans.csv',delimiter=' ')#/(10**9)
#                   skip_header=1)
data4 = pl.genfromtxt(homedir+years[3]+'/'+shed+'/catchmeans.csv',delimiter=' ')#/(10**9)
#                  skip_header=1)
data5 = pl.genfromtxt(homedir+years[4]+'/'+shed+'/catchmeans.csv',delimiter=' ')#/(10**9)
#                skip_header=1)

tot = pl.vstack((data1[:,0],data2[:,0],data3[:,0],data4[:,0],data5[:,0]))
tot_mn = pl.mean(tot,axis=0)

atl = pl.vstack((data1[:,1],data2[:,1],data3[:,1],data4[:,1],data5[:,1]))
atl_mn = pl.mean(atl,axis=0)

ind = pl.vstack((data1[:,2],data2[:,2],data3[:,2],data4[:,2],data5[:,2]))
ind_mn = pl.mean(ind,axis=0)

pac = pl.vstack((data1[:,3],data2[:,3],data3[:,3],data4[:,3],data5[:,3]))
pac_mn = pl.mean(pac,axis=0)

arc = pl.vstack((data1[:,4],data2[:,4],data3[:,4],data4[:,4],data5[:,4]))
arc_mn = pl.mean(arc,axis=0)

sou = pl.vstack((data1[:,5],data2[:,5],data3[:,5],data4[:,5],data5[:,5]))
sou_mn = pl.mean(sou,axis=0)

sts = pl.vstack((data1[:,6],data2[:,6],data3[:,6],data4[:,6],data5[:,6]))
sts_mn = pl.mean(sts,axis=0)

una = pl.vstack((data1[:,7],data2[:,7],data3[:,7],data4[:,7],data5[:,7]))
una_mn = pl.mean(una,axis=0)

#mean10 = pl.mean(data1,axis=0); mean11 = pl.mean(data2,axis=0)
#mean12 = pl.mean(data3,axis=0); mean13 = pl.mean(data4,axis=0)
#mean14 = pl.mean(data5,axis=0)

sscyc = pl.vstack((tot_mn,atl_mn,ind_mn,pac_mn,arc_mn,sou_mn,sts_mn,una_mn))
#profiles = pl.vstack((tot_mn,atl_mn,ind_mn,pac_mn,arc_mn,sou_mn))#,sts_mn,una_mn))

f = open('/home/np838619/Trajectory/fluxpart/season_cycs_'+shed+'.csv','w')
f.write('Total Atlantic Indian Pacific Arctic Southern\n')
pl.savetxt(f,sscyc.T)
f.close()

means1 = pl.mean(data1,axis=0)
means2 = pl.mean(data2,axis=0)
means3 = pl.mean(data3,axis=0)
means4 = pl.mean(data4,axis=0)
means5 = pl.mean(data5,axis=0)