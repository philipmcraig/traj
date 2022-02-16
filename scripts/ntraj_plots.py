# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 18:15:34 2018

@author: np838619
"""

import pylab as pl
from scipy import stats

def my_autopct(pct):
    return ('%1.1f%%' % pct) if pct > 3 else ''

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

pl.close('all')
clusdir = '/glusterfs/scenario/users/np838619/traj/'
homedir = '/home/np838619/Trajectory/fluxpart/'
sheddir = '/home/np838619/Watershed/shed_defs/'

sheds = ['amr','afr','eaa','ara','ari','arp','soa','soi','sop']
years = ['2010','2011','2012','2013','2014']
filenames = pl.zeros([len(years),len(sheds)],dtype='S63')

#for sd in range(len(sheds)):
#    S = sheds[sd]
for yr in range(len(years)):
    path = clusdir+years[yr]+'/'
    filenames[yr] = PrintFiles(path,'props_')
    filenames[yr] = pl.sort(filenames[yr])

NTRAJ_yrs = pl.zeros([filenames.shape[0],filenames.shape[1],7])

for yr in range(len(years)):
    for sd in range(len(sheds)):
        fname = clusdir+years[yr]+'/'+filenames[yr,sd]
        data = pl.genfromtxt(fname,skip_header=1,delimiter=',')
        data = data[:7]
        NTRAJ_yrs[yr,sd,:] = data

NTRAJ_mn = pl.mean(NTRAJ_yrs,axis=0)
res = 100 - pl.sum(NTRAJ_mn,axis=1)
NTRAJ_mn = NTRAJ_mn + res[:,None]/7

S = pl.zeros([len(sheds)])
for i in range(len(sheds)):
    for j in range(len(sheds)):
        if sheds[i] in filenames[0][j]:
                q = pl.where(filenames[0]==filenames[0][j])
                S[i] = q[0]

fig, ax = pl.subplots(3,3,figsize=(14,14))
c = ['r','b','g','deeppink','darkgoldenrod','skyblue','white']
T = ['(a) Americas','(b) Africa','(c) South-East Asia',
     '(d) Arctic Atlantic','(e) Arctic Indian','(f) Arctic Pacific',
    '(g) Southern Atlantic','(h) Southern Indian','(i) Southern Pacific']
n = ['9,625,655','13,289,614','10,246,665',
     '12,233,897','6,396,403','13,662,220',
     '11,488,685','11,426,584','13,165,412']

for i in range(ax.size):
    axx = pl.subplot(3,3,i+1)
    axx.pie(NTRAJ_mn[S[i]],autopct=my_autopct,colors=c,
            textprops={'fontsize':'16'},labeldistance=1.7,
            wedgeprops={'linewidth':'0.5'})
    axx.axis('equal')
    pl.title(T[i]+',\n $n=$'+n[i],fontsize=18)
    
    if i == 0:
        pl.text(0.65,0.8,'Atlantic',fontsize=15)
        pl.text(0.73,-0.8,'no origin',fontsize=15)
    elif i == 2:
        pl.text(-1.01,0.8,'Pacific',fontsize=15)
        pl.text(0.65,0.8,'Indian',fontsize=15)
    elif i == 5:
        pl.text(-1.2,0.6,'Arctic',fontsize=15)
        pl.text(-1.55,-0.8,'Stratosphere',fontsize=15)
    elif i == 8:
        pl.text(-1.4,0.6,'Southern',fontsize=15)

pl.tight_layout()
#pl.savefig(homedir+'annmean_plots/panels/ntraj_pies.png')