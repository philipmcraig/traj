# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 13:13:00 2017

@author: np838619
"""

import pylab as pl

homedir = '/home/np838619/Trajectory/fluxpart/'
year ='2014'#; shed = 'soi'
sheds = ['amr','afr','eaa','ara','ari','arp','soa','soi','sop']
#sheds = ['amr']

for i in range(len(sheds)):
    shed = sheds[i]
    catchints = pl.genfromtxt(homedir+year+'/'+shed+'/catchints.csv',skip_header=1,
                                                                      delimiter=' ')
    
    jan = pl.zeros([62,8]); mar = pl.zeros([62,8]); apr = pl.zeros([60,8])
    may = pl.zeros([62,8]); jun = pl.zeros([60,8]); jul = pl.zeros([62,8])
    aug = pl.zeros([62,8]); sep = pl.zeros([60,8]); ocr = pl.zeros([62,8])
    nov = pl.zeros([60,8])
    
    if year == '2012':
        feb = pl.zeros([58,8])
    else:
        feb = pl.zeros([56,8])
    
    if year == '2014':
        dec = pl.zeros([63,8])
    else:
        dec = pl.zeros([62,8])
    
    jan[:] = catchints[:len(jan)]; jan_mn = pl.mean(jan,axis=0)
    start = len(jan); end = start + len(feb)
    feb[:] = catchints[start:end]; feb_mn = pl.mean(feb,axis=0)
    start = end; end = start + len(mar)
    mar[:] = catchints[start:end]; mar_mn = pl.mean(mar,axis=0)
    start = end; end = start + len(apr)
    apr[:] = catchints[start:end]; apr_mn = pl.mean(apr,axis=0)
    start = end; end = start + len(may)
    may[:] = catchints[start:end]; may_mn = pl.mean(may,axis=0)
    start = end; end = start + len(jun)
    jun[:] = catchints[start:end]; jun_mn = pl.mean(jun,axis=0)
    start = end; end = start + len(jul)
    jul[:] = catchints[start:end]; jul_mn = pl.mean(jul,axis=0)
    start = end; end = start + len(aug)
    aug[:] = catchints[start:end]; aug_mn = pl.mean(aug,axis=0)
    start = end; end = start + len(sep)
    sep[:] = catchints[start:end]; sep_mn = pl.mean(sep,axis=0)
    start = end; end = start + len(ocr)
    ocr[:] = catchints[start:end]; ocr_mn = pl.mean(ocr,axis=0)
    start = end; end = start + len(nov)
    nov[:] = catchints[start:end]; nov_mn = pl.mean(nov,axis=0)
    start = end#; end = start + len(apr)
    dec[:] = catchints[start:]; dec_mn = pl.mean(dec,axis=0)
    
    x = (jan_mn,feb_mn,mar_mn,apr_mn,may_mn,jun_mn,jul_mn,aug_mn,sep_mn,ocr_mn,
                                                         nov_mn,dec_mn)
    means = pl.vstack(x)
    
    f = open(homedir+year+'/'+shed+'/'+'catchmeans.csv','w')
    pl.savetxt(f,means)
    f.close()