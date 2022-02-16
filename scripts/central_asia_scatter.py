# -*- coding: utf-8 -*-
"""
Created on Sat Jul 29 16:37:01 2017

@author: np838619
"""

import cartopy
import cartopy.crs as ccrs

m = Basemap(projection='robin',lon_0=90.)
#m.drawcoastlines()
#
x = pl.reshape(newlabs,(newlabs.shape[0],NCLUST,rlslabs.shape[0],newlabs.shape[-1]))
china = x[:,:,:73,:]
#seasia = newlabs[:,:61,:]
#
#china = pl.asarray(china)
china = pl.reshape(china,(china.shape[0]*china.shape[1]*china.shape[2],china.shape[3]))
#seasia = pl.asarray(seasia)
#seasia = pl.reshape(seasia,(seasia.shape[0]*seasia.shape[1],seasia.shape[2]))
##
#for i in range(china.shape[0]):
#    m.scatter(china[i,:,:,-1],china[i,:,:,-2],latlon=True)
#    m.plot(rlslabs[:73,0],rlslabs[:73,1],latlon=True,color='r',lw=2.)
#    m.drawcoastlines()
#    pl.savefig('/home/np838619/Trajectory/china_scatter/map'+str(i)+'.png')
#    pl.close()
#
lat_x = pl.reshape(lat,(32,17,155,57)); lon_x = pl.reshape(lon,(32,17,155,57))
#lon_x = lon_x[:,:,:73,:]; lat_x = lat_x[:,:,:73,:]
#
#for i in range(lon_x.shape[0]):
for i in range(lat_x.shape[3]):
    ax = pl.axes(projection=ccrs.Robinson(central_longitude=-90))
    ax.set_global()
    ax.coastlines()
    pl.plot(rlslabs[:,0],rlslabs[:,1],color='b',lw=2.,transform=ccrs.Geodetic())
    for j in range(lon_x.shape[2]):
        pl.plot(lon_x[0,0,j,:i+1],lat_x[0,0,j,:i+1],transform=ccrs.Geodetic())
    pl.tight_layout()
    #pl.savefig('/home/np838619/Trajectory/gif/0'+str(i)+'map.png')
    pl.close()