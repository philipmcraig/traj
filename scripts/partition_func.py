# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 14:02:07 2017

@author: np838619
"""

sheddir = '/home/np838619/Watershed/shed_defs/' 
endpts1 = pl.genfromtxt(sheddir+'Am_clicks.txt',skip_header=5)
endpts2 = pl.genfromtxt(sheddir+'AfMe_clicks.txt',skip_header=5)
endpts3 = pl.genfromtxt(sheddir+'EAA_clicks.txt',skip_header=5)

endpts4 = pl.genfromtxt(sheddir+'ArA_clicks.txt',skip_header=5)
endpts5 = pl.genfromtxt(sheddir+'ArI_clicks.txt',skip_header=5)
endpts6 = pl.genfromtxt(sheddir+'ArP_clicks.txt',skip_header=5)

endpts7 = pl.genfromtxt(sheddir+'SOA_clicks.txt',skip_header=5)
endpts8 = pl.genfromtxt(sheddir+'SOI_clicks.txt',skip_header=5)
endpts9 = pl.genfromtxt(sheddir+'SOP_clicks.txt',skip_header=5)

endpts10 = pl.genfromtxt(sheddir+'Ar_clicks.txt',skip_header=5)
endpts11 = pl.genfromtxt(sheddir+'SO_clicks.txt',skip_header=5)
#for i in range(endpts10.shape[0]):
#    if endpts10[i,0] < 0.:
#        endpts10[i,0] = endpts10[i,0] + 360.

# NEED 1 MAP FOR ATLANTIC/INDIAN CATCHMENTS, ANOTHER FOR PACIFIC
#m1 = Basemap(projection='robin',lon_0=0.) # Atlantic & Indian
m1 = Basemap(projection='cyl',llcrnrlon=-180,llcrnrlat=-80,urcrnrlon=180,urcrnrlat=80)
m2 = Basemap(projection='cyl',llcrnrlon=0,llcrnrlat=-80,urcrnrlon=360,urcrnrlat=90)
m3 = Basemap(projection='npaeqd',boundinglat=25,lon_0=270,round=True)
m4 = Basemap(projection='spstere',boundinglat=-13,lon_0=270,round=True)

# DEFINE POLYGONS FOR ATLANTIC, INDIAN AND PACIFIC CATCHMENT
atl_bnd = CatchLine(endpts1,endpts4,endpts2,endpts7)
ind_bnd = CatchLine(endpts2,endpts5,endpts3,endpts8)
#arW_bnd, arE_bnd = ArcticSectors(endpts4,endpts5,endpts6)
#arW_bnd[:,0], arW_bnd[:,1] = m3(arW_bnd[:,0],arW_bnd[:,1])
#arE_bnd[:,0], arE_bnd[:,1] = m3(arE_bnd[:,0],arE_bnd[:,1])
e6 = endpts6; e1 = endpts1; e1[:,0] = e1[:,0] + 360.
for i in range(e6.shape[0]):
    if e6[i,0] < 0.:
        e6[i,0] = e6[i,0] + 360.
pac_bnd = CatchLine(endpts3,e6,e1,endpts9)
arc_bnd = pl.zeros_like(endpts10); arc_bnd[:,0], arc_bnd[:,1] = m3(endpts10[:,0],endpts10[:,1])
sou_bnd = pl.zeros_like(endpts11); sou_bnd[:,0], sou_bnd[:,1] = m4(endpts11[:,0],endpts11[:,1])

atl_ply = CatchPoly(atl_bnd,m1); ind_ply = CatchPoly(ind_bnd,m1)
pac_ply = CatchPoly(pac_bnd,m2); arc_ply = CatchPoly(arc_bnd,m3)
#arW_ply = CatchPoly(arW_bnd,m1); arE_ply = CatchPoly(arE_bnd,m1)
sou_ply = CatchPoly(sou_bnd,m4)

# SET COUNTERS FOR EACH BASIN TO ZERO
atlcount = 0.; indcount = 0.; paccount = 0.; arccount = 0.; soucount = 0.
A = []; I = []; P = []; B = []
# CONVERT ORIGIN CO-ORDINATES TO MAP CO-ORDINATES
for traj in range(newlabs[0].shape[0]):
    # SHOULD ALSO CHECK IF ORIGIN TIMESTEP=0 AS THIS MEANS THE ORIGIN CO-ORDINATES
    # ARE THE RELEASE POINT AND NOT IN A PARTICULAR CATCHMENT
    if newlabs[0,traj,0] == 1. or newlabs[0,traj,0] == 2.:
        if newlabs[0,traj,1] == 0.:
            if lon[0,traj,1] > 180.:
                lonpt = lon[0,traj,1] - 360.
            else:
                lonpt = lon[0,traj,1]
            # DEAL WITH THIS BY CHECKING WHICH CATCHMENT TRAJ WAS IN AT TIMESTEP=1
            a,b = m1(lonpt,lat[0,traj,1],inverse=False)
            c,d = m2(lon[0,traj,1],lat[0,traj,1],inverse=False)
            e,f = m3(lon[0,traj,1],lat[0,traj,1],inverse=False)
            g,h = m4(lon[0,traj,1],lat[0,traj,1],inverse=False)
        else:
            if newlabs[0,traj,-1] > 180.:
                lonpt = newlabs[0,traj,-1] - 360.
            else:
                lonpt = newlabs[0,traj,-1]
            a,b = m1(lonpt,newlabs[0,traj,2],inverse=False)
            c,d = m2(newlabs[0,traj,-1],newlabs[0,traj,2],inverse=False)
            e,f = m3(newlabs[0,traj,-1],newlabs[0,traj,2],inverse=False)
            g,h = m4(newlabs[0,traj,-1],newlabs[0,traj,2],inverse=False)
        
        atlPath = mplPath.Path(list(atl_bnd)); indPath = mplPath.Path(list(ind_bnd))
        pacPath = mplPath.Path(list(pac_bnd)); arcPath = mplPath.Path(list(arc_bnd))
        arWPath = mplPath.Path(list(arW_bnd)); arEPath = mplPath.Path(list(arE_bnd))
        souPath = mplPath.Path(list(sou_bnd))
        X = atlPath.contains_point((a,b))
        Y = indPath.contains_point((a,b))
        Z = pacPath.contains_point((c,d))#; W3 = pacPath.contains_point((newlabs[0,traj,-1],newlabs[0,traj,2]))
        W = arcPath.contains_point((e,f))
        #W1 = arWPath.contains_point((e,f)); W2 = arEPath.contains_point((e,f))
        V = souPath.contains_point((g,h))
        if X == True:
            atlcount = atlcount + 1.; A.append(newlabs[0,traj])
        if Y == True:
            indcount = indcount + 1.; I.append(newlabs[0,traj])
        if Z == True:
            paccount = paccount + 1.; P.append(newlabs[0,traj])
        if W == True:# or W2 == True:# and W3 == False:
            #lt = NearestIndex(endpts10[:,0],newlabs[0,traj,-1])
            #if newlabs[0,traj,2] > endpts10[lt,1]:
            arccount = arccount + 1.; B.append(newlabs[0,traj])
            #else:
            #    pass
        if V == True:
            soucount = soucount + 1.#; print traj
print atlcount, indcount, paccount, arccount, soucount
A = pl.asarray(A); B = pl.asarray(B); I = pl.asarray(I); P = pl.asarray(P)