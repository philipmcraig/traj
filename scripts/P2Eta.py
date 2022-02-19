# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 16:37:36 2017

@author: np838619
"""

def P2Eta(psurf):
    """
    """
    p0 = 10**5
    psurf = psurf
    coeffs = pl.genfromtxt(trajdir + 'coeffs.txt')
    a = coeffs[:,1]; b = coeffs[:,2]
    
    pL = a + b*(psurf/p0)
    
    return pL