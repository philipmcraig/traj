# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 13:29:35 2017

@author: np838619
"""

import multiprocessing
import timeit

def runfile(filename):
    exec(open('/home/np838619/Trajectory/'+filename).read())

exec(open('/home/np838619/PminusE_data/ERA_Int/functions.py').read())

hometraj= '/home/np838619/Trajectory'
filenames = PrintFiles(hometraj,'readtraj_')

start_time = timeit.default_timer()

pool = multiprocessing.Pool(processes=9)
res = pool.map(runfile,filenames)
pool.close()
pool.join()

elapsed = timeit.default_timer() - start_time
print elapsed

#for i in filenames:
#        p = multiprocessing.Process(target=runfile(i))
#        p.start()