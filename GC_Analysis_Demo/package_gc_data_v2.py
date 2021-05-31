#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 10:41:23 2020

@author: jgreenhalgh2
"""

# folder= a string that matches the name of the folder to be scanned in
#           example, if opening a folder called GC_data, folder='GC_data'
#target_rts= a list of the target retention times (e.g. [1.0, 3.3, 4.3, 6.3] etc)

import sys
import gc_tools
import pickle

#Get target retention times if applicable
if len(sys.argv)==3:
    target_rts=sys.argv[2]
else:
    target_rts=[]
    
gc_data=gc_tools.read_files(sys.argv[1],target_rts)

pickle.dump(gc_data,open(sys.argv[1]+'_data.p','wb'))


#def package_gc_data(folder,target_rts=[]):
#    import gc_tools
#    import pickle
#    
#    #Call gc_tools and read the files
#    
#    gc_data=gc_tools.read_files(folder,target_rts)
#    
#    #Package data into a .p file for quick access
#     
#    pickle.dump(gc_data,open(folder+'_data.p','wb'))