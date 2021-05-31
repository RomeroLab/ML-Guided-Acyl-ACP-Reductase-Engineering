#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 15:59:37 2018

@author: jonathangreenhalgh
"""
# folder= a string that matches the name of the folder to be scanned in
#           example, if opening a folder called GC_data, folder='GC_data'
#target_rts= a list of the target retention times (e.g. [1.0, 3.3, 4.3, 6.3] etc)

def package_gc_data(folder,target_rts=[]):
    import gc_tools
    import pickle
    
    #Call gc_tools and read the files
    
    gc_data=gc_tools.read_files(folder,target_rts)
    
    #Package data into a .p file for quick access
     
    pickle.dump(gc_data,open(folder+'_data.p','wb'))