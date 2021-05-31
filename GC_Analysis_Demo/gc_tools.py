 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 12:57:30 2018

@author: jonathangreenhalgh
"""

# New GC tools functions
import os
import numpy
import matplotlib.pyplot as plt

#Define a function to pull all relevant data from all files within a folder
def read_files(folder,target_rts=[]):
    filenames=[]

    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.endswith('.TXT'):
                filenames.append(file)
    
    #Use the filenames to open each file and read out the data (leave out the standards the hexane washes)
    data=[]
    if target_rts==[]:
        for name in filenames:
            data.append(get_data(folder,name))
            print('Getting '+name)
    else:
        for name in filenames:
            data.append(get_data_fast(folder,name,target_rts))
            print('Getting '+name)
            
    
    return [filenames,data]


##############################################################################################################################
#Define a function to pull all relevant data from a file   
def get_data(folder,filename):
    
    #Test reading text file for GC outputs
    
    gc_file=folder+'/'+filename
    
    file=open(gc_file,'r')
    
    #Read all data in the textfile into a list
    file_list=[]
    for line in file:
        newline=line.split('\t')                    #Split the current line by tabs, creating a new line (newline)
        newline[-1]=newline[-1].replace('\n','')    #remove '\n' from last entry of the list
        file_list.append(newline)                   #append the file_list with the edited new line
     
    #Find line where retention times begin by locating the row where 'Intensity' occurs in the final element of the row list
    #Also, find the number of time points
    for i in range(len(file_list)):
        if file_list[i][0]=='# of Points':
            chromatograph_number_of_pts=int(file_list[i][-1])
        if file_list[i][-1]=='Intensity':
            chromatograph_start_row=i
        #Loop through file and determine start and end point of the peak table
        if file_list[i][0]=='Peak#':
            peak_table_start_row=i+1
        if file_list[i][0]=='[Compound Results (Ch1)]':
            peak_table_end_row=i-1
        
    chromatograph=[]
    for i in range(chromatograph_start_row+1, chromatograph_start_row+1+chromatograph_number_of_pts):    #The +1 is to move the start row to the row where data is rather than the row of the header
        chromatograph.append([float(file_list[i][0]),float(file_list[i][1])])
     
    chromatograph=numpy.array(chromatograph)

    #Beggining at start_row and ending at end_row in file_list, parse through the file and save the rows containing peak information
     
    rt=[]       #empty list for retention times
    area=[]     #empty list for areas
    height=[]   #empty list for height of peaks
    
    for i in range(peak_table_start_row,peak_table_end_row):
        rt.append(float(file_list[i][1]))
        area.append(float(file_list[i][4]))
        height.append(float(file_list[i][5]))
        
    return [chromatograph,rt,area, height]  

##############################################################################################################################
#Define a function to find retention times of standards (assuming the max peak is the standard)
    
def get_standard_rts(folder,max_standard_list):
    gc_data=read_files(folder)
    filenames=gc_data[0]
    data=gc_data[1]
    
    standard_rts=[]
    for i in range(len(filenames)):
        if filenames[i] in max_standard_list:
            for j in range(len(data[i][2])):
                if data[i][2][j]==numpy.max(data[i][2]):
                    standard_rts.append(data[i][1][j])
                    break
    
    return [standard_rts]

##############################################################################################################################              
#potentially faster version of get_data that generates a smaller peak table
def get_data_fast(folder,filename,target_rts):
    
    #Test reading text file for GC outputs
    
    gc_file=folder+'/'+filename
    
    file=open(gc_file,'r')
    
    #Read all data in the textfile into a list
    file_list=[]
    for line in file:
        newline=line.split('\t')                    #Split the current line by tabs, creating a new line (newline)
        newline[-1]=newline[-1].replace('\n','')    #remove '\n' from last entry of the list
        file_list.append(newline)                   #append the file_list with the edited new line
     
    #Find line where retention times begin by locating the row where 'Intensity' occurs in the final element of the row list
    #Also, find the number of time points
    for i in range(len(file_list)):
        if file_list[i][0]=='# of Points':
            chromatograph_number_of_pts=int(file_list[i][-1])
        if file_list[i][-1]=='Intensity':
            chromatograph_start_row=i
        #Loop through file and determine start and end point of the peak table
        if file_list[i][0]=='Peak#':
            peak_table_start_row=i+1
        if file_list[i][0]=='[Compound Results (Ch1)]':
            peak_table_end_row=i-1
        
    chromatograph=[]
    for i in range(chromatograph_start_row+1, chromatograph_start_row+1+chromatograph_number_of_pts):    #The +1 is to move the start row to the row where data is rather than the row of the header
        chromatograph.append([float(file_list[i][0]),float(file_list[i][1])])
     
    chromatograph=numpy.array(chromatograph)

    #Beggining at start_row and ending at end_row in file_list, parse through the file and save the rows containing peak information
     
    rt=[]       #empty list for retention times
    area=[]     #empty list for areas
    height=[]   #empty list for height of peaks
    
    window=0.1
    
    for i in range(peak_table_start_row,peak_table_end_row):
        for time in target_rts:
            if float(file_list[i][1])<time+window and float(file_list[i][1])>time-window:
                rt.append(float(file_list[i][1]))
                area.append(float(file_list[i][4]))
                height.append(float(file_list[i][5]))
        
    return [chromatograph,rt,area, height]  

##############################################################################################################################

def get_peaks(data,target_rts,windows=[]):
    #data= gc data decoupled from the filename (but in the same order)
    #Loop through each data file and find the peaks associated with each retention time
    
    if windows==[]:
        for i in range(len(target_rts)):
            window=[]
            window.append(target_rts[i]-0.1)
            window.append(target_rts[i]+0.1)
            windows.append(window)
    
    peaks=numpy.zeros([len(data),len(windows)])
    
    for i in range(len(data)):
        for j in range(len(data[i][1])):
            for k in range(len(windows)):
                if data[i][1][j]<windows[k][1] and data[i][1][j]>windows[k][0]:
                   peaks[i,k]=data[i][2][j]
                   
    return peaks
    
##############################################################################################################################
# A function for manual integration of peaks (gives peak areas)
    
def get_peak_area_manual(filename,chromatogram,target_rt,window=[]):
    
    #Chromatogram is a numpy array where the first column gives the retention time and the second column gives the signal intensity
    #target_rt is the retention time target for the peak in question
    #window is a list where the two elements correspond to times that bracket the target_rt
    
    if window==[]:
        window.append(target_rt-0.1)
        window.append(target_rt+0.1)
    #Loop through the chromatogram and pull out all times and intensities that are within the specified window

    x=[]
    y=[]
    
    for i in range(len(chromatogram)):
        if chromatogram[:,0][i]<=window[1] and chromatogram[:,0][i]>=window[0]:
            x.append(chromatogram[:,0][i])
            y.append(chromatogram[:,1][i])
    
    #Determine the slope and intercept of the line that connects the two end points. 
    line_params=numpy.polyfit([window[0],window[1]],[y[0],y[-1]],1)
    line_points=numpy.array(x)*line_params[0]+line_params[1]
    #print([filename,line_params])
    
    #Create a vector where the intensity at each point is the lower of the chromatogram and the fit line.  
    #This way if the fit line crosses the chromatogram the area below the fit line and above the chromatogram is not double counted 
    
    lower_y=[]
    
    for i in range(len(y)):
        lower_y.append(numpy.min([line_points[i],y[i]]))
    
    manual_area=numpy.trapz(y,x)-numpy.trapz(lower_y,x)
        
    return [manual_area, filename, line_params]

##############################################################################################################################
# A function for visually checking manual peak integration windows
    
def check_manual_peaks(chromatogram_list,line_parameters,filenames,data,zoom=[]):

    cmap=plt.cm.viridis
    cmaplist=[cmap(i) for i in range(cmap.N)]
    no_colors=int(len(cmaplist)/len(chromatogram_list))
    
    plt.figure()
    j=0
    legend=[]
    for i in chromatogram_list:
        plt.plot(data[i][0][:,0],data[i][0][:,1],'-',color=cmaplist[no_colors*j])
        plt.plot(data[i][0][:,0],data[i][0][:,0]*line_parameters[j][0]+line_parameters[j][1],'--',color=cmaplist[no_colors*j])
        legend.append(filenames[i])
        legend.append(filenames[i]+' Line')
        j=j+1
    plt.legend(legend,loc='lower right')
    
    if zoom!=[]:
        plt.xlim([zoom[0],zoom[1]])
    
        
    #plt.legend(legend,bbox_to_anchor=(1.001,1),borderaxespad=0,prop={'size':5})
    #plt.tight_layout(pad=4)
    
##############################################################################################################################
# A function for getting multiple peaks by manual integration.
    
def get_peaks_manual(files,filenames,data,target_rts,windows):
    #data= gc data decoupled from the filename (but in the same order). It must be the full data vector, it cannot be just the sample data
    #Loop through each data file and find the peaks associated with each retention time
    #files are the numbers of the relevant files
    peaks=numpy.zeros([len(files),len(windows)])
    
    names=[]

    j=0               
    for i in files:
        for k in range(len(windows)):
            peaks[j,k]=get_peak_area_manual(filenames[i],data[i][0],target_rts[k],windows[k])[0]
        j=j+1
        names.append(filenames[i][:-4])
        
    return [peaks, names]

##############################################################################################################################
# A function for plotting specified chromatograms

def plot_samples(sample_list,filenames,data):
    plot_legend=[]
    cmap=plt.cm.viridis
    cmaplist=[cmap(i) for i in range(cmap.N)]
    no_colors=int(len(cmaplist)/len(sample_list))
    
    j=0
    for i in range(len(filenames)):
        if filenames[i] in sample_list:
            plt.plot(data[i][0][:,0],data[i][0][:,1],'-',color=cmaplist[j*no_colors])
            j=j+1
            plot_legend.append(filenames[i][:-4])
    
    plt.legend(plot_legend)
    plt.xlabel('Retention Time (min)')
    plt.ylabel('Signal Intensity')
