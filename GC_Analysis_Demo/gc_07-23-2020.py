#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 13:24:38 2020

@author: jgreenhalgh2
"""




import pickle
import matplotlib.pyplot as plt 
import numpy
import gc_tools
import csv
plt.close('all')

#MAKE SURE TO UPDATE THE EXPERIMENT DATE FOR NEW FILES!########################
experiment_date='07-23-2020' #Must be in mm-dd-yyyy
extra_tag=''

###############################################################################

gc_data=pickle.load(open('gc_'+experiment_date+extra_tag+'_data.p','rb'))
filenames=gc_data[0]    #Extract filenames
data=gc_data[1]         #Extract gc data (excluding filenames)

## Make a dictionary of filenames:data

data_dict={}
for i in range(len(filenames)):
    data_dict.update({filenames[i]:data[i]})

#Specify Volumes in mL and stock concentrations of internal standards
media_vol=51
dodecane_vol=10
int_standard_vol=0.15
antibiotic_vol=0.05
culture_vol=0.5
IPTG_vol=0.5

total_vol=antibiotic_vol+media_vol+dodecane_vol+culture_vol+IPTG_vol    #Made a change to this after first analysis and post round modeling (had been multiplying antibiotic_vol by 3 accidentally)   

                         #C3     C5     C7       C9     C11     C13     C15     C17
stock_concs=numpy.array([10.008, 10.008, 9.998, 10.024, 10.018, 10.004, 10.27, 9.96]) #mg/mL

media_conc=int_standard_vol*stock_concs/total_vol #mg/mL

#### Specify Integration Windows and retention times 
windows=[[8.45, 8.57],          #propanol
 [13.142999999999999, 13.263],  #butanol
#windows=[   
 [14.51, 14.76],    #pentanol
 [15.9, 16.13],    #hexanol
 [17.3, 17.5],     #heptanol
 [18.53, 18.79],    #octanol
 [19.62, 19.86],     #nonanol
 [20.66, 20.88],      #decanol
 [21.63, 21.83],      #undecanol
 [22.55, 22.74],      #dodecanol
 [23.42, 23.57],    #tridecanol
 [24.27, 24.41],      #tetradecanol
 [25.07, 25.21],    #pentadecanol
 [25.85, 26.0],     #hexadecanol
 [26.59, 26.74]]     #heptadecanol

standard_rts=[8.51, #Propanol
 13.203,            #Butanol
#standard_rts=[
 16.22,             #Pentanol
 17.378,            #Hexanol
 18.778,            #Heptanol
 20.018,            #Octanol
 21.13,             #Nonanol
 22.175,            #Decanol
 23.167,            #Undecanol
 24.105,            #Dodecanol
 25.001,            #Tridecanol
 25.864,            #Tetradecanol
 26.686,            #Pentadecanol
 27.47,             #Hexadecanol
 28.363]            #Heptadecanol

###############################################################################
sample_list=[]
sample_data=[]
sample_labels=[]
for i in range(len(filenames)):
    if filenames[i][0] in ['M','A','E','B','b','m','0','1','2']:# or filenames[i][1]=='8': #or filenames[i][:2]=='c8':
        sample_list.append(filenames[i])
        sample_data.append(data[i])
        sample_labels.append(filenames[i].split('.')[0])

sample_list_numbers=[]
for i in range(len(filenames)):
    if filenames[i] in sample_list:
        sample_list_numbers.append(i)

sample_peaks_manual=gc_tools.get_peaks_manual(sample_list_numbers,filenames,data,standard_rts,windows)

###############################################################################
## Plot chromatographs

plt.figure()
for sample in sample_list:
    plt.plot(data_dict[sample][0][:,0],data_dict[sample][0][:,1])         

plt.legend(sample_labels)  

##Uncomment to show where windows cross the chromatograph
for i in range(len(windows)):
    w=windows[i]
    plt.vlines(w[0],0,200000,linestyles='dashed')
    plt.vlines(w[1],0,200000,linestyles='dashed')
    plt.text((w[0]+w[1])/2,200000,'c'+str(i+3))
    plt.text((w[0]+w[1])/2,0,'c'+str(i+3))




##################################################################################
#Convert areas to concentrations

sample_evens=numpy.zeros([len(sample_list),7])
sample_odds=numpy.zeros([len(sample_list),8])

sample_evens_manual=numpy.zeros([len(sample_list),7])
sample_odds_manual=numpy.zeros([len(sample_list),8])


for i in range(len(sample_evens[0])):
    sample_evens_manual[:,i]=sample_peaks_manual[0][:,2*i+1]
    
for i in range(len(sample_odds[0])):
    sample_odds_manual[:,i]=sample_peaks_manual[0][:,2*i]

even_concs_manual=numpy.zeros([len(sample_list),len(sample_evens[0])])

for i in range(len(sample_evens[0])):
    for j in range(len(sample_list)):
        #mg/L
        even_concs_manual[j,i]=1000*sample_evens_manual[j,i]/(numpy.average([sample_odds_manual[j,i+1],sample_odds_manual[j,i]]))*(numpy.average([media_conc[i],media_conc[i+1]]))

################################################################################
##Sort the array

s=sorted(zip(sample_list,even_concs_manual))  
sample_list,even_concs_manual=map(list,zip(*s))  
even_concs_manual=numpy.array(even_concs_manual)  
sample_labels.sort()
##############################################################################
#Plot even concentrations (manual) stacked

plt.figure()

index=numpy.arange(len(sample_list))
#index=numpy.arange(len(sample_list))
bw=1/(len(even_concs_manual[0])+1)
cmap=plt.cm.viridis
cmaplist=[cmap(i) for i in range(cmap.N)]
no_colors=int(len(cmaplist)/len(even_concs_manual[0]))

bottoms=numpy.zeros(len(sample_list))

for i in range(1,len(even_concs_manual[0])):
    plt.bar(index,even_concs_manual[:,i],color=cmaplist[i*no_colors],label='C'+str(4+2*i),bottom=bottoms)

    bottoms=bottoms+even_concs_manual[:,i]
 
plt.legend()    
plt.xticks(index,sample_labels,size=6,rotation=90)
plt.xlabel('Samples')
plt.ylabel('mg/L')      
plt.title('')

#################################################################################
#Average replicates:
#
#avg_labels=[]
#
#for label in sample_labels:
#    if label[0]=='A':
#        label=label[:6]
#    elif label[:6]=='MA-ACR':
#        label=label[:6]
#    elif label[:7] in ['MAB-ACR','MAT-ACR']:
#        label=label[:7]
#    elif label[0]=='b':
#        label=label
#      
#    if label not in avg_labels:
#        avg_labels.append(label)
##
#     
#  
##Loop through the data, append samples that are replicates to a list 
#        
#avg_concs=[]
#err_concs=[]
#for label in avg_labels:
#    samples=[]
#    for i in range(len(sample_labels)):
#        if sample_labels[i][:6]==label[:6]:
#            samples.append(even_concs_manual[i])
#    samples=numpy.array(samples)
#    avg_concs.append(numpy.average(samples,0))
#    err_concs.append(numpy.std(samples,0))
#    
#avg_concs=numpy.array(avg_concs)
#err_concs=numpy.array(err_concs)
#    
################################################################################
###Plot averaged even concentrations (manual) stacked
##
#plt.figure()
#
#index=numpy.arange(len(avg_labels))
#bw=1/(len(avg_concs[0])+1)
#cmap=plt.cm.viridis
#cmaplist=[cmap(i) for i in range(cmap.N)]
#no_colors=int(len(cmaplist)/len(avg_concs[0]))
#
#bottoms=numpy.zeros(len(avg_labels))
#
#for i in range(1,len(avg_concs[0])):
#    plt.bar(index,avg_concs[:,i],color=cmaplist[i*no_colors],label='C'+str(4+2*i),bottom=bottoms,yerr=err_concs[:,i])
#    bottoms=bottoms+avg_concs[:,i]
# 
#plt.legend()    
#plt.xticks(index,avg_labels,size=6)
#plt.xlabel('Samples')
#plt.ylabel('mg/L')      
#plt.title('')
        

#
#################################################################################
##Write the concentration data to a csv file. 
with open('gc_'+experiment_date+extra_tag+'.csv','w') as csvFile:
    csvwriter=csv.writer(csvFile)
#    for i in range(len(concentrations)):
#        csvwriter.writerow(concentrations[i])
    csvwriter.writerow(['Sample','C4','C6','C8','C10','C12','C14','C16','Date','Sum C6-C16'])
    for i in range(len(sample_labels)):
        row=[sample_labels[i]]

        for j in range(len(even_concs_manual[0])):
            row.append(even_concs_manual[i][j])     #append each species conc. to the row
        
        row.append(experiment_date)
        row.append(sum(even_concs_manual[i][1:])) #Sum titers for C6-C16
        csvwriter.writerow(row)