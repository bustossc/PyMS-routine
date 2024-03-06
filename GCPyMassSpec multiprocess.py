# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 14:12:00 2020

@author: bustosc
"""

#Routine for detecting and aligning peaks from GC-MS files in format CDF with multiproccessing using PyMassSpec library based on PyMS.
#Created by Carlos Bustos from scripts of PyMassSpec (Dominic Davis-Foster) and EasyGC (David Kainer).



# In[1]:
#Define parameters    

import sys, os
from pyms.BillerBiemann import BillerBiemann, num_ions_threshold, rel_threshold
from pyms.Experiment import Experiment
from pyms.GCMS.IO.ANDI import ANDI_reader
from pyms.IntensityMatrix import build_intensity_matrix_i
from pyms.Noise.SavitzkyGolay import savitzky_golay
from pyms.Peak.Function import peak_sum_area, peak_top_ion_areas
from pyms.TopHat import tophat

from pyms.Noise.Analysis import window_analyzer


#Indicate where your CDF files are located i.e. "C:/Users/.../" "/Users/../"
data_directory = "/Users/"

# Define the data files to process and folders, single or list (without extension)
expr_codes = ["filename1", "filename2"]


#Peak detection parameters
window = 9;     # width of window over which local ion maxima are detected
scans = 3;      # distance at which locally apexing ions can be combined into one peak
r = 5;          # min percentage of mass intensity for each ion peak relative to max peak intensity (below, the ion is removed from the peak)
n = 3;          # min number of ions with intensity above a threshold in a peak (below, the peak is removed)
noise_mult = 2  # peak ion intensity must be at least this multiple of noise level to be included in 'n'
top_ions = 5    # Number of the most abundant ions to be included in 'aligned_ions' file

#Align parameters
Dw = 1.1           # Within state (local) alignment rt modulation [s] this is the tolerance of RT shift. Cannot be an integer!
Gw = 0.35          # Within state local alignment gap penalty. Lower G is preferable as higher G favours peak-mixing
Db = 1             # Between state (global) rt modulation [s]
Gb = 0.30          # Between state (global) gap penalty

comm=10          #Minimum number of samples with a given peak to include that peak in the list

#Mass spectra range
lo_mass = 40
hi_mass = 340 

#Use the same retention time range for all experiments
lo_rt_limit = "3m"
hi_rt_limit = "10.75m"



# In[2]:
#Prefix for the folders and files of the output. No need to change this.

output_prefix = "w"+str(window)+"s"+str(scans)+"r"+str(r)+"i"+str(n)+"n"+str(noise_mult)
file_prefix = "com"+str(comm)+output_prefix+"d"+str(Dw)+"g"+str(Gw)

#Name of the folder where to store the output
output_directory = data_directory + "PyMassSpec_out_" + output_prefix+"/"


    
# In[3]:
#Functions for peak detection
#Loop over the experiments and perform the peak detaction


#Detect peaks in one file
def detect_one_code(code):

	andi_file = data_directory + code + ".CDF"

	data = ANDI_reader(andi_file)

	im = build_intensity_matrix_i(data)
	im.crop_mass(lo_mass, hi_mass)
	n_scan, n_mz = im.size

	# Preprocess the data (Savitzky-Golay smoothing and Tophat baseline detection)
	for ii in range(n_mz):
		ic = im.get_ic_at_index(ii)
		ic1 = savitzky_golay(ic)
		ic_smooth = savitzky_golay(ic1)  # Why the second pass here?
		ic_bc = tophat(ic_smooth, struct="1.5m")
		im.set_ic_at_index(ii, ic_bc)
 
	# noise level calculation
	tic = data.tic
	tic1 = savitzky_golay(tic)
	tic2 = tophat(tic1, struct="1.5m")	
	noise_level = window_analyzer(tic2)
	print("Noise level in TIC: ",noise_level)      

	# Peak detection
	pl = BillerBiemann(im, window, scans)
	print("Initial number of Peaks found:", len(pl))
    
	# Trim the mass spec of each peak by relative intensity
	apl = rel_threshold(pl, percent=r)

	# Trim the peak list by noise threshold
	peak_list = num_ions_threshold(apl, n, cutoff=noise_level*noise_mult)

	print("\t -> Number of Peaks found:", len(peak_list))

	print("\t -> Executing peak post-processing and quantification...")

	# Set the mass range, remove unwanted ions and estimate the TIC and top ions peak area
	# For peak alignment, all experiments must have the same mass range

	for peak in peak_list:
		#Crop mass range
		#peak.crop_mass(40, 340)
		#Remove mass
		#peak.null_mass(147)

		area = peak_sum_area(im, peak)
		peak.area = area
		area_dict = peak_top_ion_areas(im, peak, top_ions)
		peak.ion_areas = area_dict

	# Create an Experiment
	expr = Experiment(code, peak_list)


	print(f"\t -> Selecting retention time range between '{lo_rt_limit}' and '{hi_rt_limit}'")

	expr.sele_rt_range([lo_rt_limit, hi_rt_limit])

	# Save the experiment to disk.
	output_file = output_directory + f"{code}.expr"
	print(f"\t -> Saving the result as '{output_file}'")

	expr=expr.dump(output_file)
	return expr



import multiprocessing as mp
from multiprocessing import freeze_support    
import tqdm


#Peak detection with multiprocessing
def detect_peaksID(runs):
    if os.path.isdir(data_directory) == False:
        os.mkdir(data_directory)
	    
    #Can change the number of files per batch of multiprocess    
    pool = mp.Pool(8)

#    pool.map(detect_one_code, runs)
    for _ in tqdm.tqdm(pool.imap(detect_one_code, runs), total=len(runs)):
        pass
    


# In[4]:
#Function for aligning peaks from different samples
    
from pyms.Experiment import load_expr
from pyms.DPA.PairwiseAlignment import PairwiseAlignment, align_with_tree
from pyms.DPA.Alignment import exprl2alignment

# do the alignment
def align(runs):
    print('Aligning experiments')
    expr_list = []
    expr_dir = output_directory
    for code in runs:
        file_name = os.path.join(expr_dir, code + ".expr")
        expr = load_expr(file_name)
        expr_list.append(expr)
    F1 = exprl2alignment(expr_list)
    T1 = PairwiseAlignment(F1, Dw, Gw)
    A1 = align_with_tree(T1, min_peaks=comm)

    top_ion_list = A1.common_ion()

    A1.write_csv(expr_dir+file_prefix+'_aligned_rt.csv', expr_dir+file_prefix+'_aligned_area.csv')
    A1.write_common_ion_csv(expr_dir+file_prefix+'_area_common_ion.csv', top_ion_list)
    A1.write_ion_areas_csv(expr_dir+file_prefix+'_aligned_ions.csv')



# In[5]:
#Run peak detection and alignment
    
if __name__ == '__main__':
    freeze_support()
    detect_peaksID(expr_codes)
    align(expr_codes)



# In[6]:
#run for plot a GC file

import matplotlib.pyplot as plt
from pyms.Display import *

#Change filename
andi_file = data_directory + "filename"
data = ANDI_reader(andi_file)
data.info()

print("default intensity matrix, bin interval = 1, boundary +/- 0.5")
im = build_intensity_matrix_i(data)

tic = data.tic

tic1 = savitzky_golay(tic)
tic2 = tophat(tic1, struct="1.5m")
noise_level = window_analyzer(tic2)
print("Noise level in TIC: ",noise_level) 

fig, ax = plt.subplots(1, 1, figsize=(10,5))  

plot_ic(ax, tic1, label  = 'TIC for filename')

print('\007')    
