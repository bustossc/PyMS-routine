# PyMS-routine

This script uses the libraries PyMS and PyMassSpec to extract and align peaks from multiple GC-MS files in CDF format. Multiprocessing capabilities are added for faster peak extraction.

PyMassSpec (https://github.com/PyMassSpec/PyMassSpec) which is forked from the original PyMS Repository: https://github.com/ma-bio21/pyms. Originally by Andrew Isaac, Sean O'Callaghan and Vladimir Likić. The original publication can be found here: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-115


######################################
#INSTALLATION

Update microsoft C++ build tools
https://visualstudio.microsoft.com/visual-cpp-build-tools/

Install Anaconda for Python3

Install spyder in a virtual environment
-go to Anaconda navigator an create a new environment (don't use root environment), installing spyder

Activate virtual environment (in Anaconda console)
>activate venv_name

Install PyMassSpec (in Anaconda console)
>python -m pip install PyMassSpec

For Updating
>python -m pip install PyMassSpec -U

Install package 'tqdm' in conda terminal
>pip install tqdm

Launch spyder inside the virtual environment in navigator (or from windows menu if available)



#######################################
#USING PyMassSpec

Open script 'GCPyMassSpec multiprocess.py'
Save as script 'GCPyMassSpec multiprocess(name your experiment).py' (only after this you can modify the script without changing the original one)

Modify 'data_directory' (where your CDF files are located)
Place file names in 'expr_codes' object (without '.CDF' extension)


Change parameters to desired values

To detect and align peaks
Run script with 'detect_peaks(expr_codes)' and 'align(expr_codes)'

To only align peaks from several already extracted samples
Run script with 'detect_peaks(expr_codes)' off (# to comment), and activate 'align(expr_codes)'

If there are too few peaks detected, reduce noise multiplier or minimum ions (parameters 'noise_mult' and n'). Increase those parameters if there are many unrelevant, small peaks. Changing ion percentage('r') has a similar effect, but 5% is normally OK.

If alignment is not great change Dw and Gw, higher Dw favors aligning peaks that are further away, higher Gw favors peak mixing of peaks (it is the penalty for gaps in the alignment list).


#########################################
#Output

The script delivers four output files in .csv format. 

"..._aligned_areas.csv"   -Contains the TIC areas of each aligned peak for each chromatogram with the average retention time per peak. "..._aligned rt.csv"      -Contains the retention times for each chromatogram.
"..._aligned_ions.csv"    -Contains the main ions used for separating each peak.
"..._area_common_ion.csv" -Containes the area of one ion used for quantficiation (this ion is not necessarily the dominant ion).
