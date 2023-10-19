from numba import jit
import numpy as np
from scipy.linalg import lstsq
from matplotlib.figure import figaspect
from math import acos, pi, degrees, sin, cos, sqrt, log, log10, tanh, remainder, ceil
from scipy.interpolate import interp2d
import random
import codecs as cd
from zipfile import ZipFile
import shutil
import time 
import threading


import ipywidgets as widgets
from ipywidgets import VBox, HBox
import codecs as cd
import matplotlib.pyplot as plt
import copy as copy
import os, re
from scipy.linalg import lstsq
from numpy.polynomial.polynomial import polyfit
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib as matplotlib
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from scipy.interpolate import griddata

mu0 = 4*pi*1e-7
ms0 = 270000.0
kb = 1.3806503e-23
tau = 10e-9
roottwohffield = 2**(0.5)  
eul = 1.78107 #exp(eulers const)

#class

class ElapsedTimeThread(threading.Thread):
    """"Stoppable thread that prints the time elapsed"""
    def __init__(self):
        super(ElapsedTimeThread, self).__init__()
        self._stop_event = threading.Event()

    def stop(self):
        self._stop_event.set()

    def stopped(self):
        return self._stop_event.is_set()

    def run(self):
        thread_start = time.time()
        while not self.stopped():
            #print("\rElapsed Time {:.0f} seconds".format(time.time()-thread_start), end="")
            #include a delay here so the thread doesn't uselessly thrash the CPU
            time.sleep(1)


#find and upzip zipped file 
def open_zip_file():
    path = os.getcwd()
    try:
        os.mkdir('raw_data')
    except:
        print('This step has already been carried out and a raw_data folder exists. To add more files, restart the kernel and re-run the previous cells too.')
    for item in os.listdir(path): # loop through items in dir
        
        if item.endswith('.zip'): # check for ".zip" extension
            file_name = os.path.abspath(item) # get full path of files
            zip_ref = ZipFile(file_name) # create zipfile object
            zip_ref.extractall(os.path.join(path+ os.sep,'raw_data')) # extract file to dir
            zip_ref.close() # close file
            #os.remove(file_name) # delete zipped file
    return


def find_files(X):
    path = os.getcwd()
    path = os.path.join(path+os.sep,'raw_data')

    #sample names
    
    names = []
    for root, dir, file in os.walk(path): #loop file
        i =0
        for f in file:
            if re.match('.*.frc', f):
                names.append(f[:-4])
    print('Sample to test:', names[0])
    X['names'] = names



    forc_file = []
    for root, dir, file in os.walk(path):
        i =0
        for f in file:
            if re.match('.*.frc', f):
                forc_file.append(os.path.join(root + os.sep,f))
                i+=1
    X['fn'] = forc_file[0]

    nrm_file = []
    for root, dir, file in os.walk(path):
        i =0
        for f in file:
            if re.match('.*nrm.dat', f):
                nrm_file.append(os.path.join(root + os.sep,f))
                i+=1
    if (i==0):
        print('No NRM files found')

    X['nrm_files'] = nrm_file

    sirm_file = []
    for root, dir, file in os.walk(path):
        i =0
        for f in file:
            if re.match('.*sirm.dat', f):
                sirm_file.append(os.path.join(root + os.sep,f))
                i+=1
    if (i==0):
        print('No SIRM files found')

    X['sirm_files'] = sirm_file    
    X = field_range(X) # move function from below inside
    return X

def find_files2(X):
    path = os.getcwd()
    #path = os.path.join(path+os.sep,'raw_data')

    #sample names
    search1 = 0
    search2 = 0
    search3 = 0
    while (search1 == 0) or (search2 == 0) or (search3 == 0):
        try: 
            input_name = input('Please enter the sample name: ')
            names = []
            i =0
            for root, dir, file in os.walk(path): #loop file
                
                for f in file:
                    if re.match('{}.frc'.format(input_name), f):
                        names.append(f[:-4])
            
            X['names'] = names
            print('Sample to test:', names[0])

            forc_file = []
            i =0
            for root, dir, file in os.walk(path):
                
                for f in file:
                    if re.match('{}.frc'.format(input_name), f):
                        forc_file.append(os.path.join(root + os.sep,f))
                        search1 = 1
                        i+=1
        
                        
            X['fn'] = forc_file[0]
            if (i==0):
                print('No FORC diagram can be found for sample {}'.format(input_name))

            nrm_file = []
            i =0
            for root, dir, file in os.walk(path):
            
                for f in file:
                    if re.match('{}.nrm'.format(input_name), f):
                        nrm_file.append(os.path.join(root + os.sep,f))
                        search2 = 1
                        i+=1
            if (i==0):
                print('No NRM file can be found for sample {}'.format(input_name))

            X['nrm_files'] = nrm_file

            sirm_file = []
            i =0 #moved i=0 outside to avoid resetting when search next folder 
            for root, dir, file in os.walk(path):
                
                for f in file:
                    if re.match('{}.sirm'.format(input_name), f):
                        sirm_file.append(os.path.join(root + os.sep,f))
                        search3 = 1
                        i+=1
            if (i==0):
                print('No SIRM file can be found for sample {}'.format(input_name))

            X['sirm_files'] = sirm_file    

            #search for results file, if not exisit make it, and close and if exist see how many prev times values run
            file_exists = os.path.exists('paleo_results.dat')
            sample_copy = 1    
            if (file_exists == True):
                #see how many times
                bigf = open('paleo_results.dat', 'r')
                bigf.readline()

                for my_line in bigf:
                    line = my_line.split('\t') #tab delimnated file
                    if (line[0] == names[0]):
                        sample_copy +=1
                        #get variables from here
                        X['min_field'] = float(line[10])
                        X['max_field'] = float(line[11])
                        X['growth_rate'] = float(line[12])
                        X['curie_t'] = float(line[13])
                        X['afval'] = float(line[14])
                        X['SF'] = float(line[15])
                        X['reset_limit_hc'] = float(line[16])
                        X['reset_limit_hi'] = float(line[17])   
                        X['sf_list_correct'] = float(line[18])             
            
            else:
                #create file
                bigf = open('paleo_results.dat', 'w') 
                bigf.write('Sample name \t Repeat \t AF steps \t Range \t AF min \t AF max \t mean PI \t std mean \t median \t iqr median \t Min field \t Max field \t growth time \t curie temp \t AF value \t SF \t max hc \t max hi \t SF_factor \n')
            bigf.close() #close regardless
            X['sample_copy'] = sample_copy
            if (os.path.exists(os.path.join(os.getcwd(), 'sample_{}_V_{}'.format(names[0], str(X['sample_copy']))))):
                pass
            else:
                os.mkdir('sample_{}_V_{}'.format(names[0], str(X['sample_copy'])))
            
            print('Attempt {} of running sample {}'.format(sample_copy,names[0]))
            X = field_range(X) # move function from below inside
        except:
            time.sleep(1)
            print('No files for that sample name, re-enter sample name')

    return(X)

#field boundaries
def field_range(X):
    if (X['sample_copy'] > 1): #only if second go 
        field_1 = input('The most recent field range used for sample {} was {} \u03BCTto {} \u03BCT. If you want to re-use these variables enter K, to change them enter any other charactor'.format(X['names'][0], X['min_field'], X['max_field']))
    else:
        field_1 = 'L'
    if (field_1 != 'K') or (X['sample_copy'] < 2):
        X['min_field'] = 40
        X['max_field'] = 60
        suggest_field = input("The standard field bounds are {} \u03BCT to {} \u03BCT, if you want to keep these enter K, else enter any other charactor:".format(X['min_field'], X['max_field']))
        if (suggest_field != 'K'):
            while True:
                    
                    minf = (input("Pick the lower bound of field range to be tested in \u03BCT:" ))
                
                    try:
                        minf = float(minf)
                        if (minf > 0):# (minf <= maxSF)
                            #print('in bounds')
                            break
                    except ValueError:
                        print('Not a number')
                        True

            while True:
                    maxf = (input("Pick the upper bound of field range to be tested in \u03BCT:" ))
                
                    try:
                        maxf = float(maxf)
                        if (maxf > minf):
                            #print('in bounds')
                            break
                    except ValueError:
                        print('Not a number')
                        True
            print('Expected field range: {} \u03BCT to {} \u03BCT'.format(minf, maxf))
            X['max_field'] = maxf
            X['min_field'] = minf
        else:
            pass

    return(X)


#process FORC data for each forc to get nested results 
def proccess_all(X):

    #4
    #process data
    style = {'description_width': 'initial'}
    fn = X['fn'] 
    sample, unit, mass = sample_details(fn)

    sample_widge = widgets.Text(value=sample,description='Sample name:',style=style) 
    prop_title = widgets.HTML(value='<h3>Sample preprocessing options</h3>')
    mass_title = widgets.HTML(value='To disable mass normalization use a value of -1')
    if mass == "N/A":
        mass_widge = widgets.FloatText(value=-1, description = 'Sample mass (g):',style=style)
    else:
        mass_widge = widgets.FloatText(value=mass, description = 'Sample mass (g):',style=style)
    mass_widge1 = HBox([mass_widge,mass_title])

    X["sample"] = sample_widge
    X["mass"] = mass_widge
    X["unit"] = unit

    H, Hr, M, Fk, Fj, Ft, dH = parse_measurements(X["fn"])
    Hcal, Mcal, tcal = parse_calibration(X["fn"])
    Hc1, Hc2, Hb1, Hb2 = measurement_limts(X)

 

    X["H"] = H
    X["Hr"] = Hr
    X["M"] = M
    X["dH"] = dH
    X["Fk"] = Fk
    X["Fj"] = Fj
    X["Ft"] = Ft
    X["Hcal"] = Hcal
    X["Mcal"] = Mcal
    X["tcal"] = tcal
    X["Hc1"] = Hc1
    X["Hc2"] = Hc2
    X["Hb1"] = Hb1
    X["Hb2"] = Hb2
    slope = 70.
    X["slope"] = slope 

    if X['unit']=='Cgs': 
        X = CGS2SI(X)
    
    X = drift_correction(X) 
   
    X = slope_correction(X)
    X = remove_fpa(X) 
 
    X = remove_lpa(X)
    X = lowerbranch_subtract(X)

    return(X)
    
def prod_FORCs(X): #process and name forc from file 
    name = X['fn'].split("\\")
    name2 = name[-1].split('.')
    X['name'] = name2[0]


    maxSF = 8 
    X['maxSF1'] = maxSF
   
    X = create_arrays(X, maxSF)
    X = nan_values2(X)
   
    for sf in range(2, maxSF+1):
        X = calc_rho(X, sf)
        sf+=1   
    X = nan_values(X, maxSF) 
    X = rotate_FORC(X)
    return(X)    
    



def parse_header(file,string):

    output=-1 
    with cd.open(file,"r",encoding='latin9') as fp:
        for line in lines_that_start_with(string, fp): 
            idx = line.find('=') 
            if idx>-1.: 
                output=float(line[idx+1:]) 
            else: 
                idx = len(string) 
                output=float(line[idx+1:])  

    return output


def parse_measurements(file):


    dum=-9999.99 
    N0=int(1E6) 
    H0=np.zeros(N0)*np.nan #initialize NaN array to contain field values
    M0=np.zeros(N0)*np.nan #initialize NaN array to contain magnetization values
    H0[0]=dum #first field entry is dummy value
    M0[0]=dum #first magnetization entry is dummy value 

    count=0 #counter to place values in arrays
    with cd.open(file,"r",encoding='latin9') as fp: 
        for line in find_data_lines(fp): 
            count=count+1 #increase counter
            idx = line.find(',') #no comma indicates a blank linw
            if idx>-1: #line contains a comma
                #print(line)
                H0[count]=float(line[0:idx]) #assign field value (1st column)
                line=line[idx+1:] #remove the leading part of the line (only characters after the first comma remain)
                idx = line.find(',') #find next comman
                if idx>-1: #comma found in line
                    M0[count]=float(line[0:idx]) #read values up to next comma (assumes 2nd column is magnetizations)
                else: #comma wasn't found   
                    M0[count]=float(line) # magnetization value is just the remainder of the line 
            else:
                H0[count]=dum #line is blank, so fill with dummy value
                M0[count]=dum #line is blank, so fill with dummy value

    idx_start=np.argmax(H0!=dum) #find the first line that contains data            
    M0=M0[idx_start-1:-1] #strip out leading dummy values from magnetizations, leaving 1 dummy at start of vector           
    M0=M0[~np.isnan(M0)] #remove any NaNs at the end of the array
    H0=H0[idx_start-1:-1] #strip out leading dummy values from magnetizations, leaving 1 dummy at start of vector
    H0=H0[~np.isnan(H0)] #remove any NaNs at the end of the array

    ## determine indicies of each FORC
    idxSAT = np.array(np.where(np.isin(H0, dum))) #find start address of each blank line
    idxSAT = np.ndarray.squeeze(idxSAT) #squeeze into 1D
    idxSTART = idxSAT[1::2]+1 #find start address of each FORC
    idxEND = idxSAT[2::2]-1 ##find end address of each FORC

    
    #Extract first FORC to initialize arrays 
    M=M0[idxSTART[0]:idxEND[0]+1] #Magnetization values
    H=H0[idxSTART[0]:idxEND[0]+1] #Field values
    Hr=np.ones(idxEND[0]+1-idxSTART[0])*H0[idxSTART[0]] #Reversal field values
    Fk=np.ones(idxEND[0]+1-idxSTART[0]) #index number of FORC
    Fj=np.arange(1,1+idxEND[0]+1-idxSTART[0])# measurement index within given FORC

    #Extract remaining FORCs one by one into into a long-vector
    #print(idxSTART[0], idxEND[0]+1)
    #print(i, idxSTART, idxSTART.size)
    #print(len(M))
    for i in range(1,idxSTART.size):
        M=np.concatenate((M,M0[idxSTART[i]:idxEND[i]+1]))
        H=np.concatenate((H,H0[idxSTART[i]:idxEND[i]+1]))
        Hr=np.concatenate((Hr,np.ones(idxEND[i]+1-idxSTART[i])*H0[idxSTART[i]]))
        Fk=np.concatenate((Fk,np.ones(idxEND[i]+1-idxSTART[i])+i))
        Fj=np.concatenate((Fj,np.arange(1,1+idxEND[i]+1-idxSTART[i])))
    
    unit = parse_units(file) #Ensure use of SI units
    
    if unit=='Cgs':
        H=H/1E4 #Convert Oe into T
        Hr=Hr/1E4 #Convert Oe into T
        M=M/1E3 #Convert emu to Am^2

    dH = np.mean(np.diff(H[Fk==np.max(Fk)])) #mean field spacing

    Ft=measurement_times(file,Fk,Fj) #estimated time of each measurement point

    return H, Hr, M, Fk, Fj, Ft, dH
def parse_units(file):

    string = 'Units of measure' #header definition of units
    with cd.open(file,"r",encoding='latin9') as fp: #open the data file (latin9 encoding seems to work, UTF and ASCII don't)
        for line in lines_that_start_with(string, fp): #find the line starting with the setting name
            idxSI = line.find('Hybrid SI') #will return location if string is found, otherwise returns -1
            idxCGS = line.find('Cgs') #will return location if string is found, otherwise returns -1
    
    if idxSI>idxCGS: #determine which unit string was found in the headerline and output
        return 'SI'
    else:
        return 'Cgs'
def parse_mass(file):

    output = 'N/A'
    string = 'Mass' #header definition of units
    with cd.open(file,"r",encoding='latin9') as fp: #open the data file (latin9 encoding seems to work, UTF and ASCII don't)
        for line in lines_that_start_with(string, fp): #find the line starting with the setting name
            idx = line.find('=') #Some file formats may contain an '='
            if idx>-1.: #if '=' found
                output=(line[idx+1:]) #value taken as everything to right of '='
            else: # '=' not found
                idx = len(string) #length of the setting string 
                output=(line[idx+1:])  #value taken as everything to right of the setting name
        
            if output.find('N/A') > -1:
                output = 'N/A'
            else:
                output = float(output)

    return output
def measurement_times(file,Fk,Fj):
   
    unit=parse_units(file) #determine measurement system (CGS or SI)

    string='PauseRvrsl' #Pause at reversal field (new file format, -1 if not available)
    tr0=parse_header(file,string)
    
    string='PauseNtl' #Pause at reversal field (old file format, -1 if not available)
    tr1=parse_header(file,string)

    tr=np.max((tr0,tr1)) #select Pause value depending on file format
    
    string='Averaging time' #Measurement averaging time 
    tau=parse_header(file,string)

    string='PauseCal' #Pause at calibration point
    tcal=parse_header(file,string)

    string='PauseSat' #Pause at saturation field
    ts=parse_header(file,string)

    string='SlewRate' #Field slewrate
    alpha=parse_header(file,string)

    string='HSat' #Satuation field
    Hs=parse_header(file,string)

    string='Hb2' #upper Hb value for the FORC box
    Hb2=parse_header(file,string)

    string='Hb1' #lower Hb value for the FORC box
    Hb1=parse_header(file,string)

    string='Hc2' #upper Hc value for the FORC box (n.b. Hc1 is assumed to be 0)
    Hc2=parse_header(file,string)

    string='NForc' # Numer of measured FORCs (new file format, -1 if not available)
    N0=parse_header(file,string)

    string='NCrv'  # Numer of measured FORCs (old file format, -1 if not available)
    N1=parse_header(file,string)

    N=np.max((N0,N1)) #select Number of FORCs depending on file format

    if unit=='Cgs':
        alpha=alpha/1E4 #convert from Oe to T
        Hs=Hs/1E4 #convert from Oe to T
        Hb2=Hb2/1E4 #convert from Oe to T
        Hb1=Hb1/1E4 #convert from Oe to T

    dH = (Hc2-Hb1+Hb2)/N #estimated field spacing
    
    #now following Elgi's estimate of the measurement time
    nc2 = Hc2/dH

    Dt1 = tr + tau + tcal + ts + 2.*(Hs-Hb2-dH)/alpha
    Dt3 = Hb2/alpha

    Npts=int(Fk.size)
    Ft=np.zeros(Npts)
    
    for i in range(Npts):
        if Fk[i]<=1+nc2:
            Ft[i]=Fk[i]*Dt1+Dt3+Fj[i]*tau+dH/alpha*(Fk[i]*(Fk[i]-1))+(tau-dH/alpha)*(Fk[i]-1)**2
        else:
            Ft[i]=Fk[i]*Dt1+Dt3+Fj[i]*tau+dH/alpha*(Fk[i]*(Fk[i]-1))+(tau-dH/alpha)*((Fk[i]-1)*(1+nc2)-nc2)

    return Ft
def parse_calibration(file):


    dum=-9999.99 #dum value to indicate break in measurement seqence between FORCs and calibration points
    N0=int(1E6) #assume that any file will have less than 1E6 measurements
    H0=np.zeros(N0)*np.nan #initialize NaN array to contain field values
    M0=np.zeros(N0)*np.nan #initialize NaN array to contain magnetization values
    H0[0]=dum #first field entry is dummy value
    M0[0]=dum #first magnetization entry is dummy value 

    count=0 #counter to place values in arrays
    with cd.open(file,"r",encoding='latin9') as fp: #open the data file (latin9 encoding seems to work, UTF and ASCII don't)
        for line in find_data_lines(fp): #does the current line contain measurement data
            count=count+1 #increase counter
            idx = line.find(',') #no comma indicates a blank linw
            if idx>-1: #line contains a comma
                H0[count]=float(line[0:idx]) #assign field value (1st column)
                line=line[idx+1:] #remove the leading part of the line (only characters after the first comma remain)
                idx = line.find(',') #find next comman
                if idx>-1: #comma found in line
                    M0[count]=float(line[0:idx]) #read values up to next comma (assumes 2nd column is magnetizations)
                else: #comma wasn't found   
                    M0[count]=float(line) # magnetization value is just the remainder of the line 
            else:
                H0[count]=dum #line is blank, so fill with dummy value
                M0[count]=dum #line is blank, so fill with dummy value

    idx_start=np.argmax(H0!=dum) #find the first line that contains data            
    M0=M0[idx_start-1:-1] #strip out leading dummy values from magnetizations, leaving 1 dummy at start of vector           
    M0=M0[~np.isnan(M0)] #remove any NaNs at the end of the array
    H0=H0[idx_start-1:-1] #strip out leading dummy values from magnetizations, leaving 1 dummy at start of vector
    H0=H0[~np.isnan(H0)] #remove any NaNs at the end of the array

    ## now need to pull out the calibration points, will be after alternate -9999.99 entries
    idxSAT = np.array(np.where(np.isin(H0, dum))) #location of dummy values
    idxSAT = np.ndarray.squeeze(idxSAT) #squeeze into 1D
    idxSAT = idxSAT[0::2]+1 #every second index+1 should be calibration points

    Hcal=H0[idxSAT[0:-1]] #calibration fields
    Mcal=M0[idxSAT[0:-1]] #calibration magnetizations
    tcal=calibration_times(file,Hcal.size) #estimate the time of each calibratio measurement

    unit = parse_units(file)
    
    if unit=='Cgs': #ensure SI units
        Hcal=Hcal/1E4 #convert from Oe to T
        Mcal=Mcal/1E3 #convert from emu to Am^2

    return Hcal, Mcal, tcal
def calibration_times(file, Npts):
 
    unit=parse_units(file) #determine measurement system (CGS or SI)

    string='PauseRvrsl' #Pause at reversal field (new file format, -1 if not available)
    tr0=parse_header(file,string)
    
    string='PauseNtl' #Pause at reversal field (old file format, -1 if not available)
    tr1=parse_header(file,string)

    tr=np.max((tr0,tr1)) #select Pause value depending on file format
    
    string='Averaging time' #Measurement averaging time 
    tau=parse_header(file,string)

    string='PauseCal' #Pause at calibration point
    tcal=parse_header(file,string)

    string='PauseSat' #Pause at saturation field
    ts=parse_header(file,string)

    string='SlewRate' #Field slewrate
    alpha=parse_header(file,string)

    string='HSat' #Satuation field
    Hs=parse_header(file,string)

    string='Hb2' #upper Hb value for the FORC box
    Hb2=parse_header(file,string)

    string='Hb1' #lower Hb value for the FORC box
    Hb1=parse_header(file,string)

    string='Hc2' #upper Hc value for the FORC box (n.b. Hc1 is assumed to be 0)
    Hc2=parse_header(file,string)

    string='NForc' # Numer of measured FORCs (new file format, -1 if not available)
    N0=parse_header(file,string)

    string='NCrv'  # Numer of measured FORCs (old file format, -1 if not available)
    N1=parse_header(file,string)

    N=np.max((N0,N1)) #select Number of FORCs depending on file format

    if unit=='Cgs':
        alpha=alpha/1E4 #convert from Oe to T
        Hs=Hs/1E4 #convert from Oe to T
        Hb2=Hb2/1E4 #convert from Oe to T
        Hb1=Hb1/1E4 #convert from Oe to T
    
    dH = (Hc2-Hb1+Hb2)/N #estimated field spacing
    
    #now following Elgi's estimate of the measurement time
    nc2 = Hc2/dH
    Dt1 = tr + tau + tcal + ts + 2.*(Hs-Hb2-dH)/alpha
    Dt2 = tr + tau + (Hc2-Hb2-dH)/alpha

    Npts=int(Npts)
    tcal_k=np.zeros(Npts)
    
    for k in range(1,Npts+1):
        if k<=1+nc2:
            tcal_k[k-1]=k*Dt1-Dt2+dH/alpha*k**2+(tau-dH/alpha)*(k-1)**2
        else:
            tcal_k[k-1]=k*Dt1-Dt2+dH/alpha*k**2+(tau-dH/alpha)*((k-1)*(1+nc2)-nc2)

    return tcal_k
def sample_details(fn):

    sample = fn.split('/')[-1]
    sample = sample.split('.')
    
    if type(sample) is list:
        sample=sample[0]

    units=parse_units(fn)
    mass=parse_mass(fn)
  
    return sample, units, mass
def measurement_limts(X):
    
    
    string='Hb2' #upper Hb value for the FORC box
    Hb2=parse_header(X["fn"],string)

    string='Hb1' #lower Hb value for the FORC box
    Hb1=parse_header(X["fn"],string)

    string='Hc2' #upper Hc value for the FORC box
    Hc2=parse_header(X["fn"],string)

    string='Hc1' #lower Hc value for the FORC box
    Hc1=parse_header(X["fn"],string)

    if X['unit']=='Cgs': #convert CGS to SI
        Hc2=Hc2/1E4 #convert from Oe to T
        Hc1=Hc1/1E4 #convert from Oe to T
        Hb2=Hb2/1E4 #convert from Oe to T
        Hb1=Hb1/1E4 #convert from Oe to T  

    return Hc1, Hc2, Hb1, Hb2

#### Unit conversion ####
def CGS2SI(X):
    
    X["H"] = X["H"]/1E4 #convert Oe into T
    X["M"] = X["M"]/1E3 #convert emu to Am2
      
    return X

#### low-level IO routines
def find_data_lines(fp):

    return [line for line in fp if ((line.startswith('+')) or (line.startswith('-')) or (line.strip()=='') or line.find(',')>-1.)]
def lines_that_start_with(string, fp):
 
    return [line for line in fp if line.startswith(string)]

def remove_lpa(X):
    
    #unpack
    Fj = X["Fj"]
    H = X["H"]    
    Hr = X["Hr"]
    M = X["M"]
    Fk = X["Fk"]
    Ft = X["Ft"]
    
    #remove last point artifact
    Nforc = int(np.max(Fk))
    W = np.ones(Fk.size)
    
    for i in range(Nforc):      
        Fj_max=np.sum((Fk==i))
        idx = ((Fk==i) & (Fj==Fj_max))
        W[idx]=0.0
    
    idx = (W > 0.5)
    H=H[idx]
    Hr=Hr[idx]
    M=M[idx]
    Fk=Fk[idx]
    Fj=Fj[idx]
    Ft=Ft[idx]
    Fk=Fk-np.min(Fk)+1. #reset FORC number if required
    
    #repack
    X["Fj"] = Fj
    X["H"] = H   
    X["Hr"] = Hr
    X["M"] = M
    X["Fk"] = Fk
    X["Ft"] = Ft        
    
    return X

def remove_fpa(X):
    
    #unpack
    Fj = X["Fj"]
    H = X["H"]    
    Hr = X["Hr"]
    M = X["M"]
    Fk = X["Fk"]
    Fj = X["Fj"]
    Ft = X["Ft"]
    
    #remove first point artifact
    idx=((Fj==1.0))
    H=H[~idx]
    Hr=Hr[~idx]
    M=M[~idx]
    Fk=Fk[~idx]
    Fj=Fj[~idx]
    Ft=Ft[~idx]
    Fk=Fk-np.min(Fk)+1. #reset FORC number if required
    Fj=Fj-1.
    
    #repack
    X["Fj"] = Fj
    X["H"] = H   
    X["Hr"] = Hr
    X["M"] = M
    X["Fk"] = Fk
    X["Ft"] = Ft        
    
    return X

def drift_correction(X):
  
    #unpack
    M = X["M"]
    Mcal = X["Mcal"]    
    Ft = X["Ft"]
    tcal = X["tcal"]
  
    #perform drift correction
    M=M*Mcal[0]/np.interp(Ft,tcal,Mcal,left=np.nan) #drift correction
  
    #repack
    X["M"] = M
  
    return X

def FORC_extend(X):
    
    Ne = 20 #extend up to 20 measurement points backwards
    
    #unpack
    H = X["H"]    
    Hr = X["Hr"]
    M = X["M"]
    Fk = X["Fk"]
    Fj = X["Fj"]
    dH = X["dH"]
    
    for i in range(int(X['Fk'][-1])):
        M0 = M[Fk==i+1]
        H0 = H[Fk==i+1]
        Hr0 = Hr[Fk==i+1][0]
        
        M1 = M0[0] - (np.flip(M0)[1:]-M0[0])
        H1 = H0[0] - (np.flip(H0)[1:]-H0[0])
            
        if M1.size>Ne:
            H1 = H1[-Ne-1:-1]
            M1 = M1[-Ne-1:-1]
        
        if i==0:    
            N_new = np.concatenate((M1,M0)).size
            H_new = np.concatenate((H1,H0))
            M_new = np.concatenate((M1,M0))
            Hr_new = np.ones(N_new)*Hr0
            Fk_new = np.ones(N_new)
            Fj_new = np.arange(N_new)+1-M1.size
        else:
            N_new = np.concatenate((M1,M0)).size
            H_new = np.concatenate((H_new,H1,H0))
            M_new = np.concatenate((M_new,M1,M0))
            Hr_new = np.concatenate((Hr_new,np.ones(N_new)*Hr0))
            Fk_new = np.concatenate((Fk_new,np.ones(N_new)+i))
            Fj_new = np.concatenate((Fj_new,np.arange(N_new)+1-M1.size))
            
    #pack up variables
    X['H'] = H_new
    X['Hr'] = Hr_new
    X['M'] = M_new
    X['Fk'] = Fk_new
    X['Fj'] = Fj_new
    
    return X

def lowerbranch_subtract(X):

    
    #unpack
    H = X["H"]    
    Hr = X["Hr"]
    M = X["M"]
    Fk = X["Fk"]
    Fj = X["Fj"]
    dH = X["dH"]
    
    Hmin = np.min(H)
    Hmax = np.max(H)


    Nbar = 10
    nH = int((Hmax - Hmin)/dH)
    Hi = np.linspace(Hmin,Hmax,nH*50+1)
    Mi = np.empty(Hi.size)
    
    #perform basic loess
    for i in range(Hi.size):
        idx = (H>=Hi[i]-2.5*dH) & (H<=Hi[i]+2.5*dH)
        Mbar = M[idx]
        Hbar = H[idx]
        Fbar = Fk[idx]
        F0 = np.sort(np.unique(Fbar))
        if F0.size>Nbar:
            F0=F0[-Nbar]
        else:
            F0=np.min(F0)
        idx = Fbar>=F0
        
        p = np.polyfit(Hbar[idx],Mbar[idx],2)
        Mi[i] = np.polyval(p,Hi[i])
    
    Hlower = Hi
    Mlower = Mi
    Mcorr=M-np.interp(H,Hlower,Mlower,left=np.nan,right=np.nan) #subtracted lower branch from FORCs via interpolation

    Fk=Fk[~np.isnan(Mcorr)] #remove any nan
    Fj=Fj[~np.isnan(Mcorr)] #remove any nan
    H=H[~np.isnan(Mcorr)] #remove any nan
    Hr=Hr[~np.isnan(Mcorr)] #remove any nan
    M=M[~np.isnan(Mcorr)] #remove any nan
    Mcorr = Mcorr[~np.isnan(Mcorr)] #remove any nan
    
    #repack
    X["H"] = H    
    X["Hr"] = Hr
    X["M"] = M
    X["Fk"] = Fk
    X["Fj"] = Fj
    X["DM"] = Mcorr
    
    return X

    ###### HELPER FUNCTIONS TO READ FROM FILE

def slope_correction(X):
  
    #unpack
    H = X["H"]
    M = X["M"]
  

    Hidx = H > (X["slope"]/100) * np.max(H)
    p = np.polyfit(H[Hidx],M[Hidx],1)
    M = M - H*p[0]
  
    #repack
    X["M"]=M
  
    return X

def create_arrays(X, maxSF):

    Fk_int = (X['Fk'].astype(int)) #turn to int to use bincount
    counts = np.bincount(Fk_int) #time each no. appears in Fk (no. FORC on)
    max_FORC_len = np.max(counts) #max occurance of a FORC no. = longest FORC length = no. columns
    no_FORC = np.argmax(counts) #max FORC no.   = rows

    H_A = np.zeros((no_FORC, max_FORC_len)) #initialize arrays
    Hr_A = np.zeros((no_FORC, max_FORC_len))
    M_A = np.zeros((no_FORC, max_FORC_len))
    Fk_A = np.zeros((no_FORC, max_FORC_len))
    Rho = np.zeros((maxSF+1, no_FORC, max_FORC_len))
    #initialize zero values
    H_A[0,0] = X['H'][0]
    Hr_A[0,0] = X['Hr'][0]
    M_A[0,0] = X['M'][0]
    Fk_A[0,0] = X['Fk'][0]

    j=0 # just filled first point in first row
    i=0 # start at first row
    for cnt in range(1,len(X['Fk']+1)):
        if (X['Fk'][cnt] == X['Fk'][cnt-1]): #if Fk no is the same, stay on same row and fill data
            j +=1 #add one more to column and repeat
            H_A[i][j] = X['H'][cnt]
            Hr_A[i][j] = X['Hr'][cnt]
            M_A[i][j] = X['M'][cnt]     
        else:
            i +=1 #new row
            j = 0 #set column index back to zero
            H_A[i][j] = X['H'][cnt]
            Hr_A[i][j] = X['Hr'][cnt]
            M_A[i][j] = X['M'][cnt]            
        cnt +=1 #next point
    X['H_A'] = H_A
    X['Hr_A'] = Hr_A
    X['M_A'] = M_A
    X['rho'] = Rho
    X['no_FORC'] = no_FORC
    X['max_FORC_len'] = max_FORC_len
    return(X)


def nan_values2(X):
    H_A = X['H_A']
    Hr_A = X['Hr_A']
    M_A = X['M_A']
    for i in range(len(H_A)):
        for j in range(len(Hr_A[0])):
            if (H_A[i][j] == 0.0):
                H_A[i][j] = 'NaN'
                Hr_A[i][j] = 'NaN'
                M_A[i][j] = 'NaN'
                

    X['H_A'] = H_A
    X['Hr_A'] = Hr_A
    X['M_A'] = M_A
    return(X)

def calc_rho(X, SF):
    no_FORC = X['no_FORC']
    max_FORC_len = X['max_FORC_len']
    H_A = X['H_A']
    Hr_A = X['Hr_A']
    M_A = X['M_A']
    Rho = X['rho']

    for i in range(no_FORC): #find main points
        for j in range(max_FORC_len): #find each j indice
            #locate smoothing grids

            cnt = 0
            h1 = min(i, SF) #row from SF below and SF above
            h2 = min(SF, (no_FORC - i)) #loop over all points
            k1 = min(j, SF) #point to left, 1j if 0 etc or SF is in middle
            k2 = min(SF, (max_FORC_len-j)) #right hand side - either SF or if near edge do total - j (point at)

            A = np.zeros(((h2+h1+1)*(k1+k2+1),6))
            b = np.zeros(((h2+h1+1)*(k1+k2+1)))
            A[:,:] = np.nan
            b[:] = np.nan

            #if (M_A[i][j] != 0. and H_A[i][j] !=0 and Hr_A[i][j] != 0): 
            if (H_A[i][j] > Hr_A[i][j]):
                for h in range((-h1), (h2+1)): #loop over row in smoothing window
                    for k in range((-k1), (k2+1)): #loop over columns in smoothing window
                        if ((j+h+k) >= 0 and (j+k+h) < (max_FORC_len) and (i+h) >= 0 and (i+h) < (no_FORC)): 
                              
                            A[cnt, 0] = 1.
                            A[cnt, 1] = Hr_A[i+h][j+k+h] - Hr_A[i][j]
                            A[cnt, 2] = (Hr_A[i+h][j+k+h] - Hr_A[i][j])**2.
                            A[cnt, 3] = H_A[i+h][j+k+h] - H_A[i][j]
                            A[cnt, 4] = (H_A[i+h][j+k+h] - H_A[i][j])**2.
                            A[cnt, 5] = (Hr_A[i+h][j+k+h] - Hr_A[i][j])*(H_A[i+h][j+k+h] - H_A[i][j])
                            b[cnt] = M_A[i+h][j+k+h]

                            cnt+=1 #count number values looped over
                A = A[~np.isnan(A).any(axis=1)]
                b = b[~np.isnan(b)]
                if (len(A)>=2): #min no. points to need to smooth over
 
                    dmatrix, res, rank, s = lstsq(A,b)
                    Rho[SF][i][j] = (-1.*(dmatrix[5]))/2.

                else:
                    Rho[SF][i][j] = 0.
            else:
                Rho[SF][i][j] = 0.
            j +=1
        i += 1

    X['H_A'] = H_A #repack variables
    X['Hr_A'] = Hr_A
    X['M_A'] = M_A
    X['rho'] = Rho
    X['no_FORC'] = no_FORC
    X['max_FORC_len'] = max_FORC_len
    return(X)
    
    
def nan_values(X, maxSF):
    H_A = X['H_A']
    Hr_A = X['Hr_A']
    Rho = X['rho']
    for i in range(len(H_A)):
        for j in range(len(Hr_A[0])):
            if (H_A[i][j] == 0.0):
                H_A[i][j] = 'NaN'
                Hr_A[i][j] = 'NaN'
                
    for k in range(maxSF+1):
        for i in range(len(H_A)):
            for j in range(len(Hr_A[0])):
                if (Rho[k][i][j] == 0.0):
                    Rho[k][i][j] = 'NaN'
    X['H_A'] = H_A
    X['Hr_A'] = Hr_A
    X['rho'] = Rho
    return(X)
    
  
def rotate_FORC(X):
    H_A = X['H_A']
    Hr_A = X['Hr_A']
    Hc= (H_A - Hr_A)/2. #x axis
    Hu = (H_A + Hr_A)/2. #y axis
    X['Hc'] = Hc
    X['Hu'] = Hu
    return(X)  

#FWHM function
def half_max_test(fwHu_c, fwRho_c, ym):
    arr_L = np.where(fwRho_c == ym)[0]
    L = arr_L[0]
    half_ym = ym/2. #half max
    b = L+1

    while (b < len(fwRho_c)):

        if(fwRho_c[b] < half_ym):
            
            break
        b = b + 1
    
    top = fwRho_c[b-1] - fwRho_c[b]
    bot = fwHu_c[b-1] - fwHu_c[b]
   
    mo_test = top/bot
    r0 = fwHu_c[b] + ((half_ym - fwRho_c[b])/mo_test)
    #print(r0, mo_test, )
    u = L-1

    while (u > 0): 
       
        if (fwRho_c[u] < half_ym):
            
            break
        u = u - 1
    
    #interpolation to get half maximum for each in hu 
    m1 = (fwRho_c[u] - fwRho_c[u+1])/(fwHu_c[u] - fwHu_c[u+1])

    r1 = fwHu_c[u+1] + ((half_ym - fwRho_c[u+1])/m1)
  
    fwhm = r1 - r0
   
    return fwhm, r0, r1


def find_fwhm(X, SF, sample_name): 
    fwhmlist = X['fwhmlist']
    #print('fwhmlist', fwhmlist)
    Rho = X['rho'] 
    #print('Rho', Rho)
    #print('Hc', X['Hc'])
    Hc = X['Hc']   
    Hu = X['Hu']

    #indices = np.unravel_index(np.nanargmax(Rho[SF]),Rho[SF].shape)
    indices = np.unravel_index(np.nanargmax(Rho[SF][:][:133]),Rho[SF][:][:133].shape) #fix for model as noise at edge of FORC diagram
    #if (indices[0] > 145):

    #print(SF, indices, np.nanargmax(Rho[SF]), Rho[SF].shape)
    fwHu = []
    fwRho = []
    fwHc = []
    for i in range(len(Rho[SF])):
        fwHu.append(Hu[i][indices[1]]) 
        fwHc.append(Hc[i][indices[1]]) 
        fwRho.append(Rho[SF][i][indices[1]])
        i+=1

    fwHu = np.array(fwHu)
    fwRho = np.array(fwRho)
    fwHc = np.array(fwHc)
    fwHu = fwHu[~np.isnan(fwHu)]
    fwRho = fwRho[~np.isnan(fwRho)] 
    fwHc = fwHc[~np.isnan(fwHc)]
    r0 = 1
    r1 = -1
    #print('fwHu', fwHu)
    #print('fwRho', fwRho)
    #print('fwHc', fwHc)

    
    #if (len(fwHu > 1)):
    #    loc_o = np.argmin(abs(fwHu))
    #else:
    #    loc_o = 0

    loc_o = np.argmin(abs(fwHu))
    #print('loc_o', loc_o)
    fwHu_f = fwHu[:loc_o] 
    fwRho_f = fwRho[:loc_o]
    #print('fwRho_f', fwRho_f)
    loc_m = np.argmin(abs(fwRho_f))
    fwHu_c = fwHu[loc_m:(loc_o +(loc_o - loc_m))]
    fwRho_c = fwRho[loc_m:(loc_o +(loc_o - loc_m))]
    X['fwRho_L'].append(fwRho_c)
    X['fwHu_L'].append(fwHu_c) 
    #plt.plot(fwHu_c, fwRho_c, label = SF)
    #plt.legend()
    #plt.show

    m_rho_a = np.sort(fwRho_c)
    i = 1
    while ((r0 >0) or (r1 < 0)): #opposte to FWHM crossing 0 
        ym = m_rho_a[-i]

        # find the two crossing points
        try:
  
            fwhm, r0, r1 = half_max_test(fwHu_c, fwRho_c, ym)

        except:
            print('Error in calculating FWHM for SF %',SF)
            pass
        
        if (i >5):
            print('SF {} is too noisy'.format(SF))
            fwhm = 'Nan'
            r0 = 'NaN'
            r1 = 'NaN'
            break
        i+=1
  
    fwhmlist.append(fwhm)

    half = max(fwRho_c)/2.0
    X['Hu_rs'].append([r0, r1])
    X['Rho_half'].append([half, half])
    #plt.plot([r0, r1], [half, half], label = SF)
    #plt.xlabel('$\mathrm{h_{s}}$ (mT)')
    #plt.ylabel('FORC weighting')
    #plt.legend()
    #plt.title('Plot of the cross sections of the FWHM at each smoothing factor')

    #plt.show
    X['fwhmlist'] = fwhmlist

    #plot this and the second part all in 1 figure 
    return(X)

def find_plot_fwhm(X):
    if (X['sample_copy'] > 1): #only if second go 
        sf_1 = input('Using the FWHM the previous SF used for sample {} was {}. If you want to re-use theis variable enter K, to change them enter any other character'.format(X['names'][0], X['SF']))
    else:
        sf_1 = 'L'
    if (sf_1 != 'K') or (X['sample_copy'] < 2):
        maxSF = 8
        SFlist = []
        fwhmlist = []
        name = X['name']
        X['fwhmlist'] = fwhmlist
        X['fwRho_L'] = []
        X['fwHu_L'] = []
        X['Hu_rs'] = []
        X['Rho_half'] = []


        for SF in range(2,maxSF+1):
            SFlist.append(SF)
        X['SF_list'] = SFlist    
        col_list = ['blue', 'red', 'green', 'orange', 'purple', 'cyan', 'lime', 'hotpink', 'peru', 'orange']
        #make figure with 2 subfigures next to each other 
        #before plot this ask if they want to choose 
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        for i in range(len(SFlist)):
            X = find_fwhm(X, SFlist[i], name)
            #for each SF do the plots - add to the plots
            hu_rsP = np.array(X['Hu_rs'][i][:])
            if (hu_rsP[0] != 'NaN'):
                hu_rsP2 = hu_rsP*1000
                ax1.plot(hu_rsP2, X['Rho_half'][i], color = col_list[i])
            ax1.plot(X['fwHu_L'][i][:]*1000, X['fwRho_L'][i], color = col_list[i], label = SFlist[i])
            
            i+=1 #0,1,2,3      
        #once all calcualted we can plot them
        ax1.set_xlabel('$\mathrm{h_{s}}$ (mT)')
        ax1.set_ylabel('FORC weighting')
        ax1.legend()
        ax1.set_title('Plot of the cross sections of the FWHM at each SF')
        #part 2 of the plot

        SFlist = X['SF_list']
        fwhmlist = X['fwhmlist']
        st_line_SFlist = []
        polyfwhm = []
        polySF = []

        maxSF1 = X['maxSF1']
        for i in range(maxSF1+1):
            st_line_SFlist.append(i)
            i +=1

        st_line_SFlist= np.array(st_line_SFlist)
        SFlist = np.array(SFlist)

        for i in range(len(fwhmlist)):
            if (fwhmlist[i] != 'Nan') and (fwhmlist[i] != 'NaN'): #add in fwhmlist[i] != nan
                polyfwhm.append(float(fwhmlist[i]))
                polySF.append(float(SFlist[i]))

        poly_fwhm_plot = np.array(polyfwhm)
        b, m = polyfit(polySF, (poly_fwhm_plot), 1)
        X['b'] = b
        ax2.scatter(polySF, poly_fwhm_plot*1000)
        bplot, mplot = polyfit(polySF, (poly_fwhm_plot*1000), 1)
        ax2.set_title('Smoothing factor (SF) vs FWHM using SF 2-5')
        ax2.set_xlabel('SF')
        ax2.set_ylabel('FWHM (mT)')
        ax2.plot(st_line_SFlist, bplot + mplot * st_line_SFlist, '-')
        plt.show 
        X['fwhmlist'] = fwhmlist
        plt.pause(1)
    print('1st attempt', X['b'], X['fwhmlist'])
    X['sf_1'] = sf_1
    return(X)

def check_fwhm(X):
    if (X['sf_1'] != 'K'):
        SFlist = X['SF_list']
        answer = None
        answer2 = None
        fwhmlist = X['fwhmlist']
        maxSF1 = X['maxSF1']
        print('1st', fwhmlist, maxSF1)
        while answer not in ("Y", "N"):
            answer = input("Are any of the FWHM unreliable? Enter Y or N: ")
            if (answer == "Y"):
                sf_int = 0
                while (sf_int == 0):
                    try:
                        sf_pick = (input("Which SF is unrealiable and needs to be removed?:" ))
                        sf_pick = int(sf_pick) #ask for interger and check it is an interger, if not ask again
                        sf_int = 1
                        print('pick', sf_pick)
                        if (sf_pick >= 2) and (sf_pick <= maxSF1):
                            sf_int = 1

                        else:
                            sf_int = 0
                            print('Not an interger between 2 and 5. Please input an interger between 2 and 5.')

                    except ValueError:
                        print('Not an interger. Please input an interger between 2 and 5.')


                while answer2 not in ("Y", "N"):
                    answer2 = input("Are any other FWHM unreliable? Enter yes or no: ")
            
                    if (answer2 == "Y"):
                        sf_int2 = 0
                        while (sf_int2 == 0):
                            try:
                                sf_pick2 = (input("Which other SF is unrealiable and needs to be removed?:" ))
                                sf_pick2 = int(sf_pick2) #ask for interger and check it is an interger, if not ask again
                                sf_int2 = 1

                                if (sf_pick2 >= 2) and (sf_pick2 <= maxSF1):
                                    sf_int2 = 1

                                else:
                                    sf_int2 = 0
                                    print('Not an interger between 2 and 5. Please input an interger between 2 and 5.')

                            except ValueError:
                                print('Not an interger. Please input an interger between 2 and 5.')
                            


                    elif (answer2 == "N"):
                        print(answer2)

                    
                    elif (isinstance(answer2, str)):
                        print("Please enter Y or N.")
                        
                    
                

                fwhmlist[sf_pick-2] = 'che' 
                            
                
                if (answer2 == "Y"):

                    fwhmlist[sf_pick2-2] = 'che' 
                X['fwhmlist'] = fwhmlist

                X = plot_fwhm(X) 
                

            elif answer == "N":
                print('1', X['fwhmlist'])
                X = plot_fwhm(X) 
            
            elif (isinstance(answer, str)):
                print("Please enter yes or no.")
    
        fwhmlist = np.array(fwhmlist)
        X['fwhmlist'] = fwhmlist

        print(X['fwhmlist'])
    return(X)

def plot_fwhm(X):
    SFlist = X['SF_list']
    fwhmlist = X['fwhmlist']
    st_line_SFlist = []
    polyfwhm = []
    polySF = []
    print('redoing the calcualtion')
    maxSF1 = X['maxSF1']
    for i in range(maxSF1+1):
        st_line_SFlist.append(i)
        i +=1

    st_line_SFlist= np.array(st_line_SFlist)

    SFlist = np.array(SFlist)

    print('looking for che', fwhmlist)
    for i in range(len(SFlist)):
        if (fwhmlist[i] != 'che') and (fwhmlist[i] != 'Nan'): #remove value 

            polyfwhm.append(float(fwhmlist[i]))
            polySF.append(float(SFlist[i]))
    
    poly_fwhm_plot = np.array(polyfwhm)
    print('removed the nan and che', poly_fwhm_plot)
    b, m = polyfit(polySF, (poly_fwhm_plot), 1)
    X['b'] = b
    plt.scatter(polySF, poly_fwhm_plot*1000)
    bplot, mplot = polyfit(polySF, (poly_fwhm_plot*1000), 1)

    plt.scatter(polySF, poly_fwhm_plot*1000)
    plt.plot(st_line_SFlist, bplot + mplot * st_line_SFlist, '-')
    plt.xlabel('SF')
    plt.ylabel('FWHM (mT)')
    plt.title('Smoothing factor (SF) vs FWHM using SF 2-5')
    plt.show 
    print('redoing calcualtion the results', X['b'])
    Hu = X['Hu']

    i=0
    for i in range(len(SFlist)):

        if (fwhmlist[i] == 'Nan') or (fwhmlist[i] == 'che'):

            fwhmlist[i] = float(m*SFlist[i] + b)


    X['fwhmlist'] = fwhmlist
    print('new fwhmlist', X['fwhmlist'])
    return(X)

"""
def finding_fwhm(X):
    maxSF = 5
    SFlist = []
    fwhmlist = []
    name = X['name']
    X['fwhmlist'] = fwhmlist
    for SF in range(2,maxSF+1):
        SFlist.append(SF)
    X['SF_list'] = SFlist    
    for i in range(len(SFlist)):
        X = find_fwhm(X, SFlist[i], name)
        i+=1 #0,1,2,3    
    return (X)
    
#FWHM function
def half_max_test(fwHu_c, fwRho_c, ym):
    arr_L = np.where(fwRho_c == ym)[0]
    L = arr_L[0]
    half_ym = ym/2. #half max
    b = L+1

    while (b < len(fwRho_c)):

        if(fwRho_c[b] < half_ym):
            
            break
        b = b + 1
    
    top = fwRho_c[b-1] - fwRho_c[b]
    bot = fwHu_c[b-1] - fwHu_c[b]
   
    mo_test = top/bot
    r0 = fwHu_c[b] + ((half_ym - fwRho_c[b])/mo_test)
   
    u = L-1

    while (u > 0): 
       
        if (fwRho_c[u] < half_ym):
            
            break
        u = u - 1
    

    m1 = (fwRho_c[u] - fwRho_c[u+1])/(fwHu_c[u] - fwHu_c[u+1])

    r1 = fwHu_c[u+1] + ((half_ym - fwRho_c[u+1])/m1)
  
    fwhm = r1 - r0
   
    return fwhm, r0, r1


def find_fwhm(X, SF, sample_name): 
    fwhmlist = X['fwhmlist']
    Rho = X['rho'] 

    Hu = X['Hu']

    indices = np.unravel_index(np.nanargmax(Rho[SF]),Rho[SF].shape)

    fwHu = []
    fwRho = []
    for i in range(len(Rho[SF])):
        fwHu.append(Hu[i][indices[1]]) 
        fwRho.append(Rho[SF][i][indices[1]])
        i+=1

    fwHu = np.array(fwHu)
    fwRho = np.array(fwRho)
    fwHu = fwHu[~np.isnan(fwHu)]
    fwRho = fwRho[~np.isnan(fwRho)] 
    r0 = 1
    r1 = -1

    
    loc_o = np.argmin(abs(fwHu))
    fwHu_f = fwHu[:loc_o] 
    fwRho_f = fwRho[:loc_o]

    loc_m = np.argmin(abs(fwRho_f))
    fwHu_c = fwHu[loc_m:(loc_o +(loc_o - loc_m))]
    fwRho_c = fwRho[loc_m:(loc_o +(loc_o - loc_m))]

    plt.plot(fwHu_c, fwRho_c)
    plt.show
    m_rho_a = np.sort(fwRho_c)
    i = 1
    while ((r0 >0) or (r1 < 0)): #opposte to FWHM crossing 0 
        ym = m_rho_a[-i]

        # find the two crossing points
        try:
  
            fwhm, r0, r1 = half_max_test(fwHu_c, fwRho_c, ym)

        except:
            print('Error in calculating FWHM for SF %',SF)
            pass
        
        if (i >5):
            print('SF {} is too noisy'.format(SF))
            fwhm = 'Nan'
            r0 = 'NaN'
            r1 = 'NaN'
            break
        i+=1
  
    fwhmlist.append(fwhm)

    half = max(fwRho_c)/2.0

    plt.plot([r0, r1], [half, half], label = SF)
    plt.xlabel('$\mathrm{h_{s}}$ (mT)')
    plt.ylabel('FORC weighting')
    plt.legend()
    plt.title('Plot of the cross sections of the FWHM at each smoothing factor')

    plt.show
    X['fwhmlist'] = fwhmlist
    return(X)
    

def plot_fwhm_1(X):
    SFlist = X['SF_list']
    fwhmlist = X['fwhmlist']
    st_line_SFlist = []
    polyfwhm = []
    polySF = []

    maxSF1 = X['maxSF1']
    for i in range(maxSF1+1):
        st_line_SFlist.append(i)
        i +=1

    st_line_SFlist= np.array(st_line_SFlist)
    SFlist = np.array(SFlist)

    for i in range(len(fwhmlist)):
       if (fwhmlist[i] != 'Nan') and (fwhmlist[i] != 'NaN'): #add in fwhmlist[i] != nan
            polyfwhm.append(float(fwhmlist[i]))
            polySF.append(float(SFlist[i]))


    plt.scatter(polySF, polyfwhm)

    b, m = polyfit(polySF, polyfwhm, 1)
    X['b'] = b
    plt.title('Smoothing factor (SF) vs FWHM using SF 2-5')
    plt.xlabel('SF')
    plt.ylabel('FWHM')
    plt.plot(st_line_SFlist, b + m * st_line_SFlist, '-')

    plt.show 

    Hu = X['Hu']

    i=0

    X['fwhmlist'] = fwhmlist

    return(X)
    

def plot_fwhm(X):
    SFlist = X['SF_list']
    fwhmlist = X['fwhmlist']
    st_line_SFlist = []
    polyfwhm = []
    polySF = []

    maxSF1 = X['maxSF1']
    for i in range(maxSF1+1):
        st_line_SFlist.append(i)
        i +=1

    st_line_SFlist= np.array(st_line_SFlist)

    SFlist = np.array(SFlist)


    for i in range(len(SFlist)):
        if (fwhmlist[i] != 'che') and (fwhmlist[i] != 'Nan'): #remove value 

            polyfwhm.append(float(fwhmlist[i]))
            polySF.append(float(SFlist[i]))

    plt.scatter(polySF, polyfwhm)

    b, m = polyfit(polySF, polyfwhm, 1)
    X['b'] = b
    plt.xlabel('SF')
    plt.ylabel('FWHM')
    plt.title('SF versus FWHM plot with the accepted SFs')
    plt.plot(st_line_SFlist, b + m * st_line_SFlist, '-')

    plt.show 

    Hu = X['Hu']

    i=0
    for i in range(len(SFlist)):

        if (fwhmlist[i] == 'Nan') or (fwhmlist[i] == 'che'):

            fwhmlist[i] = float(m*SFlist[i] + b)


    X['fwhmlist'] = fwhmlist

    return(X)
    
    
def check_fwhm(X):
    SFlist = X['SF_list']
    answer = None
    answer2 = None
    fwhmlist = X['fwhmlist']
    maxSF1 = X['maxSF1']
    while answer not in ("yes", "no"):
        answer = input("Are any of the FWHM unreliable? Enter yes or no: ")
        if (answer == "yes"):
            sf_int = 0
            while (sf_int == 0):
                try:
                    sf_pick = (input("Which SF is unrealiable and needs to be removed?:" ))
                    sf_pick = int(sf_pick) #ask for interger and check it is an interger, if not ask again
                    sf_int = 1

                    if (sf_pick >= 2) and (sf_pick <= maxSF1):
                        sf_int = 1

                    else:
                        sf_int = 0
                        print('Not an interger between 2 and 5. Please input an interger between 2 and 5.')

                except ValueError:
                    print('Not an interger. Please input an interger between 2 and 5.')


            while answer2 not in ("yes", "no"):
                answer2 = input("Are any other FWHM unreliable? Enter yes or no: ")
        
                if (answer2 == "yes"):
                    sf_int2 = 0
                    while (sf_int2 == 0):
                        try:
                            sf_pick2 = (input("Which other SF is unrealiable and needs to be removed?:" ))
                            sf_pick2 = int(sf_pick2) #ask for interger and check it is an interger, if not ask again
                            sf_int2 = 1
 
                            if (sf_pick2 >= 2) and (sf_pick2 <= maxSF1):
                                sf_int2 = 1

                            else:
                                sf_int2 = 0
                                print('Not an interger between 2 and 5. Please input an interger between 2 and 5.')

                        except ValueError:
                            print('Not an interger. Please input an interger between 2 and 5.')
                        


                elif (answer2 == "no"):
                    print(answer2)

                
                elif (isinstance(answer2, str)):
                    print("Please enter yes or no.")
                    
                
            

            fwhmlist[sf_pick-2] = 'che' 
                        
              
            if (answer2 == "yes"):

                fwhmlist[sf_pick2-2] = 'che' 
            X['fwhmlist'] = fwhmlist

            X = plot_fwhm(X) 
            

        elif answer == "no":
        
            X = plot_fwhm(X) 
           
        elif (isinstance(answer, str)):
            print("Please enter yes or no.")
   
    fwhmlist = np.array(fwhmlist)
    X['fwhmlist'] = fwhmlist

    
    return(X)
    """

def norm_rho_all(X):
    Rho = X['rho']
    Rho_n = np.copy(Rho)
    a, x, y = np.shape(Rho)
    X['max_Rho'] = np.zeros((a))


    k=0
    for k in range(2,a):
        i = 0
        j = 0
        max_Rho = np.nanmax(Rho_n[k])

        for i in range(x):
            for j in range(y):
                Rho_n[k][i][j] = Rho[k][i][j]/max_Rho
        X['rho_n'] = Rho_n
        X['max_Rho'][k] = max_Rho
        k+=1

    return(X)    

  
def plot_sample_FORC(x, y, z, SF, sample_name):
    z = z[SF]
    zn = np.copy(z)

    xp = x*1000
    yp = y*1000

    con = np.linspace(0.1, 1, 9)



    cmap, vmin, vmax = FORCinel_colormap(zn) #runs FORCinel colormap

    plt.contourf(xp, yp, zn, 50, cmap= cmap, vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(ticks=[0,0.2, 0.4, 0.6, 0.8, 1], format = '%.1f')
    plt.contour(xp, yp, zn, con, colors = 'k')
    plt.xlabel('$\mathrm{h_{c}}$ (mT)', fontsize=14)
    plt.ylabel('$\mathrm{h_{s}}$ (mT)', fontsize=14)

    plt.xlim(0, np.nanmax(xp))
    plt.ylim(np.nanmin(yp), np.nanmax(yp))
    plt.tick_params(axis='both', which='major', labelsize=12)

    cbar.ax.tick_params(labelsize=12)
    plt.title('{} FORC diagram, SF = {}'.format(sample_name, SF))

    plt.tight_layout()
    plt.savefig('FORC_diagram_SF_{}_{}.pdf'.format(SF, sample_name))
    #path = os.getcwd
    #plt.savefig(path + os.sep+'sample_{}_V_{}'.format(sample_name, sample_copy) + os.sep+'FORC_diagram_sample_{}_SF_{}.pdf'.format(sample_name,SF))
    plt.show
    plt.pause(1)
    return     
    
def plot_sample_FORC2(x, y, z, SF, sample_name, xm, ym2):
    path = os.getcwd() #current directory
    z = z[SF]
    zn = np.copy(z)
    
    xp = x*1000
    yp = y*1000

    con = np.linspace(0.1, 1, 9)
   
    cmap, vmin, vmax = FORCinel_colormap(zn) #runs FORCinel colormap
    
    plt.contourf(xp, yp, zn, 50, cmap= cmap, vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(ticks=[0,0.2, 0.4, 0.6, 0.8, 1], format = '%.1f')
    plt.contour(xp, yp, zn, con, colors = 'k')
    plt.xlabel('$\mathrm{h_{c}}$ (mT)', fontsize=14)
    plt.ylabel('$\mathrm{h_{s}}$ (mT)', fontsize=14)
    plt.xlim(0, xm)
    plt.ylim(-ym2, ym2)

    plt.tick_params(axis='both', which='major', labelsize=12)

    cbar.ax.tick_params(labelsize=12)
    plt.title('FORC diagram for sample {}, SF = {}'.format(sample_name, SF))
    plt.savefig(os.path.join(path+ os.sep,'processed_data','{}.pdf'.format(sample_name,SF)), bbox_inches='tight')
    plt.close()

    #zip up file 
    
    shutil.make_archive('processed_data', 'zip', 'processed_data')
    shutil.rmtree('processed_data') 

    return      

def divide_mu0(X):
    mu0 = 4.*pi*1e-7
    X['Hc_mu'] = X['Hc']/mu0

    X['Hu_mu'] = X['Hu']/mu0
    
    return(X)

def FORCinel_colormap(Z):

    #setup initial colormap assuming that negative range does not require extension
    cdict = {'red':     ((0.0,  127/255, 127/255),
                         (0.1387,  255/255, 255/255),
                         (0.1597,  255/255, 255/255),
                         (0.1807,  255/255, 255/255),
                         (0.3193,  102/255, 102/255),
                       (0.563,  204/255, 204/255),
                       (0.6975,  204/255, 204/255),
                       (0.8319,  153/255, 153/255),
                       (0.9748,  76/255, 76/255),
                       (1.0, 76/255, 76/255)),

            'green':   ((0.0,  127/255, 127/255),
                         (0.1387,  255/255, 255/255),
                         (0.1597,  255/255, 255/255),
                         (0.1807,  255/255, 255/255),
                       (0.3193,  178/255, 178/255),
                        (0.563,  204/255, 204/255),
                       (0.6975,  76/255, 76/255),
                       (0.8319,  102/255, 102/255),
                       (0.9748,  25/255, 25/255),
                       (1.0, 25/255, 25/255)),

             'blue':   ((0.0,  255/255, 255/255),
                         (0.1387,  255/255, 255/255),
                         (0.1597,  255/255, 255/255),
                         (0.1807,  255/255, 255/255),
                       (0.3193,  102/255, 102/255),
                        (0.563,  76/255, 76/255),
                       (0.6975,  76/255, 76/255),
                       (0.8319,  153/255, 153/255),
                       (0.9748,  76/255, 76/255),
                       (1.0, 76/255, 76/255))}
    
    if np.abs(np.min(Z))<=np.nanmax(Z):#*0.19: #negative extension is not required
        vmin = -np.nanmax(Z)#*0.19
        vmax = np.nanmax(Z)
 
    else: #negative extension is required

        vmin=np.nanmin(Z)
        vmax=np.nanmax(Z)
      
    anchors = np.zeros(10)    
    
    anchors[1]=(0.0005*vmax-vmin)/(vmax-vmin)
    anchors[2]=(0.005*vmax-vmin)/(vmax-vmin)   
    anchors[3]=(0.025*vmax-vmin)/(vmax-vmin)
    anchors[4]=(0.19*vmax-vmin)/(vmax-vmin)
    anchors[5]=(0.48*vmax-vmin)/(vmax-vmin)
    anchors[6]=(0.64*vmax-vmin)/(vmax-vmin)
    anchors[7]=(0.80*vmax-vmin)/(vmax-vmin)
    anchors[8]=(0.97*vmax-vmin)/(vmax-vmin)
    anchors[9]=1.0

    
    anchors = abs(anchors)
 

    Rlst = list(cdict['red'])
    Glst = list(cdict['green'])
    Blst = list(cdict['blue'])


    for i in range(9):
        Rlst[i] = tuple((anchors[i],Rlst[i][1],Rlst[i][2]))
        Glst[i] = tuple((anchors[i],Glst[i][1],Glst[i][2]))
        Blst[i] = tuple((anchors[i],Blst[i][1],Blst[i][2]))
        
    cdict['red'] = tuple(Rlst)
    cdict['green'] = tuple(Glst)
    cdict['blue'] = tuple(Blst)
   
    
    cmap = matplotlib.colors.LinearSegmentedColormap('forc_cmap', cdict)
    
    return cmap, vmin, vmax
    
def inter_FORC(X, SF):
    Hu_f = X['Hu_mu'].flatten() 
    Hc_f = X['Hc_mu'].flatten()
    Rho_f = X['rho'][SF].flatten() 
    
    #remove nan
    Hu_f = Hu_f[~np.isnan(Hu_f)]
    Hc_f = Hc_f[~np.isnan(Hc_f)]
    Rho_f = Rho_f[~np.isnan(Rho_f)]
    
    step_xi = np.nanmax(X['Hc_mu'])/181.
    step_yi = (np.nanmax(X['Hu_mu']) - np.nanmin(X['Hu_mu']))/146.

    # target grid to interpolate to
    xi = np.arange(0,np.nanmax(X['Hc_mu']),step_xi) 
    yi = np.arange(np.nanmin(X['Hu_mu']),np.nanmax(X['Hu_mu']),step_yi) 
    xi1,yi1 = np.meshgrid(xi,yi) 



    zi = griddata((Hc_f,Hu_f),Rho_f,(xi1,yi1),method='cubic') 

    X['xi1'] = xi1
    X['yi1'] = yi1
    X['zi_{}'.format(SF)] = zi
    return (X)
    
def inter_rho(xi_s_f, yi_s_f, zi_s_f, hys, i): # when call use xi_s etc, i is no. hysteron to test
    xi1_row = xi_s_f[0,:] 

    up_hc = xi1_row[xi1_row > hys[i,0]].min()   
    lo_hc = xi1_row[xi1_row < hys[i,0]].max() 

    up_hc_idx = list(xi1_row).index(up_hc) 
    lo_hc_idx = list(xi1_row).index(lo_hc)
    
    yi1_col = yi_s_f[:,0] 

    up_hi = yi1_col[yi1_col > hys[i,1]].min()   
    lo_hi = yi1_col[yi1_col < hys[i,1]].max() 

    up_hi_idx = list(yi1_col).index(up_hi) 
    lo_hi_idx = list(yi1_col).index(lo_hi)

    x_arr = np.array([xi_s_f[lo_hi_idx,lo_hc_idx], xi_s_f[up_hi_idx, lo_hc_idx], xi_s_f[up_hi_idx, up_hc_idx], xi_s_f[lo_hi_idx, up_hc_idx]])
    y_arr = np.array([yi_s_f[lo_hi_idx,lo_hc_idx], yi_s_f[up_hi_idx, lo_hc_idx], yi_s_f[up_hi_idx, up_hc_idx], yi_s_f[lo_hi_idx, up_hc_idx]])
    z_arr = np.array([zi_s_f[lo_hi_idx,lo_hc_idx], zi_s_f[up_hi_idx, lo_hc_idx], zi_s_f[up_hi_idx, up_hc_idx], zi_s_f[lo_hi_idx, up_hc_idx]])

   
    xarr_sum = np.sum(x_arr)
    xarr_has_nan = np.isnan(xarr_sum)
    yarr_sum = np.sum(y_arr)
    yarr_has_nan = np.isnan(yarr_sum)
    zarr_sum = np.sum(z_arr)
    zarr_has_nan = np.isnan(zarr_sum)    

    
    if (xarr_has_nan != True) and (yarr_has_nan != True) and (zarr_has_nan != True):
        f = interp2d(x_arr, y_arr, z_arr, kind='linear')
        hys[i,3] = f(hys[i,0], hys[i,1]) 
    else:
        hys[i,3] = -0.001
    


    return hys
    
    
def sym_FORC(X, SF):
    xi1 = X['xi1']
    yi1 = X['yi1']
    zi = X['zi_{}'.format(SF)] 
    yi_axis = np.copy(yi1)

    yi_axis = abs(yi_axis) - 0

    
    indices = np.unravel_index(np.nanargmin(yi_axis),yi_axis.shape)
   

    xi_s = np.copy(xi1)
    yi_s = np.copy(yi1)
    zi_s = np.copy(zi)
  
    x=1
    j=0

    while x < (len(xi1) - indices[0]): 

        j=0
        for j in range(len(xi_s[0])):


            find_mean = np.array([zi_s[indices[0]+x][j],zi_s[indices[0]-x][j]])
          
            find_mean = find_mean[~np.isnan(find_mean)]
            if (len(find_mean) > 0):
                zi_s[indices[0]+x][j] = np.mean(find_mean)
             
                zi_s[indices[0]-x][j] = zi_s[indices[0]+x][j]


            j+=1
        x+=1        
    
    lower_x_bound = indices[0] - (len(xi_s) - indices[0])
    upper_x_bound = len(xi_s)
    

    xi_s_cut = xi_s[lower_x_bound:upper_x_bound,:]
    yi_s_cut = yi_s[lower_x_bound:upper_x_bound,:]
    zi_s_cut = zi_s[lower_x_bound:upper_x_bound,:]


    X['xis'] = xi_s_cut
    X['yis'] = yi_s_cut
    X['zis_{}'.format(SF)] = zi_s_cut
    
    
    return(X)
    

def sym_norm_forcs(X):
    if (X['sf_1'] != 'K'):
        SFlist = X['SF_list']
        for i in range(len(SFlist)): #loop over each SF 2,3,4,5
            X = inter_FORC(X, SFlist[i])    
            X = sym_FORC(X, SFlist[i])
            print('check picking SF doe sthis correct', SFlist[i])
        for i in range(len(SFlist)):
            X = norm_z2(X, SFlist[i])
        fwhmlist = X['fwhmlist'] #correct FHWM values 
        X['sf_list_correct'] = (X['b']/fwhmlist) 
    else:
        X = inter_FORC(X, X['SF'])    
        X = sym_FORC(X, X['SF'])    
        X = norm_z2(X, X['SF'])    

        

    return(X)
        

def norm_z2(X, SF):

    z_pre_norm = X['zis_{}'.format(SF)]

    maxval_z = np.nanmax(z_pre_norm)

    minval_z = np.nanmin(z_pre_norm)
    
 
    z_norm = (z_pre_norm)/(maxval_z)

    X['zis_{}'.format(SF)] = z_norm
    return(X)    

def user_input(X):

    #FORC limits
    #if done previously
    if (X['sample_copy'] > 1): #only if second go 
        lim_1 = input('The most recent maximum hc and maximum hi used for sample {} were {} mT and {} mT respectively. If you want to re-use these variables enter K, to change them enter any other charactor'.format(X['names'][0], X['reset_limit_hc'], X['reset_limit_hi']))
    else:
        lim_1 = 'L'
    if (lim_1 != 'K') or (X['sample_copy'] < 2):    

        quest_res = 0.
        while (quest_res == 0):
            quest = input('Do you want to change the maximum hc and/or hi value from max hc = {} and max hi = {} from thr FORC file? (Y or N)'.format(X['Hc2']*1000, X['Hb2']*1000))
            if (quest == 'Y'):
                quest_res = 1.
                hc_l = 0
                while (hc_l == 0):                
                    try:
                        max_input_hc = (input("Set maximum hc (mT) using the above FORC diagram, if unchanged enter 0:" ))
                        max_input_hc = float(max_input_hc) #ask for interger and check it is an interger, if not ask again
                        if (max_input_hc == 0):
                            max_input_hc = X['Hc2']*1000
                        hc_l = 1
            
            
                    except ValueError:
                        print('Not a number. Please input an number.') #excpet for question 
                        hc_l = 0

                hi_l = 0
                while (hi_l == 0):
                    try:
                        max_input_hi = (input("Set absolute maximum hi (mT) using the above FORC diagram, if unchanged enter 0:" ))
                        max_input_hi = float(max_input_hi) #ask for interger and check it is an interger, if not ask again
                        if (max_input_hi == 0):
                            max_input_hi = X['Hb2']*1000                
                        hi_l = 1
                    except ValueError:
                        print('Not a number. Please input an number.')
                        hi_l = 0
            if (quest == 'N'):
                quest_res = 1.
                max_input_hc = X['Hc2']*1000
                max_input_hi = X['Hb2']*1000
            else:
                print('Enter Y or N')


    X['reset_limit_hc'] = float(max_input_hc)
    X['reset_limit_hi'] = float(max_input_hi)
    print('FORC limits: Max hc = {}, max hi = {}'.format(X['reset_limit_hc'], X['reset_limit_hi']))  
    return X
  

#demag data
def demag_data(X):
    af_step = []
    af_nrm_n = []
    af_sirm_n = []
    af_step_irm = []

    af_nrm_data = open(X['nrm_files'][0], "r") #write file
    for myline1 in af_nrm_data:
        line = myline1.split(' ') #split into elements by spaces

        af_step.append(float(line[1])) 
       # print('af_step', af_step)
        af_nrm_n.append(float(line[3])) 
        #print('af_nrm_n', af_nrm_n)

    af_nrm_data.close()

    af_arm_data = open(X['sirm_files'][0], "r") #write file
    for myline in af_arm_data:
        line = myline.split(' ') #split into elements by spaces
    
        af_arm_n.append(float(line[3])) 
        #print(af_nrm)

    af_arm_data.close()

    #no af steps
    cntfield = len(af_step)

    af_step = np.array(af_step)
    af_nrm_n = np.array(af_nrm_n)
    af_arm_n = np.array(af_arm_n)

    af_arm_n = af_arm_n[:len(af_nrm_n)]

    afnorm = af_nrm_n[0]/af_arm_n[0] #only used for write out if want to use need to make in X
    #af_step = af_step/10.
    X['af_step'] = af_step
    X['af_nrm'] = af_nrm_n
    X['af_arm'] = af_arm_n
    X['cntfield'] = cntfield
    return(X)
    
def demag_data_mur(X):
    af_step = []
    af_nrm_n = []
    af_sirm_n = []
    af_step_irm = []

    af_nrm_data = open(X['nrm_files'][0], "r") #write file
    for myline1 in af_nrm_data:
        line = myline1.split('\t') #split into elements by tabs

        af_step.append(float(line[2])) 
       # print('af_step', af_step)
        af_nrm_n.append(float(line[1])) 
        #print('af_nrm_n', af_nrm_n)

    af_nrm_data.close()

    af_irm_data = open(X['sirm_files'][0], "r") #write file
    for myline in af_irm_data:
        line = myline.split('\t') #split into elements by tab
    
        af_sirm_n.append(float(line[1])) 
        #print(af_nrm)

    af_irm_data.close()

    #no af steps
    cntfield = len(af_step)

    af_step = np.array(af_step)
    af_nrm_n = np.array(af_nrm_n)
    af_sirm_n = np.array(af_sirm_n)

    af_sirm_n = af_sirm_n[:len(af_nrm_n)]

    afnorm = af_nrm_n[0]/af_sirm_n[0] #only used for write out if want to use need to make in X
    #af_step = af_step/10.
    X['af_step'] = af_step
    X['af_step'] = X['af_step']/10.
    X['af_nrm'] = af_nrm_n
    X['af_arm'] = af_sirm_n
    X['cntfield'] = cntfield
    return(X)

def demag_data_generic(X):
    af_step = []

    af_nrm_n = []
    af_sirm_n = []
    af_step_sirm = []
    af_nrm_dec = []
    af_nrm_inc = []
    af_sirm_dec = []
    af_sirm_inc = []


    af_nrm_data = open(X['nrm_files'][0], "r") #write file
    af_nrm_data.readline()
    for myline1 in af_nrm_data:
        line = myline1.split() #split into elements by tabs

        af_step.append(float(line[1])) 
        af_nrm_n.append(float(line[2]))
        af_nrm_dec.append(float(line[3]))
        af_nrm_inc.append(float(line[4]))


    af_nrm_data.close()

    af_irm_data = open(X['sirm_files'][0], "r") #write file
    af_irm_data.readline()
    for myline in af_irm_data:
        line = myline.split('\t') #split into elements by tab

        af_step_sirm.append(float(line[1])) 
        af_sirm_n.append(float(line[2]))
        af_sirm_dec.append(float(line[3]))
        af_sirm_inc.append(float(line[4]))

    af_irm_data.close()


    if (af_step == af_step_sirm):
        pass
    else:
        print('AF demagnetisation steps are different for NRM and SIRM. Re-import data.')

    #no af steps
    cntfield = len(af_step)

    af_step = np.array(af_step)
    af_nrm_n = np.array(af_nrm_n)
    af_sirm_n = np.array(af_sirm_n)
    af_nrm_dec = np.array(af_nrm_dec)
    af_sirm_dec = np.array(af_sirm_dec)
    af_nrm_inc = np.array(af_nrm_inc)
    af_sirm_inc = np.array(af_sirm_inc)

    af_sirm_n = af_sirm_n[:len(af_nrm_n)]

    afnorm = af_nrm_n[0]/af_sirm_n[0] #only used for write out if want to use need to make in X

    X['af_step'] = af_step
    X['af_nrm'] = af_nrm_n
    X['af_sirm'] = af_sirm_n
    X['af_nrm_dec'] = af_nrm_dec
    X['af_irm_dec'] = af_sirm_dec
    X['af_nrm_inc'] = af_nrm_inc
    X['af_irm_inc'] = af_sirm_inc
    X['cntfield'] = cntfield
    return(X)

def demag_data_generic2(X):
    af_step = []

    af_nrm_n = []
    af_sirm_n = []
    af_step_sirm = []
    af_nrm_dec = []
    af_nrm_inc = []
    af_sirm_dec = []
    af_sirm_inc = []


    af_nrm_data = open(X['nrm_files'][0], "r") #read file
    for myline1 in af_nrm_data:
        line = myline1.split() #split into elements by tabs

        af_step.append(float(line[1])) 
        af_nrm_n.append(float(line[3]))
        af_nrm_dec.append(float(line[4]))
        af_nrm_inc.append(float(line[5]))


    af_nrm_data.close()

    af_irm_data = open(X['sirm_files'][0], "r") #write file
    for myline in af_irm_data:
        line = myline.split() #split into elements by tab

        af_step_sirm.append(float(line[1])) 
        af_sirm_n.append(float(line[3]))
        af_sirm_dec.append(float(line[4]))
        af_sirm_inc.append(float(line[5]))

    af_irm_data.close()


    if (af_step == af_step_sirm):
        pass
    else:
        print('AF demagnetisation steps are different for NRM and SIRM. Re-import data.')

    #no af steps
    cntfield = len(af_step)

    af_step = np.array(af_step)
    af_nrm_n = np.array(af_nrm_n)
    af_sirm_n = np.array(af_sirm_n)
    af_nrm_dec = np.array(af_nrm_dec)
    af_sirm_dec = np.array(af_sirm_dec)
    af_nrm_inc = np.array(af_nrm_inc)
    af_sirm_inc = np.array(af_sirm_inc)
    #af_sirm_n = af_sirm_n*18.5
    af_sirm_n = af_sirm_n[:len(af_nrm_n)]

    afnorm = af_nrm_n[0]/af_sirm_n[0] #only used for write out if want to use need to make in X

    X['af_step'] = af_step
    X['af_nrm'] = af_nrm_n
    X['af_sirm'] = af_sirm_n
    X['af_nrm_dec'] = af_nrm_dec
    X['af_irm_dec'] = af_sirm_dec
    X['af_nrm_inc'] = af_nrm_inc
    X['af_irm_inc'] = af_sirm_inc
    X['cntfield'] = cntfield
    return(X)
   

def hys_angles():
    
    angle = random.random()

    phi = acos(2*angle - 1)
    if(phi > (pi/2.)): 
        phi = pi-phi
    angle2 = random.random()
    phistatic = acos(2*angle2 - 1)
    if(phistatic > (pi/2.)):
        phistatic = pi-phistatic
    
    angle3 = random.random()
    thetastatic = 2*pi*angle3
    return phi, phistatic, thetastatic
    
def calc_hk_arrays(hys, num, V): 
    tm = V['tm']
    ms = V['ms']
    hc = np.copy(hys[:,0]) 
    tempt = 300.
    #print('inside calc_hk_arrays testing')
    hf = ((hc/(sqrt(2)))**(0.54))*(10**(-0.52)) 
    #hf = ((hc)**(0.54))*(10**(-0.52)) 
    
    phi = np.copy(hys[:,5])
    
    phitemp = (((np.sin(phi))**(2./3.))+((np.cos(phi))**(2./3.)))**(-3./2.) 

    gphitemp = (0.86 + (1.14*phitemp)) 
  
    hatmp = hc/(sqrt(2)) 

    ht = hf*(log(tm/tau)) 
    
    hktmp = hatmp +ht + (2*hatmp*ht+ht**2)**0.5 
    hktmpmax = hktmp
    
    
    hktmpstore = hktmp
    i=0
    for i in range(int(num)):
        factor = 0.5
        searchme = 1000000.0     
        hstep = (hktmp[i]-hatmp[i])/5.

        while (abs(searchme)> 5):
            searchme = hktmp[i] - hatmp[i] - hktmp[i]*((2*ht[i]*phitemp[i]/hktmp[i])**(1/gphitemp[i])) #
            hktmpstore[i] = hktmp[i]

            if (searchme > 0):
                hktmp[i] = hktmp[i] - hstep
            else:
                hktmp[i] = hktmp[i] + hstep
                hstep = hstep*factor 

    hkphi = hktmpstore 
    hk = hkphi[:int(num)]/phitemp[:int(num)]

    hys[:int(num),9] = hk 
    kb=1.3806503e-23 
    
    hys[:int(num),11] = (kb*tempt)/(mu0*ms) #v_coeff
    v_act = (kb*tempt)/(hf*mu0*ms)
    vol = v_act[:int(num)]/(1 - (((hys[:int(num),0])/sqrt(2))/hk))
    hys[:int(num),10] = vol
    #for k in range(5):
    #    print(k, 'hk[k]', hk[k], 'v_act',v_act[k], 'v_coeff', hys[k,11], 'vol', vol[k])
    #print('end of testing inside hk array')
    """
    plt.scatter(hys[:2000,0],vol[:2000])
    plt.xlabel('hc (hys[0]')
    plt.ylabel('volume')
    plt.title('hc vs volume for 2000 hysterons')
    plt.savefig('hc_vol.png')
    plt.close()

    plt.scatter(hys[:2000,9],vol[:2000])
    plt.xlabel('hk (hys[9]')
    plt.ylabel('volume')
    plt.title('hk vs volume for 2000 hysterons')
    plt.savefig('hk_vol.png')
    plt.close()
    """
    return(hys)


def blocked_hys(V):
    hys_pre = np.copy(V['hys'])
    hys_new = np.zeros_like(hys_pre)
    hys_block = np.zeros_like(hys_pre)
	
    ms = V['ms']
    t_vol_a = hys_pre[:,10]*hys_pre[:,3]
    t_vol = sum(t_vol_a)
    #print('vol of all hys', t_vol, 'tot hys', len(hys_pre))
    i=0
    j=0
    k=0
    beta = (1-(300.-273)/(float(V['curie_t'])))**0.43
    temp = 300
    #how many end up in remanance region  
    for i in range(len(hys_new)):
        phitemp=((sin(hys_pre[i,5])**(2./3.))+(cos(hys_pre[i,5])**(2./3.)))**(-3./2.) 
        hc=(sqrt(2))*(hys_pre[i,9]*beta)
        hi = hys_pre[i,1]*beta

        g=0.86+1.14*phitemp

        v_co = (1 - ((abs(hi/sqrt(2)))/hys_pre[i,9]*beta))
        v_act_c = hys_pre[i,10]*v_co

            #hf = hys[i,11]*(1/(v_act_c))   
       
        hf = (kb*temp)/(ms*mu0*v_act_c)

        tm = 60.0
        ht = (roottwohffield)*hf*(log(tm/tau))
        field = 0.0
        bracket = 1-(2*ht*phitemp/hc)**(1/g)
        hiflipplus = hc*bracket
        hiflipminus=-hc*bracket
		
        if (hc >= (2*ht*phitemp)):
            if ((hi > hiflipminus) and (hi < hiflipplus)):
                hys_new[j,:] = hys_pre[i,:]
                j+=1
    i+=1


    hys_new = hys_new[hys_new[:,0] > 0.00000000001]
    
    vol_rem_a = hys_new[:,10]*hys_new[:,3]
    vol_rem = sum(vol_rem_a)
    #vol_rem = np.sum(hys_new[:,3])
    V['vol_rem'] = vol_rem
    #print('all sum= ', t_vol_a, ' blocked sum= ', vol_rem)
    #I['hys'] = hys_pre
    return(V)
    
    
def pop_hys(num_hys, X, V, SF): 
    #populate hysterons
    ms = V['ms']
    corr = X['sf_list_correct'][SF - 2]
    if (corr > 1.):
        corr = 1.
    hys = np.zeros((num_hys,12)) 
    
    print('corr', corr)
    num_pop = num_hys/2
    
    xi_s_cut = X['xis']
    yi_s_cut = X['yis']
    zi_s_cut = X['zis_{}'.format(SF)] 
    maxHc = np.nanmax(xi_s_cut) 
    maxHi = np.nanmax(yi_s_cut)

    xlim_res = X['reset_limit_hc']/1000.
    ylim_res = X['reset_limit_hi']/1000.
    maxHc = min(xlim_res,X['Hc2'])
    maxHi = min(ylim_res,X['Hb2'])
    #print(xlim_res,X['Hc2'], X['reset_limit_hc'])
    #print(ylim_res,X['Hb2'],X['reset_limit_hi'])
    #print(np.nanmax(xi_s_cut),np.nanmax(yi_s_cut))
    #print('B4 /mu0', 'maxHc', maxHc, 'maxHi', maxHi)
   
    i=0

    while (i < int(num_pop)):
        z1 = random.random()
        z2 = random.random()
        z3 = random.random()
  
        hys[i,0] = (z2*maxHc)/mu0

        hys[i,1] = (z3*maxHi)/mu0
        
        
        hys = inter_rho(xi_s_cut, yi_s_cut, zi_s_cut, hys, i) 
        hys[i,1] = hys[i,1]*corr
        hys[i,5], hys[i,6], hys[i,7] = hys_angles() #calc for half hys
        
        if (hys[i,3] >= 0) and (hys[i,3] <= 1): # and (hys[i,3] >= z1): #pick if FORC value greater than random  ((hys[i,1]) <= (hys[i,0])) and 
            i +=1 

    hys = calc_hk_arrays(hys, int(num_pop), V) # calc hk using arrrys 
    #print('populated')
    hys[:,4] = 1
 
    hys[:,8] = hys[:,5]*hys[:,4] 
    num_pop = int(num_pop)
    j=0
    for j in range(num_pop):
        hys[(j+num_pop),:] = hys[j,:]
        hys[j+num_pop,1] = -hys[j,1]
        j+=1
    #print('check 5 example hysterons')
    #for print_test in range(5):
    #    print(print_test, hys[print_test,:])
    #for print_test in range(num_pop, (num_pop+5), 1):
    #    print(print_test, hys[print_test,:])
    #print('end of test of hys population inside pop_hys')
    return hys, num_pop


#TRM acquisition on cooling 





"""
@jit(nopython = True)
def block_loc(var_1, hys, blockg, boltz):

    #unpack var_1
    num_hyss = int(var_1[0])
    beta = float(var_1[1])
    blockper = float(var_1[2])
    temp = float(var_1[3])
    aconst = float(var_1[4])
    curie_t = float(var_1[5])
    rate = float(var_1[6])    
    tempmin = float(var_1[7])
    field = float(var_1[8])
    tm = float(var_1[9])

    tau = 10e-9
    roottwohffield = 2**(0.5)
    hfstore = np.zeros((num_hyss))
    histore = np.zeros((num_hyss))    
    hcstore = np.zeros((num_hyss))   
    blocktemp = np.zeros((num_hyss))
    i=0
 
    for i in range(num_hyss): 
        phitemp=((sin(hys[i,5])**(2./3.))+(cos(hys[i,5])**(2./3.)))**(-3./2.) 
        hc=(sqrt(2))*(hys[i,9])*beta
       
        hcstore[i] = hc/(sqrt(2)) 

        hi = hys[i,1]*beta*blockper

        histore[i] = hi/(sqrt(2)) 
 


        g=0.86+1.14*phitemp

        hf=((hys[i,9]**0.54))*(10**(-0.52))*temp/(300*beta) 
   
        hfstore[i] = hf 

        if (rate == 1): 

            r = (1./aconst)*(temp-tempmin)
  
            tm = (temp/r)*(1 - (temp/(curie_t+273)))/log((2*temp)/(tau*r)*((1. - (temp/(curie_t+273)))))    
            #tm=(field*temp*aconst)/((temp-tempmin)*hf)   

            
            if (tm == 0.0): 
                tm = 60.0

        ht = (roottwohffield)*hf*(log(tm/tau))  

        bracket = 1-(2*ht*phitemp/hc)**(1/g)
        
        hiflipplus = hc*bracket+field*(roottwohffield) 

        hiflipminus=-hc*bracket+field*(roottwohffield) 
   
        if (hc >= (2*ht*phitemp)): 

            if ((hi > hiflipminus) and (hi < hiflipplus)):
              

                if ((blockg[i] == 0) or (blockg[i] == 2) or (blockg[i] == -2)): 
                    
                    if (hi >= (field*roottwohffield)): 
                        
                        blocktemp[i] = -1
                       
                    else:

                        blocktemp[i] = 1 
                      
                elif (blockg[i] == -1): 
                 
                    blocktemp[i] = -1
            

                elif (blockg[i] == 1):
  
                    blocktemp[i] = 1
                   
                else:
                    #write to screen, blockg[i]
                    print(blockg[i], blocktemp[i]) 
                    print('----', i)
                
            elif (hi >= hiflipplus):
               
                blocktemp[i] = -2
              
            else:
               
                blocktemp[i] = 2
        else: 

            if ((hi < hiflipminus) and (hi > hiflipplus)): 
               
                blocktemp[i] = 0
                boltz[i] = 0.0

            else: 
                if (hi >= hiflipminus):

                    blocktemp[i] = -2
                else:

                    blocktemp[i] = 2

    return hfstore, histore, boltz, blocktemp  


@jit(nopython = True)
def block_val(hys, histore, hfstore, blocktemp, beta, num_hyss, boltz, blockg, field, afswitch, afstore):
   
    i=0
    totalm = 0.0
    blockper = 0   
    for i in range(num_hyss):
        x = blocktemp[i]
        blockg[i] = x
        absblockg = abs(blockg[i])
        
        if (absblockg == 1): 
            if (boltz[i] < 0.00000001) and (boltz[i] > -0.000000001): #only zero

                boltz[i] = tanh((field - histore[i])/hfstore[i])
                if (afswitch == 1):

                    boltz[i] = 0.
        if (blockg[i] == -2):
         
            moment = -0 
        elif (blockg[i] == 2):

            moment = 0
        else:

            moment = blockg[i] 

        if (afswitch == 1): #afswtich given in function, 1 or 0

            hi = histore[i]*(sqrt(2))
            hc = (sqrt(2))*(hys[i,9]*beta)*(((sin(hys[i,6])**(2./3.))+(cos(hys[i,6])**(2./3.)))**(-3./2.))
    
            af = afstore/(1000*mu0) #define later. afstore number
       
            if (hi >=0) and (hi > (hc-af)): #whats this bit mean
                moment = 0
               
                boltz[i] = 0
                blockg[i] = -2
            
        
            if (hi <= 0) and (hi < (-hc + af)):
                moment = 0
                blockg[i] = 2
                boltz[i] = 0

       

        totalm = totalm + abs(moment)*abs(cos(hys[i,5]))*hys[i,3]*beta*(boltz[i]) 


        blockper=blockper+abs(moment)*1.0*hys[i,3] #remove hys[3] 
        i+=1

        
    blockper=blockper/(1.0*np.sum(hys[:,3]))

    return blockper, totalm, boltz, blockg
    

def blockfind(temp, field, afswitch, V, X): 

    hys = V['hys']
    num_hyss = V['num_hyss']
    hcstore = V['hcstore']
    histore = V['histore']
    beta = V['beta']
    rate = V['rate']
    aconst = V['aconst']
    tempmin = V['tempmin']

    
    tm = V['tm'] 

    hfstore = np.zeros(num_hyss)

    blockper = V['blockper'] 

    blocktemp = V['blocktemp'] 
    boltz = V['boltz']
    blockg = V['blockg']

    curie_t = X['curie_t']

    totalm = V['totalm']

    afstore = V['afstore']


    var_1 = np.array((num_hyss, beta, blockper, temp, aconst, curie_t, rate, tempmin, field, tm))

    hfstore, histore, boltz, blocktemp = block_loc(var_1, hys, blockg, boltz)  
    
    blockper, totalm, boltz, blockg = block_val(hys, histore, hfstore, blocktemp, beta, num_hyss, boltz, blockg, field, afswitch, afstore)
    
    V['blockper'] = blockper
    V['blocktemp'] = blocktemp
    V['boltz'] = boltz
    V['blockg'] = blockg
    V['totalm'] = totalm
    V['tm'] = tm

    return(V) 
"""

def blockfind_SC(t, field, afswitch, V, X): #try removing counttime as input/output

    hys = V['hys']
    prop = V['prop']
   
    #max_total_m_a = hys[:,10]*hys[:,3]
    #max_total = sum(max_total_m_a)
	
    vol_rem = V['vol_rem']
	
    num_hyss = V['num_hyss']

    histore = V['histore']
    beta = V['beta']
    rate = V['rate']

    tm = V['tm'] #try setting this for use in AF demag
    tstep = V['tstep']
    
    temp = V['temp']
    #temp = tempmin + tempstep
   
    hfstore = np.zeros(num_hyss)


    blockper = V['blockper'] 

    blocktemp = V['blocktemp'] 
    boltz = V['boltz']
    blockg = V['blockg']

    #end_mag = V['end_mag']

    totalm = V['totalm']
    afstore = V['afstore']


    max_total = np.sum(hys[:,3]) # sum of FORC distributions 
    var_1 = np.array((num_hyss, beta, blockper, temp, t, tstep, rate, field, tm))

    hfstore, histore, boltz, blocktemp = block_loc_C(var_1, hys, prop, blockg, boltz)
  
    blockper, totalm, boltz, blockg = block_val_C(hys, prop, histore, hfstore, blocktemp, beta, num_hyss, boltz, blockg, field, afswitch, afstore, max_total, vol_rem, t, tstep)
 
    V['aftotalm'] = totalm
    V['sir'] = totalm

    V['blockper'] = blockper
    V['blocktemp'] = blocktemp
    V['boltz'] = boltz
    V['blockg'] = blockg
    V['totalm'] = totalm

    V['tm'] = tm

    return #totalm

#CRM acquisition functions
#CRM acquisition and AF demag 
@jit(nopython = True)
def block_loc_C(var_1, hys, prop, blockg, boltz):
     #unpack variables from input array
    #unpack var_1
    num_hyss = int(var_1[0])
    beta = float(var_1[1])
    blockper = float(var_1[2])
    temp = float(var_1[3])
    t = float(var_1[4])
    tstep = float(var_1[5])
    rate = float(var_1[6])    
    field = float(var_1[7])
    tm = float(var_1[8])


    tau = 10e-9
 
    roottwohffield = 2**(0.5)
    hfstore = np.zeros((num_hyss))
    histore = np.zeros((num_hyss))    
    hcstore = np.zeros((num_hyss))   
    blocktemp = np.zeros((num_hyss))
    i=0 
    for i in range(num_hyss): #dynamic contributions - to hc dormann 1988 used?
        if (t >= (prop[i,3] + (0.5*tstep))):

            phitemp=((sin(hys[i,5])**(2./3.))+(cos(hys[i,5])**(2./3.)))**(-3./2.) #phitemp from hys[5] -> phi               
            hc=(sqrt(2))*(hys[i,9]*beta) #test and time HCbysqrt2 to get hc and matfhes similar order mag to hys[i,0/mu0 -> hc in same units using in hk - seems correct]
            
            hcstore[i] = hc/(sqrt(2.)) #different to hys[0]
        
            
            hi = hys[i,1]*beta*blockper #divide here by mu0
            
    
            histore[i] = hi/(sqrt(2.)) #should all be zero, without blockper, all look milar to hys[1,0] in mg

            g_ang=0.86+1.14*phitemp
            #print(hys[i,9], beta)
            h_hs = field - histore[i]
            v_co = (1 - ((abs(h_hs))/(hys[i,9]*beta)))
            if (v_co <= 0):
                v_co = 0.00001               
            
            v_act_c = prop[i,5]*v_co

            hf = hys[i,11]*(1/(v_act_c))
        
            hfstore[i] = hf #values at expected as hc/mu0

            if (rate == 1): #rate first set to 1 so shold do this
                
                tm= prop[i,10]
                
                if (tm == 0.0): #unsure
                    tm = 60.0
    
            #print('tm', tm)
            
            #if (end_mag == 1.0):
            #    tm = 60.0
            ht = (roottwohffield)*hf*(log(tm/tau)) #new tm 
        
            
            bracket = 1-(2*ht*phitemp/hc)**(1/g_ang)
            
            hiflipplus = hc*bracket+field*(roottwohffield) # using hs=-hi then field is +ve not -ve, 
    
            hiflipminus=-hc*bracket+field*(roottwohffield) #see which way fields flips
            trialtemp=0

            if (hc >= (2*ht*phitemp)): #still loop over each hys, +1

                if ((hi > hiflipminus) and (hi < hiflipplus)): #+2 blocked

                    if (blockg[i] == 0) or (blockg[i] == 2) or (blockg[i] == -2): #+3 prev blocked until this point
                        
                        if (hi >= (field*roottwohffield)): #+4
                            
                            blocktemp[i] = -1
                            
                        else:
                            #print('blockg = 0,2,-2 and hi < field*rootwo', 'blockg', blockg[i], 'hi', hi, 'field*roottwo', field*roottwohffield)
                            blocktemp[i] = 1 #end +3 unsure if sholud ended both or just one
                            
                    elif (blockg[i] == -1): #this line, already block stay , not need
                        
                        blocktemp[i] = -1
                

                    elif (blockg[i] == 1):
                        blocktemp[i] = 1

                elif (hi >= hiflipplus):#else field blocking above ht, this is hi test hiflip etc
                    # if ( abs(blockg[i]) != 1) and (rate == 1): #this occurs if blocked a hysteron during cooling
                    blocktemp[i] = -2
                
                    
                else:
                    blocktemp[i] = 2

            else: #hc < 2*ht*phitemp. this is correctm meaning else above isnt

        
                if ((hi < hiflipminus) and (hi > hiflipplus)): 
                    blocktemp[i] = 0
                else: #field blocking - below hc
                    if (hi >= hiflipminus):
                        blocktemp[i] = -2

                    else:

                        blocktemp[i] = 2

            if (temp < trialtemp):
                #need print to screen
                print('blocktemp', blocktemp[i])
                
    return hfstore, histore, boltz, blocktemp

@jit(nopython = True)
def find_tm(prop, num_hyss, t, g, Ro, eul, tau, tstep):  
    l=0
    for l in range(int(num_hyss)): #just do for half points  - (prop[i,2] + (0.5*tstep))

        if (t >= (prop[l,3] + (0.5*tstep))) : #time grater than nucleation time, less than final growth time

            prop[l,4] = Ro + g*(t - prop[l,3]) #current radius

            prop[l,11] = prop[l,5] #prev volume before change it
            prop[l,5] = ((4*pi)/3.)*(prop[l,4])**3  #volume current
            Vdot = (4*pi*g*(prop[l,4])**2) #rate volume growth
            prop[l,10] = ((eul*prop[l,5])/Vdot)/(log((eul*prop[l,5])/(tau*Vdot))) # use Vdot here and overwrite for next time
        l+=1
    return prop



@jit(nopython = True)
def block_val_C(hys, prop, histore, hfstore, blocktemp, beta, num_hyss, boltz, blockg, field, afswitch, afstore, max_total, vol_rem, t, tstep):

    totalm = 0.0
    totalmoment = 0

    i=0
    num_blocking = 0
    blockper = 0.
    #print('field ofr CRM acq inside', field)
    for i in range(num_hyss):

        
        blockg[i] = blocktemp[i]
        absblockg = abs(blockg[i])
        if (t >= (prop[i,3] + (0.5*tstep))):
           
            if (absblockg == 1): #if blockg 1 or -1
                if (boltz[i] < 0.00000001) and (boltz[i] > -0.000000001): #only zero
                    #print('setting boltz as its 0 (check)', boltz[i])
                    boltz[i] = tanh((field - histore[i])/hfstore[i])
                    #print('boltz', boltz[i])
                    if (afswitch == 0):
                        num_blocking+=1
                    if (afswitch == 1):
                   # print('inside tanh reset during af demag', (((50E-6/mu0) - histore[i])/hfstore[i]))
                        boltz[i] = 0.

            if (blockg[i] == -2):
                #print('blockg = -2')
                moment = 0 #changed from zero to minues zero
                #boltz[i] = -1
            elif (blockg[i] == 2):
                #print('blockg = +2')
                moment = 0
                #boltz[i] = 1 #but need to overwrite this
            else:

                moment = blockg[i] #where does moment come from? - set to blockg if blockg not equal to 1


            
            
            if (afswitch == 1): #afswtich given in function, 1 or 0
               # afstore =  CT['afstore'] # move where relevant as no values unti demag - poss give dummy name to keep with rest of unpack of V
                #break
                #print('should not see this message, afswtich wrong')
                hi = histore[i]*(sqrt(2))

                #hys9 change due to field?
                #hc is just hk and not hf 

                hc = (sqrt(2))*(hys[i,9]*beta)*(((sin(hys[i,6])**(2./3.))+(cos(hys[i,6])**(2./3.)))**(-3./2.))
        
                af = afstore/(1000*mu0) #define later. afstore number
                if (hi >= 0) and (hi > (hc-af)): #whats this bit mean
                    moment = 0
                    blockg[i] = -2

                    boltz[i] = 0.
                if (hi <= 0) and (hi < (-hc + af)):
                    moment = 0
                    blockg[i] = 2

                    boltz[i] = 0.0
            #totalm = totalm + (abs(moment)*abs(cos(hys[i,5]))*beta*(boltz[i])*((hys[i,3]*prop[i,5])/max_total)) #add in forcdist (hys[i,3])
            totalm = totalm + (abs(moment)*abs(cos(hys[i,5]))*beta*(boltz[i])*(hys[i,3])) #*prop[i,5])/max_total)) 
            #totalm = totalm + (abs(moment)*abs(cos(hys[i,5]))*beta*(boltz[i])*((prop[i,5])/(max_total)))
            #once done totalm return field blocked hys boltz to zero so can reblock 
            if (blockg[i] == -1) or (blockg[i] == 1):
                moment_block = 1
            else:
                moment_block = 0
		
            blockper=blockper+ (abs(moment_block))*hys[i,3] #*prop[i,5]
            totalmoment=totalmoment+moment

        

    if (blockper != 0.):
        blockper= (blockper)/(vol_rem)
    else:
        #blockper = 0.
        blockper = 1E-15
    #print(num_blocking)
    return blockper, totalm, boltz, blockg

#main gCRM acquisition 

def CRM_acq(X, V):
    
    SF = X['SF']
    #print('SF', SF)
    curie_t = X['curie_t']
          
    mu0 = 4*pi*1e-7
    kb = 1.3806503e-23
    tau = 10e-9
    roottwohffield = 2**(0.5)  

    afone = 1
    afzero = 0
    num_hyss =  1000 #300000 #200000 #600000
    V['num_hyss'] = num_hyss

    tempmax = float(X['curie_t'] + 273)
    tempmin=300 
    V['tempmin'] = tempmin
    V['tempmax'] = tempmax
    tempstep=1
    #for each field
    #print('here1')
    #arrays whole thing
    V['cntfield'] = len(X['af_step'])
    cntfield = V['cntfield']
    V['sirm'] = np.zeros((9, cntfield)) #if sirm is global it should be updated and used in function ok
    V['afmag'] = np.zeros((9, cntfield)) #is this the same cntfield size
    V['arm'] = np.zeros((9, cntfield)) #if sirm is global it should be updated and used in function ok
    V['af_step'] = np.copy(X['af_step'])



    V['crm_a'] = np.zeros((9, 5000))
    V['time_a'] = np.zeros((9, 5000))
    #track field used
    fields = np.zeros(10)    

    field = (X['min_field']*1E-6)/mu0
    fieldmax = (X['max_field']*1E-6)/mu0
    fieldstep = (fieldmax-field)/3.
    fields_list_muT = np.array([25., 50., 100., 150.])
    field_arr = (fields_list_muT*1E-6)/mu0
    ct = 0
    #arm_acq_steps_list = np.array([3001, 101])

    for ct in range(len(field_arr)): #(field < (fieldmax)): 
        #print('inside field loop')
        field = field_arr[ct]
        #set CRM field and v quick CRM acquistion
        
        #print('field start CRM acq', field)
        #arm_acq_step = arm_acq_steps_list[count_no]
        num_hyss = 70000 #300000 #hysteron_list[count_no]
        V['num_hyss'] = num_hyss       
        print('number hys', num_hyss) 
        #print('arm_acq_number_steps', arm_acq_step)
        blockper = 0.0
        V['blockper'] = blockper
        hcstore = np.zeros(num_hyss)
        V['hcstore'] = hcstore
        histore = np.zeros(num_hyss)
        V['histore'] = histore
        V['totalm'] = 0
        V['afstore'] = 0. #when TRM acq

        i=0 
        temp=300
        tempt=300

        tm = 0.2 
        V['tm'] = tm
        V['beta'] = (1-((temp-273)/578.0))**0.43 #MATCH TO INPUTTED CURIE T
        V['ms'] = ms0*V['beta']
        ms = V['ms']
        hys1, num_pop = pop_hys(num_hyss, X, V, SF)   
        #print('pop hys')



        
        ifield = ct
        g = X['growth_rate']

        blockg = np.zeros(num_hyss) 
        boltz = np.zeros(num_hyss)
        blocktemp = np.zeros(num_hyss) 
        temp = 300.
        V['blocktemp'] = blocktemp
        V['boltz'] = boltz
        V['blockg'] = blockg
        V['totalm'] = 0

        V['temp'] = temp

        rate = 1.
        t = 0
        V['rate'] = rate


        i=0 
        temp=300
        tm = 0.2 
        V['tm'] = tm

        afswitch = afzero

        V['hys'] = copy.deepcopy(hys1)


        V = blocked_hys(V)

              

        hys = V['hys']
        #print('ms', ms)
        Ro = 1E-09 #nucleation radius - same for all particles 
        prop = np.zeros((num_hyss,12))
        prop[:,9] = ms*mu0*hys[:,9] #eo for constant de/dv but unsure about hk being constant lose mu0
        prop[:,0] = hys[:,10] #Vf
        prop[:,1] = (3*prop[:,0]/(4.*pi))**(1./3.) #Rf (radius from forc arond 40 nm)
        final_R = prop[:,1]
        prop[:,2] = (prop[:,1] - Ro)/g #growth time in seconds
        maxt = np.max(prop[:,2]) #max growth time 
        prop[:,3] = (maxt- prop[:,2])#nucleation time for each 

        V['prop'] = prop
        V['totalvol'] = sum(prop[:,0])
        
        tm = 0

        V['tm'] = tm
        track = 0

        #CRM acquition

        t = 0
        tstep = maxt/8000. #chane from 0.5 (change back )
        V['tstep'] = tstep

        tot_time_steps = (np.max(prop[:,2])+0.5*tstep)/tstep
        tot_time_steps = ceil(tot_time_steps)

        tm = 0

        V['aftotalm'] = 0
        V['tm'] = tm #to use for af demag bit delcare here
        x=0 #loop over time steps
        #print('field', field)
        #temperature for CRM acquisition
        temp  = 300. #350. + 273. 
        V['temp'] = 350. + 273. 
        V['beta'] = (1-((temp-273)/578.0))**0.43 #MATCH TO INPUTTED CURIE T
        while (t <= (np.max(prop[:,2]) + 0.5*tstep)): #big main model for 1 field #remove  + 0.5*tstep 
            #print('crm acq')
            #find volumes

            prop = find_tm(prop, num_hyss, t, g, Ro, eul, tau, tstep)
            #print('here2')
            V['track'] = x
            V['prop'] = prop
            blockfind_SC(t, field, afzero, V, X)
            #print('field during CRM acq', field)
            #print(t, field, afzero, V, X)
            #print('here6')
            #print('blockper', V['blockper'])
            V['crm_a'][ifield,x] = V['totalm'] # /total_vol_step
            V['time_a'][ifield,x] = t

            x+=1
            
            if (t <= (0.4*np.max(prop[:,2]) + 0.5*tstep)):
                t+= 120*tstep
            else:
                t+= tstep #time steps of 2 hours - run about 6 times

        t = t - tstep

        rate = 0
        V['rate'] = rate
        tm = 60
        V['tm'] = tm
        temp = 300.
        beta = (1-(temp-273)/578.0)**0.43
        V['beta'] = beta
        fieldzero = 0.0

        #print(t, )

        #cool it down to room temperature and measure the CRM 
        #first approx nothing else block so just times by the new Ms 
        
        blockfind_SC(t, fieldzero, afzero, V, X) #check if totalm, spcounti work still

        #plot CRM acquisition
        plt.scatter((V['time_a'][ifield,:x]/(60*60)), V['crm_a'][ifield,:x])
        plt.xlabel('time in hours')
        plt.ylabel('crm no units')
        plt.title('CRM acquisition for field = {} muT'.format((field*mu0)*1E6))
        plt.savefig('CRM_acq_field_{}.png'.format(ifield))
        plt.close()
        track+=1
       
        #AF demag

        #calc AF demag of CRM at 273 fK
        V['afmag'][ct, 0] = V['totalm'] #/(sum(hys[:,10]))#/num_hyss #set value here #afmag is 2 thing? and check right totalm
        for i in range(cntfield): #try 1,cntfield
            V['afstore'] = V['af_step'][i]
            blockfind_SC(t, fieldzero, afone, V, X) #unsure about aftotalm - output
            
            #print('CRM demag', V['totalm'])
            V['afmag'][ct, i] = V['totalm'] #*4. #/(sum(hys[:,10]))#/num_hyss #af mag for each field and cnt step - 3 lines of total m for af demag
            #print('CRM', V['totalm']) #*4.)
        track2 = 0
        temp = 300.1



        #ARM acquisitin 300 K after NRM demag

        tm = 60.0 
        prop[:,10] = 60.0 #lab timescales for ARM
        V['prop'] = prop        
        #all blocked at the start but should unblock all first with big AF field
        V['arm_steps'] = np.linspace(100,0,3001)
        #V['blocktemp'] = np.ones(num_hyss) #is it all saturation before ARM acqusiton - no there is no magnetisation
        #V['boltz'] = np.ones(num_hyss)
        #V['blockg'] = np.ones(num_hyss)
        V['arm_blocked'] = np.zeros((num_hyss)) #treat asif none blocked
        #V['blockper'] = 1. # start with non blocked

        #ARM acquisition simulate smooth decrease in intensity with 100 - 0 mT AF field in 1 mT steps, once blocked say blocked ?
        #expect result to be same as if do 100 mT and 0 mT

        

        V['afstore'] = 100. #V['arm_steps'][0] #AF field used for ARM acquisition
        blockfind_SC_arm(t, fieldzero, afone, V, X)
        V['arm'][ifield,i]
        #initial ARM array for demag, same size as SIRM and afmag



        for i in range(len(V['arm_steps'])):
            V['afstore'] = V['arm_steps'][i] #AF field used for ARM acquisition
            blockfind_SC_arm(t, fieldzero, afone, V, X)

            #plt.tricontourf(V['hys'][:,0], V['hys'][:,1], V['boltz'][:], levels=np.linspace(-1.1,1.1,100))
            #plt.colorbar()
            #plt.savefig('contour_final_distribution_{}.png'.format(V['afstore']))
            #plt.close()
            #print('blockper outside function', V['blockper'])
            #current_hi = V['hys'][:,1]*V['blockper']
            #print('original hi', V['hys'][:,1])
            #print('current blockper', V['blockper'])
            #print('current hi', current_hi)

            #if (V['blockper'] > 0.2):
            #    plt.tricontourf(V['hys'][:,0], current_hi[:], V['boltz'][:], levels=np.linspace(-1.1,1.1,100))
            #    plt.colorbar()
            #    plt.savefig('contour_current_distribution_{}.png'.format(V['afstore']))
            #    plt.close()

                #plt.tricontourf(V['hys'][:,0], V['hys'][:,1], V['boltz'][:], levels=np.linspace(-1.1,1.1,100))
                #plt.colorbar()
                #plt.title('Preisach space hc, hi, boltz arm 0 dc acq step {}'.format(V['afstore']))
                #plt.ylabel('hi')
                #plt.xlabel('hc')
                #plt.savefig('preisach_arm_acq_boltz_0dc_{}.png'.format(V['afstore']))
                #plt.close()                

            #print('arm acq', V['totalm'])

        #V['arm'] =  #track totalm dueing ARM acquisition
        tm = 60.0 
        prop[:,10] = 60.0 #lab timescales for ARM
        V['prop'] = prop




        #plt.tricontourf(V['hys'][:,0], V['hys'][:,1], V['boltz'][:], levels=np.linspace(-1.1,1.1,100))
        #plt.colorbar()
        #plt.title('Preisach space hc, hi, boltz arm 0 dc acq step {}'.format(V['afstore']))
        #plt.ylabel('hi')
        #plt.xlabel('hc')
        #plt.savefig('preisach_arm_acq_boltz_0dc_{}.png'.format(V['num_hyss']))
        #plt.close()           
        #acquire ARM
        #V['afstore'] = 100. #AF field used for ARM acquisition
        #blockfind_SC_arm(t, fieldzero, afone, V, X)

        #V['afstore'] = 0. #AF field used for ARM acquisition
        #blockfind_SC_arm(t, fieldzero, afone, V, X)

        #V['arm_track'] = np.zeros(len(V['arm_steps']))
        V['arm'][ifield,0] = V['totalm']  #just running once, same idea as afmagstore ?? 
        print('arm', V['totalm'])

        #AF demag ARM

        for i in range(cntfield): 
            afstore = V['af_step'][i]
            V['afstore'] = afstore            
            blockfind_SC(t, fieldzero, afone, V, X) #af demag arm as normal like did for afmag and sirm
            V['arm'][ifield,i] = V['totalm'] #*4. #rewrite the first point?
            #print('arm demag', V['arm'][ifield,i])

            #plt.tricontourf(V['hys'][:,0], V['hys'][:,1], V['boltz'][:], levels=np.linspace(-1.1,1.1,100))
            #plt.colorbar()
            #plt.title('Preisach space hc, hi, boltz arm (0 dc) demag step {}'.format(V['afstore']))
            #plt.ylabel('hi')
            #plt.xlabel('hc')
            #plt.savefig('preisach_arm_demag_boltz_0dc_{}.png'.format(V['afstore']))
            #plt.close()



        V['arm_mag']= V['arm']
        #print('ams', field, V['arm'][ifield,0])        
        #print('nrm',V['afmag'][ct, 0])

        #calc SIRM at 300 K
        blockg = np.ones(num_hyss)
        boltz = np.ones(num_hyss)
        blocktemp = np.ones(num_hyss) 
        V['blocktemp'] = blocktemp
        V['boltz'] = boltz
        V['blockg'] = blockg
        #print('blockg match one above', V['blockg'])
        fieldzero = 0.0
        #aftotalm referred to as sir
        V['afstore'] = 0. #for first value
        blockfind_SC(t, fieldzero, afone, V, X)

        
        V['sirm'][ifield,0] = V['totalm'] #just running once, same idea as afmagstore ?? 


        for i in range(cntfield): 
            afstore = V['af_step'][i]
            V['afstore'] = afstore

            blockfind_SC(t, fieldzero, afone, V, X)
            V['sirm'][ifield,i] = V['totalm'] #rewrite the first point?


        #impart ARM aswell 
        #V['arm'] =  #track totalm dueing ARM acquisition

        fields[ct] = field
        #field = field + fieldstep 

    
        #negs = 0.0 #negative sirm??
        
        #negsirm[ct] = -negs

        #V['aconst'] = aconst
        #V['afmag'] = afmag
        #V['sirm'] = sirm #use to plot rom dictionary
        #V['ifield'] = ifield
        V['fields'] = fields
   
        ifield = ifield +1 
        #V['aconst'] = aconst

        V['ifield'] = ifield
        #V['fields'] = fields
        #V['track'] = track
        #V['track2'] = track2        
       
        ct +=1


    return X,V

def plot_file_arm(V):
    arm_int = np.copy(V['arm'][:,0]) # gives me the maximum value
    crm_int = np.copy(V['afmag'][:,0]) # first af demag step for crm i think - check it from 
    #hysteron_list_copy = np.copy(V['hysteron_list'])
    arm_acq_steps_list_copy = np.copy(V['arm_acq_steps_list'])
    plt.scatter(arm_acq_steps_list_copy, arm_int[:len(arm_acq_steps_list_copy)])
    plt.xlabel('number arm acq steps')
    plt.ylabel('max ARM intensity')
    plt.title('max ARM intensity vs number arm acq steps')
    plt.savefig('arm_steps_vs_arm.png')
    plt.close()
    arm_file = open("arm_steps_vs_arm.txt", "w") #'w' will overwrite file is already exists
    arm_file.write('arm_steps' + '\t' + 'arm_intensity' + '\t' + 'crm intensity' '\n') 
    for write_list in range(len(arm_acq_steps_list_copy)):
        arm_file.write(str(arm_acq_steps_list_copy[write_list]) + '\t' + str(arm_int[write_list]) + '\t' + str(crm_int[write_list]) + '\n')
    arm_file.close()



def sirm_test(V, X):
    sirm = V['sirm']
    cntfield = X['cntfield']
    name = X['name']
    ifield = V['ifield']
    w, h = figaspect(1)
    fig, ax = plt.subplots(figsize=(w,h))
    demagstep = X['af_step']
    #set first value as v v low just so remove zeros to plot easily
    demagstep2 = demagstep
    demagstep2[0] = 0.0001

    demagstep2 = demagstep2[demagstep2 != 0]

    sirm2 = sirm
    sirmn = np.copy(sirm2)

    for i in range(V['ifield']):
        for j in range(cntfield): #change back to 23 
            sirmn[i][j] = (sirm2[i][j]/(np.mean(sirm2[i][0:3])))


    af_sirm_n = X['af_irm']

    norm = np.mean(af_sirm_n[0:3])
    

    af_sirm_n_n = np.copy(af_sirm_n)
    for i in range(len(af_sirm_n)):
        af_sirm_n_n[i] = af_sirm_n[i]/norm

    sirm_p = sirmn[:ifield,:cntfield] #change 4 to no fields

    plt.plot(demagstep2, sirm_p[2], marker = 'o', color = 'r') #x = 22, y = 100
    
    plt.plot(demagstep2, af_sirm_n_n, marker = 'o', color = 'b')        
    plt.ylabel('SIRM (normalized %)')
    plt.xlabel('af demag step')
    plt.title("Demagnetisation of the SIRM for sample '{0}'".format(name))

    plt.text(30, 0.9, r'measured ', fontsize=12) #poss change to percentage across line (max/50 etc)
    plt.plot(25,0.9, marker = 'o', color='b')
    plt.text(30, 0.8, r'simulated', fontsize=12)
    plt.plot(25, 0.8, marker = 'o', color='r')
    plt.savefig('sirm_test_{}.pdf'.format(name))
    plt.show
    return
    
def calc_pal(X, V):
    cntfield = X['cntfield']
    ifield = V['ifield']
    afmag = V['afmag']
    sirm = V['sirm'] #sirm seems global so shoulld have worked anwyay
    af_nrm_n = X['af_nrm']
    af_sirm_n = X['af_irm']
    name = X['name']
    aconst = V['aconst']
    fields = V['fields']
    af_step = X['af_step']
    tempmin = V['tempmin']
    tempmax = V['tempmax']
    afnegsirm = 0.0         
    averageme = 0.0
    averageunweight = 0.0
    averagemeadj = 0.0
    averagecount = 0
    flatcount = np.zeros((10))
    noflatsub = 0
    averageflat = np.zeros((10))
    trmsirmratio = 0.0
    sigmtot = 0.0
    selrms = 0.0
    sumegli = 0

    ###################
    afpick = X['af_pick'] #6 #pick from user - Z-plot

    #######################
    sit = ifield
    std = np.zeros((100))
    ys = np.zeros((100))
    used = np.zeros((50))
    flatused = np.zeros((50))
    ystore = np.zeros((100))
    sigmstore = np.zeros((100))
    flatstore = np.zeros((10,100))
    flatdi = np.zeros((10))
    sigmflat = np.zeros((10))
    flatmedianme = np.zeros((10))
    sigmtotflat = np.zeros((10))
    cntfieldaf = cntfield
    actfield = 1 #random number but need to see what this is
    dispersion = 0
    vtwo = 0
    fortyone = 0

    shortdivlist = []
    shortminuslist = []
    ilist = []
    af_step_list = []


    flat = 0 #set here as well as earlier to save time when testing
    afratio = 0
    #do the shortdiv etc seperate

    for i in range(cntfieldaf): #cntfieldsaf

        xlist = []
        ylist = []
        sumx = 0.
        sumy = 0.
        sumxy = 0.
        sumxx = 0.

        for j in range(sit): #for each field - a
            
            
            y = fields[j]*mu0*1e6 #field time something, x, y ust values
            x = afmag[j,i]/sirm[j,i] #norm af mag to sirm - sirm after that stage of demagetisation
            #plotting points
            xlist.append(x)
            ylist.append(y)

            if (sirm[j,i] == 0):

                for i in range(afpick, averagecount): #this all moves to diff point
                    #afpick given by file called sirmms.dat
                    dispersion = dispersion + (((ystore[i] - averageme/sigmtot))**2)/sigmstore[i]
                    vtwo = vtwo + (1/sigmstore[i])**2
                    fortyone = 1

            sumx = sumx + x        
            sumxx = sumxx + x**2
            sumxy = sumxy + x*y
            sumy = sumy + y

        mfit=(((sit+000.0)*sumxy)-sumx*sumy)/(((sit+000.0)*sumxx)-sumx*sumx)
        cfit=(sumy*sumxx-sumx*sumxy)/((sit+000.0)*sumxx-sumx*sumx)

        xtest = np.linspace(min(xlist), max(xlist), 10)
        #print('xtest', xtest)
        xtest = np.array(xtest)
        ytest = np.copy(xtest)
        ytest = mfit*xtest + cfit
        sumdi = 0.
        xlist2 = []
        ylist2 = []
        dilist = []
        for j in range(sit): #x and y values set for each field again
            y= fields[j]*mu0*1e6
            x=(afmag[j,i]/sirm[j,i]) #same as above, line 1024
            di = y - (mfit*x+cfit) #see how closely x and y fit with line equation
            dilist.append(di)
            sumdi = sumdi + di**2
            xlist.append(x)
            ylist.append(y)
        sumdi = sumdi/(sit-2) 
        sigm = sit*sumdi/(sit*sumxx-sumx*sumx) #variance
        x = af_nrm_n[i]/af_sirm_n[i]   
        y = x*mfit +cfit #field assoacted witht his step for demag  
        xlist2.append(x)
        ylist2.append(y)

        ys[i] = y

        std[i] = sqrt(sigm)*x
        used[i] = 0 #set used for each cntfield to zero
        flatused[i] = 0
        
        #here
        shortdiv=abs(1-((sirm[j-1,i]/np.mean(sirm[j-1,0:3]))/(af_sirm_n[i]/np.mean(af_sirm_n[0:3]))))*100 #do for 3rd field
        shortminus=abs(((sirm[j-1,i]/np.mean(sirm[j-1,0:3]))-(af_sirm_n[i]/np.mean(af_sirm_n[0:3]))))*100
        shortdivm=abs(1-((sirm[j-2,i]/np.mean(sirm[j-2,0:3]))/(af_sirm_n[i-1]/np.mean(af_sirm_n[0:3]))))*100
        shortminusm=abs(((sirm[j-2,i]/np.mean(sirm[j-2,0:3]))-(af_sirm_n[i-1]/np.mean(af_sirm_n[0:3]))))*100
        shortdivlist.append(shortdiv)
        shortminuslist.append(shortminus)
        af_step_list.append(af_step[i])
        if (i >= afpick) and (sumx != 0.0): 

            sigm = sigm*(x**2)


            afratio = afmag[2,i]/afmag[2,0]

            if (i >= 1): #greater than 1 as 0 counts here
                afratiom = (afmag[2, i-1]/afmag[2,0]) #also do ratio for prev point - poss compares afratio and afratiom

            shortdiv=abs(1-((sirm[j-1,i]/np.mean(sirm[j-1,0:3]))/(af_sirm_n[i]/np.mean(af_sirm_n[0:3]))))*100 #do for 3rd field
            shortminus=abs(((sirm[j-1,i]/np.mean(sirm[j-1,0:3]))-(af_sirm_n[i]/np.mean(af_sirm_n[0:3]))))*100
            shortdivm=abs(1-((sirm[j-2,i]/np.mean(sirm[j-2,0:3]))/(af_sirm_n[i-1]/np.mean(af_sirm_n[0:3]))))*100
            shortminusm=abs(((sirm[j-2,i]/np.mean(sirm[j-2,0:3]))-(af_sirm_n[i-1]/np.mean(af_sirm_n[0:3]))))*100
            sigm=abs(1-sigm)

            #shortdivlist.append(shortdiv)
            #shortminuslist.append(shortminus)
            ilist.append(i)
            #af_step_list.append(af_step[i])
            if (i > 0) and (i < 30): #30 seems arbituary

                if (af_nrm_n[i]/af_nrm_n[0] > 0.01) and (afratio > 0.01): #stage requires af ratio been given value but not needed

                    if (y > 0.0):
                        #print('y>0', y)

                        if (shortdiv < 100): #acceptable ranges of S - seen in data displayed
                            if (shortminus < 20):
                                #print('conds are met')
                                selrms=selrms+abs((af_sirm_n[i]/af_sirm_n[0])- (sirm[j-1,i]/sirm[j-1,0]))               

                                sigmtot=sigmtot+(1/sigm) # sum of 1/variance

                                sumegli=sumegli+(y/actfield) #sum of fields used in predied

                                averageme=averageme+(y)/sigm 

                                averageunweight=averageunweight+y

                                averagecount=averagecount+1

                                ystore[averagecount]=y #not want averagecount to skip points?

                                sigmstore[averagecount]=sigm

                                used[i]=1 #array, this means use tihs point in calc? therefore add to used array
                                trmsirmratio=trmsirmratio+x #x is true trm/sirm measred ratio

                                if (i > 1):
                                    flatdiff = abs(y - ys[i-1])/max(y,ys[i-1])

                                    if (flatdiff < 0.2) or (y < 3.0 and ys[i-1] < 3.0): #if flatdiff < 0.2 or both <3

                                        if (noflatsub == 0): # move onto new section

                                            flat = flat +1

                                            if (i-1 < afpick) or (shortdivm > 100) or (shortminusm > 20) or (ys[i-1] < 0.0) or (af_nrm_n[i-1]/af_nrm_n[0] < 0.01) or (afratiom < 0.01):
                                                #print('rejecting prev point as out of bounds or not close enough, uses lookin 2 points back?')
                                                pass
                                            else:

                                                flatcount[flat] = flatcount[flat] + 1

                                                flatstore[flat][int(flatcount[flat])] = ys[i-1] #index erroe - change to int

                                                averageflat[flat] = averageflat[flat] + ys[i-1] #page 1121

                                                sigmflat[flat] = sigmflat[flat] + (1/sigm) #sum of 1/variance

                                                flatused[i-1] = flat

                                                flatdi[flat] = flatdi[flat]+(y-ys[i-1])/max(y,ys[i-1])

                                        noflatsub = 1

                                        flatcount[flat] = flatcount[flat] + 1

                                        flatstore[flat][int(flatcount[flat])] = y

                                        averageflat[flat] = averageflat[flat] + y #add one this point to list, some overwrite what above some diff
                                        sigmflat[flat] = sigmflat[flat] + (1/sigm) # sum of 1/variance
                                        flatdi[flat] = flatdi[flat] + (y - ys[i-1])/max(y,ys[i-1])
                                        flatused[i] = flat
                                    else: #one indent because one end if ending the else which not needed in python

                                        noflatsub = 0

                                else:
                                    noflatsub = 0

                            else:
                                noflatsub = 0


                        else:
                            noflatsub = 0

                    else:
                        noflatsub =0

                else:
                    noflatsub = 0


            else:
                noflatsub = 0

        else: 
            noflatsub = 0
    if (averagecount != 0):
        selrms=100*(selrms)/(1.0*averagecount) #line 1179

    #calc median
    temptwo = np.empty((100))
    temptwo[:] = np.nan

    for i in range(averagecount):
        temptwo[i] = ystore[i]

    temptwo = temptwo[~np.isnan(temptwo)] #want to avoid extra zeros
    medianme = np.median(temptwo)

    for jf in range(flat):

        flatc = 0

        for i in range(int(flatcount[jf])): #line 1207 - called calc flat median

            if (i >= 1): #multiple points to ave over - shoudl it be greater than just 2?

                if (flatstore[jf,i] > 1.0*flatstore[jf,i-1]): #unsure if shold be i-1 or i
                    flatc = flatc +1

                if (flatstore[jf,i] < 1.0*flatstore[jf,i-1]): #looking at point behind it
                    flatc = flatc - 1


        if (abs(flatc) == (flatcount[jf] -1)) and (flatcount[jf] > 3): #unsure if should be greater than 2 or greater than 3

            if (abs(flatstore[jf,0] - flatstore[jf,flatcount[jf]]) > 10):
                flatcount[jf] = 0
                #go to 911

                break #unsure if this is correct
        else:
            pass
            #print('conds not met') #1st point do this as zero

        if (flatcount[jf] <= 1): #assumes just needs 2 points included

            flatcount[jf] = 0
            #go to 911
            break #cut of whole loop - cut of loop for jf in range flat - seems correct

        #calc flat median #for i in range(flatcount[jf] - 1): #unsure if this is the correct limit
        flatmedianme[jf] = np.median(temptwo[:, jf])   #temptwo[i] = flatsotre range for i in length flatcount add in values , try just do upto certain amount


    #set dispersion etc - line 1257
    if (fortyone == 0):
        dispersion = 0.0
        vtwo = 0.0 #but not if reached 41 - in 41 section - set variable to one 

    for i in range(afpick, averagecount): #this all moves to diff point
        #afpick given by file called sirmms.dat
        dispersion = dispersion + (((ystore[i] - averageme/sigmtot))**2)/sigmstore[i]
        vtwo = vtwo + (1/sigmstore[i])**2


    #flat section normal districution standard distribution
    for i in range(flat):
        sigmtotflat[i] = 0
        for k in range(int(flatcount[flat])):

            sigmtotflat[i] = sigmtotflat[i] + ((flatstore[i,k] - (averageflat[i]/flatcount[i]))**2)/flatcount[i]

        sigmtotflat[i] = sqrt(sigmtotflat[i])
        #end of flat normal

    sirmrms = 0
    print("All demag and palaeointensity data is saved in afdemag-all-{0}.dat".format(name))
    print("All demag and data regarding the SIRM data is saved in afdemag-sirm-{0}.dat".format(name))
    fall = open("afdemag-all-{0}.dat".format(name), "w") #'w' will overwrite file is already exists
    fall.write('afield' + '\t' + 'step' + '\t' + 'stepPI' + '\t' + 'std' + '\t' + 'select' + '\t' + 'flatno'+ '\t' + 'shortminus' + '\t' + 'SIRMAFS%' + '\t' + 'SIRMAFM%' + '\t' + 'AFNRM/SIRM-M%' + '\t' + 'AFNRM/SIRM-S%' + '\t' +  'AFNRM/NRM-M%' + '\t' + 'AFNRM/NRM-S%' + 'shortdiv%')
    fsirm = open("afdemag-sirm-{0}.dat".format(name), "w")

    fsirm.write('afield' + '\t' + 'measured' + '\t' + 'simulated')
    for i in range(cntfieldaf): #0 to 21
        sirmrms = sirmrms + abs((af_sirm_n[i]/af_sirm_n[0]) - sirm[j-1, i]/sirm[j-1][0])
        #print(sirmrms) #increases 0.004 -> 2.6
        fall.write('\n')
        #f.write(str(af_step[i]))
        fall.write(str(af_step[i]) + '\t' + str(i) + '\t' + str(ys[i]) + '\t' + str(std[i]) + '\t' + str(used[i]) + '\t' + str(flatused[i]) + '\t' + str((((sirm[j-1,i]/sirm[j-1,0])-(af_sirm_n[i]/af_sirm_n[0])))*100) + '\t' + str(sirm[j-1,i]/sirm[j-1,0]*100) + '\t' + str((af_sirm_n[i]/af_sirm_n[0])*100) + '\t' + str((af_nrm_n[i]/af_sirm_n[i])*100) + '\t' + str((afmag[j-1,i]/sirm[j-1,i])*100) + '\t' + str((af_nrm_n[i]/af_nrm_n[0])*100) + '\t' + str((afmag[j-1,i]/afmag[j-1,0])*100) + '\t' + str(abs(1-((sirm[j-1,i]/sirm[j-1,0])/(af_sirm_n[i]/af_sirm_n[0])))*100))
        fsirm.write('\n')
        fsirm.write(str(af_step[i]) + '\t' + str(af_sirm_n[i]/af_sirm_n[0]) + '\t' + str(sirm[j-1, i]/sirm[j-1,0]))

    fsirm.close()

    fall.write('\n')

    sirmrms = sirmrms/cntfieldaf

    fall.write('SIRM MEAN DIFFERENCE % =' + '\t '+ str(100*sirmrms))
    fall.write('\n')
    if (averagecount !=0) and (averagecount !=1) and (sigmtot !=0): #only do both is average at least 1 
        sampleunbias=dispersion*(sigmtot/((sigmtot**2)-vtwo))
        dispersion=dispersion/(averagecount-1)


    if (averagecount == 1): #here equals 11
        dispersion = 1

    if (sigmtot == 0):
        sigmtot = 1
        
        print('sigm tot = 0 weighted average not possible')
        
    if (averagecount == 0):
        averagecount = 1
        print('avercount = 0 weighted average not possible')   
        
    #fall.write('weighted average = ' + '\t' + str(averageme/(sigmtot)) + '\t' + str(sqrt(dispersion/sigmtot)))
    #fall.write('\n')
    #fall.write('unweighted average = ' + '\t' + str(averageunweight/averagecount))
    #fall.write('\n')
    #fall.write('unweighted median = ' + '\t' + str(medianme))
    #fall.write('\n')
    #fall.write('cntfields = (cntfield, cntfieldaf, afpick) ' + '\t' + str(cntfield) + '\t' + str(cntfieldaf) + '\t' + str(afpick))


    #determining which is best flat 
    maxjf = 0
    minsig = 10000
    jfs = 0

    for jf in range(flat):

        if ((flatcount[jf] >= maxjf) and (flatcount[jf] > 0)):

            if ((flatcount[jf] == maxjf) or (maxjf == 1)): #change maxjf == 0

                if (sigmtotflat[jf] < minsig): #seem to do same thing regardnless of these if statements

                    minsig = sigmtotflat[jf]
                    maxjf = flatcount[jf]
                    jfs = jf

            else:
                minsig = sigmtotflat[jf]
                maxjf = flatcount[jf]
                jfs = jf


        if (flatcount[jf] > 0):

            fall.write('\n')
            fall.write('unweighted average flat (jf, averageflat(jf)/flatcount(jf), sigmtotflat(jf, flatdi(jf)) =' + '\t' + str(jf) + '\t' + str(averageflat[jf]/flatcount[jf]) + '\t' + str(sigmtotflat[jf]) + '\t' + str(flatdi[jf]) + '\t' + str(flatcount[jf]) + '\t' + str(flatdi[jf]/flatcount[jf]))
            fall.write('\n')
            fall.write('unweighted flat median (jf, flatmedianme(jf) = ' + '\t' + str(jf) + '\t' + str(flatmedianme[jf]))
        else:

            fall.write('no points selected for flat section (jf)' + '\t' + str(jf))


    fall.close()        
    aconst = -aconst*log(0.01*(tempmin)/(tempmax-tempmin))/3600
    print('Output results to averageout_{0}.dat'.format(name))
    fave = open('averageout_{}.dat'.format(name), 'w') #117 file
    if (averageme > 0):

        if (jfs != 0):

            fave.write(str(averageme/(sigmtot)) + '\t' + str(sqrt(dispersion/sigmtot)) + '\t' + str(averagecount) + '\t' + str(averageunweight/averagecount) + '\t' + str(medianme) + '\t' + str(flatcount[jfs]) + '\t' +  str(averageflat[jfs]/(1.0*flatcount[jfs])) + '\t '+  str(sigmtotflat[jfs]) + '\t' + str(flatmedianme[jfs]) + '\t' + str(aconst) + '\t' + '147' + '\t' + str(100*af_nrm_n[0]/af_sirm_n[0]) + '\t' + str(sirmrms) + '\t' + str(selrms))
        else:
            fave.write(str(averageme/sigmtot) + '\t' + str(sqrt(dispersion/sigmtot)) + '\t' + str(averagecount) + '\t' + str(averageunweight/averagecount) + '\t' + str(medianme) + '\t' + '0.0' + '\t' + '0.0' + '\t' + '0.0' + '\t' + str(aconst) + '\t' + '147' + str(100*(af_nrm_n[0]/af_sirm_n[0])) + '\t' + str(sirmrms) + '\t' + str(selrms))
    else:
        fave.write('no points selected' + '\t' + str(sirmrms))


    fave.close()
    V['shortdivlist'] = shortdivlist
    V['shortminuslist'] = shortminuslist
    V['ys'] = ys
    
    return(X,V)
    
    
def plot_sirm_check(X,V):
    af_step_list = X['af_step']
    shortdivlist = V['shortdivlist']
    shortminuslist = V['shortminuslist']
    name = X['name']
    shortdivlistplot = np.array(shortdivlist[:(len(af_step_list))])
    shortminuslistplot = np.array(shortminuslist[:(len(af_step_list))])
    shortdivlistplot[np.isnan(shortdivlistplot)] = 0
    shortminuslistplot[np.isnan(shortminuslistplot)] = 0
    twenty = []
    hundred = []

    for i in range(len(af_step_list)):
        twenty.append(20)
        hundred.append(100)
    w, h = figaspect(1) #from https://stackoverflow.com/questions/48190628/aspect-ratio-of-a-plot
    fig, ax = plt.subplots(figsize=(w,h))

    plt.plot(af_step_list, twenty, 'b') #af step list is just af_step?
    plt.plot(af_step_list, hundred, 'r')
    plt.plot([af_step_list[0], af_step_list[-1]], [-20, -20], 'b') #af step list is just af_step?
    plt.plot([af_step_list[0], af_step_list[-1]], [-100, -100], 'r')
    plt.plot(af_step_list, shortdivlistplot,  marker='o', color= 'r')
    plt.plot(af_step_list, shortminuslistplot,  marker='o', color = 'b')
    plt.title('SIRM checks looking at difference between measured and simulated. Sample {0}'.format(name))

    plt.xlabel('AF peak (mT)')
    plt.ylabel('S_diff or S_ratio (%)')
    plt.text(25, 85, r'S diff ', fontsize=12)
    plt.plot(22,87, marker = 'o', color='b')
    plt.text(25, 75, r'S ratio', fontsize=12)
    plt.plot(22, 77, marker = 'o', color='r')
    print('Figure saved as SIRM-checks-{}.svg'.format(name))
    plt.savefig('SIRM-checks-{}.pdf'.format(name))
    plt.show
    return
    
def plot_pal(V,X):
    
    ys = V['ys']
    af_step = X['af_step']
    name = X['name']
    afpick = X['af_pick']
    
    w, h = figaspect(1) #from https://stackoverflow.com/questions/48190628/aspect-ratio-of-a-plot
    fig, ax = plt.subplots(figsize=(w,h))

    #plot with red colour and calc average
    #selected_mean = np.mean(ys[8:20])
    #mean_dev = np.std(ys[8:20])
    #selected_med = np.median(ys[8:20])
    plt.plot(af_step,ys[:len(af_step)], 'b')
    plt.plot(af_step,ys[:len(af_step)],  marker='o', color= 'b')
    #plot vertical line
    
    plt.plot([af_step[afpick], af_step[afpick]], [0,np.max(ys[:len(af_step)])], color = 'green')
    #plt.plot(af_step[8:20],ys[8:20],  marker='o', color= 'r')
    
    #plt.text(20, 6, r'selected median: %.2f $\mu$T'%selected_med, fontsize=11)
    #plt.text(20, 7, r'rejected mean: %.2f $\pm$ %.2f $\mu$T'%(selected_mean ,mean_dev), fontsize=11)
    #plt.text(20, 6, r'selected median: %.2f $\mu$ T'%selected_med, fontsize=12)
    plt.xlabel('AF peak (mT)')
    plt.ylabel('paleointensity (\u03BCT)')

    plt.title('TRM PI est {}'.format(name))

    plt.savefig('ESA_MF-PI_ESTS_colour_{}.pdf'.format(name))
    plt.show
    
    ind = []
    for i in range(len(af_step)):
        ind.append(i)
        i+=1


    ind = np.array(ind)

    print("AF step index = AF step")
    for n, v in zip(ind, af_step):
        print("{} = {}".format(n, v))


    



def plot_zplot(X):
    #normalise intensity between 0 and 1 to plot nicer

    nat_g = 0
    if (X['sample_copy'] > 1): #only if second go 
        g_1 = input('The most recent growth rate used for sample {} was {} hours. If you want to re-use these variables enter K, to change them enter any other charactor'.format(X['names'][0], X['growth_rate']))
    else:
        g_1 = 'L'
    if (g_1 != 'K') or (X['sample_copy'] < 2):    
        while (nat_g == 0):
            try:
                nature_growth_rate = (input("Input the natural growth rate in hours:" )) 
                nature_growth_rate = float(nature_growth_rate) #ask for interger and check it is an interger, if not ask again
                nat_g = 1
            except ValueError:
                print('Not a number. Please input an number.')
                nat_g = 0
        X['growth_rate'] = nature_growth_rate

 
    c_l = 0

    if (X['sample_copy'] > 1) : #only if second go 

        cur_1 = input('The most recent Curie temperature used for sample {} was {}\xb0C. If you want to re-use this value enter K, to change them enter any other charactor'.format(X['names'][0], X['curie_t']))
    else:
        cur_1 = 'L'
    if (cur_1 != 'K') or (X['sample_copy'] < 2):  
        while (c_l == 0):
            try:
                curie_t = (input("Input Curie temperature in \xb0C:" ))
                curie_t = float(curie_t) #ask for interger and check it is an interger, if not ask again
                c_l = 1
            except ValueError:
                print('Not a number. Please input an number.')
                c_l = 0
        X['curie_t'] = curie_t #if no Ms-T curve - set Curie T here    

    if (X['sample_copy'] > 1): #only if second go 

        afval_1 = input('The most recent AF step used for sample {} was {} mT. If you want to re-use this value enter K, to change them enter any other charactor and the Zplot will show'.format(X['names'][0], X['afval']))
    else:
        afval_1 = 'L'
    if (afval_1 != 'K') or (X['sample_copy'] < 2):      
        norm1 = np.copy(X['af_nrm'])
        norm_int = norm1/X['af_nrm'][0]

        n_nrm=norm_int*np.cos(X['af_nrm_dec']*pi/180.0)*np.cos(X['af_nrm_inc']*pi/180.0)
        e_nrm=norm_int*np.sin(X['af_nrm_dec']*pi/180.0)*np.cos(X['af_nrm_inc']*pi/180.0)
        u_nrm=-1*norm_int*np.sin(X['af_nrm_inc']*pi/180.0)
        af_step = np.copy(X['af_step'])
        h_nrm = np.sqrt((n_nrm)**2 + (e_nrm)**2)


        #plot horizontal
        x_p = np.copy(e_nrm)
        y_p = np.copy(n_nrm)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 9))
        ax1.plot(x_p, y_p, 's-', linewidth=1, markersize=4, label='E, N', color='blue')

        #plot vertical
        x_p2 = np.copy(n_nrm)
        y_p2 = np.copy(u_nrm)
        ax1.plot(x_p2, y_p2, 's-', linewidth=1, markersize=4, markerfacecolor='none', label='N, U', color = 'red')


        #ax.scatter(x_p, y_p, marker)
        xmin, xmax, ymin, ymax = min(np.min(x_p), np.min(x_p2), 0), max(np.max(x_p), np.max(x_p2), 0), min(np.min(y_p), np.min(y_p2), 0), max(np.max(y_p), np.max(y_p2), 0)
        i=0
        for x,y in zip(x_p, y_p):
            #label = "{:.0f}".format(X['af_step'][i])
            #if (i % 2 == 0) or (i == 0):
            ax1.annotate(af_step[i], (x,y), textcoords="offset points", xytext=(0,2), ha='left', fontsize = 9) # horizontal alignment can be left, right or center
            i+=1
        i=0
        for x,y in zip(x_p2, y_p2):
            #label = "{:.0f}".format(X['af_step'][i])
            #if (i % 2 == 0) or (i == 0):        
            ax1.annotate(af_step[i], (x,y), textcoords="offset points", xytext=(0,2), ha='left', fontsize = 9) # horizontal alignment can be left, right or center
            i+=1

        # Set identical scales for both axes
        ax1.set(xlim=(xmin, xmax), ylim=(ymin, ymax), aspect='equal')
        # Set bottom and left spines as x and y axes of coordinate system
        ax1.spines['bottom'].set_position('zero')
        ax1.spines['left'].set_position('zero')
        ax1.legend(loc='best')
        # Remove top and right spines
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        #ax1.title('Zijderveld plot for sample {}. Label also highlighting corresponding (X,Y) axes'.format(X['name']))
        #ax1.text((xmax), 0, 'E,N', fontsize = 16, horizontalalignment = 'right', verticalalignment = 'center')
        #af_step[i], (x,y), textcoords="offset points", xytext=(0,2), ha='left', fontsize = 9
        ax1.text((xmax), 0, 'E,N', fontsize = 16, ha='left', va='center')
        ax1.text((xmin), 0, 'W,S', fontsize = 16, ha='right', va='center')
        ax1.text(0, (ymax), 'N,U', fontsize = 16, ha='center', va='bottom')
        ax1.text(0, (ymin), 'S,D', fontsize = 16, ha='center', va='top')   

        #qax.text(right, top, 'right bottom',
        #horizontalalignment='right',
        #verticalalignment='bottom',
        #transform=ax.transAxes)



        #ax1.set_xlabel('E,N', fontsize = 16, loc ='right')
        #ax1.set_ylabel('N,U', fontsize = 16, loc ='top')
        #ax3 = ax1.twinx()
        #ax3.set_xlabel('W,S', fontsize = 16, loc ='left')
        #ax3.set_ylabel('S,D', fontsize = 16, loc ='bottom')
        #ax2 for 2nd projection 
        #plot E-U
        x_p3 = np.copy(e_nrm)
        y_p3 = np.copy(u_nrm)

        x_p4 = np.copy(h_nrm)
        y_p4 = np.copy(u_nrm)

        ax2.plot(x_p3, y_p3, 's-', linewidth=1, markersize=4, label='E,U', color='green')
        ax2.plot(x_p4, y_p4, 's-', linewidth=1, markersize=4, label='H,U', color='orange')
        #ax.scatter(x_p, y_p, marker)
        xmin, xmax, ymin, ymax = min(np.min(x_p3), np.min(x_p4), 0), max(np.max(x_p3), np.max(x_p4), 0), min(np.min(y_p3), np.min(y_p4), 0), max(np.max(y_p3), np.max(y_p4), 0)        
        i=0
        for x,y in zip(x_p3, y_p3):
            #label = "{:.0f}".format(X['af_step'][i])
            #if (i % 2 == 0) or (i == 0):    
            ax2.annotate(af_step[i], (x,y), textcoords="offset points", xytext=(0,2), ha='left', fontsize = 9) # horizontal alignment can be left, right or center
            i+=1
        i=0
        for x,y in zip(x_p4, y_p4):
            #label = "{:.0f}".format(X['af_step'][i])
            #if (i % 2 == 0) or (i == 0):    
            ax2.annotate(af_step[i], (x,y), textcoords="offset points", xytext=(0,2), ha='left', fontsize = 9) # horizontal alignment can be left, right or center
            i+=1


        # Set identical scales for both axes
        ax2.set(xlim=(xmin, xmax), ylim=(ymin, ymax), aspect='equal')
        # Set bottom and left spines as x and y axes of coordinate system
        ax2.spines['bottom'].set_position('zero')
        ax2.spines['left'].set_position('zero')
        ax2.legend(loc='best')
        # Remove top and right spines
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.set_xlabel('E, H', fontsize = 16)
        ax2.set_ylabel('U', fontsize = 16)    

        #ax2.set_xlabel('E,H', fontsize = 16, loc = 'right')
        #ax2.set_ylabel('U', fontsize = 16, loc = 'top')
        #ax4 = ax2.twinx()
        #ax4.set_xlabel('W,H', fontsize = 16, loc = 'left')
        #ax4.set_ylabel('D', fontsize = 16, loc = 'bottom')
        ax2.text((xmax), 0, 'E,H', fontsize = 16, ha='left', va='center')
        ax2.text((xmin), 0, 'W,H', fontsize = 16, ha='right', va='center')
        ax2.text(0, (ymax), 'U', fontsize = 16, ha='center', va='bottom')
        ax2.text(0, (ymin), 'D', fontsize = 16, ha='center', va='top')   

        plt.suptitle('Normalised Zijderveld plots for sample {}'.format(X['name']))
        #plot table
        #create data
    
        plt.savefig('zijderveld_plot_{}.pdf'.format(X['name']))
        plt.show()
        cntfield = X['cntfield']

        #ask if want to plot a zoomed in version
        replot = 1
        
        while (replot == 1):
            zoom = input('Would you like to zoom in the zplot to make your decison (you can do this multiple times)? (Y or N)')
            while (zoom == 'Y'):
                index_type = 0
                while (index_type == 0):
                    step = input('What AF field would you like to look from?')
                    #replot  
                    try:
                        new_index1 = np.where(af_step == float(step))
                        index_type = 1.
                        new_index = int(new_index1[0])
                        plt.close()


                        #norm1 = X['af_nrm'][new_index:]
                        #norm_int = norm1/X['af_nrm'][new_index]

                        n_nrm=norm_int[new_index:]*np.cos(X['af_nrm_dec'][new_index:]*pi/180.0)*np.cos(X['af_nrm_inc'][new_index:]*pi/180.0)
                        e_nrm=norm_int[new_index:]*np.sin(X['af_nrm_dec'][new_index:]*pi/180.0)*np.cos(X['af_nrm_inc'][new_index:]*pi/180.0)
                        u_nrm=-1*norm_int[new_index:]*np.sin(X['af_nrm_inc'][new_index:]*pi/180.0)
                        h_nrm = np.sqrt((n_nrm)**2 + (e_nrm)**2)
                        af_step_plot = np.copy(X['af_step'][new_index:])

                        #plot horizontal
                        x_p = np.copy(e_nrm)
                        y_p = np.copy(n_nrm)
                        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 9))
                        ax1.plot(x_p, y_p, 's-', linewidth=1, markersize=4, label='horizontal (E, N)', color='blue')

                        #plot vertical
                        x_p2 = np.copy(n_nrm)
                        y_p2 = np.copy(u_nrm)
                        ax1.plot(x_p2, y_p2, 's-', linewidth=1, markersize=4, markerfacecolor='none', label='vertical (U, N)', color = 'red')


                        xmin, xmax, ymin, ymax = min(np.min(x_p), np.min(x_p2), 0), max(np.max(x_p), np.max(x_p2), 0), min(np.min(y_p), np.min(y_p2), 0), max(np.max(y_p), max(y_p2), 0)
                        i=0
    
                        for x,y in zip(x_p, y_p):
                            #label = "{:.0f}".format(X['af_step'][i])
                            #if (i % 2 == 0) or (i == 0):
                            ax1.annotate(af_step_plot[i], (x,y), textcoords="offset points", xytext=(0,2), ha='left', fontsize = 9) # horizontal alignment can be left, right or center
                            i+=1
                        i=0
                        for x,y in zip(x_p2, y_p2):
                            #label = "{:.0f}".format(X['af_step'][i])
                            #if (i % 2 == 0) or (i == 0):        
                            ax1.annotate(af_step_plot[i], (x,y), textcoords="offset points", xytext=(0,2), ha='left', fontsize = 9) # horizontal alignment can be left, right or center
                            i+=1

                        # Set identical scales for both axes
                        ax1.set(xlim=(xmin, xmax), ylim=(ymin, ymax), aspect='equal')
                        # Set bottom and left spines as x and y axes of coordinate system
                        ax1.spines['bottom'].set_position('zero')
                        ax1.spines['left'].set_position('zero')
                        ax1.legend(loc='best')
                        # Remove top and right spines
                        ax1.spines['top'].set_visible(False)
                        ax1.spines['right'].set_visible(False)
                        #ax1.title('Zijderveld plot for sample {}. Label also highlighting corresponding (X,Y) axes'.format(X['name']))
                        ax1.text((xmax), 0, 'E,N', fontsize = 16, ha='left', va='center')
                        ax1.text((xmin), 0, 'W,S', fontsize = 16, ha='right', va='center')
                        ax1.text(0, (ymax), 'N,U', fontsize = 16, ha='center', va='bottom')
                        ax1.text(0, (ymin), 'S,D', fontsize = 16, ha='center', va='top')   
                        
                        #ax2 for 2nd projection 
                        #plot E-U
                        x_p3 = np.copy(e_nrm)
                        y_p3 = np.copy(u_nrm)

                        x_p4 = np.copy(h_nrm)
                        y_p4 = np.copy(u_nrm)

                        ax2.plot(x_p3, y_p3, 's-', linewidth=1, markersize=4, label='vertical (U, E)', color='green')
                        ax2.plot(x_p4, y_p4, 's-', linewidth=1, markersize=4, label='H,U', color='orange')
                        #ax.scatter(x_p, y_p, marker)
                        xmin, xmax, ymin, ymax = min(np.min(x_p3), np.min(x_p4), 0), max(np.max(x_p3), np.max(x_p4), 0), min(np.min(y_p3), np.min(y_p4), 0), max(np.max(y_p3), np.max(y_p4), 0)    
                        i=0
                        for x,y in zip(x_p3, y_p3):
                            #label = "{:.0f}".format(X['af_step'][i])
                            #if (i % 2 == 0) or (i == 0):    
                            ax2.annotate(af_step_plot[i], (x,y), textcoords="offset points", xytext=(0,2), ha='left', fontsize = 9) # horizontal alignment can be left, right or center
                            i+=1

                        i=0
                        for x,y in zip(x_p4, y_p4):
                            #label = "{:.0f}".format(X['af_step'][i])
                            #if (i % 2 == 0) or (i == 0):    
                            ax2.annotate(af_step[i], (x,y), textcoords="offset points", xytext=(0,2), ha='left', fontsize = 9) # horizontal alignment can be left, right or center
                            i+=1


                        # Set identical scales for both axes
                        ax2.set(xlim=(xmin, xmax), ylim=(ymin, ymax), aspect='equal')
                        # Set bottom and left spines as x and y axes of coordinate system
                        ax2.spines['bottom'].set_position('zero')
                        ax2.spines['left'].set_position('zero')
                        ax2.legend(loc='best')
                        # Remove top and right spines
                        ax2.spines['top'].set_visible(False)
                        ax2.spines['right'].set_visible(False)
                        ax2.text((xmax), 0, 'E,H', fontsize = 16, ha='left', va='center')
                        ax2.text((xmin), 0, 'W,H', fontsize = 16, ha='right', va='center')
                        ax2.text(0, (ymax), 'U', fontsize = 16, ha='center', va='bottom')
                        ax2.text(0, (ymin), 'D', fontsize = 16, ha='center', va='top')  
                        
                        plt.suptitle('Zoomed in normalised Zijderveld plots for sample {}'.format(X['name']))
                        #plot table
                        #create data

                        plt.savefig('zijderveld_plot_zoomed_{}.pdf'.format(X['name']))
                        plt.show()
                        zoom2 = 0.
                        while (zoom2 == 0.):
                            try: 
                                zoom = input('Would you like to zoom in again to help your decision? (Y or N)')
                                if (zoom == 'Y'):
                                    zoom2 = 1.
                                    index_type = 0
                                elif (zoom == 'N'):
                                    index_type = 1.
                                    zoom2 =1.
                            
                            except:
                                print('Not a value')
                                index_type == 0.  
                    except:
                            print('Not a value')
                            index_type == 0.                          
                    
                    #worked and ask if want another 
                


            if (zoom == 'N'):
                replot = 0
                while True:
                    afval = (input("Pick the AF demag step where the primary componant is identified:" ))
                    try:
                        afpick_l = np.where(af_step == float(afval))
                        afpick = afpick_l[0]

                        if (afpick >= 0) and (afpick <= cntfield): #within range of AF demag steps - may break if pick too high
                            #print('in bounds')
                            break
                    except ValueError:
                        print('Not an interger')
                        True
                    # if (isinstance(sf_choose, int)):
                        print('int')
            else:
                print('Enter Y or N.')
                
                
            


        X['af_pick'] = int(afpick)
        X['afval'] = afval
        print('Selected AF step:', afval)
    else:
        afpick_l = np.where(X['af_step'] == float(X['afval']))
        afpick = afpick_l[0]
        X['af_pick'] = int(afpick)
        print('Selected AF step:', X['afval'])
    return(X)

def plot_zplot_old(X):
    n_nrm=X['af_nrm']*np.cos(X['af_nrm_dec']*pi/180.0)*np.cos(X['af_nrm_inc']*pi/180.0)
    e_nrm=X['af_nrm']*np.sin(X['af_nrm_dec']*pi/180.0)*np.cos(X['af_nrm_inc']*pi/180.0)
    u_nrm=-1*X['af_nrm']*np.sin(X['af_nrm_inc']*pi/180.0)
    af_step = np.copy(X['af_step'])

    #plot horizontal
    x_p = np.copy(e_nrm)
    y_p = np.copy(n_nrm)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 10))
    ax1.plot(x_p, y_p, 's-', linewidth=1, markersize=4, label='horizontal (E, N)', color='blue')

    #plot vertical
    x_p2 = np.copy(n_nrm)
    y_p2 = np.copy(u_nrm)
    ax1.plot(x_p2, y_p2, 's-', linewidth=1, markersize=4, markerfacecolor='none', label='vertical (U, N)', color = 'red')


    #ax.scatter(x_p, y_p, marker)
    xmin, xmax, ymin, ymax = min(np.min(x_p), np.min(x_p2)), max(np.max(x_p), np.max(x_p2)), min(np.min(y_p), np.min(y_p2)), max(np.max(y_p), max(y_p2))
    i=0
    for x,y in zip(x_p, y_p):
        #label = "{:.0f}".format(X['af_step'][i])
        if (i % 3 == 0) or (i == 0):
            ax1.annotate(af_step[i], (x,y), textcoords="offset points", xytext=(0,2), ha='left', fontsize = 12) # horizontal alignment can be left, right or center
        i+=1
    i=0
    for x,y in zip(x_p2, y_p2):
        #label = "{:.0f}".format(X['af_step'][i])
        if (i % 3 == 0) or (i == 0):        
            ax1.annotate(af_step[i], (x,y), textcoords="offset points", xytext=(0,2), ha='left', fontsize = 12) # horizontal alignment can be left, right or center
        i+=1

    # Set identical scales for both axes
    ax1.set(xlim=(xmin, xmax), ylim=(ymin, ymax), aspect='equal')
    # Set bottom and left spines as x and y axes of coordinate system
    ax1.spines['bottom'].set_position('zero')
    ax1.spines['left'].set_position('zero')
    ax1.legend()
    # Remove top and right spines
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.title('Zijderveld plot for sample {}. Label also highlighting corresponding (X,Y) axes'.format(X['name']))
    ax1.xlabel('E,N', fontsize = 16)
    ax1.ylabel('N,U', fontsize = 16)
    
    #ax2 for 2nd projection 
    #plot E-U
    x_p3 = np.copy(e_nrm)
    y_p3 = np.copy(u_nrm)

    ax2.plot(x_p3, y_p3, 's-', linewidth=1, markersize=4, label='vertical (U, E)', color='green')

    #ax.scatter(x_p, y_p, marker)
    xmin, xmax, ymin, ymax = np.min(x_p3), np.max(x_p3), np.min(y_p3), np.max(y_p3)
    i=0
    for x,y in zip(x_p3, y_p3):
        #label = "{:.0f}".format(X['af_step'][i])
        if (i % 3 == 0) or (i == 0):    
            ax2.annotate(af_step[i], (x,y), textcoords="offset points", xytext=(0,2), ha='left', fontsize = 12) # horizontal alignment can be left, right or center
        i+=1

    # Set identical scales for both axes
    ax2.set(xlim=(xmin, xmax), ylim=(ymin, ymax), aspect='equal')
    # Set bottom and left spines as x and y axes of coordinate system
    ax2.spines['bottom'].set_position('zero')
    ax2.spines['left'].set_position('zero')
    ax2.legend()
    # Remove top and right spines
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.xlabel('E', fontsize = 16)
    ax2.ylabel('U', fontsize = 16)    
    
    plt.title('Zijderveld plots for sample {}'.format(X['name']))
    plt.savefig('zijderveld_plot_{}.pdf'.format(X['name']))
    plt.show()
    cntfield = X['cntfield']

    #ask if want to plot a zoomed in version


    while True:
        afval = (input("Pick the AF demag step where the primary componant is identified:" ))
        try:
            #search for value 
            afpick = af_step.index(afval)
            if (afpick >= 0) and (afpick <= cntfield): #within range of AF demag steps - may break if pick too high
                #print('in bounds')
                break
        except ValueError:
            print('Not an interger')
            True
        # if (isinstance(sf_choose, int)):
            print('int')
    X['af_pick'] = afpick
    print(afpick)

    return(X)


def pick_SF(X):
    if (X['sf_1'] != 'K'):
        sf_2 = 0.
        if (X['sample_copy'] > 1): #only if second go 
            sf_1 = input('The most recent SF used for sample {} was {}. If you want to re-use these variables enter K, to change them enter any other charactor'.format(X['names'][0], X['SF']))
        else:
            sf_1 = 'L'
        if (sf_1 != 'K') or (X['sample_copy'] < 2):   
            while (sf_2 == 0):
                try:
                    SF = (input("Input the lowest reliable SF from the FWHM graph above:" )) 
                    X['SF'] = int(SF)
                    if (int(SF) <= 9) and (int(SF) >= 2):
                        sf_2 = 1.
                    else:
                        print('Not between 2 and 5, input an interger between 2 and 5.')
                except ValueError:
                    print('Not an interger. Please input an interger.')
    return(X)

def calc_PI_checks(V,X):
    mu0 = 4.*pi*1e-7
    sirm = V['sirm']
    print('sirm', V['sirm'])

    arm = V['arm']
    print('arm', arm)
    cntfield = X['cntfield']
    name = X['name']
    ifield = V['ifield']
    demagstep = X['af_step']
    #set first value as v v low just so remove zeros to plot easily
    demagstep2 = demagstep
    demagstep2[0] = 0.0001
    demagstep2 = demagstep2[demagstep2 != 0]
    sirm2 = sirm
    sirmn = np.copy(sirm2) #this is normalising the sirm demag
    for i in range(V['ifield']):
        for j in range(cntfield): #change back to 23 
            sirmn[i][j] = (sirm2[i][j]/(np.mean(sirm2[i][0])))

    af_sirm_n = X['af_arm']
    norm = np.mean(af_sirm_n[0])
    af_sirm_n_n = np.copy(af_sirm_n)
    for i in range(len(af_sirm_n)):
        af_sirm_n_n[i] = af_sirm_n[i]/norm

    V['af_sirm_n_n'] = af_sirm_n_n

    sirm_p = sirmn[:ifield,:cntfield] #change 4 to no fields
    #write to file at some point
    V['sirm_plot'] = sirm_p[2]

    #normalise ARM demag

    arm2 = arm
    armn = np.copy(arm2) #this is normalising the sirm demag
    for i in range(V['ifield']):
        for j in range(cntfield): #change back to 23 
            armn[i][j] = (arm2[i][j]/(np.mean(arm2[i][0])))

    af_arm_n = X['af_arm']
    norm = np.mean(af_arm_n[0])
    af_arm_n_n = np.copy(af_arm_n)
    for i in range(len(af_arm_n)):
        af_arm_n_n[i] = af_arm_n[i]/norm

    V['af_arm_n_n'] = af_arm_n_n

    arm_p = armn[:ifield,:cntfield] #change 4 to no fields
    #write to file at some point
    V['arm_plot'] = arm_p[2]

    #prepare CRM demag 
    afmag = V['afmag'] #simulated af demag
    afmag2 = afmag
    afmagn = np.copy(afmag2)
    for i in range(V['ifield']):
        for j in range(cntfield): #change back to 23 
            afmagn[i][j] = (afmag2[i][j]/(afmag2[i][0])) # norm to max not first 3 

    af_nrm_n = X['af_nrm']
    norm_af = af_nrm_n[0]
    af_nrm_n_n = np.copy(af_nrm_n)
    for i in range(len(af_nrm_n_n)):
        af_nrm_n_n[i] = af_nrm_n[i]/norm_af
    
    V['af_nrm_n_n'] = af_nrm_n_n
    V['afmagn'] = afmagn

    V['af_arm_n_n'] = af_arm_n_n

    arm_p = armn[:ifield,:cntfield] #change 4 to no fields
    #write to file at some point
    V['arm_plot'] = arm_p[2]



    #initiate figure
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize = (20,5))
    #w, h = figaspect(1)
    #fig, ax = plt.subplots(figsize=(w,h))
    #ax1.set(aspect=1)
    ax1.plot(demagstep2, arm_p[2], marker = 'o', color = 'r', label = 'simulated') #x = 22, y = 100 one plotted is largest field? - all same
    ax1.plot(demagstep2, af_arm_n_n, marker = 'o', color = 'b', label = 'measured')        
    ax1.set_ylim(0,1.3)
    ax1.set_xlim(0,np.nanmax(demagstep2))
    ax1.set_ylabel('ARM (normalized)')
    ax1.set_xlabel('AF peak (mT)')
    ax1.legend(loc='upper right')
    ax1.set_title('ARM demagnetization spectra for sample {}'.format(name))
    ax1.grid()

    # get rest ready 
    #afmag = V['afmag'] #simulated af demag
    #measured data
    #af_nrm_n = X['af_nrm'] 
    af_arm_n = X['af_arm']
    name = X['name']
    fields = V['fields']
    af_step = X['af_step']


    ###################
    afpick = X['af_pick'] #6 #pick from user - Z-plot

    #######################
    sit = ifield

    ys = np.zeros((100))

    cntfieldaf = cntfield

    af_step_list = []
    shortdiv = np.zeros((cntfield))
    shortminus = np.zeros((cntfield))

    res = []
    flat = 0 #set here as well as earlier to save time when testing
    afratio = 0
    #do the shortdiv etc seperate

    for i in range(cntfieldaf): #cntfieldsaf
        #likely remove this -----
        xlist = []
        ylist = []
        sumx = 0.
        sumy = 0.
        sumxy = 0.
        sumxx = 0.
        print('afmag', afmag)
        print('sit', sit) #sit is field
        print('i', i)
        print('arm', arm)
        af_arm = afmag[:sit,i]/arm[:sit,i] #changed to divide by max value
        field_t = fields[:sit]*mu0*1E6
        af_arm_new = np.zeros((int(sit+1)))
        field_t_new = np.zeros((int(sit+1)))
        af_arm_new[1:] = af_arm
        field_t_new[1:] = field_t
        print(af_arm_new)
        print(field_t_new)
        print('af_arm_new', af_arm_new)
        print('field_t_new', field_t_new)
        #p,residuals, _, _, _ = np.polyfit((afmag[:sit,i]/sirm[:sit,i]), (fields[:sit]*mu0*1E6), 1, full = True) # not need loop- all fields for each point  
        p,residuals, _, _, _ = np.polyfit((af_arm_new), (field_t_new), 1, full = True)
        mfit = p[0]
        cfit = p[1]
    
        res.append(residuals)

        # ----------- to here 
        x = af_nrm_n[i]/af_arm_n[i]   #change this to norm by max measured value
        y = x*mfit +cfit #field assoacted witht his step for demag  

        ys[i] = y

        j = sit-1
        #include these divides for plotting at each point - rewrite as arrays 
        shortdiv[i]=(1-((arm[j-1,i]/np.mean(arm[j-1,0]))/(af_arm_n[i]/np.mean(af_arm_n[0]))))*100 #do for 3rd field
        shortminus[i]=(((arm[j-1,i]/np.mean(arm[j-1,0]))-(af_arm_n[i]/np.mean(af_arm_n[0]))))*100
        af_step_list.append(af_step[i])

    #do 2 plots and finish figure and write the data to files then can add in rest

    #middle fig
    af_step_list = X['af_step']
    
    twenty = []
    hundred = []

    for i in range(len(af_step_list)):
        twenty.append(20)
        hundred.append(100)

    ax2.plot(af_step_list, twenty, 'b') #af step list is just af_step?
    ax2.plot(af_step_list, hundred, 'r')
    ax2.set_ylim(0,130)
    ax2.plot([af_step_list[0], af_step_list[-1]], [-20, -20], 'b') #af step list is just af_step?
    ax2.plot([af_step_list[0], af_step_list[-1]], [-100, -100], 'r')
    ax2.plot(af_step_list, shortdiv,  marker='o', color= 'r', label='S$_{ratio}$')
    ax2.plot(af_step_list, shortminus,  marker='o', color = 'b', label = 'S$_{diff}$')
    ax2.set_title('ARM checks for sample {}'.format(name))
    ax2.set_xlabel('AF peak (mT)')
    ax2.set_ylabel('S$_{diff}$ or S$_{ratio}$ (%)')
    ax2.set_xlim(0,np.nanmax(af_step_list))
    #ax2.set_ylim(min(np.min(shortdiv), np.min(shortminus), -100),max(np.max(shortdiv), np.max(shortminus), 100))
    ax2.legend(loc='upper right')
    ax2.grid()

    #plot PI estimates - mark where SIRM checks fit + in file write the values and also sigm for each point

    ax3.plot(af_step_list,ys[:len(af_step_list)], 'b', label = 'All')
    ax3.plot(af_step_list,ys[:len(af_step_list)],  marker='o', color= 'b')
    #plot ones which are accepted by SIRM checks in diff colour and after af pick here
    #ax3.plot(af_step_list,ys[:len(af_step_list)], 'b', label = 'Pass')
    #ax3.plot(af_step_list,ys[:len(af_step_list)],  marker='o', color= 'b')
    #plot vertical line
    
    ax3.plot([af_step_list[afpick], af_step_list[afpick]], [0,np.max(ys[:len(af_step)])], color = 'green')
    #plt.plot(af_step[8:20],ys[8:20],  marker='o', color= 'r')
    ax3.set_xlim(0,np.nanmax(af_step_list))
    #ax3.set_ylim(0,np.nanmax(ys[:len(af_step)]))
    #plt.text(20, 6, r'selected median: %.2f $\mu$T'%selected_med, fontsize=11)
    #plt.text(20, 7, r'rejected mean: %.2f $\pm$ %.2f $\mu$T'%(selected_mean ,mean_dev), fontsize=11)
    #plt.text(20, 6, r'selected median: %.2f $\mu$ T'%selected_med, fontsize=12)
    ax3.set_xlabel('AF peak (mT)')
    ax3.set_ylabel('paleointensity (\u03BCT)')
    ax3.grid()
    ax3.set_title('PI for each AF step for sample {}'.format(name))

    ax4.plot(demagstep2, afmagn[0], marker = 'o', color = 'r', label = 'simulated CRM 1 ') #x = 22, y = 100 one plotted is largest field? - all same
    ax4.plot(demagstep2, afmagn[1], marker = 'o', color = 'r', label = 'simulated CRM 2 ') #x = 22, y = 100 one plotted is largest field? - all same
    ax4.plot(demagstep2, afmagn[2], marker = 'o', color = 'r', label = 'simulated CRM 3 ') #x = 22, y = 100 one plotted is largest field? - all same
    ax4.plot(demagstep2, afmagn[3], marker = 'o', color = 'r', label = 'simulated CRM 4 ') #x = 22, y = 100 one plotted is largest field? - all same
    ax4.plot(demagstep2, af_nrm_n_n, marker = 'o', color = 'b', label = 'measured NRM')        
    ax4.set_ylim(0,1.3)
    ax4.set_xlim(0,np.nanmax(demagstep2))
    ax4.set_ylabel('NRM or CRM (normalized to max NRM or CRM)')
    ax4.set_xlabel('AF peak (mT)')
    ax4.legend(loc='upper right')
    ax4.set_title('NRM and CRM demagnetization spectra for sample {}'.format(name))
    ax4.grid()

    # get rest ready 
    afmag = V['afmag'] #simulated af demag
    #measured data
    af_nrm_n = X['af_nrm'] 
    af_arm_n = X['af_arm']
    name = X['name']
    fields = V['fields']
    af_step = X['af_step']


    plt.show()
    #------------ to here 

    #write to the file 

    #write lots of data files 
    #make directory for folder and save all these ine tc 

    path = os.getcwd()
    fall = open(path + os.sep+'sample_{}_V_{}'.format(name,(X['sample_copy'])) + os.sep+'full_data_sample_{}_run_{}.dat'.format(name,(X['sample_copy'])), 'w')
    fall.write('Sample \t Run \t AF step \t PI est \SRR \t meaured NRM/SIRM \t meausured SIRM/SIRM_) = \tsimulated SIRM/SIRM_0 \t SIRM ratio \t SIRM diff')
    for i in range(cntfield):
        fall.write(str(name) + '\t' + str(X['sample_copy']) + '\t' + str(af_step_list[i]) + '\t' + str(ys[i]) + '\t' + str(res[i]) + '\t' + str(af_nrm_n[i]/af_sirm_n[i]) + '\t' + str(af_sirm_n[i]) + '/t' + str(sirm_p[2]) + '\t' + str(shortminus[i]) + '\t' + str(shortdiv[i]) + '\n')
    fall.close()

    V['shortdiv'] = shortdiv
    V['shortminus'] = shortminus
    V['ys'] = ys  

    return(X,V)

def fin_pal(X,V):
    ys = V['ys']
    af_step = X['af_step']
    name = X['name']
    cntfield = X['cntfield']
    
    while True:
            low_b_val = (input("Pick the AF step for the lower bound of the platau palaeointensity region to calcualte the median palaeointensity from:"))
            
            #print(type(afpick))
            try:
                low_b1 = np.where(af_step == float(low_b_val))
                low_b = int(low_b1[0])
                if (low_b >= 0) and (low_b <= cntfield): #within range of AF demag steps - may break if pick too high
                    #print('in bounds')
                    break
            except:
                print('Not an AF step')
                True
           # if (isinstance(sf_choose, int)):
                #print('int')
            
    while True:
            up_b_val = (input("Pick the AF step for the upper bound of the platau palaeointensity region to calcualte the median palaeointensity from:"))
            
            try:
                up_b1 = np.where(af_step == float(up_b_val))
                up_b = int(up_b1[0])
                if (up_b >= low_b) and (up_b <= cntfield): #within range of AF demag steps - may break if pick too high
                    
                    break
                else:
                    print('out of bounds, must be above the lower bound')
            except :
                print('Not an AF step')
                True
           # if (isinstance(sf_choose, int)):
                #print('int')
            
            #plot restulting  same graoh with labelled in red and wiht median etc 
            
    ys = V['ys']
    af_step = X['af_step']

    #plot final figure 
 
    name = X['name']

    #initiate figure
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize = (20,5))
    #w, h = figaspect(1)
    #fig, ax = plt.subplots(figsize=(w,h))
    #ax1.set(aspect=1)
    ax1.plot(af_step, V['arm_plot'], marker = 'o', color = 'r', label = 'simulated') #x = 22, y = 100 one plotted is largest field? - all same
    ax1.plot(af_step, V['af_arm_n_n'], marker = 'o', color = 'b', label = 'measured')     
    ax1.set_ylim(0,1.3)
    ax1.set_xlim(0,np.nanmax(af_step))
    ax1.set_ylabel('ARM (normalized)')
    ax1.set_xlabel('AF peak (mT)')
    ax1.legend(loc='upper right')
    ax1.set_title('ARM demagnetization spectra for sample {}'.format(name))
    ax1.grid()

    # get rest ready 
    afmag = V['afmag'] #simulated af demag
    #measured data
    af_nrm_n = X['af_nrm'] 
    af_arm_n = X['af_arm']
    name = X['name']
    fields = V['fields']
    af_step = X['af_step']


    ###################
    afpick = X['af_pick'] #6 #pick from user - Z-plot

    
    twenty = []
    hundred = []

    for i in range(len(af_step)):
        twenty.append(20)
        hundred.append(100)


    ax2.plot(af_step, twenty, 'b') #af step list is just af_step?
    ax2.plot(af_step, hundred, 'r')
    ax2.set_ylim(-110,130)
    ax2.plot([af_step[0], af_step[-1]], [-20, -20], 'b') #af step list is just af_step?
    ax2.plot([af_step[0], af_step[-1]], [-100, -100], 'r')
    ax2.plot(af_step, V['shortdiv'],  marker='o', color= 'r', label='S$_{ratio}$')
    ax2.plot(af_step, V['shortminus'],  marker='o', color = 'b', label = 'S$_{diff}$')
    ax2.set_title('ARM checks for sample {}'.format(name))
    ax2.set_xlabel('AF peak (mT)')
    ax2.set_ylabel('S$_{diff}$ or S$_{ratio}$ (%)')
    ax2.set_xlim(0,np.nanmax(af_step))
    #ax2.set_ylim(min(np.min(shortdiv), np.min(shortminus), -100),max(np.max(shortdiv), np.max(shortminus), 100))
    ax2.legend(loc='upper right')
    ax2.grid()



    #plot PI estimates - mark where SIRM checks fit + in file write the values and also sigm for each point

    ax3.plot(af_step,ys[:len(af_step)], 'b', label = 'All')
    ax3.plot(af_step,ys[:len(af_step)],  marker='o', color= 'b')
    #plot ones which are accepted by SIRM checks in diff colour and after af pick here
    #ax3.plot(af_step_list,ys[:len(af_step_list)], 'b', label = 'Pass')
    #ax3.plot(af_step_list,ys[:len(af_step_list)],  marker='o', color= 'b')
    #plot vertical line
    
    ax3.set_xlim(0,np.nanmax(af_step))
    ax3.set_ylim(0,(np.max(ys) + 0.1*np.max(ys)))
    #plot with red colour and calc average
    selected_mean = np.mean(ys[low_b:up_b+1]) #check these include these values
    mean_dev = np.std(ys[low_b:up_b+1])
    selected_med = np.median(ys[low_b:up_b+1])
    q3, q1 = np.percentile(ys[low_b:up_b+1], [75, 25])
    iqr = q3-q1
    ax3.plot(af_step,ys[:len(af_step)], 'b')
    ax3.plot(af_step,ys[:len(af_step)],  marker='o', color= 'b')
    ax3.plot(af_step[low_b:up_b+1],ys[low_b:up_b+1],  marker='o', color= 'r')
    ax3.plot([af_step[afpick], af_step[afpick]], [0,np.max(ys[:len(af_step)])], color = 'green')   

    ax3.text(max(af_step)/2, -(0.18*np.max(ys)), r'median: %.2f $\mu$T'%selected_med)
    ax3.text(max(af_step)/2, -(0.25*np.max(ys)), r'mean: %.2f $\pm$ %.2f $\mu$T'%(selected_mean ,mean_dev))
    ax3.grid()
    ax3.set_title('PI for each AF step for sample {}'.format(name))
       
    ax3.set_xlabel('AF peak (mT)')
    ax3.set_ylabel('paleointensity (\u03BCT)')


    af_nrm_n_n = V['af_nrm_n_n']
    agmagn = V['afmagn']

    #add ax4
    ax4.plot(af_step, V['afmagn'][0], marker = 'o', color = 'r', label = 'simulated CRM 1 ') #x = 22, y = 100 one plotted is largest field? - all same
    ax4.plot(af_step, V['afmagn'][1], marker = 'o', color = 'r', label = 'simulated CRM 2 ') #x = 22, y = 100 one plotted is largest field? - all same
    ax4.plot(af_step, V['afmagn'][2], marker = 'o', color = 'r', label = 'simulated CRM 3 ') #x = 22, y = 100 one plotted is largest field? - all same
    ax4.plot(af_step, V['afmagn'][3], marker = 'o', color = 'r', label = 'simulated CRM 4 ') #x = 22, y = 100 one plotted is largest field? - all same
    ax4.plot(af_step, af_nrm_n_n, marker = 'o', color = 'b', label = 'measured NRM')        
    ax4.set_ylim(0,1.3)
    ax4.set_xlim(0,np.nanmax(af_step))
    ax4.set_ylabel('NRM or CRM (normalized to max NRM or CRM)')
    ax4.set_xlabel('AF peak (mT)')
    ax4.legend(loc='upper right')
    ax4.set_title('NRM and CRM demagnetization spectra for sample {}'.format(name))
    ax4.grid()


    plt.suptitle('Restuls for sample {}'.format(name))
    #print('Figure saved as PI_est_{0}.pdf'.format(name))
    path = os.getcwd()
    plt.savefig(path + os.sep+'sample_{}_V_{}'.format(name, X['sample_copy']) + os.sep+'PI_results_{}.pdf'.format(name))
    plt.show

    plt.pause(1)
    #open big file
    bigf = open('paleo_results.dat', 'a')
    bigf.write(str(name) + '\t' + str(X['sample_copy']) + '\t' + str(up_b-low_b) + '\t' + str(af_step[up_b]-af_step[low_b]) +'\t' + str(af_step[low_b]) + '\t' + str(af_step[up_b]) + '\t' + '{:.6f}'.format(selected_mean) + '\t' + '{:.6f}'.format(mean_dev) + '\t' + '{:.6f}'.format(selected_med) + '\t' + '{:.6f}'.format(iqr) + '\t' + str(X['min_field']) + '\t' + str(X['max_field']) + '\t' + str(X['growth_rate']) + '\t' + str(X['curie_t']) + '\t' +  str(X['afval']) + '\t' + str(X['SF']) + '\t' + str(X['reset_limit_hc']) + '\t' + str(X['reset_limit_hi']) + '\n')
    bigf.close()
    return

####### FOR SIRM ###############

def calc_PI_checks_SIRM(V,X):
    mu0 = 4.*pi*1e-7
    sirm = V['sirm']
    print('sirm', V['sirm'])

    #arm = V['arm']
    #print('arm', arm)
    cntfield = X['cntfield']
    name = X['name']
    ifield = V['ifield']
    demagstep = X['af_step']
    #set first value as v v low just so remove zeros to plot easily
    demagstep2 = demagstep
    demagstep2[0] = 0.0001
    demagstep2 = demagstep2[demagstep2 != 0]
    sirm2 = sirm
    sirmn = np.copy(sirm2) #this is normalising the sirm demag
    for i in range(V['ifield']):
        for j in range(cntfield): #change back to 23 
            sirmn[i][j] = (sirm2[i][j]/(np.mean(sirm2[i][0:3]))) #norm by 3

    af_sirm_n = X['af_sirm']
    norm = np.mean(af_sirm_n[0:3])
    af_sirm_n_n = np.copy(af_sirm_n)
    for i in range(len(af_sirm_n)):
        af_sirm_n_n[i] = af_sirm_n[i]/norm

    V['af_sirm_n_n'] = af_sirm_n_n

    sirm_p = sirmn[:ifield,:cntfield] #change 4 to no fields
    #write to file at some point
    V['sirm_plot'] = sirm_p[2]

    #normalise ARM demag

    #arm2 = arm
    #armn = np.copy(arm2) #this is normalising the sirm demag
    #for i in range(V['ifield']):
    #    for j in range(cntfield): #change back to 23 
    #        armn[i][j] = (arm2[i][j]/(np.mean(arm2[i][0])))

    #af_arm_n = X['af_arm']
    #norm = np.mean(af_arm_n[0])
    #af_arm_n_n = np.copy(af_arm_n)
    #for i in range(len(af_arm_n)):
    #    af_arm_n_n[i] = af_arm_n[i]/norm

    #V['af_arm_n_n'] = af_arm_n_n

    #arm_p = armn[:ifield,:cntfield] #change 4 to no fields
    #write to file at some point
    #V['arm_plot'] = arm_p[2]
    ###################
    afpick = X['af_pick'] #6 #pick from user - Z-plot

    #prepare CRM demag 
    afmag = V['afmag'] #simulated af demag
    afmag2 = afmag
    afmagn = np.copy(afmag2)
    #norm_afmag = afmag2[i][afpick]
    for i in range(V['ifield']):
        for j in range(cntfield): #change back to 23 
            afmagn[i][j] = (afmag2[i][j]/(afmag2[i][afpick])) # norm to max not first 3 

    af_nrm_n = X['af_nrm']
    norm_af = af_nrm_n[afpick]
    af_nrm_n_n = np.copy(af_nrm_n)
    for i in range(len(af_nrm_n_n)):
        af_nrm_n_n[i] = af_nrm_n[i]/norm_af
    
    V['af_nrm_n_n'] = af_nrm_n_n
    V['afmagn'] = afmagn

    #V['af_arm_n_n'] = af_arm_n_n

    #arm_p = armn[:ifield,:cntfield] #change 4 to no fields
    #write to file at some point
    #V['arm_plot'] = arm_p[2]



    #initiate figure
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize = (20,5))
    #w, h = figaspect(1)
    #fig, ax = plt.subplots(figsize=(w,h))
    #ax1.set(aspect=1)
    ax1.plot(demagstep2, sirm_p[2], marker = 'o', color = 'r', label = 'simulated') #x = 22, y = 100 one plotted is largest field? - all same
    ax1.plot(demagstep2, af_sirm_n_n, marker = 'o', color = 'b', label = 'measured')        
    ax1.set_ylim(0,1.3)
    ax1.set_xlim(0,np.nanmax(demagstep2))
    ax1.set_ylabel('SIRM (normalised)')
    ax1.set_xlabel('AF peak (mT)')
    ax1.legend(loc='upper right')
    ax1.set_title('SIRM demagnetization spectra for sample {}'.format(name))
    ax1.grid()

    # get rest ready 
    #afmag = V['afmag'] #simulated af demag
    #measured data
    #af_nrm_n = X['af_nrm'] 
    #af_arm_n = X['af_arm']
    name = X['name']
    fields = V['fields']
    af_step = X['af_step']




    #######################
    sit = ifield

    ys = np.zeros((100))

    cntfieldaf = cntfield

    af_step_list = []
    shortdiv = np.zeros((cntfield))
    shortminus = np.zeros((cntfield))

    res = []
    flat = 0 #set here as well as earlier to save time when testing
    afratio = 0
    #do the shortdiv etc seperate

    for i in range(cntfieldaf): #cntfieldsaf
        #likely remove this -----
        xlist = []
        ylist = []
        sumx = 0.
        sumy = 0.
        sumxy = 0.
        sumxx = 0.
        #print('afmag', afmag)
        #print('sit', sit) #sit is field
        #print('i', i)
        #print('arm', arm)
        af_sirm = afmag[:sit,i]/sirm[:sit,i] #changed to divide by max value
        field_t = fields[:sit]*mu0*1E6
        af_sirm_new = np.zeros((int(sit+1)))
        field_t_new = np.zeros((int(sit+1)))
        af_sirm_new[1:] = af_sirm
        field_t_new[1:] = field_t
        #print(af_arm_new)
        #print(field_t_new)
        #print('af_arm_new', af_arm_new)
        #print('field_t_new', field_t_new)
        #p,residuals, _, _, _ = np.polyfit((afmag[:sit,i]/sirm[:sit,i]), (fields[:sit]*mu0*1E6), 1, full = True) # not need loop- all fields for each point  
        p,residuals, _, _, _ = np.polyfit((af_sirm_new), (field_t_new), 1, full = True)
        mfit = p[0]
        cfit = p[1]
    
        res.append(residuals)

        # ----------- to here 
        x = af_nrm_n[i]/af_sirm_n[i]   #change this to norm by max measured value
        y = x*mfit +cfit #field assoacted witht his step for demag  

        ys[i] = y

        j = sit-1
        #include these divides for plotting at each point - rewrite as arrays 
        shortdiv[i]=(1-((sirm[j-1,i]/np.mean(sirm[j-1,0]))/(af_sirm_n[i]/np.mean(af_sirm_n[0]))))*100 #do for 3rd field
        shortminus[i]=(((sirm[j-1,i]/np.mean(sirm[j-1,0]))-(af_sirm_n[i]/np.mean(af_sirm_n[0]))))*100
        af_step_list.append(af_step[i])

    #do 2 plots and finish figure and write the data to files then can add in rest

    #middle fig
    af_step_list = X['af_step']
    
    twenty = []
    hundred = []

    for i in range(len(af_step_list)):
        twenty.append(20)
        hundred.append(100)

    ax2.plot(af_step_list, twenty, 'b') #af step list is just af_step?
    ax2.plot(af_step_list, hundred, 'r')
    ax2.set_ylim(-130,130)
    ax2.plot([af_step_list[0], af_step_list[-1]], [-20, -20], 'b') #af step list is just af_step?
    ax2.plot([af_step_list[0], af_step_list[-1]], [-100, -100], 'r')
    ax2.plot(af_step_list, shortdiv,  marker='o', color= 'r', label='S$_{ratio}$')
    ax2.plot(af_step_list, shortminus,  marker='o', color = 'b', label = 'S$_{diff}$')
    ax2.set_title('SIRM checks for sample {}'.format(name))
    ax2.set_xlabel('AF peak (mT)')
    ax2.set_ylabel('S$_{diff}$ or S$_{ratio}$ (%)')
    ax2.set_xlim(0,np.nanmax(af_step_list))
    #ax2.set_ylim(min(np.min(shortdiv), np.min(shortminus), -100),max(np.max(shortdiv), np.max(shortminus), 100))
    ax2.legend(loc='upper right')
    ax2.grid()

    #plot PI estimates - mark where SIRM checks fit + in file write the values and also sigm for each point

    ax3.plot(af_step_list,ys[:len(af_step_list)], 'b', label = 'All')
    ax3.plot(af_step_list,ys[:len(af_step_list)],  marker='o', color= 'b')
    #plot ones which are accepted by SIRM checks in diff colour and after af pick here
    #ax3.plot(af_step_list,ys[:len(af_step_list)], 'b', label = 'Pass')
    #ax3.plot(af_step_list,ys[:len(af_step_list)],  marker='o', color= 'b')
    #plot vertical line
    
    ax3.plot([af_step_list[afpick], af_step_list[afpick]], [0,np.max(ys[:len(af_step)])], color = 'green')
    #plt.plot(af_step[8:20],ys[8:20],  marker='o', color= 'r')
    ax3.set_xlim(0,np.nanmax(af_step_list))
    #ax3.set_ylim(0,np.nanmax(ys[:len(af_step)]))
    #plt.text(20, 6, r'selected median: %.2f $\mu$T'%selected_med, fontsize=11)
    #plt.text(20, 7, r'rejected mean: %.2f $\pm$ %.2f $\mu$T'%(selected_mean ,mean_dev), fontsize=11)
    #plt.text(20, 6, r'selected median: %.2f $\mu$ T'%selected_med, fontsize=12)
    ax3.set_xlabel('AF peak (mT)')
    ax3.set_ylabel('paleointensity (\u03BCT)')
    ax3.grid()
    ax3.set_title('PI for each AF step for sample {}'.format(name))

    ax4.plot(demagstep2, afmagn[0], marker = 'o', color = 'r', label = 'simulated CRM 1 ') #x = 22, y = 100 one plotted is largest field? - all same
    ax4.plot(demagstep2, afmagn[1], marker = 'o', color = 'r', label = 'simulated CRM 2 ') #x = 22, y = 100 one plotted is largest field? - all same
    ax4.plot(demagstep2, afmagn[2], marker = 'o', color = 'r', label = 'simulated CRM 3 ') #x = 22, y = 100 one plotted is largest field? - all same
    ax4.plot(demagstep2, afmagn[3], marker = 'o', color = 'r', label = 'simulated CRM 4 ') #x = 22, y = 100 one plotted is largest field? - all same
    ax4.plot(demagstep2, af_nrm_n_n, marker = 'o', color = 'b', label = 'measured NRM')        
    ax4.set_ylim(0,1.3)
    ax4.set_xlim(0,np.nanmax(demagstep2))
    ax4.set_ylabel('NRM or CRM (normalized to max NRM or CRM)')
    ax4.set_xlabel('AF peak (mT)')
    ax4.legend(loc='upper right')
    ax4.set_title('NRM and CRM demagnetization spectra for sample {}'.format(name))
    ax4.grid()

    # get rest ready 
    afmag = V['afmag'] #simulated af demag
    #measured data
    af_nrm_n = X['af_nrm'] 
    #af_arm_n = X['af_arm']
    name = X['name']
    fields = V['fields']
    af_step = X['af_step']


    plt.show()
    #------------ to here 

    #write to the file 

    #write lots of data files 
    #make directory for folder and save all these ine tc 

    #path = os.getcwd()
    #fall = open(path + os.sep+'sample_{}_V_{}'.format(name,(X['sample_copy'])) + os.sep+'full_data_sample_{}_run_{}.dat'.format(name,(X['sample_copy'])), 'w')
    #fall.write('Sample \t Run \t AF step \t PI est \SRR \t meaured NRM/SIRM \t meausured SIRM/SIRM_) = \tsimulated SIRM/SIRM_0 \t SIRM ratio \t SIRM diff')
    #for i in range(cntfield):
    #    fall.write(str(name) + '\t' + str(X['sample_copy']) + '\t' + str(af_step_list[i]) + '\t' + str(ys[i]) + '\t' + str(res[i]) + '\t' + str(af_nrm_n[i]/af_sirm_n[i]) + '\t' + str(af_sirm_n[i]) + '/t' + str(sirm_p[2]) + '\t' + str(shortminus[i]) + '\t' + str(shortdiv[i]) + '\n')
    #fall.close()

    V['shortdiv'] = shortdiv
    V['shortminus'] = shortminus
    V['ys'] = ys  

    return(X,V)

def fin_pal_SIRM(X,V):
    ys = V['ys']
    af_step = X['af_step']
    name = X['name']
    cntfield = X['cntfield']
    
    while True:
            low_b_val = (input("Pick the AF step for the lower bound of the platau palaeointensity region to calcualte the median palaeointensity from:"))
            
            #print(type(afpick))
            try:
                low_b1 = np.where(af_step == float(low_b_val))
                low_b = int(low_b1[0])
                if (low_b >= 0) and (low_b <= cntfield): #within range of AF demag steps - may break if pick too high
                    #print('in bounds')
                    break
            except:
                print('Not an AF step')
                True
           # if (isinstance(sf_choose, int)):
                #print('int')
            
    while True:
            up_b_val = (input("Pick the AF step for the upper bound of the platau palaeointensity region to calcualte the median palaeointensity from:"))
            
            try:
                up_b1 = np.where(af_step == float(up_b_val))
                up_b = int(up_b1[0])
                if (up_b >= low_b) and (up_b <= cntfield): #within range of AF demag steps - may break if pick too high
                    
                    break
                else:
                    print('out of bounds, must be above the lower bound')
            except :
                print('Not an AF step')
                True
           # if (isinstance(sf_choose, int)):
                #print('int')
            
            #plot restulting  same graoh with labelled in red and wiht median etc 
            
    ys = V['ys']
    af_step = X['af_step']

    #plot final figure 
 
    name = X['name']

    #initiate figure
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize = (20,5))
    #w, h = figaspect(1)
    #fig, ax = plt.subplots(figsize=(w,h))
    #ax1.set(aspect=1)
    ax1.plot(af_step, V['sirm_plot'], marker = 'o', color = 'r', label = 'simulated') #x = 22, y = 100 one plotted is largest field? - all same
    ax1.plot(af_step, V['af_sirm_n_n'], marker = 'o', color = 'b', label = 'measured')     
    ax1.set_ylim(0,1.3)
    ax1.set_xlim(0,np.nanmax(af_step))
    ax1.set_ylabel('SIRM (normalised)')
    ax1.set_xlabel('AF peak (mT)')
    ax1.legend(loc='upper right')
    ax1.set_title('SIRM demagnetization spectra for sample {}'.format(name))
    ax1.grid()

    # get rest ready 
    afmag = V['afmag'] #simulated af demag
    #measured data
    af_nrm_n = X['af_nrm'] 
    af_sirm_n = X['af_sirm']
    name = X['name']
    fields = V['fields']
    af_step = X['af_step']


    ###################
    afpick = X['af_pick'] #6 #pick from user - Z-plot

    
    twenty = []
    hundred = []

    for i in range(len(af_step)):
        twenty.append(20)
        hundred.append(100)


    ax2.plot(af_step, twenty, 'b') #af step list is just af_step?
    ax2.plot(af_step, hundred, 'r')
    ax2.set_ylim(-130,130)
    ax2.plot([af_step[0], af_step[-1]], [-20, -20], 'b') #af step list is just af_step?
    ax2.plot([af_step[0], af_step[-1]], [-100, -100], 'r')
    ax2.plot(af_step, V['shortdiv'],  marker='o', color= 'r', label='S$_{ratio}$')
    ax2.plot(af_step, V['shortminus'],  marker='o', color = 'b', label = 'S$_{diff}$')
    ax2.set_title('SIRM checks for sample {}'.format(name))
    ax2.set_xlabel('AF peak (mT)')
    ax2.set_ylabel('S$_{diff}$ or S$_{ratio}$ (%)')
    ax2.set_xlim(0,np.nanmax(af_step))
    #ax2.set_ylim(min(np.min(shortdiv), np.min(shortminus), -100),max(np.max(shortdiv), np.max(shortminus), 100))
    ax2.legend(loc='upper right')
    ax2.grid()



    #plot PI estimates - mark where SIRM checks fit + in file write the values and also sigm for each point

    ax3.plot(af_step,ys[:len(af_step)], 'b', label = 'All')
    ax3.plot(af_step,ys[:len(af_step)],  marker='o', color= 'b')
    #plot ones which are accepted by SIRM checks in diff colour and after af pick here
    #ax3.plot(af_step_list,ys[:len(af_step_list)], 'b', label = 'Pass')
    #ax3.plot(af_step_list,ys[:len(af_step_list)],  marker='o', color= 'b')
    #plot vertical line
    
    ax3.set_xlim(0,np.nanmax(af_step))
    ax3.set_ylim(0,(np.max(ys) + 0.1*np.max(ys)))
    #plot with red colour and calc average
    selected_mean = np.mean(ys[low_b:up_b+1]) #check these include these values
    mean_dev = np.std(ys[low_b:up_b+1])
    selected_med = np.median(ys[low_b:up_b+1])
    q3, q1 = np.percentile(ys[low_b:up_b+1], [75, 25])
    iqr = q3-q1
    ax3.plot(af_step,ys[:len(af_step)], 'b')
    ax3.plot(af_step,ys[:len(af_step)],  marker='o', color= 'b')
    ax3.plot(af_step[low_b:up_b+1],ys[low_b:up_b+1],  marker='o', color= 'r')
    ax3.plot([af_step[afpick], af_step[afpick]], [0,np.max(ys[:len(af_step)])], color = 'green')   

    ax3.text(max(af_step)/2, -(0.18*np.max(ys)), r'median: %.2f $\mu$T'%selected_med)
    ax3.text(max(af_step)/2, -(0.25*np.max(ys)), r'mean: %.2f $\pm$ %.2f $\mu$T'%(selected_mean ,mean_dev))
    ax3.grid()
    ax3.set_title('PI for each AF step for sample {}'.format(name))
       
    ax3.set_xlabel('AF peak (mT)')
    ax3.set_ylabel('paleointensity (\u03BCT)')


    af_nrm_n_n = V['af_nrm_n_n']
    agmagn = V['afmagn']

    #add ax4
    ax4.plot(af_step, V['afmagn'][0], marker = 'o', color = 'r', label = 'simulated CRM 1 ') #x = 22, y = 100 one plotted is largest field? - all same
    ax4.plot(af_step, V['afmagn'][1], marker = 'o', color = 'r', label = 'simulated CRM 2 ') #x = 22, y = 100 one plotted is largest field? - all same
    ax4.plot(af_step, V['afmagn'][2], marker = 'o', color = 'r', label = 'simulated CRM 3 ') #x = 22, y = 100 one plotted is largest field? - all same
    ax4.plot(af_step, V['afmagn'][3], marker = 'o', color = 'r', label = 'simulated CRM 4 ') #x = 22, y = 100 one plotted is largest field? - all same
    ax4.plot(af_step, af_nrm_n_n, marker = 'o', color = 'b', label = 'measured NRM')        
    #ax4.set_ylim(0,1.3)
    ax4.set_xlim(0,np.nanmax(af_step))
    ax4.set_ylabel('NRM or CRM (normalized to max NRM or CRM)')
    ax4.set_xlabel('AF peak (mT)')
    ax4.legend(loc='upper right')
    ax4.set_title('NRM and CRM demagnetization spectra for sample {}'.format(name))
    ax4.grid()


    plt.suptitle('PI results for sample {}'.format(name))
    #print('Figure saved as PI_est_{0}.pdf'.format(name))
    path = os.getcwd()
    plt.savefig(path + os.sep+'sample_{}_V_{}'.format(name, X['sample_copy']) + os.sep+'PI_results_{}.pdf'.format(name))
    plt.show

    plt.pause(1)
    #open big file
    #bigf = open('paleo_results.dat', 'a')
    #bigf.write(str(name) + '\t' + str(X['sample_copy']) + '\t' + str(up_b-low_b) + '\t' + str(af_step[up_b]-af_step[low_b]) +'\t' + str(af_step[low_b]) + '\t' + str(af_step[up_b]) + '\t' + '{:.6f}'.format(selected_mean) + '\t' + '{:.6f}'.format(mean_dev) + '\t' + '{:.6f}'.format(selected_med) + '\t' + '{:.6f}'.format(iqr) + '\t' + str(X['min_field']) + '\t' + str(X['max_field']) + '\t' + str(X['growth_rate']) + '\t' + str(X['curie_t']) + '\t' +  str(X['afval']) + '\t' + str(X['SF']) + '\t' + str(X['reset_limit_hc']) + '\t' + str(X['reset_limit_hi']) + '\n')
    #bigf.close()
    return



def sirm_test(V, X):
    sirm = V['sirm']
    cntfield = X['cntfield']
    name = X['name']
    ifield = V['ifield']
    w, h = figaspect(1)
    fig, ax = plt.subplots(figsize=(w,h))
    demagstep = X['af_step']
    #set first value as v v low just so remove zeros to plot easily
    demagstep2 = demagstep
    demagstep2[0] = 0.0001

    demagstep2 = demagstep2[demagstep2 != 0]

    sirm2 = sirm
    sirmn = np.copy(sirm2)

    for i in range(V['ifield']):
        for j in range(cntfield): #change back to 23 
            sirmn[i][j] = (sirm2[i][j]/(np.mean(sirm2[i][0:3])))


    af_sirm_n = X['af_irm']

    norm = np.mean(af_sirm_n[0:3])
    

    af_sirm_n_n = np.copy(af_sirm_n)
    for i in range(len(af_sirm_n)):
        af_sirm_n_n[i] = af_sirm_n[i]/norm

    sirm_p = sirmn[:ifield,:cntfield] #change 4 to no fields

    plt.plot(demagstep2, sirm_p[2], marker = 'o', color = 'r') #x = 22, y = 100
    
    plt.plot(demagstep2, af_sirm_n_n, marker = 'o', color = 'b')        
    plt.ylabel('SIRM (normalized %)')
    plt.xlabel('af demag step')
    plt.title("Demagnetisation of the SIRM for sample '{0}'".format(name))

    plt.text(30, 0.9, r'measured ', fontsize=12) #poss change to percentage across line (max/50 etc)
    plt.plot(25,0.9, marker = 'o', color='b')
    plt.text(30, 0.8, r'simulated', fontsize=12)
    plt.plot(25, 0.8, marker = 'o', color='r')
    plt.savefig('sirm_test_{}.svg'.format(name))
    plt.show
    return
    
def calc_pal(X, V):
    cntfield = X['cntfield']
    ifield = V['ifield']
    afmag = V['afmag']
    sirm = V['sirm'] #sirm seems global so shoulld have worked anwyay
    af_nrm_n = X['af_nrm']
    af_sirm_n = X['af_irm']
    name = X['name']
    aconst = V['aconst']
    fields = V['fields']
    af_step = X['af_step']
    tempmin = V['tempmin']
    tempmax = V['tempmax']
    afnegsirm = 0.0         
    averageme = 0.0
    averageunweight = 0.0
    averagemeadj = 0.0
    averagecount = 0
    flatcount = np.zeros((10))
    noflatsub = 0
    averageflat = np.zeros((10))
    trmsirmratio = 0.0
    sigmtot = 0.0
    selrms = 0.0
    sumegli = 0

    ###################
    afpick = X['af_pick'] #6 #pick from user - Z-plot

    #######################
    sit = ifield
    std = np.zeros((100))
    ys = np.zeros((100))
    used = np.zeros((50))
    flatused = np.zeros((50))
    ystore = np.zeros((100))
    sigmstore = np.zeros((100))
    flatstore = np.zeros((10,100))
    flatdi = np.zeros((10))
    sigmflat = np.zeros((10))
    flatmedianme = np.zeros((10))
    sigmtotflat = np.zeros((10))
    cntfieldaf = cntfield
    actfield = 1 #random number but need to see what this is
    dispersion = 0
    vtwo = 0
    fortyone = 0

    shortdivlist = []
    shortminuslist = []
    ilist = []
    af_step_list = []


    flat = 0 #set here as well as earlier to save time when testing
    afratio = 0
    #do the shortdiv etc seperate

    for i in range(cntfieldaf): #cntfieldsaf

        xlist = []
        ylist = []
        sumx = 0.
        sumy = 0.
        sumxy = 0.
        sumxx = 0.

        for j in range(sit): #for each field - a
            
            
            y = fields[j]*mu0*1e6 #field time something, x, y ust values
            x = afmag[j,i]/sirm[j,i] #norm af mag to sirm - sirm after that stage of demagetisation
            #plotting points
            xlist.append(x)
            ylist.append(y)

            if (sirm[j,i] == 0):

                for i in range(afpick, averagecount): #this all moves to diff point
                    #afpick given by file called sirmms.dat
                    dispersion = dispersion + (((ystore[i] - averageme/sigmtot))**2)/sigmstore[i]
                    vtwo = vtwo + (1/sigmstore[i])**2
                    fortyone = 1

            sumx = sumx + x        
            sumxx = sumxx + x**2
            sumxy = sumxy + x*y
            sumy = sumy + y

        mfit=(((sit+000.0)*sumxy)-sumx*sumy)/(((sit+000.0)*sumxx)-sumx*sumx)
        cfit=(sumy*sumxx-sumx*sumxy)/((sit+000.0)*sumxx-sumx*sumx)

        xtest = np.linspace(min(xlist), max(xlist), 10)
        #print('xtest', xtest)
        xtest = np.array(xtest)
        ytest = np.copy(xtest)
        ytest = mfit*xtest + cfit
        sumdi = 0.
        xlist2 = []
        ylist2 = []
        dilist = []
        for j in range(sit): #x and y values set for each field again
            y= fields[j]*mu0*1e6
            x=(afmag[j,i]/sirm[j,i]) #same as above, line 1024
            di = y - (mfit*x+cfit) #see how closely x and y fit with line equation
            dilist.append(di)
            sumdi = sumdi + di**2
            xlist.append(x)
            ylist.append(y)
        sumdi = sumdi/(sit-2) 
        sigm = sit*sumdi/(sit*sumxx-sumx*sumx) #variance
        x = af_nrm_n[i]/af_sirm_n[i]   
        y = x*mfit +cfit #field assoacted witht his step for demag  
        xlist2.append(x)
        ylist2.append(y)

        ys[i] = y

        std[i] = sqrt(sigm)*x
        used[i] = 0 #set used for each cntfield to zero
        flatused[i] = 0
        
        #here
        shortdiv=abs(1-((sirm[j-1,i]/np.mean(sirm[j-1,0:3]))/(af_sirm_n[i]/np.mean(af_sirm_n[0:3]))))*100 #do for 3rd field
        shortminus=abs(((sirm[j-1,i]/np.mean(sirm[j-1,0:3]))-(af_sirm_n[i]/np.mean(af_sirm_n[0:3]))))*100
        shortdivm=abs(1-((sirm[j-2,i]/np.mean(sirm[j-2,0:3]))/(af_sirm_n[i-1]/np.mean(af_sirm_n[0:3]))))*100
        shortminusm=abs(((sirm[j-2,i]/np.mean(sirm[j-2,0:3]))-(af_sirm_n[i-1]/np.mean(af_sirm_n[0:3]))))*100
        shortdivlist.append(shortdiv)
        shortminuslist.append(shortminus)
        af_step_list.append(af_step[i])
        if (i >= afpick) and (sumx != 0.0): 

            sigm = sigm*(x**2)


            afratio = afmag[2,i]/afmag[2,0]

            if (i >= 1): #greater than 1 as 0 counts here
                afratiom = (afmag[2, i-1]/afmag[2,0]) #also do ratio for prev point - poss compares afratio and afratiom

            shortdiv=abs(1-((sirm[j-1,i]/np.mean(sirm[j-1,0:3]))/(af_sirm_n[i]/np.mean(af_sirm_n[0:3]))))*100 #do for 3rd field
            shortminus=abs(((sirm[j-1,i]/np.mean(sirm[j-1,0:3]))-(af_sirm_n[i]/np.mean(af_sirm_n[0:3]))))*100
            shortdivm=abs(1-((sirm[j-2,i]/np.mean(sirm[j-2,0:3]))/(af_sirm_n[i-1]/np.mean(af_sirm_n[0:3]))))*100
            shortminusm=abs(((sirm[j-2,i]/np.mean(sirm[j-2,0:3]))-(af_sirm_n[i-1]/np.mean(af_sirm_n[0:3]))))*100
            sigm=abs(1-sigm)

            #shortdivlist.append(shortdiv)
            #shortminuslist.append(shortminus)
            ilist.append(i)
            #af_step_list.append(af_step[i])
            if (i > 0) and (i < 30): #30 seems arbituary

                if (af_nrm_n[i]/af_nrm_n[0] > 0.01) and (afratio > 0.01): #stage requires af ratio been given value but not needed

                    if (y > 0.0):
                        #print('y>0', y)

                        if (shortdiv < 100): #acceptable ranges of S - seen in data displayed
                            if (shortminus < 20):
                                #print('conds are met')
                                selrms=selrms+abs((af_sirm_n[i]/af_sirm_n[0])- (sirm[j-1,i]/sirm[j-1,0]))               

                                sigmtot=sigmtot+(1/sigm) # sum of 1/variance

                                sumegli=sumegli+(y/actfield) #sum of fields used in predied

                                averageme=averageme+(y)/sigm 

                                averageunweight=averageunweight+y

                                averagecount=averagecount+1

                                ystore[averagecount]=y #not want averagecount to skip points?

                                sigmstore[averagecount]=sigm

                                used[i]=1 #array, this means use tihs point in calc? therefore add to used array
                                trmsirmratio=trmsirmratio+x #x is true trm/sirm measred ratio

                                if (i > 1):
                                    flatdiff = abs(y - ys[i-1])/max(y,ys[i-1])

                                    if (flatdiff < 0.2) or (y < 3.0 and ys[i-1] < 3.0): #if flatdiff < 0.2 or both <3

                                        if (noflatsub == 0): # move onto new section

                                            flat = flat +1

                                            if (i-1 < afpick) or (shortdivm > 100) or (shortminusm > 20) or (ys[i-1] < 0.0) or (af_nrm_n[i-1]/af_nrm_n[0] < 0.01) or (afratiom < 0.01):
                                                #print('rejecting prev point as out of bounds or not close enough, uses lookin 2 points back?')
                                                pass
                                            else:

                                                flatcount[flat] = flatcount[flat] + 1

                                                flatstore[flat][int(flatcount[flat])] = ys[i-1] #index erroe - change to int

                                                averageflat[flat] = averageflat[flat] + ys[i-1] #page 1121

                                                sigmflat[flat] = sigmflat[flat] + (1/sigm) #sum of 1/variance

                                                flatused[i-1] = flat

                                                flatdi[flat] = flatdi[flat]+(y-ys[i-1])/max(y,ys[i-1])

                                        noflatsub = 1

                                        flatcount[flat] = flatcount[flat] + 1

                                        flatstore[flat][int(flatcount[flat])] = y

                                        averageflat[flat] = averageflat[flat] + y #add one this point to list, some overwrite what above some diff
                                        sigmflat[flat] = sigmflat[flat] + (1/sigm) # sum of 1/variance
                                        flatdi[flat] = flatdi[flat] + (y - ys[i-1])/max(y,ys[i-1])
                                        flatused[i] = flat
                                    else: #one indent because one end if ending the else which not needed in python

                                        noflatsub = 0

                                else:
                                    noflatsub = 0

                            else:
                                noflatsub = 0


                        else:
                            noflatsub = 0

                    else:
                        noflatsub =0

                else:
                    noflatsub = 0


            else:
                noflatsub = 0

        else: 
            noflatsub = 0
    if (averagecount != 0):
        selrms=100*(selrms)/(1.0*averagecount) #line 1179

    #calc median
    temptwo = np.empty((100))
    temptwo[:] = np.nan

    for i in range(averagecount):
        temptwo[i] = ystore[i]

    temptwo = temptwo[~np.isnan(temptwo)] #want to avoid extra zeros
    medianme = np.median(temptwo)

    for jf in range(flat):

        flatc = 0

        for i in range(int(flatcount[jf])): #line 1207 - called calc flat median

            if (i >= 1): #multiple points to ave over - shoudl it be greater than just 2?

                if (flatstore[jf,i] > 1.0*flatstore[jf,i-1]): #unsure if shold be i-1 or i
                    flatc = flatc +1

                if (flatstore[jf,i] < 1.0*flatstore[jf,i-1]): #looking at point behind it
                    flatc = flatc - 1


        if (abs(flatc) == (flatcount[jf] -1)) and (flatcount[jf] > 3): #unsure if should be greater than 2 or greater than 3

            if (abs(flatstore[jf,0] - flatstore[jf,flatcount[jf]]) > 10):
                flatcount[jf] = 0
                #go to 911

                break #unsure if this is correct
        else:
            pass
            #print('conds not met') #1st point do this as zero

        if (flatcount[jf] <= 1): #assumes just needs 2 points included

            flatcount[jf] = 0
            #go to 911
            break #cut of whole loop - cut of loop for jf in range flat - seems correct

        #calc flat median #for i in range(flatcount[jf] - 1): #unsure if this is the correct limit
        flatmedianme[jf] = np.median(temptwo[:, jf])   #temptwo[i] = flatsotre range for i in length flatcount add in values , try just do upto certain amount


    #set dispersion etc - line 1257
    if (fortyone == 0):
        dispersion = 0.0
        vtwo = 0.0 #but not if reached 41 - in 41 section - set variable to one 

    for i in range(afpick, averagecount): #this all moves to diff point
        #afpick given by file called sirmms.dat
        dispersion = dispersion + (((ystore[i] - averageme/sigmtot))**2)/sigmstore[i]
        vtwo = vtwo + (1/sigmstore[i])**2


    #flat section normal districution standard distribution
    for i in range(flat):
        sigmtotflat[i] = 0
        for k in range(int(flatcount[flat])):

            sigmtotflat[i] = sigmtotflat[i] + ((flatstore[i,k] - (averageflat[i]/flatcount[i]))**2)/flatcount[i]

        sigmtotflat[i] = sqrt(sigmtotflat[i])
        #end of flat normal

    sirmrms = 0
    print("All demag and palaeointensity data is saved in afdemag-all-{0}.dat".format(name))
    print("All demag and data regarding the SIRM data is saved in afdemag-sirm-{0}.dat".format(name))
    fall = open("afdemag-all-{0}.dat".format(name), "w") #'w' will overwrite file is already exists
    fall.write('afield' + '\t' + 'step' + '\t' + 'stepPI' + '\t' + 'std' + '\t' + 'select' + '\t' + 'flatno'+ '\t' + 'shortminus' + '\t' + 'SIRMAFS%' + '\t' + 'SIRMAFM%' + '\t' + 'AFNRM/SIRM-M%' + '\t' + 'AFNRM/SIRM-S%' + '\t' +  'AFNRM/NRM-M%' + '\t' + 'AFNRM/NRM-S%' + 'shortdiv%')
    fsirm = open("afdemag-sirm-{0}.dat".format(name), "w")

    fsirm.write('afield' + '\t' + 'measured' + '\t' + 'simulated')
    for i in range(cntfieldaf): #0 to 21
        sirmrms = sirmrms + abs((af_sirm_n[i]/af_sirm_n[0]) - sirm[j-1, i]/sirm[j-1][0])
        #print(sirmrms) #increases 0.004 -> 2.6
        fall.write('\n')
        #f.write(str(af_step[i]))
        fall.write(str(af_step[i]) + '\t' + str(i) + '\t' + str(ys[i]) + '\t' + str(std[i]) + '\t' + str(used[i]) + '\t' + str(flatused[i]) + '\t' + str((((sirm[j-1,i]/sirm[j-1,0])-(af_sirm_n[i]/af_sirm_n[0])))*100) + '\t' + str(sirm[j-1,i]/sirm[j-1,0]*100) + '\t' + str((af_sirm_n[i]/af_sirm_n[0])*100) + '\t' + str((af_nrm_n[i]/af_sirm_n[i])*100) + '\t' + str((afmag[j-1,i]/sirm[j-1,i])*100) + '\t' + str((af_nrm_n[i]/af_nrm_n[0])*100) + '\t' + str((afmag[j-1,i]/afmag[j-1,0])*100) + '\t' + str(abs(1-((sirm[j-1,i]/sirm[j-1,0])/(af_sirm_n[i]/af_sirm_n[0])))*100))
        fsirm.write('\n')
        fsirm.write(str(af_step[i]) + '\t' + str(af_sirm_n[i]/af_sirm_n[0]) + '\t' + str(sirm[j-1, i]/sirm[j-1,0]))

    fsirm.close()

    fall.write('\n')

    sirmrms = sirmrms/cntfieldaf

    fall.write('SIRM MEAN DIFFERENCE % =' + '\t '+ str(100*sirmrms))
    fall.write('\n')
    if (averagecount !=0) and (averagecount !=1) and (sigmtot !=0): #only do both is average at least 1 
        sampleunbias=dispersion*(sigmtot/((sigmtot**2)-vtwo))
        dispersion=dispersion/(averagecount-1)


    if (averagecount == 1): #here equals 11
        dispersion = 1

    if (sigmtot == 0):
        sigmtot = 1
        
        print('sigm tot = 0 weighted average not possible')
        
    if (averagecount == 0):
        averagecount = 1
        print('avercount = 0 weighted average not possible')   
        
    #fall.write('weighted average = ' + '\t' + str(averageme/(sigmtot)) + '\t' + str(sqrt(dispersion/sigmtot)))
    #fall.write('\n')
    #fall.write('unweighted average = ' + '\t' + str(averageunweight/averagecount))
    #fall.write('\n')
    #fall.write('unweighted median = ' + '\t' + str(medianme))
    #fall.write('\n')
    #fall.write('cntfields = (cntfield, cntfieldaf, afpick) ' + '\t' + str(cntfield) + '\t' + str(cntfieldaf) + '\t' + str(afpick))


    #determining which is best flat 
    maxjf = 0
    minsig = 10000
    jfs = 0

    for jf in range(flat):

        if ((flatcount[jf] >= maxjf) and (flatcount[jf] > 0)):

            if ((flatcount[jf] == maxjf) or (maxjf == 1)): #change maxjf == 0

                if (sigmtotflat[jf] < minsig): #seem to do same thing regardnless of these if statements

                    minsig = sigmtotflat[jf]
                    maxjf = flatcount[jf]
                    jfs = jf

            else:
                minsig = sigmtotflat[jf]
                maxjf = flatcount[jf]
                jfs = jf


        if (flatcount[jf] > 0):

            fall.write('\n')
            fall.write('unweighted average flat (jf, averageflat(jf)/flatcount(jf), sigmtotflat(jf, flatdi(jf)) =' + '\t' + str(jf) + '\t' + str(averageflat[jf]/flatcount[jf]) + '\t' + str(sigmtotflat[jf]) + '\t' + str(flatdi[jf]) + '\t' + str(flatcount[jf]) + '\t' + str(flatdi[jf]/flatcount[jf]))
            fall.write('\n')
            fall.write('unweighted flat median (jf, flatmedianme(jf) = ' + '\t' + str(jf) + '\t' + str(flatmedianme[jf]))
        else:

            fall.write('no points selected for flat section (jf)' + '\t' + str(jf))


    fall.close()        
    aconst = -aconst*log(0.01*(tempmin)/(tempmax-tempmin))/3600
    print('Output results to averageout_{0}.dat'.format(name))
    fave = open('averageout_{}.dat'.format(name), 'w') #117 file
    if (averageme > 0):

        if (jfs != 0):

            fave.write(str(averageme/(sigmtot)) + '\t' + str(sqrt(dispersion/sigmtot)) + '\t' + str(averagecount) + '\t' + str(averageunweight/averagecount) + '\t' + str(medianme) + '\t' + str(flatcount[jfs]) + '\t' +  str(averageflat[jfs]/(1.0*flatcount[jfs])) + '\t '+  str(sigmtotflat[jfs]) + '\t' + str(flatmedianme[jfs]) + '\t' + str(aconst) + '\t' + '147' + '\t' + str(100*af_nrm_n[0]/af_sirm_n[0]) + '\t' + str(sirmrms) + '\t' + str(selrms))
        else:
            fave.write(str(averageme/sigmtot) + '\t' + str(sqrt(dispersion/sigmtot)) + '\t' + str(averagecount) + '\t' + str(averageunweight/averagecount) + '\t' + str(medianme) + '\t' + '0.0' + '\t' + '0.0' + '\t' + '0.0' + '\t' + str(aconst) + '\t' + '147' + str(100*(af_nrm_n[0]/af_sirm_n[0])) + '\t' + str(sirmrms) + '\t' + str(selrms))
    else:
        fave.write('no points selected' + '\t' + str(sirmrms))


    fave.close()
    V['shortdivlist'] = shortdivlist
    V['shortminuslist'] = shortminuslist
    V['ys'] = ys
    
    return(X,V)
    
    
def plot_sirm_check(X,V):
    af_step_list = X['af_step']
    shortdivlist = V['shortdivlist']
    shortminuslist = V['shortminuslist']
    name = X['name']
    shortdivlistplot = np.array(shortdivlist[:(len(af_step_list))])
    shortminuslistplot = np.array(shortminuslist[:(len(af_step_list))])
    shortdivlistplot[np.isnan(shortdivlistplot)] = 0
    shortminuslistplot[np.isnan(shortminuslistplot)] = 0
    twenty = []
    hundred = []

    for i in range(len(af_step_list)):
        twenty.append(20)
        hundred.append(100)
    w, h = figaspect(1) #from https://stackoverflow.com/questions/48190628/aspect-ratio-of-a-plot
    fig, ax = plt.subplots(figsize=(w,h))

    plt.plot(af_step_list, twenty, 'b') #af step list is just af_step?
    plt.plot(af_step_list, hundred, 'r')
    plt.plot([af_step_list[0], af_step_list[-1]], [-20, -20], 'b') #af step list is just af_step?
    plt.plot([af_step_list[0], af_step_list[-1]], [-100, -100], 'r')
    plt.plot(af_step_list, shortdivlistplot,  marker='o', color= 'r')
    plt.plot(af_step_list, shortminuslistplot,  marker='o', color = 'b')
    plt.title('SIRM checks looking at difference between measured and simulated. Sample {0}'.format(name))

    plt.xlabel('AF peak (mT)')
    plt.ylabel('S_diff or S_ratio (%)')
    plt.text(25, 85, r'S diff ', fontsize=12)
    plt.plot(22,87, marker = 'o', color='b')
    plt.text(25, 75, r'S ratio', fontsize=12)
    plt.plot(22, 77, marker = 'o', color='r')
    print('Figure saved as SIRM-checks-{}.svg'.format(name))
    plt.savefig('SIRM-checks-{}.svg'.format(name))
    plt.show
    return
    
def plot_pal(V,X):
    
    ys = V['ys']
    af_step = X['af_step']
    name = X['name']
    afpick = X['af_pick']
    
    w, h = figaspect(1) #from https://stackoverflow.com/questions/48190628/aspect-ratio-of-a-plot
    fig, ax = plt.subplots(figsize=(w,h))

    #plot with red colour and calc average
    #selected_mean = np.mean(ys[8:20])
    #mean_dev = np.std(ys[8:20])
    #selected_med = np.median(ys[8:20])
    plt.plot(af_step,ys[:len(af_step)], 'b')
    plt.plot(af_step,ys[:len(af_step)],  marker='o', color= 'b')
    #plot vertical line
    
    plt.plot([af_step[afpick], af_step[afpick]], [0,np.max(ys[:len(af_step)])], color = 'green')
    #plt.plot(af_step[8:20],ys[8:20],  marker='o', color= 'r')
    
    #plt.text(20, 6, r'selected median: %.2f $\mu$T'%selected_med, fontsize=11)
    #plt.text(20, 7, r'rejected mean: %.2f $\pm$ %.2f $\mu$T'%(selected_mean ,mean_dev), fontsize=11)
    #plt.text(20, 6, r'selected median: %.2f $\mu$ T'%selected_med, fontsize=12)
    plt.xlabel('AF peak (mT)')
    plt.ylabel('paleointensity (\u03BCT)')

    plt.title('TRM PI est {}'.format(name))

    plt.savefig('ESA_MF-PI_ESTS_colour_{}.svg'.format(name))
    plt.show
    
    ind = []
    for i in range(len(af_step)):
        ind.append(i)
        i+=1


    ind = np.array(ind)

    print("AF step index = AF step")
    for n, v in zip(ind, af_step):
        print("{} = {}".format(n, v))



#ARM blockfind function

def blockfind_SC_arm(t, field, afswitch, V, X): #try removing counttime as input/output

    hys = V['hys']
    prop = V['prop']
   
    max_total_m_a = hys[:,10]*hys[:,3]
    max_total = sum(max_total_m_a)
	
    vol_rem = V['vol_rem']
	
    num_hyss = V['num_hyss']

    histore = V['histore']
    beta = V['beta']
    rate = V['rate']

    tm = V['tm'] #try setting this for use in AF demag
    tstep = V['tstep']
    
    temp = V['temp']
    #temp = tempmin + tempstep
   
    hfstore = np.zeros(num_hyss)


    blockper = V['blockper'] 

    blocktemp = V['blocktemp'] 
    boltz = V['boltz']
    blockg = V['blockg']

    #end_mag = V['end_mag']

    totalm = V['totalm']
    afstore = V['afstore']
    arm_blocked = V['arm_blocked']

    af_field = V['afstore']
    #print(af_field)


    var_1 = np.array((num_hyss, beta, blockper, temp, t, tstep, rate, field, tm))
    #print('var_1', var_1)
    hfstore, histore, boltz, blocktemp = block_loc_C_arm(var_1, hys, prop, blockg, boltz, af_field)
  
    blockper, totalm, boltz, blockg, arm_blocked = block_val_C_arm(hys, prop, histore, hfstore, blocktemp, beta, num_hyss, boltz, blockg, field, afswitch, afstore, max_total, vol_rem, t, tstep, arm_blocked, blockper)
 
    V['aftotalm'] = totalm
    V['sir'] = totalm

    V['blockper'] = blockper
    V['blocktemp'] = blocktemp
    V['boltz'] = boltz
    V['blockg'] = blockg
    V['totalm'] = totalm
    V['arm_blocked'] = arm_blocked
    V['tm'] = tm

    return #totalm

#CRM acquisition functions
#CRM acquisition and AF demag 
@jit(nopython = True)
def block_loc_C_arm(var_1, hys, prop, blockg, boltz, af_field):
     #unpack variables from input array
    #unpack var_1
    num_hyss = int(var_1[0])
    beta = float(var_1[1])
    blockper = float(var_1[2])
    temp = float(var_1[3])
    t = float(var_1[4])
    tstep = float(var_1[5])
    rate = float(var_1[6])    
    field = float(var_1[7])
    tm = float(var_1[8])
    #print('blockper for calc histore in ARM', blockper)

    tau = 1e-9
 
    roottwohffield = 2**(0.5)
    hfstore = np.zeros((num_hyss))
    histore = np.zeros((num_hyss))    
    hcstore = np.zeros((num_hyss))   
    blocktemp = np.zeros((num_hyss))
    i=0 
    for i in range(num_hyss): #dynamic contributions - to hc dormann 1988 used?
        if (t >= (prop[i,3] + (0.5*tstep))):

            phitemp=((sin(hys[i,6])**(2./3.))+(cos(hys[i,6])**(2./3.)))**(-3./2.) #phitemp from hys[5] -> phi               
            hc=(sqrt(2))*(hys[i,9]*beta) #test and time HCbysqrt2 to get hc and matfhes similar order mag to hys[i,0/mu0 -> hc in same units using in hk - seems correct]
            
            hcstore[i] = hc/(sqrt(2.)) #different to hys[0]
        
            
            hi = hys[i,1]*beta*blockper #divide here by mu0
            
    
            histore[i] = hi/(sqrt(2.)) #should all be zero, without blockper, all look milar to hys[1,0] in mg

            g_ang=0.86+1.14*phitemp
            #print(hys[i,9], beta
            af_field_cor = (af_field*1E-6)/mu0 # 100 mT
            h_hs = af_field_cor - histore[i] #100 applied field 
            v_co = (1 - ((abs(h_hs))/(hys[i,9]*beta)))
            if (v_co <= 0):
                v_co = 0.00001               
            
            v_act_c = prop[i,5]*v_co

            hf = hys[i,11]*(1/(v_act_c))
        
            hfstore[i] = hf #values at expected as hc/mu0
            #this needs to be the correct hf for arm acquisition 

            if (rate == 1): #rate first set to 1 so shold do this
                
                tm= 60. #prop[i,10]
                #60. 
                if (tm == 0.0): #unsure
                    tm = 60.0
    
            #print('tm', tm)
            
            #if (end_mag == 1.0):
            #    tm = 60.0
            ht = (roottwohffield)*hf*(log(tm/tau)) #new tm 
        
            
            bracket = 1-(2*ht*phitemp/hc)**(1/g_ang)
            
            hiflipplus = hc*bracket+field*(roottwohffield) # using hs=-hi then field is +ve not -ve, 
    
            hiflipminus=-hc*bracket+field*(roottwohffield) #see which way fields flips
            trialtemp=0

            if (hc >= (2*ht*phitemp)): #still loop over each hys, +1
                
                if ((hi > hiflipminus) and (hi < hiflipplus)): #+2 blocked

                    if (blockg[i] == 0) or (blockg[i] == 2) or (blockg[i] == -2): #+3 prev blocked until this point
                        
                        if (hi >= (field*roottwohffield)): #+4
                            
                            blocktemp[i] = -1
                            
                        else:
                            #print('blockg = 0,2,-2 and hi < field*rootwo', 'blockg', blockg[i], 'hi', hi, 'field*roottwo', field*roottwohffield)
                            blocktemp[i] = 1 #end +3 unsure if sholud ended both or just one
                            
                    elif (blockg[i] == -1): #this line, already block stay , not need
                        
                        blocktemp[i] = -1
                

                    elif (blockg[i] == 1):
                        blocktemp[i] = 1

                elif (hi >= hiflipplus):#else field blocking above ht, this is hi test hiflip etc
                    # if ( abs(blockg[i]) != 1) and (rate == 1): #this occurs if blocked a hysteron during cooling
                    blocktemp[i] = -2
                
                        
                else:
                    blocktemp[i] = 2

            else: #hc < 2*ht*phitemp. this is correctm meaning else above isnt
                #print('SP', hc, 2*ht*phitemp)
        
                if ((hi < hiflipminus) and (hi > hiflipplus)): 
                    blocktemp[i] = 0
                else: #field blocking - below hc
                    if (hi >= hiflipminus):
                        blocktemp[i] = -2

                    else:

                        blocktemp[i] = 2

            if (temp < trialtemp):
                #need print to screen
                print('blocktemp', blocktemp[i])
                
    return hfstore, histore, boltz, blocktemp

@jit(nopython = True)
def block_val_C_arm(hys, prop, histore, hfstore, blocktemp, beta, num_hyss, boltz, blockg, field, afswitch, afstore, max_total, vol_rem, t, tstep, arm_blocked, blockper_current):

    totalm = 0.0
    totalmoment = 0

    i=0
    num_blocking = 0
    blockper = 0.
    #print('blockper current', blockper_current)
    #use old blockper here
    field_block_pos = 0.
    field_block_neg = 0.
    dc_blocked_pos = 0.
    dc_blocked_neg = 0.
    dc = (50.*1E-6)/mu0 #0.05/(1000.*(mu0))
    af = afstore/(1000.*(mu0))
    #print('af field', af)
    #print('dc field used boundaries', dc)
    for i in range(num_hyss):

        
        blockg[i] = blocktemp[i]
        absblockg = abs(blockg[i])
        if (t >= (prop[i,3] + (0.5*tstep))):
           
            if (absblockg == 1): #if blockg 1 or -1
                if (boltz[i] < 0.00000001) and (boltz[i] > -0.000000001): #only zero
                    #print('setting boltz as its 0 (check)', boltz[i])
                    boltz[i] = tanh((field - histore[i])/hfstore[i])
                    if (afswitch == 0):
                        num_blocking+=1
                    if (afswitch == 1):
                   # print('inside tanh reset during af demag', (((50E-6/mu0) - histore[i])/hfstore[i]))
                        boltz[i] = 0.

            if (blockg[i] == -2):
                #print('blockg = -2')
                moment = -0 #changed from zero to minues zero
            elif (blockg[i] == 2):
                #print('blockg = +2')
                moment = 0
            else:

                moment = blockg[i] #where does moment come from? - set to blockg if blockg not equal to 1


            
            
            if (afswitch == 1): #afswtich given in function, 1 or 0
                # afstore =  CT['afstore'] # move where relevant as no values unti demag - poss give dummy name to keep with rest of unpack of V
                #break
                #print('should see this message during ARM acq, afswtich 1')
                hi = histore[i]*(sqrt(2))*blockper_current #hi depend on number blocked - occur slowly 
                #print('blockper current', blockper_current)
                #print('blockper', blockper) 
                hc = (sqrt(2))*(hys[i,9]*beta)*(((sin(hys[i,5])**(2./3.))+(cos(hys[i,5])**(2./3.)))**(-3./2.))
                #uses hf
                #what is hf during ARM acquisition 

                af = afstore/(1000.*(mu0)) #define later. afstore number do in 100 mTfield 100 mT, 0.05 mT
                #dc = (50.*1E-6)/mu0 #0.05/(1000.*(mu0))


                
                if (hi > (hc - af + dc)):#  and (hi >= 0): remove need for hi > 0 so it does SP all in 1 step 
                    boltz[i] = 0.
                    arm_blocked[i] = 0. #FB so no blockper contribution
                    moment = -1.
                elif (hi < (-hc + af + dc)): #remove need to be hi less than zero so it does SP too
                    boltz[i] = 0.
                    arm_blocked[i] = 0. #FB so no blockper contribution
                    moment = 1.
                #most should overright 
                elif (hi <= dc) and (boltz[i] == 0.): #if not in these bounds its remanance
                    boltz[i] = 1.
                    arm_blocked[i] = 1.
                    moment = 1. 
                elif (boltz[i] == 0.): #this should be the bit above - which is slightly less
                    boltz[i] = -1.                    
                    arm_blocked[i] = -1.
                    moment = -1.
                else:
                    #everything stay the same
                    boltz[i] = boltz[i]
                    arm_blocked[i] = arm_blocked[i]
                    moment = moment
            #totalm = totalm + (abs(moment)*abs(cos(hys[i,5]))*beta*(boltz[i])*((hys[i,3]*prop[i,5])/max_total)) #add in forcdist (hys[i,3])
            totalm = totalm + (abs(moment)*abs(cos(hys[i,5]))*beta*(boltz[i])*(hys[i,3])) #*prop[i,5])/max_total))
            #totalm = totalm + (abs(moment)*abs(cos(hys[i,5]))*beta*(boltz[i])*hys[i,3]) #*((prop[i,5])/(max_total)))
            #changed angle here to hys[6]
            #each hysteron
            #if (i<3) or (i>9997):
            #    print(i, '\t', hys[i,1], '\t', hi, '\t', blockper_current, '\t', boltz[i], '\t', arm_blocked[i], '\t', (hc-af+dc), '\t', (-hc+af+dc))
            #    print('blockper current', blockper_current)
                #print('blockper', blockper) 
            #if (blockg[i] == -1) or (blockg[i] == 1):
            #blockper come from arm blocked
            if (arm_blocked[i] == -1) or (arm_blocked[i] == 1): #(blockg[i] == -1) or (blockg[i] == 1):
                moment_block = 1
            else:
                moment_block = 0
		
            blockper=blockper+ ((abs(moment_block))*prop[i,5]*hys[i,3]) #*hys[i,10])
            totalmoment=totalmoment+moment

    #just do each step 
    #print('sum armblocked', np.sum(arm_blocked))
    #print('af step', afstore, 'sum boltz', np.sum(boltz[:]))


    
        
    #print('printed all numbers each time')
    #print('af at end of this af step')
    if (blockper != 0.):
        blockper= (blockper)/(vol_rem)
    else:
        #blockper = 0.
        blockper = 1E-15
    #print(num_blocking)
    return blockper, totalm, boltz, blockg, arm_blocked