import numpy as np              
import matplotlib.pyplot as plt 
from scipy.signal import find_peaks
import os
#parameters and function defintion
n=1 # air
c=3e8 # speed of light in air 
lambda_m =np.linspace(400,800,10000000) * 1e-9 #10,000 points between 400nm and 800nm 
wavelength_step_nm = (lambda_m[1] - lambda_m[0]) * 1e9
T_by_R = {}
Peaks_by_R = {}
def compute_transmission(R,L,n,wavelength):
    delta_m= 4*np.pi*L*n/wavelength
    F=(4*R)/(1-R)**2
    T= 1/(1+F*(np.sin(delta_m/2)**2))
    return F,T
L=10e-6 # cavity length
Reflectivety_values=[0.5,0.8,0.95]
plt.figure(figsize=(15,4))
for R_1 in Reflectivety_values:
    F_loop,T_Loop= compute_transmission(R=R_1, L=L,n=n , wavelength=lambda_m)
    plt.plot(lambda_m*1e9,T_Loop,label=f'R={R_1}')
    #Part 3 , To get finesse from graph 
    T_by_R[R_1] = T_Loop 
# ok now i have 3 vectors representing each R , to find FSR and fwhm for each 
    Peaks_Loop,_= find_peaks(T_by_R[R_1])
    Peaks_by_R[R_1]=Peaks_Loop
    mid_idx = len(Peaks_Loop) // 2
    middle_peak_idx = Peaks_Loop[mid_idx]
    next_peak_idx = Peaks_Loop[mid_idx + 1]
    FSR = (next_peak_idx - middle_peak_idx) * wavelength_step_nm
    
    # Slice a local window 
    FSR_indices = next_peak_idx - middle_peak_idx
    window_radius = FSR_indices // 2
    window_start = middle_peak_idx - window_radius
    window_end = middle_peak_idx + window_radius
    T_window = T_Loop[window_start:window_end]
    
    indices_above_half = np.where(T_window >= 0.5)[0]
    
    # index width 
    fwhm_index_width = indices_above_half[-1] - indices_above_half[0]
    
    FWHM = fwhm_index_width * wavelength_step_nm
    
    # Calcuate Finesse 
    measured_finesse = FSR / FWHM
    print(f"R={R_1} ,F={measured_finesse}")
plt.title('Transmittance Plot at different Reflectivities')
plt.xlabel('Wavelength(nm)')
plt.ylabel('Transmittance')
plt.grid()
plt.legend(loc="upper right")
#Part 2 
L_values=np.array([5,10,20])*1e-6
plt.figure(figsize=(15,4))
plt.gca().set_prop_cycle(color=['#5E3C99', '#FDB863', '#2CA02C'])
for L_1 in L_values: 
    _,T_Loop2= compute_transmission(R=0.9,L=L_1,n=n,wavelength=lambda_m)
    plt.plot(lambda_m*1e9,T_Loop2,label=f'L={np.ceil(L_1/1e-6)}μm')
plt.title('FSR comparison with changing cavity Lengths')
plt.xlabel('Wavelength(nm)')
plt.ylabel('Transmittance')
plt.grid()
plt.legend(loc="upper right")
plt.xlim(400,500)
    #PART 3 
# Calculate F from each R 
R_array= np.array(Reflectivety_values)
Finesse_values_equation = (np.pi * np.sqrt(R_array)) / (1 - R_array)
print(f"R= 0.5, F={Finesse_values_equation[0]}")
print(f"R=0.8 , F={Finesse_values_equation[1]}")
print(f"R=0.95 , F={Finesse_values_equation[2]}")
plt.show()



