# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#import matplotlib.pyplot as plt
#import matplotlib.ticker as ticker
import numpy as np
import struct
import sys
from array import array
import itertools
import os
import sys
from scipy import signal
from scipy.optimize import curve_fit
import subprocess
import shutil

pi = np.pi
clight = 299792458
working_folder = 'eSR/1RF_1type/Transient_V/'
home = os.getcwd()
cwd = os.path.join(home,working_folder)
inputfile = 'input0.txt'
inputfile = os.path.join(cwd,inputfile)
tempinput = {}
with open(inputfile) as inputfile:
    for line in inputfile:
        if len(line.split())>1:
            tempinput[line.split()[0]] = line.split()[1:]
for i in tempinput:
    for j in range(len(tempinput[i])):
        tempinput[i][j] = float(tempinput[i][j])
    
nRF = int(tempinput['nRF'][0])
nRF1 = int(tempinput['nRF1'][0])
nRF2 = int(tempinput['nRF2'][0])
nRFc = int(tempinput['nRFc'][0])

E0Au = 196.9665687*931.5e6
E0Elec = 0.51099895000e6
nTurns = int(tempinput['n_turns'][0])
nfill = int(tempinput['n_fill'][0])
n_q_ramp = int(tempinput['n_q_ramp'][0])
NpRF = int(tempinput['N_bins'][0])
h = [int(i) for i in tempinput['h']]
detune_ini = np.array([i for i in tempinput['detune_ini']])
detune_final = np.array([i for i in tempinput['detune_final']])

step = int(tempinput['step_store'][0])
fill_step = int(tempinput['fill_step'][0])
nBeam = int(tempinput['nBeam'][0])
beam_shift = int(tempinput['beam_shift'][0])
nBunch = int(tempinput['n_bunches'][0])
nPar = int(tempinput['Npar'][0])
NperBunch = int(tempinput['NperBunch'][0])
nTot = nBunch*nPar*nBeam
Gamma0 = tempinput['Gamma'][0]
Rring = tempinput['R'][0]
n_record = nTurns/step
clight = 299792458
beta = np.sqrt(1-1/Gamma0**2)
T0 = 2*np.pi*Rring/(clight*beta)
f0 = 1/T0
V0 = [i for i in tempinput['Vref_I']]
V0Q = [i for i in tempinput['Vref_Q']]
II = [i for i in tempinput['Iref_I']]
IQ = [i for i in tempinput['Iref_Q']]
mainRF = int(tempinput['mainRF'][0])
if int(tempinput['type'][0]==2):
    atomicZ = 79
    Ek = Gamma0*E0Au
else:
    atomicZ =1
if int(tempinput['type'][0]==1):  
    Ek = Gamma0*E0Elec
    
GMTSQ = tempinput['GMTSQ'][0]
Ek_damp = tempinput['Ek_damp'][0]

eta = 1/GMTSQ-1/Gamma0**2
Qs = np.sqrt(h[mainRF]*atomicZ*np.abs(V0[mainRF])*eta/(2*np.pi*Ek))

omegarf = 2*np.pi*(np.array(h)*f0)
omegac = 2*np.pi*(np.array(h)*f0+detune_final)
Trf = 2*np.pi/omegarf
RoQ = [i for i in tempinput['RoQ']]
QL = [i for i in tempinput['QL']]
R = [RoQ[i]*QL[i] for i in range(nRF)]

Th = 2*np.pi/omegarf[0]
dthat =Th/NpRF

pattern = 'd'+'dd'*nBeam+3*nRF*'d'
n_stride = 1+2*nBeam+3*nRF
stride = len(pattern)*8
test = array('d')
bucket_height = 2*Qs/(h[mainRF]*eta)*Gamma0

print(bucket_height)
print(Ek)
print(Qs)

iMin = 2.5
iMax = 2.5
nParMin = iMin/f0/nBunch/1.6e-19
nParMax = iMax/f0/nBunch/1.6e-19


N_samples = 1 
N_thetaL = 1 
ThetaL_min = np.zeros(nRF)#17.5
ThetaL_max = np.zeros(nRF)#17.5

ThetaL_min[0] = 0
ThetaL_max[0] = 0

dnPar = (nParMax-nParMin)/N_samples
dThetaL = (ThetaL_max-ThetaL_min)/N_thetaL

for charge_factor in range(N_samples):
    for thetaL_factor in range(N_thetaL):
        # arguments
        ParType = 1 # 0 means proton, 1 means electron, 2 means gold

        mainRF = 0
        main_detune = 0
        detune_slow_factor = 1

        nTurn = 10000
        step_store = 1000
        n_record = nTurns/step_store

        n_dynamic = 30000

        n_fill = 1000
        n_q_ramp = 2000
        n_detune_start = 1000
        n_detune_ramp = 3000
        n_detune_ramp_tot = 3000 # last turn of detuning process
        n_I_ramp_start = 1000
        n_I_ramp_end = 3000

        R_ring = 610.1754 
        GMTSQ = 961.0 
        Gamma0 = 19600.0 


        t_rad_long = 0.03555 
        rms_poverp = 5.8e-14 # quantum excitation caused equalibrium dp/p
        siglong = rms_poverp*Gamma0# 30.76 # rms d_gamma
        Ek_damp = 5e8 # artificial dampping

        nRF = 1
        nRF1 = 1.0 
        nRFc = 0.0 
        nRF2 = 0.0 

        nCav = np.array([14,0,0]) # number of fundamental cavities
        h = np.array([7560,7560,7560*3]) 
        RoQ = np.array([73,73,15])# RoQ per cavity, circuit defination
        gII = np.zeros(nRF)
        gQQ = np.zeros(nRF)
        for i in range(nRF):
            gII[i] = 3
            gQQ[i] = 3
        delay = 350 # number of RF turns 
        nBunch = 580
        nPar0 = 17.2e10*2
        Prad0 = 10e6
        fill_step = 12
        nPar = nParMin+dnPar*charge_factor #nPar0/N_samples*(charge_factor+1)
        Prad = Prad0*nPar/nPar0

        N_macro = 1 
        nBins = 33 
        NC = nCav[0]+nCav[1]
        NF = nCav[0]
        ND = nCav[1]
        
        thetaL = np.zeros(nRF)
        Vs = np.zeros(nRF)
        Vq = np.zeros(nRF)
        Phis = np.zeros(nRF)
        PhisPhasor = np.zeros(nRF)
        PhiIQ = np.zeros(nRF)
        PhiIQIg = np.zeros(nRF)
        VrefI = np.zeros(nRF)
        VrefQ = np.zeros(nRF)
        Vreftot = np.zeros(nRF)
        IrefI = np.zeros(nRF)
        IrefQ = np.zeros(nRF)
        QL = np.zeros(nRF)
        Rsh = np.zeros(nRF)
        
    # Need to calculate the required voltage and phase
    # then calculate the inputs for the code, namely VrefI, VrefQ, IrefI, IrefQ.
    
    # for fundamental, 
        Vtot = 23.7e6 # total voltage 
        Urad0 = Prad/(nBunch*nPar*1.6e-19*f0) # radiation caused Voltage total
        U_loss = Urad0/NC #+2*pi*h*f0*RoQ*NC/4*nPar*1.6e-19 # loss per cavity for fundamental
        V0 = Vtot/NC
        Vsynch_need = U_loss
       # Vsynch_need = 3.7e6/NC
        Vquard_need = V0*np.sin(np.arccos(U_loss/V0))

        
        
        Vnew = np.sqrt(Vsynch_need**2+(NC/(NF-ND)*Vquard_need)**2) # new cavity voltage per cavity, assuming the new phiSynch are the same (different sign) between two types of cavity
        Vs[0] = Vsynch_need
        Vq[0] = NC/(NF-ND)*Vquard_need
        
        # new synchronous phase if we change the number of focusing and defocusing cavity.
        
        PhisPhasor[0] = np.arctan(Vq[0]/Vs[0])
        
        
        
        Pbeam0 = Prad0/NC # beam power per fundamental cavity, 

        Pbeam = Prad0*nPar/nPar0/NC # beam power per fundamental cavity
        RoQacc = RoQ*2
        IbDC = nBunch*nPar*1.6e-19*f0
        f = h*f0

        Qbeam0 = Vnew**2/(RoQacc*Pbeam)
        
        Qbeam = Qbeam0
        QL[0] = Qbeam[0]
        
        Rsh[0] = RoQ[0]*QL[0]
        
        
        Vreftot = np.sqrt(Vs**2+Vq**2)
        
    # Now calculate the inputs
    
        thetaL = (ThetaL_min+dThetaL*thetaL_factor)/180.0*pi  # angle between Ig and Vc
        Vbr = 2*IbDC*Rsh
        # special case for V_trans calculation
        #Vbr = 2*IbDC*Rsh*630/580
        Vgr = Vreftot/np.cos(thetaL)*(1+Vbr/Vreftot*np.cos(PhisPhasor))
        
        tgPhi = -(Vbr*np.sin(PhisPhasor)/Vreftot+(1+Vbr*np.cos(PhisPhasor)/Vreftot)*np.tan(thetaL))
        tgPhi_ini = -np.tan(thetaL)
        delta_f_ini = f*(tgPhi_ini/2/QL+np.sqrt((tgPhi_ini/2/QL)**2+1))-f
        delta_f = f*(tgPhi/2/QL+np.sqrt((tgPhi/2/QL)**2+1))-f#-f*tgPhi/2/Qbeam
        #delta_f_ini = -f*tgPhi_ini/2/QL
        #delta_f = -f*tgPhi/2/QL #-f*tngPhi/2/Qbeam

        
        VrefI = Vreftot*np.sin(PhisPhasor)
        VrefQ = -Vreftot*np.cos(PhisPhasor)
        
        I_I = 2*IbDC*(tgPhi+np.tan(PhisPhasor))/tgPhi*np.cos(PhisPhasor)/np.cos(thetaL)*np.sin(PhisPhasor+thetaL) # becareful not to forget the factor of '2'
        I_Q = -2*IbDC*(tgPhi+np.tan(PhisPhasor))/tgPhi*np.cos(PhisPhasor)/np.cos(thetaL)*np.cos(PhisPhasor+thetaL)
        
        I_I = Vgr/Rsh*np.sin(PhisPhasor+thetaL) # becareful not to forget the factor of '2'
        I_Q = -Vgr/Rsh*np.cos(PhisPhasor+thetaL)
        
        I_I_ini = Vreftot/(Rsh)/np.cos(thetaL)*np.sin(PhisPhasor+thetaL)
        I_Q_ini = -Vreftot/(Rsh)/np.cos(thetaL)*np.cos(PhisPhasor+thetaL)
        
        print("Vnew : ",Vreftot)
        print("QL : ",QL)
        print("ThetaL : ",thetaL/pi*180, " [degree]")
        
        print("Tan(PhisPhasor) : ",np.tan(PhisPhasor))
        print("PhisPhasor : ",PhisPhasor/pi*180)

        print("detune tan : ", tgPhi)
        print("detune angle : ", np.arctan(tgPhi)/pi*180, " [degree]")
        print("delta_f : ",delta_f, " [Hz]")
        print("VrefI : ",VrefI)
        print("VrefQ : ",VrefQ)
        
        print("II : ",I_I)
        print("IQ : ",I_Q)
        print("II_ini : ",I_I_ini)
        print("IQ_ini : ",I_Q_ini)

        print("VrefTot : ",Vreftot)
        print("IbDC : ", IbDC)
        

        tempinput['n_turns'][0] = nTurn
        tempinput['step_store'][0] = step_store
        tempinput['n_fill'][0] = n_fill
        tempinput['n_q_ramp'][0] = n_q_ramp
        tempinput['n_detune_start'][0] = n_detune_start
        tempinput['n_detune_ramp'][0] = n_detune_ramp
        tempinput['n_detune_ramp_tot'][0] = n_detune_ramp_tot
        tempinput['n_I_ramp_start'][0] = n_I_ramp_start
        tempinput['n_I_ramp_end'][0] = n_I_ramp_end
        tempinput['n_dynamicOn'][0] = n_dynamic
        tempinput['Npar'][0] = N_macro
        tempinput['NperBunch'][0] = nPar
        tempinput['N_bins'][0] = nBins

        tempinput['n_bunches'][0] = nBunch
        tempinput['Prad'][0] = Pbeam*NC
        for i in range(nRF):
            tempinput['QL'][i] = QL[i]
            tempinput['nCav'][i] = nCav[i]
            tempinput['RoQ'][i] = RoQ[i]*nCav[i]

            tempinput['Vref_I'][i] = VrefI[i]*nCav[i]
            tempinput['Vref_Q'][i] = VrefQ[i]*nCav[i]
            tempinput['Iref_I'][i] = I_I[i]
            tempinput['Iref_Q'][i] = I_Q[i]
            tempinput['I_I_ref_ini'][i] = I_I_ini[i]
            tempinput['I_I_ref_final'][i] = I_I[i]
            tempinput['I_Q_ref_ini'][i] = I_Q_ini[i]
            tempinput['I_Q_ref_final'][i] = I_Q[i]
            tempinput['gII'][i] = gII[i]
            tempinput['gQQ'][i] = gQQ[i]
            tempinput['detune'][i] = delta_f_ini[i]
            tempinput['detune_ini'][i] = delta_f_ini[i]
            tempinput['detune_mid'][i] = (delta_f[i]-delta_f_ini[i])/2
            tempinput['detune_final'][i] = delta_f[i]

        tempinput['siglong'][0] = siglong
        tempinput['Ek_damp'][0] = Ek_damp
        tempinput['t_rad_long'][0] = t_rad_long
        tempinput['fill_step'][0] = fill_step
        fn1 = 'input.txt'
        inputfile1 = os.path.join(cwd,fn1)
        with open(inputfile1,'w') as wrt_to_input:
            for i in tempinput:
                wrt_to_input.write(str(i)+' ')
                #print(i)
                for j in range(len(tempinput[i])):
                    wrt_to_input.write(str(tempinput[i][j])+' ')
                    #print(tempinput[i][j])
                wrt_to_input.write('\n')

        args = ("/home/txin/Dropbox/code/Cpp/APES_pack/APES8.2/APES")
        popen = subprocess.Popen(args, stdout=subprocess.PIPE,cwd=cwd)
        popen.wait()
        output = popen.stdout.read()
        print(output)
        outfile = 'out'
        

        path = os.path.join(cwd,"nmacro{0:.0f}".format(N_macro)+"_nBin{0:.0f}".format(nBins)+"_Idc{0:.2f}A".format(nBunch*nPar*f0*1.6e-19)+"_ThetaL{0:.1f}degree".format(ThetaL_max[0]/N_thetaL*thetaL_factor))
        
        try:
            os.mkdir(path)
        except OSError:
            print ("Creation of the directory %s failed" % path)
        else:
            print ("Successfully created the directory %s" % path)
        files = os.listdir(cwd)
        jpgFn = [i for i in files if i[-3:]=='jpg' or i[-3:]=='bin']
        for i in jpgFn:
            path_jpg = os.path.join(cwd,i)
            subprocess.call(["cp",path_jpg,path])
        path_out = os.path.join(cwd,"out")
        subprocess.call(["cp",path_out,path])
        path_in = os.path.join(cwd,"input.txt")
        subprocess.call(["cp",path_in,path])

  
