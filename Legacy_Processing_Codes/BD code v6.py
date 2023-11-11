import numpy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import math
from math import pi
from math import sqrt
import statistics
import scipy.integrate
from scipy.signal import find_peaks
from scipy.signal import argrelmin
from numpy import trapz
from scipy.integrate import cumtrapz

# SETUP VARIABLES - USER INPUTS
BD = 2 #Bluedrop file is from 
fileNum = '02B7' # write the bin file number you want to analyze (do not include 'bLog' or '.bin')
soiltype = "s" #s = sand, c=clay, m=mixed, u=unknown
atype = 'm'  # m = mantle area (best for sands), p = projected area (best for clays)
tiptype = 'c'  # c = cone, p = parabolic, b = blunt
offset = 1 # this value is subtracted from the accelerometer readings
droptype = 'w' #w = water, #a = air
sign = "uk" #enter an effective unit weight value in kg/m^3 or "uk" if unknown 
# paste the filepath to the folder where the BD data is stored
binFilepath = Path("H:\My Drive\CEE 5904 - Project & Report/2021 FRF Data\Bluedrop\October/14 October 2021 AM\BD2 - Pier\Pier")
#paste the filepath to an excel file that the analysis results will be printed in
outputPath = Path("H:\My Drive\CEE 5904 - Project & Report/2021 FRF Data\Bluedrop\October/14 October 2021 AM\Analysis Results-14 October 2021 AM.xlsx") #  Path to pre-existing Excel File
plotPath = Path("H:\My Drive\CEE 5904 - Project & Report/2021 FRF Data\Bluedrop\October/14 October 2021 AM\Analysis Figures- 14 October 2021 AM")
#if applicable, paste the filepath to an excel file that troubleshooting data will be printed in
troubleshootingPath = Path("H:\My Drive\CEE 5904 - Project & Report/2021 FRF Data\Bluedrop\October/14 October 2021 AM\Troubleshooting.xlsx")

# FUNCTIONS FOR USE IN ANALYSIS
def masslength(tiptype): #sets the mass and length of the pentrometer based on the tip
    global mass
    global length
    if tiptype == 'c':
        mass = 7.71 #kg
        length = 8.48-0.93 #originally 7.87, 7.57 for perfect 60 deg, 8.48 measured - .93 measured 90 deg
    elif tiptype == 'p':
        mass = 9.15
        length = 8.26
    elif tiptype == 'b':
        mass = 10.30
        length = 8.57

def dropstartend(peak): #after locating the peaks, this function chops the minute long file into a smaller segment immediately before and after the peak
    global dropstart
    global dropend
    if peak <= 1500:
        dropstart = 1
        dropend = peak + 500
    elif peak > 119500:
        dropstart = peak - 1500
        dropend = 120000
    else:
        dropstart = peak - 1500
        dropend = peak + 500

def accPick(dg, d): #this function picks the smallest accelerometer that's not maxed out to perform the integration on
    maxAcc = dg["250g (g)"].max()
    global acc
    global accName
    global accNameg
    global accg
    if maxAcc < 5 - offset:
        if dg["2g (g)"].max() < 1.8 - offset:  # does an extra check for the 2g because of noise
            acc = d["2g (m/s^2)"]
            accg = dg["2g (g)"]
            accName = "2g (m/s^2)"
            accNameg = "2g (g)"
        else:
            acc = d["18g (m/s^2)"]
            accg = dg["18g (g)"]
            accName = "18g (m/s^2)"
            accNameg = "18g (g)"
    elif maxAcc < 18 - offset:
        acc = d["18g (m/s^2)"]
        accg = dg["18g (g)"]
        accName = "18g (m/s^2)"
        accNameg = "18g (g)"
    elif maxAcc < 50 - offset:
        acc = d["50g (m/s^2)"]
        accg = dg["50g (g)"]
        accName = "50g (m/s^2)"
        accNameg = "50g (g)"
    else:
        acc = d["250g (m/s^2)"]
        accg = dg["50g (g)"]
        accName = "250g (m/s^2)"
        accNameg = "250g (g)"

def findchangepts(): #This function picks the moment that the BD impacts the ground
    global drop
    jlist = list()
    global jindex
    print("finding start of drop...")
    for i in range(4,len(accg)-4):
        p1 = 1
        #print(p1)
        p2 = i
        #print(p2)
        p3 = len(accg)
        #print(p3)
        sample1 = list(accg[p1:p2-1])
        #print(sample1)
        sample2 = list(accg[p2:p3])
        #print(sample2)
        stat1 = math.log(statistics.variance(sample1))
        stat2 = math.log(statistics.variance(sample2))
        #print(stat1)

        j1 = (i-1)*stat1
        j2 = ((len(accg)-1)-i+1)*stat2
        j = j1+j2
        #print(j)
        jlist.append(j)

    drop = min(jlist)
    #print("drop is")
    #print(drop)
    jlist = np.array(jlist)
    #print(jlist)
    #print(jlist.size)
    jlist = np.insert(jlist, 0, (0,0,0,0)) #reshape to match up with dataframe d
    jlist = np.append(jlist, (0,0,0,0)) #reshape to match up with dataframe d
    #print(jlist.size) #should be 2000
    jindex = np.where(jlist==drop) #finds the index of the drop start
    jindex = int(jindex[0]) #converts the index into a number from a tuple

def finddropend(n): #finds the location where the deceleration is 1-offset after the peak
    global num1
    global num2
    below0list = list()
    #for i in range(dropstart+jindex, dropend, 1):
    for i in range(peaksArray[n], dropend, 1):
        if accg[i] < 1 - offset:
            num1 = i - dropstart
            #num2 = i-jindex-1
            below0list = np.append(below0list, num1)
    num1=int(min(below0list))

def integration(d): #integrates the deceleration data to solve for velocity and penetration depth
    global vel
    global maxvel
    global dep
    global maxdep
    accint = acc[jindex:num1]
    print(len(d))
    print(len(accint))
    vel = scipy.integrate.cumtrapz(accint, x=d["Time (s)"])  # solves for velocity
    vel = np.array(vel)
    vel = numpy.insert(vel, 0, 0) #not sure what to insert here, but it makes it the right size
    vel = np.flip(vel)
    maxvel = vel.max()
    maxvel = round(maxvel,1)
    dep = scipy.integrate.cumtrapz(vel, x=d["Time (s)"])  # solves for penetration depth
    dep = numpy.insert(dep, 0, 0)  # not sure what to insert here, but it makes it the right size
    maxdep = dep.max()
    maxdep = round(maxdep,4)
    d.insert(9, "Velocity (m/s)", vel)
    d.insert(10, "Penetration Depth (m)", dep)

def areafind(): #finds the embedded area based on the penetration depth, area type, and the tip
    global area
    a1 = list() #placeholder for the penetrated area at that specific depth
    d1 = dep*100 #penetration depth array, in cm
    #print(len(d1))
    r = list() #placeholder for the radius at that specific depth
    if tiptype == 'c':
        if atype == 'm':
            for k in range(0,len(d1)):
                if d1[k]<length:
                    r.append(.22+d1[k]*((4.375-.22)/7.55))
                    #r.append(d1[k]*np.tan(30*pi/180)) this is from the original matlabcode
                    a1.append(pi*r[k]*(sqrt((r[k]**2)+(d1[k]**2))))
                elif d1[k]>=length:
                    r.append(4.375)
                    a1.append(pi*r[k]*(sqrt((r[k]**2)+(length**2))))
                a1[k] = a1[k]/10000
            area = a1
            #print(r)
            print(np.column_stack((d1,r, area)))
        elif atype == 'p':
            for k in range(0,len(d1)):
                if d1[k]<length:
                    r.append(.22+d1[k]*((4.375-.22)/7.55))
                    # r.append(d1[k]*np.tan(30*pi/180))
                    a1.append(pi*r[k]**2)
                elif d1[k]>=length:
                    r.append(4.375)
                    a1.append(pi*(r[k])**2)
                a1[k] = a1[k]/10000
            area = a1

    elif tiptype == 'b':
        if atype =='m':
            for k in range(0,len(d1)):
                if d1[k]<length:
                    r.append(4.375)
                    a1.append(pi*r[k]**2 + 2*pi*r[k]*d1[k])
                if d1(k)>=length:
                    r.append(4.375)
                    a1.append(pi*r[k]**2 + 2*pi*r[k]*length)
                a1[k]=a1[k]/10000
            area = a1
        elif atype == 'p':
            for k in range(0,len(d1)):
                a1.append(pi*4.375**2)
                a1[k]=a1[k]/10000
            area = a1

    '''elif tiptype == "p":
        if atype == 'm':
            for k in range(1,len(d)):
                if d1[k]<length:
                     r[k]=sqrt(2.4184*d1[k])
                           polarfun = @(theta,r) r.*sqrt(0.745*r.^2 + 1);
                           A1(k)= integral2(polarfun,0,2*pi,0,r(k));
            
        elif atype == 'p':
            for k in range(1,len(d)):'''

def qdynfun(acc): # calculates dynamic bearing capacity
    global qdyn
    global maxbcdep
    global maxqdyn
    global qdyntable

    buoy = 1020*0.002473
    if droptype == "w": #water drops
        Fbe = (mass-buoy)*acc #drop force

    elif droptype =="a": #air drops
        Fbe = mass*acc

    qdyn = (Fbe/area)/1000 #Dynamic bearing capacity (kPa)
    qdyntable = pd.DataFrame()
    qdyntable["Qdyn"] = qdyn
    qdyntable["Depth"] = dep
    #print(qdyntable)

    maxqdyn = qdyn.max() #gives maximum qdyn
    maxbcdep = qdyntable.loc[qdyntable["Qdyn"]==maxqdyn, 'Depth'].item()

    maxqdyn = round(maxqdyn,1)
    maxbcdep = round(maxbcdep,4)
    #print("maxbcdep = ", maxbcdep) #finds the penetration depth corresponding to the maximum qdyn


    """bcmin = argrelmin(qdyn)
    peaksArray = np.array(bcmin)
    print("Peaks Array = ", peaksArray)"""

def qsbcfun(acc): #calculates quasi-static bearing capacity
    global bctable
    global qsbc
    global maxHqsbc
    global maxLqsbc
    global maxAqsbc

    if droptype == "w":
        srcv = np.log10(vel/0.02) #Velocity portion of the strain rate correction. 
        srfk = [0.2, 0.4, 1, 1.5] #list of strain rate factors to run; must be in increasing order

        bctable = pd.DataFrame()
        for i in range(0,len(srfk)):
            fsr = 1+srfk[i]*srcv
            qsbc = qdyn/fsr
            bctable.insert(i, "qsbc for srf = "+str(srfk[i]), qsbc)

        bctable_avgs = pd.DataFrame()
        bctable_avgs["High Average qsbc"]= (bctable.iloc[:,0]+bctable.iloc[:,1])/2 #Average of the lowest 2 strain rate factors
        bctable_avgs["Low Average qsbc"] = (bctable.iloc[:,int(len(srfk))-1]+bctable.iloc[:,int(len(srfk)-2)])/2 #average of the highest 2 strain rate factors
        maxHqsbc = bctable_avgs["High Average qsbc"].max()
        maxHqsbc = round(maxHqsbc, 1)
        maxLqsbc = bctable_avgs["Low Average qsbc"].max()
        maxLqsbc = round(maxLqsbc,1)
        maxAqsbc = "-"

        print(bctable)

    elif droptype == "a":
        srcv = np.log10(vel/0.02) #Velocity portion of the strain rate correction. 
        srfk = 0.31*maxvel
        fsr = 1+srfk*srcv
        qsbc = qdyn/fsr
        maxHqsbc = "-"
        maxLqsbc = "-"
        maxAqsbc = qsbc.max()
        maxAqsbc = round(maxAqsbc,1)

def dr(acc): #Albatal 2019 method for calculating relative density (Duck, NC specific!)
    global Dr
    maxacc = accg.max() #maximum deceleration, in g
    #print("Max deceleration = ", maxacc)
    Dr = -2.18*10**-4*maxacc**3+1.29*10**-2*maxacc**2+1.61*maxacc-13.09 #Albatal 2019
    Dr = round(Dr, 1)
    #print("Relative density is ", dr)

def duncan_correlation(): #duncan correlation for caluclating friction angle
    global phi
    global sign
    #coefficients A, B, C, and D used in the Duncan Correlation. Duncan[34,10,3,2], Albatal[34,10,2,5]
    A = 34
    B = 10
    C = 3
    D = 2
    pa = 101.325 #kPa
    print("max penetration depth = ", maxdep)
    if sign == "uk":
        sign = 1120*maxdep #1120 = 70pcf in kg/m^3

    phi = A + B*Dr/100-(C+D*Dr/100)*math.log(sign/pa)
    phi = round(phi, 1)
    print("phi =", phi)

def firmnessfactor(acc):
    global maxacc
    global ff
    maxacc = acc.max()
    g = 9.80665
    tp = len(drop1)/2000

    ff = maxacc/(g*tp*maxvel)
    ff = round(ff,1)
    print("Firmness factor = ", ff)

#Plots for exporting
def overviewplot(): #Plot showing all accellerometers and pore pressure readings
    fig, (ax1, ax2, ax3) = plt.subplots(3)
    plt.tight_layout()
    ax1.plot(time, g250g, label="250g")
    ax1.plot(time, g50g, label="50g")
    ax1.plot(time, g18g, label="18g")
    ax1.plot(time, g2g, label="2g")
    ax1.legend()
    ax1.set(ylabel="Deceleration (g)")
    ax1.set(xlabel="Time (s)")
    ax1.set_title("BD file "+fileNum)
    ax1.grid()

    ax2.plot(time, ppm, label="Pore Pressure")
    ax2.set(ylabel="Pore Pressure (kPa)")
    ax2.set(xlabel="Time (s)")
    ax2.grid()

    ax3.plot(time, gX55g, label="X 55g")
    ax3.plot(time, gY55g, label="Y 55g")
    ax3.legend()
    ax3.set(ylabel="Deceleration (g)")
    ax3.set(xlabel="Time (s)")
    ax3.grid()
    plotName = fileNum+" Overview.png"
    plt.savefig(plotPath / plotName)
    plt.show()

def integplot(drop, accName): #Deceleration,Velocity,and penetration depth vs time plots
    fig, (ax1, ax2, ax3) = plt.subplots(3)
    plt.tight_layout()
    timereset = drop["Time (s)"]- drop["Time (s)"].iloc[0]
    ax1.plot(timereset, drop[accName], color = "k" ) #,marker = 11
    ax1.set(ylabel="Deceleration (m/s^2)", xlabel="Time(s)")
    ax1.set_title("BD file "+fileNum+ " Drop #"+str(n))
    ax1.grid()

    ax2.plot(timereset, drop['Velocity (m/s)'] , color = "k" )
    ax2.set(ylabel="Velocity (m/s)", xlabel="Time(s)")
    ax2.grid()

    ax3.plot(timereset, drop["Penetration Depth (m)"] , color = "k")
    ax3.set(ylabel="Penetration Depth (m)", xlabel="Time(s)")
    ax3.grid()

    plotName = fileNum+" "+str(n)+" Integ"
    plt.savefig(plotPath / plotName)
    plt.show()

def dep_qsbc_comboplot(dropg, drop, accNameg, qdyn):
    fig, (ax1, ax2) = plt.subplots(1,2)
    plt.tight_layout()
    ax1.plot(dropg[accNameg], drop["Penetration Depth (m)"]*100, color = "k", linestyle = "solid") #marker = 11
    ax1.plot(drop["Velocity (m/s)"], drop["Penetration Depth (m)"]*100, color = "k", linestyle = "dashed")
    ax1.plot(drop["'X55g (m/s^2)'"]), drop["Penetration Depth (m)"]*100, color = "blue", linestyle = "dashed"
    ax1.plot(drop["'Y55g (m/s^2)'"]), drop["Penetration Depth (m)"]*100, color = "orange", linestyle = "dashed"
    ax1.set(xlabel="Deceleration (g) and Velocity (m/s)", ylabel="Penetration Depth (cm)")
    ax1.invert_yaxis()
    ax1.legend(["Deceleration (g)", "Velocity (m/s)"])
    ax1.set_title("BD file "+fileNum)
    ax1.grid()

    if droptype == "a": #for partial saturation, only one qsbc curve
        ax2.plot(qdyn, drop["Penetration Depth (m)"]*100, label="Qdyn") #color = "k", marker = 11 #Plots Qdyn
        ax2.plot(qsbc, drop["Penetration Depth (m)"]*100, label="QSBC") #plots QSBC
        ax2.set(xlabel="Bearing Capacity (kPa)")
        ax2.set_xlim(0,)
        ax2.invert_yaxis()
        ax2.legend(["Qdyn", "QSBC"])
        ax2.set_title("Drop #"+str(n))
        ax2.grid()

    else: #for saturated drops, use a range of strain-rate factors
        ax2.plot(qdyn, drop["Penetration Depth (m)"]*100, label="Qdyn") #color = "k", marker = 11 #Plots Qdyn
        ax2.plot(bctable.iloc[:,0], drop["Penetration Depth (m)"]*100, label=str(bctable.columns[0])) #plots QSBC with smallest srf
        ax2.plot(bctable.iloc[:,len(bctable.columns)-1], drop["Penetration Depth (m)"]*100, label=str(bctable.columns[len(bctable.columns)-1])) #plots QSBC with largest srf
        for i in range(1,len(bctable.columns)-1): #plots qsbc for all other srfs
            ax2.plot(bctable.iloc[:,i], drop["Penetration Depth (m)"]*100, label=str(bctable.columns[i]), color = '.6')
        ax2.set(xlabel="Bearing Capacity (kPa)")
        ax2.set_xlim(0,)
        ax2.invert_yaxis()
        ax2.legend(["Qdyn", str(bctable.columns[0]), str(bctable.columns[len(bctable.columns)-1])])
        ax2.set_title("Drop #"+str(n))
        ax2.grid()


    plotName = fileNum+" "+str(n)+" BC"
    plt.savefig(plotPath / plotName)
    plt.show()

#plots for troubleshooting 
def depthplot(dropg, drop, accNameg): #Velocity and develeration vs. penetration depth 
    fig, (ax1) = plt.subplots(1)
    ax1.plot(dropg[accNameg], drop["Penetration Depth (m)"]*100, color = "k", linestyle = "solid") #marker = 11
    ax1.plot(drop["Velocity (m/s)"], drop["Penetration Depth (m)"]*100, color = "k", linestyle = "dashed")
    ax1.set(xlabel="Deceleration (g) and Velocity (m/s)", ylabel="Penetration Depth (cm)")
    ax1.invert_yaxis()
    ax1.legend(["Deceleration (g)", "Velocity (m/s)"])
    ax1.set_title("BD file "+fileNum)
    plt.show()

def peakplot(): # Plot showing peak deceleration
    peakplot = plt.plot(x)
    peakplot = plt.plot(peaks, x[peaks], "x")
    plt.show()

def qdynplot(drop, qdyn): #Plot showing dynamic bearing capacity vs depth
    #print(len(drop))
    #print(len(qdyn))
    #print((bctable.iloc[:,0]))
    #print(len(bctable.iloc[:,len(bctable.columns)-1]))
    fig, (ax1) = plt.subplots(1)

    if droptype == "a": #for partial saturation, only one qsbc curve
        ax1.plot(qdyn, drop["Penetration Depth (m)"]*100, label="Qdyn") #color = "k", marker = 11 #Plots Qdyn
        ax1.plot(qsbc, drop["Penetration Depth (m)"]*100, label="QSBC") #plots QSBC
        ax1.set(xlabel="Bearing Capacity (kPa)")
        ax1.set_xlim(0,)
        ax1.invert_yaxis()
        ax1.legend(["Qdyn", "QSBC"])
        ax1.set_title("Bearing Capacity- "+fileNum+ " "+str(n))
        ax1.grid()

    else:
        ax1.plot(qdyn, drop["Penetration Depth (m)"]*100, label="Qdyn") #color = "k", marker = 11 #Plots Qdyn
        ax1.plot(bctable.iloc[:,0], drop["Penetration Depth (m)"]*100, label=str(bctable.columns[0])) #plots QSBC with smallest srf
        ax1.plot(bctable.iloc[:,len(bctable.columns)-1], drop["Penetration Depth (m)"]*100, label=str(bctable.columns[len(bctable.columns)-1])) #plots QSBC with largest srf
        for i in range(1,len(bctable.columns)-1): #plots qsbc for all other srfs
            ax1.plot(bctable.iloc[:,i], drop["Penetration Depth (m)"]*100, label=str(bctable.columns[i]), color = "k")
        ax1.set(xlabel="Bearing Capacity (kPa)", ylabel="Penetration Depth (cm)")
        ax1.set_xlim(0,)
        ax1.invert_yaxis()
        ax1.legend(["Qdyn", str(bctable.columns[0]), str(bctable.columns[len(bctable.columns)-1])])
        ax1.set_title("Bearing Capacity- "+fileNum+ " "+str(n))
    
    plt.show()

def tiltplot():
    fig, (ax1) = plt.subplots(1)
    plt.plot(time, gX55g, label="2g")
    plt.plot(time, gY55g, label="2g")
    ax1.legend()
    ax1.set(ylabel="Deceleration (g)")
    ax1.set(xlabel="Time (s)")
    ax1.set_title("BD file "+fileNum)

#Exporting Functions
def exporttoexcel():
    with pd.ExcelWriter(outputPath, mode="a", if_sheet_exists='replace') as writer:
        output_table.to_excel(writer, sheet_name = fileNum, index=False)

def troubleshooting_export():
        with pd.ExcelWriter(troubleshootingPath, mode="a", if_sheet_exists = 'new') as writer:
            drop1.to_excel(writer, sheet_name = str(n), index = False)
            qdyntable.to_excel(writer, sheet_name = "qdyn", index = False)

#Set the penetrometer mass and length
masslength(tiptype)

# READ BD DATA IN
data_array = []  # creates an empty array for us to fill with bd data
fileName = 'Copy of bLog'+fileNum+".bin"
# print(fileName)
newPath = binFilepath / fileName
#print(newPath)
file = open(newPath, 'rb')  # read file
element = file.read(3)  # create a byte list with each element having 3 bytes

while element:
    # Convert to signed integer before adding to data array
    iVAl = int.from_bytes(element, byteorder='big', signed=True)
    data_array.append(iVAl)  # adds the reshaped data from the bd file to the data frame
    element = file.read(3)

np_array = np.array(data_array)  # create numpy array from the list
np_array = np.reshape(np_array, (-1, 10))  # convert the 1d array to 2d array with 10 cols

#print(np_array.shape)
# print(np_array)

df = pd.DataFrame(np_array) # Creates a Dataframe in pandas from the bd data
df.columns = ['Count', 'no clue', 'g2g', 'g18g', 'g50g', 'ppm', 'g200g', 'gX55g', 'gY55g', 'g250g'] # names columns
# print(dfCal)

# APPLY CALIBRATION FACTORS
if BD == 3:  # calibration factors from July 2019
    g2g = (df['g2g']-38285.6)/1615800.9 - offset# accelerometers are in g
    g18g = (df['g18g']+13738)/163516.8 - offset
    g50g = (df['g50g']-238520.6)/63666 - offset
    ppm = ((df['ppm']-139040.1)/20705) * 6.89475729 # converts to kPa
    g200g = ((df['g200g'] +12142.6)/27751.9) - offset
    gX55g = (df['gX55g']-90237)/65351.5  
    gY55g = (df['gY55g']-57464.2)/65545.
    g250g = (df['g250g']-40420.3)/13636.9 - offset

if BD == 2: # calibration factors from Aug 26, 2021
    g2g = (df['g2g']+37242.2)/1639250.2 - offset# accelerometers are in g
    g18g = (df['g18g']-26867.0)/160460.5 - offset
    g50g = (df['g50g']-213923.3)/64080.7- offset
    ppm = ((df['ppm']+55518.9)/18981.7) * 6.89475729 # converts to kPa
    g200g = (df['g200g']-171448.6)/30334.2 - offset
    gX55g = (df['gX55g']-54242.6)/64767.7 
    gY55g = (df['gY55g']-40574.2)/66343.1 
    g250g = (df['g250g']-40614.9)/13654.6 - offset

if BD == 1: # calibration factors from July 2020
    g2g = (df['g2g']-42590.9)/1626361.1 - offset # accelerometers are in g
    g18g = (df['g18g']-44492.9)/161125.5 - offset
    g50g = (df['g50g']-171656.1)/64020.3 - offset
    ppm = ((df['ppm']+31776.1)/20679.7) * 6.89475729 # this is kPa
    g200g = (df['g200g'] -723404.8)/32209.7  - offset
    gX55g = (df['gX55g'] -54881.1)/64858.6 
    gY55g = (df['gY55g']-28735.5)/63839.9 
    g250g = (df['g250g']+13299.7)/13697.1 - offset

time = (df['Count']-df['Count'].iloc[0]+1)/2000 # gives time in s
count = df["Count"]

# make a new dataframe of the calibrated values in units of g
dfCalg = pd.DataFrame([time, g2g, g18g, g50g, g200g, g250g, gX55g, gY55g, ppm])
dfCalg = dfCalg.T
dfCalg.columns = ['Time (s)', '2g (g)', '18g (g)', '50g (g)', '200g (g)', '250g (g)', 'X55g (g)', 'Y55g (g)', 'Pore Pressure (kPa)']  # names columns
#print(dfCalg)

#make a new dataframe of the calibrated values in units of m/s^2
dfCal = pd.DataFrame([time, g2g, g18g, g50g, g200g, g250g, gX55g, gY55g, ppm])
dfCal = dfCal.T
dfCal.columns = ['Time (s)', '2g (m/s^2)', '18g (m/s^2)', '50g (m/s^2)', '200g (m/s^2)', '250g (m/s^2)', 'X55g (m/s^2)', 'Y55g (m/s^2)', 'Pore Pressure (kPa)']  # names columns
dfCal['2g (m/s^2)'] = dfCal['2g (m/s^2)'] * 9.80665
dfCal['18g (m/s^2)'] = dfCal['18g (m/s^2)'] * 9.80665
dfCal['50g (m/s^2)'] = dfCal['50g (m/s^2)'] * 9.80665
dfCal['200g (m/s^2)'] = dfCal['200g (m/s^2)'] * 9.80665
dfCal['250g (m/s^2)'] = dfCal['250g (m/s^2)'] * 9.80665
dfCal['X55g (m/s^2)'] = dfCal['X55g (m/s^2)'] * 9.80665
dfCal['Y55g (m/s^2)'] = dfCal['Y55g (m/s^2)'] * 9.80665
#print(dfCal)

#Locate the drops
x = np.array(g250g)  # what accelerometer to get the peaks from - use 250 because it will never be maxed out
max250 = max(g250g)

peaks, _ = find_peaks(x, height = 2, distance=10000, prominence=3)  # finds the largest peaks more than 2g spaced at least 10000 counts apart

peaksArray = np.array(peaks)  # prints a list of the count where the peaks occur
#peaksArray = filter(filterpeak, peaksArray)

#print(peaksArray)
q = (peaksArray.shape) #gives number of peaks
nDrops = int(q[0]) #number of drops in the file
#print(nDrops)

# For each drop, find the start and end points and integrate to solve for velocity and acceleration
#peakplot()

a = 0 #index variable for drop analysis
n = 1 #index variable for drop analysis

overviewplot()

while n <= nDrops :
    peak = int(peaksArray[a])  # count at the ath drop
    dropstartend(peak) #zooms in the drop file to only consider 500 counts before and 1500 counts after the peak deceleration
    #print(dropstart, dropend)

    if n == 1 :
        drop1 = dfCal[dropstart:dropend]  # start and end points of the drop in m/s^2
        drop1g = dfCalg[dropstart:dropend]  # start and end points of the drop in g
        drop1 = pd.DataFrame(drop1)  # makes dataframe including all data within the start and end points of the drop
        drop1g = pd.DataFrame(drop1g)
        dg = drop1g 
        d = drop1
        accPick(dg, d)  # chooses what accelerometer to use
        acc1 = acc
        acc1Name = accName
        acc1Nameg = accNameg
        #print(num1, num2)
        #print(acc)
        #print(acc.iloc[1])
        findchangepts() #automatically identifies the moment of penetration
        finddropend(n-1) #automatically find the end of the drop
        #print(drop)
        d = d[jindex:num1] #resizes the dataframe to only include data during penetration
        dg = dg[jindex:num1]
        drop1 = d
        drop1g = dg
        integration(d) #double integration to solve for velocity and penetration depth
        drop1 = d
        areafind() #area calculations for bearing capacity calculations
        acc1 = acc1[jindex:num1]
        qdynfun(acc1) #calculation of dynamic bearing capacity
        qsbcfun(acc1) #calculation of quasi-static bearing capacities
        qdyn1 = qdyn
        #qdynplot(drop1, qdyn1)
        print("Max penetration depth is: ", maxdep)
        firmnessfactor(acc1) #calculates the firmness factor
        if soiltype == "s":
            dr(acc1) #calculates relative density
            duncan_correlation() #calculates friction angle
        elif soiltype == "u":
            if maxdep <= .2:
                dr(acc1)
                duncan_correlation()
        pentime=(num1-jindex)/2 #ms #calculates the time of penetration
        integplot(drop1,acc1Name)
        dep_qsbc_comboplot(drop1g,drop1,acc1Nameg, qdyn1)

        print("pentime= ", pentime)

        output_table1 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype, maxdep*100, maxacc, maxvel, maxbcdep*100, maxqdyn, maxHqsbc, maxLqsbc, maxAqsbc, ff, pentime, Dr, phi])
        output_table1 = output_table1.T
        output_table = output_table1
        output_table.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data

        #troubleshooting_export()
 
    if n == 2 :
        drop2 = dfCal[dropstart:dropend]  # start and end points of the drop in m/s^2
        drop2g = dfCalg[dropstart:dropend]  # start and end points of the drop in g
        drop2 = pd.DataFrame(drop2)  # makes dataframe including all data within the start and end points of the drop
        drop2g = pd.DataFrame(drop2g)
        dg = drop2g # chooses what accelerometer to use based on the max g
        d = drop2
        accPick(dg, d)  # chooses what accelerometer to use
        acc2 = acc
        acc2Name = accName
        acc2Nameg = accNameg
        #print(num1, num2)
        #print(acc)
        #print(acc.iloc[1])
        findchangepts()
        finddropend(n-1)
        #print(drop)
        d = d[jindex:num1]
        dg = dg[jindex:num1]
        drop2 = d
        drop2g = dg
       # drop1plot = drop1.plot(y=accName, ylabel="Deceleration (g)", title="drop 1")
        #drop1plot = plt.plot(acc1Name, acc1Name[num1], "x")
        integration(d)
        drop2 = d
        areafind()
        acc2 = acc2[jindex:num1]
        qdynfun(acc2)
        qsbcfun(acc2)
        qdyn2 = qdyn
        print("Max penetration depth is: ", maxdep)
        firmnessfactor(acc2)
        if soiltype == "s":
            dr(acc2)
            duncan_correlation()
        elif soiltype == "u":
            if maxdep <= .2:
                dr(acc2)
                duncan_correlation()
        integplot(drop2,acc2Name)
        dep_qsbc_comboplot(drop2g,drop2,acc2Nameg, qdyn2)
        pentime=(num1-jindex)/2 #ms

        output_table2 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype, maxdep*100, maxacc, maxvel, maxbcdep*100, maxqdyn, maxHqsbc, maxLqsbc, maxAqsbc, ff, pentime, Dr, phi])
        output_table2 = output_table2.T
        output_table2.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data

        output_table = pd.concat([output_table, output_table2])
        
    if n == 3 :
        drop3 = dfCal[dropstart:dropend]  # start and end points of the drop in m/s^2
        drop3g = dfCalg[dropstart:dropend]  # start and end points of the drop in g
        drop3 = pd.DataFrame(drop3)  # makes dataframe including all data within the start and end points of the drop
        drop3g = pd.DataFrame(drop3g)
        dg = drop3g # chooses what accelerometer to use based on the max g
        d = drop3
        accPick(dg, d)  # chooses what accelerometer to use
        acc3 = acc
        acc3Name = accName
        acc3Nameg = accNameg
        findchangepts()
        finddropend(n-1)
        #print(drop)
        d = d[jindex:num1]
        dg = dg[jindex:num1]
        drop3 = d
        drop3g = dg
        integration(d)
        drop3 = d
        areafind()
        acc3 = acc3[jindex:num1]
        qdynfun(acc3)
        qsbcfun(acc3)
        qdyn3 = qdyn
        print("Max penetration depth is: ", maxdep)
        firmnessfactor(acc3)
        if soiltype == "s":
            dr(acc3)
            duncan_correlation()
        elif soiltype == "u":
            if maxdep <= .2:
                dr(acc3)
                duncan_correlation()
        integplot(drop3,acc3Name)
        dep_qsbc_comboplot(drop3g,drop3,acc3Nameg, qdyn3)
        pentime=(num1-jindex)/2 #ms

        output_table3 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype, maxdep*100, maxacc, maxvel, maxbcdep*100, maxqdyn, maxHqsbc, maxLqsbc, maxAqsbc, ff, pentime, Dr, phi])
        output_table3 = output_table3.T
        output_table3.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data


        output_table = pd.concat([output_table, output_table3])
        
    if n == 4 :
        drop4 = dfCal[dropstart:dropend]  # start and end points of the drop in m/s^2
        drop4g = dfCalg[dropstart:dropend]  # start and end points of the drop in g
        drop4 = pd.DataFrame(drop4)  # makes dataframe including all data within the start and end points of the drop
        drop4g = pd.DataFrame(drop4g)
        dg = drop4g # chooses what accelerometer to use based on the max g
        d = drop4
        accPick(dg, d)  # chooses what accelerometer to use
        acc4 = acc
        acc4Name = accName
        acc4Nameg = accNameg
        findchangepts()
        finddropend(n-1)
        d = d[jindex:num1]
        dg = dg[jindex:num1]
        drop4 = d
        drop4g = dg
        integration(d)
        drop4 = d
        areafind()
        acc4 = acc4[jindex:num1]
        qdynfun(acc4)
        qsbcfun(acc4)
        qdyn4 = qdyn
        print("Max penetration depth is: ", maxdep)
        firmnessfactor(acc4)
        if soiltype == "s":
            dr(acc4)
            duncan_correlation()
        elif soiltype == "u":
            if maxdep <= .2:
                dr(acc4)
                duncan_correlation()
        integplot(drop4,acc4Name)
        dep_qsbc_comboplot(drop4g,drop4,acc4Nameg, qdyn4)
        pentime=(num1-jindex)/2 #ms

        output_table4 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype, maxdep*100, maxacc, maxvel, maxbcdep*100, maxqdyn, maxHqsbc, maxLqsbc, maxAqsbc, ff, pentime, Dr, phi])
        output_table4 = output_table4.T
        output_table4.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data

        output_table = pd.concat([output_table, output_table4])
        
    if n == 5 :
        drop5 = dfCal[dropstart:dropend]  # start and end points of the drop in m/s^2
        drop5g = dfCalg[dropstart:dropend]  # start and end points of the drop in g
        drop5 = pd.DataFrame(drop5)  # makes dataframe including all data within the start and end points of the drop
        drop5g = pd.DataFrame(drop5g)
        dg = drop5g # chooses what accelerometer to use based on the max g
        d = drop5
        accPick(dg, d)  # chooses what accelerometer to use
        acc5 = acc
        acc5Name = accName
        acc5Nameg = accNameg
        findchangepts()
        finddropend(n-1)
        d = d[jindex:num1]
        dg = dg[jindex:num1]
        drop5 = d
        drop5g = dg
        integration(d)
        drop5 = d
        areafind()
        acc5 = acc5[jindex:num1]
        qdynfun(acc5)
        qsbcfun(acc5)
        qdyn5 = qdyn
        print("Max penetration depth is: ", maxdep)
        firmnessfactor(acc5)
        if soiltype == "s":
            dr(acc5)
            duncan_correlation()
        elif soiltype == "u":
            if maxdep <= .2:
                dr(acc5)
                duncan_correlation()
        integplot(drop5,acc5Name)
        dep_qsbc_comboplot(drop5g,drop5,acc5Nameg, qdyn5)
        pentime=(num1-jindex)/2 #ms

        output_table5 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype, maxdep*100, maxacc, maxvel, maxbcdep*100, maxqdyn, maxHqsbc, maxLqsbc, maxAqsbc, ff, pentime, Dr, phi])
        output_table5 = output_table5.T
        output_table5.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data

        output_table = pd.concat([output_table, output_table5])
        
    if n == 6 :
        drop6 = dfCal[dropstart:dropend]  # start and end points of the drop in m/s^2
        drop6g = dfCalg[dropstart:dropend]  # start and end points of the drop in g
        drop6 = pd.DataFrame(drop6)  # makes dataframe including all data within the start and end points of the drop
        drop6g = pd.DataFrame(drop6g)
        dg = drop6g # chooses what accelerometer to use based on the max g
        d = drop6
        accPick(dg, d)  # chooses what accelerometer to use
        acc6 = acc
        acc6Name = accName
        acc6Nameg = accNameg
        findchangepts()
        finddropend(n-1)
        d = d[jindex:num1]
        dg = dg[jindex:num1]
        drop6 = d
        drop6g = dg
        integration(d)
        drop6 = d
        areafind()
        acc6 = acc6[jindex:num1]
        qdynfun(acc6)
        qsbcfun(acc6)
        qdyn6 = qdyn
        print("Max penetration depth is: ", maxdep)
        firmnessfactor(acc6)
        if soiltype == "s":
            dr(acc6)
            duncan_correlation()
        elif soiltype == "u":
            if maxdep <= .2:
                dr(acc6)
                duncan_correlation()
        integplot(drop6,acc6Name)
        dep_qsbc_comboplot(drop6g,drop6,acc6Nameg, qdyn6)
        pentime=(num1-jindex)/2 #ms     

        output_table6 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype, maxdep*100, maxacc, maxvel, maxbcdep*100, maxqdyn, maxHqsbc, maxLqsbc, maxAqsbc, ff, pentime, Dr, phi])
        output_table6 = output_table6.T
        output_table6.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data

        output_table = pd.concat([output_table, output_table6])
     
    n = n + 1
    a = a + 1

print(output_table)

exporttoexcel()

print("Run Complete")

# GENERATE PLOTS

#overviewplot()
#peakplot()

'''def printplots():
    if nDrops == 1:
        integplot(drop1,acc1Name)
        depthplot(drop1g,drop1,acc1Nameg)
        qdynplot(drop1, qdyn1)

    elif nDrops ==2:
        integplot(drop1,acc1Name)
        depthplot(drop1g,drop1,acc1Nameg)
        qdynplot(drop1, qdyn1)

        integplot(drop2,acc2Name)
        depthplot(drop2g,drop2,acc2Nameg)
        qdynplot(drop2, qdyn2)

    elif nDrops ==3:
        integplot(drop1,acc1Name)
        depthplot(drop1g,drop1,acc1Nameg)
        qdynplot(drop1, qdyn1)

        integplot(drop2,acc2Name)
        depthplot(drop2g,drop2,acc2Nameg)
        qdynplot(drop2, qdyn2)

        integplot(drop3,acc3Name)
        depthplot(drop3g,drop3,acc3Nameg)
        qdynplot(drop3, qdyn3)

    elif nDrops ==4:
        integplot(drop1,acc1Name)
        depthplot(drop1g,drop1,acc1Nameg)
        qdynplot(drop1, qdyn1)

        integplot(drop2,acc2Name)
        depthplot(drop2g,drop2,acc2Nameg)
        qdynplot(drop2, qdyn2)

        integplot(drop3,acc3Name)
        depthplot(drop3g,drop3,acc3Nameg)
        qdynplot(drop3, qdyn3)

        integplot(drop4,acc4Name)
        depthplot(drop4g,drop4,acc4Nameg)
        qdynplot(drop4, qdyn4)

    elif nDrops ==5:
        integplot(drop1,acc1Name)
        depthplot(drop1g,drop1,acc1Nameg)
        qdynplot(drop1, qdyn1)

        integplot(drop2,acc2Name)
        depthplot(drop2g,drop2,acc2Nameg)
        qdynplot(drop2, qdyn2)

        integplot(drop3,acc3Name)
        depthplot(drop3g,drop3,acc3Nameg)
        qdynplot(drop3, qdyn3)

        integplot(drop4,acc4Name)
        depthplot(drop4g,drop4,acc4Nameg)
        qdynplot(drop4, qdyn4)

        integplot(drop5,acc5Name)
        depthplot(drop5g,drop5,acc5Nameg)
        qdynplot(drop5, qdyn5)

    else:
        integplot(drop1,acc1Name)
        depthplot(drop1g,drop1,acc1Nameg)
        qdynplot(drop1, qdyn1)

        integplot(drop2,acc2Name)
        depthplot(drop2g,drop2,acc2Nameg)
        qdynplot(drop2, qdyn2)

        integplot(drop3,acc3Name)
        depthplot(drop3g,drop3,acc3Nameg)
        qdynplot(drop3, qdyn3)

        integplot(drop4,acc4Name)
        depthplot(drop4g,drop4,acc4Nameg)
        qdynplot(drop4, qdyn4)

        integplot(drop5,acc5Name)
        depthplot(drop5g,drop5,acc5Nameg)
        qdynplot(drop5, qdyn5)

        integplot(drop6,acc6Name)
        depthplot(drop6g,drop6,acc6Nameg)
        qdynplot(drop6, qdyn6)

printplots()'''






