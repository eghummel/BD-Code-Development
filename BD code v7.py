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
BD = 3 #Bluedrop file is from 
fileNum = '009B' # write the bin file number you want to analyze (do not include 'bLog' or '.bin')
soiltype = "u" #s = sand, c=clay, m=mixed, u=unknown
atype = 'm'  # m = mantle area (best for sands), p = projected area (best for clays), b = both 
tiptype = 'c'  # c = cone, p = parabolic, b = blunt
#offset = 1 # this value is subtracted from the accelerometer readings
droptype = 'a' #w = water, #a = air
sign = "uk" #enter an effective unit weight value in kg/m^3 or "uk" if unknown 
# paste the filepath to the folder where the BD data is stored
binFilepath = Path("H:/.shortcut-targets-by-id/1-AJMoJ7sju0bbgBxYYfsq9uMyeECS_MX/Quinhagak June 2023/Demo data collection/PFFP/June 19, 2023/Drops Only")
#paste the filepath to an excel file that the analysis results will be printed in
outputPath = Path("H:/.shortcut-targets-by-id/1-AJMoJ7sju0bbgBxYYfsq9uMyeECS_MX/Quinhagak June 2023/Demo data collection/PFFP/June 19, 2023/Data Analysis 6-19-23.xlsx") #  Path to pre-existing Excel File
plotPath = Path("H:/.shortcut-targets-by-id/1-AJMoJ7sju0bbgBxYYfsq9uMyeECS_MX/Quinhagak June 2023/Demo data collection/PFFP/June 19, 2023/Analysis Figures")
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
    if maxAcc < 5: #- offset:
        if dg["2g (g)"].max() < 1.8: #offset:  # does an extra check for the 2g because of noise
            acc = dfCal["2g (m/s^2)"]
            accg = dfCalg["2g (g)"]
            accName = "2g (m/s^2)"
            accNameg = "2g (g)"
        else:
            acc = dfCal["18g (m/s^2)"]
            accg = dfCalg["18g (g)"]
            accName = "18g (m/s^2)"
            accNameg = "18g (g)"
    elif maxAcc < 18: #- offset
        acc = dfCal["18g (m/s^2)"]
        accg = dfCalg["18g (g)"]
        accName = "18g (m/s^2)"
        accNameg = "18g (g)"
    elif maxAcc < 50: #- offset:
        acc = dfCal["50g (m/s^2)"]
        accg = dfCalg["50g (g)"]
        accName = "50g (m/s^2)"
        accNameg = "50g (g)"
    else:
        acc = dfCal["250g (m/s^2)"]
        accg = dfCalg["50g (g)"]
        accName = "250g (m/s^2)"
        accNameg = "250g (g)"

def findchangepts(): #This function picks the moment that the BD impacts the ground (not in use currently)
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
    print("Jlist: ", jlist)

def finddropend(n): #finds the location after the peak where the deceleration is 1 
    global penend
    below0list = list()
    #for i in range(dropstart+jindex, dropend, 1):
    for i in range(peak, dropend, 1):
        if accg[i] < 1: #the offset is applied later in the code
            penend = i 
            below0list = np.append(below0list, penend)
    penend=int(min(below0list))

def integration(d, accint): #integrates the deceleration data to solve for velocity and penetration depth
    global vel
    global maxvel
    global dep
    global maxdep
    #accint = acc[penstart:penend]
    #print(len(d))
    #print(len(accint))
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

def peakpick():
    global penstart
    fig, (ax1) = plt.subplots(1)
    plt.plot(dfCalg['Count'], dfCalg[accNameg], label=accNameg)
    plt.plot(dfCalg['Count'], dfCalg['2g (g)'], label='2g')

    ax1.set_xlim(left=dropstart, right=dropend)
    ax1.legend()
    ax1.set(ylabel="Deceleration (g)")
    ax1.set(xlabel="Time")
    ax1.set_title('Double-click start of drop #' + str(n))
    startendpt = plt.ginput(2, 0)
    pentimes = []
    for t in startendpt:
        for x in t:
            pentimes.append(x)

    penstart = int(pentimes[0])
    print("start of penetration: ", penstart)

def CheckingFunction(): # (to be used with auto peak finding) This function is added to check the values of the start and end of the drop if the autofinding function doesn't work

    global jindex
    global num1
    
    xvalues=np.linspace(0,len(dg)-1,len(dg))
    plt.figure(num=1)
    temp_dg=dg['250g (g)']
    plt.plot(xvalues,temp_dg)
    plt.grid(visible=True,which='both')
    plt.show()
    print ('The start of the drop is', jindex)
    print ('The end of the drop is', num1)
    print('Are the start and end of the drop correct?(y/n)')
    drop_input=input()
    if drop_input=='n':   
        peakpick()
        jindex = penstart
        num1 = penend

def areafind(): #finds the embedded area based on the penetration depth, area type, and the tip
    global area
    global trunc_index
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
                    a1.append(pi*(r[k])**2)
                elif d1[k]>=length:
                    r.append(4.375)
                    '''if r[k-1]<r[k]:
                        trunc_index = r.index(r[k-1])+jindex+dropstart #performs analysis only on readings beyond depth of cone height'''
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
    qdyntablet = qdyntable#.truncate(before=trunc_index)
    qdynt = qdyntablet["Qdyn"]
    #print(qdyntable)

    maxqdyn = qdynt.max() #gives maximum qdyn
    maxbcdep = qdyntablet.loc[qdyntable["Qdyn"]==maxqdyn, 'Depth'].item()

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

        #print(bctable)

    elif droptype == "a":
        srcv = np.log10(vel/0.02) #Velocity portion of the strain rate correction. 
        srfk = 0.31*maxvel
        fsr = 1+srfk*srcv
        qsbc = qdyn/fsr
        maxHqsbc = "-"
        maxLqsbc = "-"
        maxAqsbc = qsbc.max()
        maxAqsbc = round(maxAqsbc,1)

def dr(accg): #Albatal 2019 method for calculating relative density (Duck, NC specific!)
    global Dr
    maxacc = accg.max() #maximum deceleration, in g
    #print("Max deceleration = ", maxacc)
    Dr = -2.18*10**-4*maxacc**3+1.29*10**-2*maxacc**2+1.61*maxacc-13.09 #Albatal 2019
    Dr = round(Dr, 1)
    #print("Relative density is ", dr)

def ignore_drop():
    global keep_drop
    print('Keep Drop? (y/n)')
    keep_drop=input()

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

def firmnessfactor(acc, drop):
    global maxacc
    global ff
    maxacc = acc.max()
    g = 9.80665
    tp = len(drop)/2000

    ff = maxacc/(g*tp*maxvel)
    ff = round(ff,1)
    print("Firmness factor = ", ff)

def bctruncate(dropg, drop, acc, qdyn):
    dropg = pd.DataFrame(dropg)
    drop = pd.DataFrame(drop)
    qdyn = pd.DataFrame(qdyn)
    acc = pd.DataFrame(acc)
    dropg = dropg.truncate(before=trunc_index)
    drop = drop.truncate(before=trunc_index)
    acc = acc.truncate(before=trunc_index)
    qdyn = qdyn.truncate(before=trunc_index)

#Plots for exporting
def overviewplot(): #Plot showing all accellerometers and pore pressure readings
    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig.set_size_inches(14, 7)
    #plt.tight_layout()
    ax1.plot(time, g250g, label="250g")
    ax1.plot(time, g50g, label="50g")
    ax1.plot(time, g18g, label="18g")
    ax1.plot(time, g2g, label="2g")
    ax1.legend(loc = "upper right")
    ax1.set(ylabel="Deceleration (g)")
    ax1.set_title("BD file "+fileNum)
    ax1.grid()

    ax2.plot(time, ppm, label="Pore Pressure")
    ax2.set(ylabel="Pore Pressure (kPa)")
    ax2.grid()

    ax3.plot(time, gX55g, label="X 55g")
    ax3.plot(time, gY55g, label="Y 55g")
    ax3.legend(loc = "upper right")
    ax3.set(ylabel="Deceleration (g)")
    ax3.set(xlabel="Time (s)")
    ax3.grid()

    fig.subplots_adjust(bottom=.1, left = .1)
    plotName = fileNum+" Overview.png"
    #plt.savefig(plotPath / plotName)
    plt.show()

def integplot(drop, accName): #Deceleration,Velocity,and penetration depth vs time plots
    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig.set_size_inches(7, 7)
    #resize to 7,7, 7,14
    #plt.tight_layout()
    timereset = drop["Time (s)"]- drop["Time (s)"].iloc[0]
    ax1.plot(timereset, drop[accName], color = "k" ) #,marker = 11
    ax1.set(ylabel="Deceleration (m/s^2)")
    ax1.set_title("BD file "+fileNum+ " Drop #"+str(n))
    ax1.grid()

    ax2.plot(timereset, drop['Velocity (m/s)'] , color = "k" )
    ax2.set(ylabel="Velocity (m/s)")
    ax2.grid()

    ax3.plot(timereset, drop["Penetration Depth (m)"] , color = "k")
    ax3.set(ylabel="Penetration Depth (m)", xlabel="Time(s)")
    ax3.grid()

    fig.subplots_adjust(bottom=.1, left = .1)
    plotName = fileNum+" "+str(n)+" Integ"
    plt.savefig(plotPath / plotName)
    plt.show()

def dep_qsbc_comboplot(dropg, drop, accNameg, qdyn):
    fig, (ax1, ax2) = plt.subplots(1,2)
    fig.set_size_inches(7, 7)
    ax1.plot(dropg[accNameg], drop["Penetration Depth (m)"]*100, color = "k", linestyle = "solid") #marker = 11
    ax1.plot(drop["Velocity (m/s)"], drop["Penetration Depth (m)"]*100, color = "k", linestyle = "dashed")
    ax1.plot(drop["X55g (m/s^2)"], drop["Penetration Depth (m)"]*100, color = "b", linestyle = "dashed")
    ax1.plot(drop["Y55g (m/s^2)"], drop["Penetration Depth (m)"]*100, color = "r", linestyle = "dashed")
    ax1.set(xlabel="Deceleration (g) and Velocity (m/s)", ylabel="Penetration Depth (cm)")
    ax1.invert_yaxis()
    ax1.legend(["Deceleration (g)", "Velocity (m/s)", "X tilt (g)", "Y tilt (g)"], loc = "upper right")
    ax1.set_title("BD file "+fileNum)
    ax1.grid()

    if droptype == "a": #for partial saturation, only one qsbc curve
        ax2.plot(qdyn, drop["Penetration Depth (m)"]*100, label="Qdyn") #color = "k", marker = 11 #Plots Qdyn
        ax2.plot(qsbc, drop["Penetration Depth (m)"]*100, label="QSBC") #plots QSBC
        ax2.set(xlabel="Bearing Capacity (kPa)")
        ax2.set_xlim(0,800)
        ax2.invert_yaxis()
        ax2.legend(["Qdyn", "QSBC"], loc = "upper right")
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
        ax2.legend(["Qdyn", str(bctable.columns[0]), str(bctable.columns[len(bctable.columns)-1])], loc = "upper right")
        ax2.set_title("Drop #"+str(n))
        ax2.grid()
    
    fig.subplots_adjust(bottom=.1, left = .1)
    
def duck_dep_qsbc_comboplot(dropg, drop, dropt, accNameg, qdynt):
    fig, (ax1, ax2) = plt.subplots(1,2)
    fig.set_size_inches(7, 7)
    ax1.plot(dropg[accNameg], drop["Penetration Depth (m)"]*100, color = "k", linestyle = "solid") #marker = 11
    ax1.plot(drop["Velocity (m/s)"], drop["Penetration Depth (m)"]*100, color = "k", linestyle = "dashed")
    #ax1.plot(drop["X55g (m/s^2)"], drop["Penetration Depth (m)"]*100, color = "b", linestyle = "dashed")
    #ax1.plot(drop["Y55g (m/s^2)"], drop["Penetration Depth (m)"]*100, color = "r", linestyle = "dashed")
    ax1.set(xlabel="Deceleration (g) and Velocity (m/s)", ylabel="Penetration Depth (cm)")
    ax1.invert_yaxis()
    ax1.legend(["Deceleration (g)", "Velocity (m/s)"])
    ax1.set_title("BD file "+fileNum)
    ax1.grid()

    if droptype == "a": #for partial saturation, only one qsbc curve
        y = dropt["Penetration Depth (m)"]*100
        ax2.plot(qdynt[10:], y[10:], label="Qdyn") #color = "k", marker = 11 #Plots Qdyn, chops off first 40 readings
        ax2.plot(qsbc[10:-20], y[10:-20], label="QSBC") #plots QSBC; chops off last 10 points which go off towards infinity (can adjust as needed)
        ax2.set(xlabel="Bearing Capacity (kPa)",)
        ax2.set_xlim(0,qdyn.max())
        ax2.invert_yaxis()
        ax2.legend(["Qdyn", "QSBC"], loc = "upper right")
        ax2.set_title("Drop #"+str(n))
        ax2.grid()

    else: #for saturated drops, use a range of strain-rate factors
        y = list(dropt["Penetration Depth (m)"]*100)
        ax2.plot(qdynt[10:-20], y[10:-20], label="Qdyn") #color = "k", marker = 11 #Plots Qdyn
        ax2.plot(bctablet.iloc[10:-20,0], y[10:-20], label=str(bctablet.columns[0])) #plots QSBC with smallest srf
        ax2.plot(bctablet.iloc[10:-20,len(bctable.columns)-1], y[10:-20], label=str(bctable.columns[len(bctable.columns)-1])) #plots QSBC with largest srf
        for i in range(1,len(bctablet.columns)-1): #plots qsbc for all other srfs
            ax2.plot(bctablet.iloc[10:-20,i], y[10:-20], label=str(bctablet.columns[i]), color = '.6')
            #bctablet.iloc[:,i].pop()
        ax2.set(xlabel="Bearing Capacity (kPa)")
        ax2.set_xlim(0,)
        ax2.invert_yaxis()
        ax2.legend(["Qdyn", str(bctablet.columns[0]), str(bctablet.columns[len(bctablet.columns)-1])], loc = "upper right")
        ax2.set_title("Drop #"+str(n))
        ax2.grid()
        

        

    fig.subplots_adjust(bottom=.1, left = .1)
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
fileName = 'bLog'+fileNum+".bin"
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
if BD == 3:  # calibration factors from March 2023
    g2g = (df['g2g']-33570.1)/1614577.9 #- offset# accelerometers are in g
    g18g = (df['g18g']+13495)/163387.2 #- offset
    g50g = (df['g50g']-238817.4)/63779.5 #- offset
    ppm = ((df['ppm']-139040.1)/20705) * 6.89475729 # converts to kPa
    g200g = (df['g200g'] -262332.4)/38888.7 #- offset
    gX55g = (df['gX55g']-70406)/59754.3
    gY55g = (df['gY55g']-69421.1)/141871.5
    g250g = (df['g250g']-39077.4)/13746.9 #- offset

if BD == 2: # calibration factors from March 2023
    g2g = (df['g2g']+36597.6)/1637627.3 #- offset# accelerometers are in g
    g18g = (df['g18g']-26185.3)/160297.2 #- offset
    g50g = (df['g50g']-212256.9)/63968.7 #- offset
    ppm = ((df['ppm']+55518.9)/18981.7) * 6.89475729 # converts to kPa
    g200g = (df['g200g']-175499.4)/30583.8 #- offset
    gX55g = (df['gX55g']-53629.9)/68590.9
    gY55g = (df['gY55g']-43694.3)/68280.3
    g250g = (df['g250g']-39619.9)/13545.8 #- offset

if BD == 1: # calibration factors from March 2023
    g2g = (df['g2g']-43727.6)/1625064 #- offset # accelerometers are in g
    g18g = (df['g18g']-45085.6)/160925.7 #- offset
    g50g = (df['g50g']-173493.4)/63944.6 #- offset
    ppm = ((df['ppm']+31776.1)/20679.7) * 6.89475729 # this is kPa
    g200g = (df['g200g'] -731734.3)/32192.4  #- offset
    gX55g = (df['gX55g'] -68837.9)/52137.3
    gY55g = (df['gY55g']-68015.9)/28074.9
    g250g = (df['g250g']+10580)/13688.7 #- offset

time = (df['Count']-df['Count'].iloc[0]+1)/2000 # gives time in s
count = df["Count"]-df['Count'].iloc[0]

# make a new dataframe of the calibrated values in units of g
dfCalg = pd.DataFrame([count, time, g2g, g18g, g50g, g200g, g250g, gX55g, gY55g, ppm])
dfCalg = dfCalg.T
dfCalg.columns = ['Count', 'Time (s)', '2g (g)', '18g (g)', '50g (g)', '200g (g)', '250g (g)', 'X55g (g)', 'Y55g (g)', 'Pore Pressure (kPa)']  # names columns
#print(dfCalg)

#make a new dataframe of the calibrated values in units of m/s^2
dfCal = pd.DataFrame([count, time, g2g, g18g, g50g, g200g, g250g, gX55g, gY55g, ppm])
dfCal = dfCal.T
dfCal.columns = ['Count','Time (s)', '2g (m/s^2)', '18g (m/s^2)', '50g (m/s^2)', '200g (m/s^2)', '250g (m/s^2)', 'X55g (m/s^2)', 'Y55g (m/s^2)', 'Pore Pressure (kPa)']  # names columns
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

peaks, _ = find_peaks(x, height = 2, distance=1000, prominence=3)  # finds the largest peaks more than 2g spaced at least 10000 counts apart

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

peakplot()
overviewplot()

while n <= nDrops :
    peak = int(peaksArray[a])  # count at the ath drop
    print("peak index:" , peak)
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
        acc1g = accg
        acc1Name = accName
        acc1Nameg = accNameg
        #findchangepts() #automatically identifies the moment of penetration
        peakpick() #interactive graph for user to choose start of penetration
        print("penstart: ", penstart)
        finddropend(n-1) #automatically find the end of the drop
        print("penend: ", penend)
        offset = acc1g[penstart] #sets the offset to be the acceleration at the moment of impact
        print("offset: ", offset)
        #print(drop)
        dg = dfCalg[penstart:penend] #chops the file to only include penetration 
        d = dfCal[penstart:penend]
        acc1g = accg[penstart:penend]
        acc1g = acc1g-offset
        acc1 = acc[penstart:penend]
        offsetms = offset * 9.80665
        acc1 = acc1-offsetms


        drop1 = d
        drop1g = dg
        integration(d, acc1) #double integration to solve for velocity and penetration depth
        drop1 = d
        areafind() #area calculations for bearing capacity calculations
        qdynfun(acc1) #calculation of dynamic bearing capacity
        qsbcfun(acc1) #calculation of quasi-static bearing capacities
        qdyn1 = qdyn
        integplot(drop1,acc1Name)
        #qdynplot(drop1, qdyn1)
        print("Max penetration depth is: ", maxdep)
        print("Pre-truncation d1", drop1)
        drop1 = pd.DataFrame(drop1)
        drop1g = pd.DataFrame(drop1g)
        acc1 = pd.DataFrame(acc1)
        qdyn1 = pd.DataFrame(qdyn1)
        if droptype == "w":
            bctable = pd.DataFrame(bctable)

        drop1t = drop1 #.truncate(before=trunc_index) #this truncates each df such that only depths below cone height are considered
        drop1gt = drop1g #.truncate(before=trunc_index)
        acc1t = acc1 #.truncate(before=trunc_index)
        acc1t = np.array(acc1t)
        acc1gt = acc1g #.truncate(before=trunc_index)
        acc1gt = np.array(acc1gt)
        drop1t = drop1 #truncate(before=trunc_index)
        qdyn1t = qdyn1#.truncate(before=trunc_index)
        if droptype == "w":
            bctablet = bctable#.truncate(before=trunc_index)

        #bctruncate(drop1g, drop1, acc1, qdyn1)
        firmnessfactor(acc1t, drop1) #calculates the firmness factor
        if soiltype == "s":
            dr(acc1gt) #calculates relative density
            duncan_correlation() #calculates friction angle
        elif soiltype == "u":
            if maxdep <= .2:
                dr(acc1gt)
                duncan_correlation()
        pentime=(penend - penstart)/2 #ms #calculates the time of penetration

        dep_qsbc_comboplot(drop1g,drop1,acc1Nameg, qdyn1)
        plt.show()
        ignore_drop()
        if keep_drop == "y":
            duck_dep_qsbc_comboplot(drop1g, drop1, drop1t, acc1Nameg, qdyn1t)
            print("pentime= ", pentime)

            output_table1 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype, maxdep*100, maxacc, maxvel, maxbcdep*100, maxqdyn, maxHqsbc, maxLqsbc, maxAqsbc, ff, pentime, Dr, phi])
            output_table1 = output_table1.T
            output_table = output_table1
            output_table.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data
        else:
            dep_qsbc_comboplot(drop1g,drop1,acc1Nameg, qdyn1)
            plotName = fileNum+" "+str(n)+" BC - REJECTED"
            plt.savefig(plotPath / plotName)
            output_table1 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype,'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'])
            output_table1 = output_table1.T
            output_table = output_table1
            output_table.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data

        #troubleshooting_export()
 
    if n == 2 :
        drop2 = dfCal[dropstart:dropend]  # start and end points of the drop in m/s^2
        drop2g = dfCalg[dropstart:dropend]  # start and end points of the drop in g
        drop2 = pd.DataFrame(drop2)  # makes dataframe including all data within the start and end points of the drop
        drop2g = pd.DataFrame(drop2g)
        dg = drop2g 
        d = drop2
        accPick(dg, d)  # chooses what accelerometer to use
        acc2 = acc
        acc2g = accg
        acc2Name = accName
        acc2Nameg = accNameg
        #findchangepts() #automatically identifies the moment of penetration
        peakpick() #interactive graph for user to choose start of penetration
        print("penstart: ", penstart)
        finddropend(n-1) #automatically find the end of the drop
        print("penend: ", penend)
        offset = acc2g[penstart] #sets the offset to be the acceleration at the moment of impact
        print("offset: ", offset)
        #print(drop)
        dg = dfCalg[penstart:penend] #chops the file to only include penetration 
        d = dfCal[penstart:penend]
        acc2g = accg[penstart:penend]
        acc2g = acc2g-offset
        acc2 = acc[penstart:penend]
        offsetms = offset * 9.80665
        acc2 = acc2-offsetms

        drop2 = d
        drop2g = dg
        integration(d, acc2) #double integration to solve for velocity and penetration depth
        drop2 = d
        areafind() #area calculations for bearing capacity calculations
        qdynfun(acc2) #calculation of dynamic bearing capacity
        qsbcfun(acc2) #calculation of quasi-static bearing capacities
        qdyn2 = qdyn
        integplot(drop2,acc2Name)
        #qdynplot(drop1, qdyn1)
        print("Max penetration depth is: ", maxdep)
        drop2 = pd.DataFrame(drop2)
        drop2g = pd.DataFrame(drop2g)
        acc2 = pd.DataFrame(acc2)
        qdyn2 = pd.DataFrame(qdyn2)
        if droptype == "w":
            bctable = pd.DataFrame(bctable)

        drop2t = drop2 #.truncate(before=trunc_index) #this truncates each df such that only depths below cone height are considered
        drop2gt = drop2g #.truncate(before=trunc_index)
        acc2t = acc2 #.truncate(before=trunc_index)
        acc2t = np.array(acc2t)
        acc2gt = acc2g #.truncate(before=trunc_index)
        acc2gt = np.array(acc2gt)
        drop2t = drop2 #truncate(before=trunc_index)
        qdyn2t = qdyn2#.truncate(before=trunc_index)
        if droptype == "w":
            bctablet = bctable#.truncate(before=trunc_index)

        #bctruncate(drop1g, drop1, acc1, qdyn1)
        firmnessfactor(acc2t, drop2) #calculates the firmness factor
        if soiltype == "s":
            dr(acc2gt) #calculates relative density
            duncan_correlation() #calculates friction angle
        elif soiltype == "u":
            if maxdep <= .2:
                dr(acc2gt)
                duncan_correlation()
        pentime=(penend - penstart)/2 #ms #calculates the time of penetration

        dep_qsbc_comboplot(drop2g,drop2,acc2Nameg, qdyn2)
        plt.show()
        ignore_drop()
        if keep_drop == "y":
            duck_dep_qsbc_comboplot(drop2g, drop2, drop2t, acc2Nameg, qdyn2t)
            print("pentime= ", pentime)

            output_table2 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype, maxdep*100, maxacc, maxvel, maxbcdep*100, maxqdyn, maxHqsbc, maxLqsbc, maxAqsbc, ff, pentime, Dr, phi])
            output_table2 = output_table2.T
            output_table2.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data
        else:
            dep_qsbc_comboplot(drop2g,drop2,acc2Nameg, qdyn2)
            plotName = fileNum+" "+str(n)+" BC - REJECTED"
            plt.savefig(plotPath / plotName)
            output_table2 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype,'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'])
            output_table2 = output_table2.T
            output_table2.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data

        output_table = pd.concat([output_table, output_table2])
        
    if n == 3 :
        drop3 = dfCal[dropstart:dropend]  # start and end points of the drop in m/s^2
        drop3g = dfCalg[dropstart:dropend]  # start and end points of the drop in g
        drop3 = pd.DataFrame(drop3)  # makes dataframe including all data within the start and end points of the drop
        drop3g = pd.DataFrame(drop3g)
        dg = drop3g 
        d = drop3
        accPick(dg, d)  # chooses what accelerometer to use
        acc3 = acc
        acc3g = accg
        acc3Name = accName
        acc3Nameg = accNameg
        #findchangepts() #automatically identifies the moment of penetration
        peakpick() #interactive graph for user to choose start of penetration
        print("penstart: ", penstart)
        finddropend(n-1) #automatically find the end of the drop
        print("penend: ", penend)
        offset = acc3g[penstart] #sets the offset to be the acceleration at the moment of impact
        print("offset: ", offset)
        #print(drop)
        dg = dfCalg[penstart:penend] #chops the file to only include penetration 
        d = dfCal[penstart:penend]
        acc3g = accg[penstart:penend]
        acc3g = acc3g-offset
        acc3 = acc[penstart:penend]
        offsetms = offset * 9.80665
        acc3 = acc3-offsetms

        drop3 = d
        drop3g = dg
        integration(d, acc3) #double integration to solve for velocity and penetration depth
        drop3 = d
        areafind() #area calculations for bearing capacity calculations
        qdynfun(acc3) #calculation of dynamic bearing capacity
        qsbcfun(acc3) #calculation of quasi-static bearing capacities
        qdyn3 = qdyn
        integplot(drop3,acc3Name)
        print("Max penetration depth is: ", maxdep)
        drop3 = pd.DataFrame(drop3)
        drop3g = pd.DataFrame(drop3g)
        acc3 = pd.DataFrame(acc3)
        qdyn3 = pd.DataFrame(qdyn3)
        if droptype == "w":
            bctable = pd.DataFrame(bctable)

        drop3t = drop3 #.truncate(before=trunc_index) #this truncates each df such that only depths below cone height are considered
        drop3gt = drop3g #.truncate(before=trunc_index)
        acc3t = acc3 #.truncate(before=trunc_index)
        acc3t = np.array(acc3t)
        acc3gt = acc3g #.truncate(before=trunc_index)
        acc3gt = np.array(acc3gt)
        drop3t = drop3 #truncate(before=trunc_index)
        qdyn3t = qdyn3#.truncate(before=trunc_index)
        if droptype == "w":
            bctablet = bctable#.truncate(before=trunc_index)

        #bctruncate(drop1g, drop1, acc1, qdyn1)
        firmnessfactor(acc3t, drop3) #calculates the firmness factor
        if soiltype == "s":
            dr(acc3gt) #calculates relative density
            duncan_correlation() #calculates friction angle
        elif soiltype == "u":
            if maxdep <= .2:
                dr(acc3gt)
                duncan_correlation()
        pentime=(penend - penstart)/2 #ms #calculates the time of penetration

        dep_qsbc_comboplot(drop3g,drop3,acc3Nameg, qdyn3)
        plt.show()
        ignore_drop()
        if keep_drop == "y":
            duck_dep_qsbc_comboplot(drop3g, drop3, drop3t, acc3Nameg, qdyn3t)
            print("pentime= ", pentime)

            output_table3 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype, maxdep*100, maxacc, maxvel, maxbcdep*100, maxqdyn, maxHqsbc, maxLqsbc, maxAqsbc, ff, pentime, Dr, phi])
            output_table3 = output_table3.T
            output_table3.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data
        else:
            dep_qsbc_comboplot(drop3g,drop3,acc3Nameg, qdyn3)
            plotName = fileNum+" "+str(n)+" BC - REJECTED"
            plt.savefig(plotPath / plotName)
            output_table3 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype,'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'])
            output_table3 = output_table3.T
            output_table3.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data

        output_table = pd.concat([output_table, output_table3])
        
    if n == 4 :
        drop4 = dfCal[dropstart:dropend]  # start and end points of the drop in m/s^2
        drop4g = dfCalg[dropstart:dropend]  # start and end points of the drop in g
        drop4 = pd.DataFrame(drop4)  # makes dataframe including all data within the start and end points of the drop
        drop4g = pd.DataFrame(drop4g)
        dg = drop4g 
        d = drop4
        accPick(dg, d)  # chooses what accelerometer to use
        acc4 = acc
        acc4g = accg
        acc4Name = accName
        acc4Nameg = accNameg
        #findchangepts() #automatically identifies the moment of penetration
        peakpick() #interactive graph for user to choose start of penetration
        print("penstart: ", penstart)
        finddropend(n-1) #automatically find the end of the drop
        print("penend: ", penend)
        offset = acc4g[penstart] #sets the offset to be the acceleration at the moment of impact
        print("offset: ", offset)
        #print(drop)
        dg = dfCalg[penstart:penend] #chops the file to only include penetration 
        d = dfCal[penstart:penend]
        acc4g = accg[penstart:penend]
        acc4g = acc4g-offset
        acc4 = acc[penstart:penend]
        offsetms = offset * 9.80665
        acc4 = acc4-offsetms

        drop4 = d
        drop4g = dg
        integration(d, acc4) #double integration to solve for velocity and penetration depth
        drop4 = d
        areafind() #area calculations for bearing capacity calculations
        qdynfun(acc4) #calculation of dynamic bearing capacity
        qsbcfun(acc4) #calculation of quasi-static bearing capacities
        qdyn4 = qdyn
        integplot(drop4,acc4Name)
        print("Max penetration depth is: ", maxdep)
        drop4 = pd.DataFrame(drop4)
        drop4g = pd.DataFrame(drop4g)
        acc4 = pd.DataFrame(acc4)
        qdyn4 = pd.DataFrame(qdyn4)
        if droptype == "w":
            bctable = pd.DataFrame(bctable)

        drop4t = drop4 #.truncate(before=trunc_index) #this truncates each df such that only depths below cone height are considered
        drop4gt = drop4g #.truncate(before=trunc_index)
        acc4t = acc4 #.truncate(before=trunc_index)
        acc4t = np.array(acc4t)
        acc4gt = acc4g #.truncate(before=trunc_index)
        acc4gt = np.array(acc4gt)
        drop4t = drop4 #truncate(before=trunc_index)
        qdyn4t = qdyn4#.truncate(before=trunc_index)
        if droptype == "w":
            bctablet = bctable#.truncate(before=trunc_index)

        #bctruncate(drop1g, drop1, acc1, qdyn1)
        firmnessfactor(acc4t, drop4) #calculates the firmness factor
        if soiltype == "s":
            dr(acc4gt) #calculates relative density
            duncan_correlation() #calculates friction angle
        elif soiltype == "u":
            if maxdep <= .2:
                dr(acc4gt)
                duncan_correlation()
        pentime=(penend - penstart)/2 #ms #calculates the time of penetration

        dep_qsbc_comboplot(drop4g,drop4,acc4Nameg, qdyn4)
        plt.show()
        ignore_drop()
        if keep_drop == "y":
            duck_dep_qsbc_comboplot(drop4g, drop4, drop4t, acc4Nameg, qdyn4t)
            print("pentime= ", pentime)

            output_table4 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype, maxdep*100, maxacc, maxvel, maxbcdep*100, maxqdyn, maxHqsbc, maxLqsbc, maxAqsbc, ff, pentime, Dr, phi])
            output_table4 = output_table4.T
            output_table4.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data
        else:
            dep_qsbc_comboplot(drop4g,drop4,acc4Nameg, qdyn4)
            plotName = fileNum+" "+str(n)+" BC - REJECTED"
            plt.savefig(plotPath / plotName)
            output_table4 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype,'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'])
            output_table4 = output_table4.T
            output_table4.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data

        output_table = pd.concat([output_table, output_table4])
        
    if n == 5 :
        drop5 = dfCal[dropstart:dropend]  # start and end points of the drop in m/s^2
        drop5g = dfCalg[dropstart:dropend]  # start and end points of the drop in g
        drop5 = pd.DataFrame(drop5)  # makes dataframe including all data within the start and end points of the drop
        drop5g = pd.DataFrame(drop5g)
        dg = drop5g 
        d = drop5
        accPick(dg, d)  # chooses what accelerometer to use
        acc5 = acc
        acc5g = accg
        acc5Name = accName
        acc5Nameg = accNameg
        #findchangepts() #automatically identifies the moment of penetration
        peakpick() #interactive graph for user to choose start of penetration
        print("penstart: ", penstart)
        finddropend(n-1) #automatically find the end of the drop
        print("penend: ", penend)
        offset = acc5g[penstart] #sets the offset to be the acceleration at the moment of impact
        print("offset: ", offset)
        #print(drop)
        dg = dfCalg[penstart:penend] #chops the file to only include penetration 
        d = dfCal[penstart:penend]
        acc5g = accg[penstart:penend]
        acc5g = acc5g-offset
        acc5 = acc[penstart:penend]
        offsetms = offset * 9.80665
        acc5 = acc5-offsetms

        drop5 = d
        drop5g = dg
        integration(d, acc5) #double integration to solve for velocity and penetration depth
        drop5 = d
        areafind() #area calculations for bearing capacity calculations
        qdynfun(acc5) #calculation of dynamic bearing capacity
        qsbcfun(acc5) #calculation of quasi-static bearing capacities
        qdyn5 = qdyn
        integplot(drop5,acc5Name)
        print("Max penetration depth is: ", maxdep)
        drop5 = pd.DataFrame(drop5)
        drop5g = pd.DataFrame(drop5g)
        acc5 = pd.DataFrame(acc5)
        qdyn5 = pd.DataFrame(qdyn5)
        if droptype == "w":
            bctable = pd.DataFrame(bctable)

        drop5t = drop5 #.truncate(before=trunc_index) #this truncates each df such that only depths below cone height are considered
        drop5gt = drop5g #.truncate(before=trunc_index)
        acc5t = acc5 #.truncate(before=trunc_index)
        acc5t = np.array(acc5t)
        acc5gt = acc5g #.truncate(before=trunc_index)
        acc5gt = np.array(acc5gt)
        drop5t = drop5 #truncate(before=trunc_index)
        qdyn5t = qdyn5#.truncate(before=trunc_index)
        if droptype == "w":
            bctablet = bctable#.truncate(before=trunc_index)

        #bctruncate(drop1g, drop1, acc1, qdyn1)
        firmnessfactor(acc5t, drop5) #calculates the firmness factor
        if soiltype == "s":
            dr(acc5gt) #calculates relative density
            duncan_correlation() #calculates friction angle
        elif soiltype == "u":
            if maxdep <= .2:
                dr(acc5gt)
                duncan_correlation()
        pentime=(penend - penstart)/2 #ms #calculates the time of penetration

        dep_qsbc_comboplot(drop5g,drop5,acc5Nameg, qdyn5)
        plt.show()
        ignore_drop()
        if keep_drop == "y":
            duck_dep_qsbc_comboplot(drop5g, drop5, drop5t, acc5Nameg, qdyn5t)
            print("pentime= ", pentime)

            output_table5 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype, maxdep*100, maxacc, maxvel, maxbcdep*100, maxqdyn, maxHqsbc, maxLqsbc, maxAqsbc, ff, pentime, Dr, phi])
            output_table5 = output_table5.T
            output_table5.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data
        else:
            dep_qsbc_comboplot(drop5g,drop5,acc5Nameg, qdyn5)
            plotName = fileNum+" "+str(n)+" BC - REJECTED"
            plt.savefig(plotPath / plotName)
            output_table5 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype,'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'])
            output_table5 = output_table5.T
            output_table5.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data

        output_table = pd.concat([output_table, output_table5])
        
    if n == 6 :
        drop6 = dfCal[dropstart:dropend]  # start and end points of the drop in m/s^2
        drop6g = dfCalg[dropstart:dropend]  # start and end points of the drop in g
        drop6 = pd.DataFrame(drop6)  # makes dataframe including all data within the start and end points of the drop
        drop6g = pd.DataFrame(drop6g)
        dg = drop6g 
        d = drop6
        accPick(dg, d)  # chooses what accelerometer to use
        acc6 = acc
        acc6g = accg
        acc6Name = accName
        acc6Nameg = accNameg
        #findchangepts() #automatically identifies the moment of penetration
        peakpick() #interactive graph for user to choose start of penetration
        print("penstart: ", penstart)
        finddropend(n-1) #automatically find the end of the drop
        print("penend: ", penend)
        offset = acc6g[penstart] #sets the offset to be the acceleration at the moment of impact
        print("offset: ", offset)
        #print(drop)
        dg = dfCalg[penstart:penend] #chops the file to only include penetration 
        d = dfCal[penstart:penend]
        acc6g = accg[penstart:penend]
        acc6g = acc6g-offset
        acc6 = acc[penstart:penend]
        offsetms = offset * 9.80665
        acc6 = acc6-offsetms

        drop6 = d
        drop6g = dg
        integration(d, acc6) #double integration to solve for velocity and penetration depth
        drop6 = d
        areafind() #area calculations for bearing capacity calculations
        qdynfun(acc6) #calculation of dynamic bearing capacity
        qsbcfun(acc6) #calculation of quasi-static bearing capacities
        qdyn6 = qdyn
        integplot(drop6,acc6Name)
        print("Max penetration depth is: ", maxdep)
        drop6 = pd.DataFrame(drop6)
        drop6g = pd.DataFrame(drop6g)
        acc6 = pd.DataFrame(acc6)
        qdyn6 = pd.DataFrame(qdyn6)
        if droptype == "w":
            bctable = pd.DataFrame(bctable)

        drop6t = drop6 #.truncate(before=trunc_index) #this truncates each df such that only depths below cone height are considered
        drop6gt = drop6g #.truncate(before=trunc_index)
        acc6t = acc6 #.truncate(before=trunc_index)
        acc6t = np.array(acc6t)
        acc6gt = acc6g #.truncate(before=trunc_index)
        acc6gt = np.array(acc6gt)
        drop6t = drop6 #truncate(before=trunc_index)
        qdyn6t = qdyn6#.truncate(before=trunc_index)
        if droptype == "w":
            bctablet = bctable#.truncate(before=trunc_index)

        #bctruncate(drop1g, drop1, acc1, qdyn1)
        firmnessfactor(acc6t, drop6) #calculates the firmness factor
        if soiltype == "s":
            dr(acc6gt) #calculates relative density
            duncan_correlation() #calculates friction angle
        elif soiltype == "u":
            if maxdep <= .2:
                dr(acc6gt)
                duncan_correlation()
        pentime=(penend - penstart)/2 #ms #calculates the time of penetration

        dep_qsbc_comboplot(drop6g,drop6,acc6Nameg, qdyn6)
        plt.show()
        ignore_drop()
        if keep_drop == "y":
            duck_dep_qsbc_comboplot(drop6g, drop6, drop6t, acc6Nameg, qdyn6t)
            print("pentime= ", pentime)

            output_table6 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype, maxdep*100, maxacc, maxvel, maxbcdep*100, maxqdyn, maxHqsbc, maxLqsbc, maxAqsbc, ff, pentime, Dr, phi])
            output_table6 = output_table6.T
            output_table6.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data
        else:
            dep_qsbc_comboplot(drop6g,drop6,acc6Nameg, qdyn6)
            plotName = fileNum+" "+str(n)+" BC - REJECTED"
            plt.savefig(plotPath / plotName)
            output_table6 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype,'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'])
            output_table6 = output_table6.T
            output_table6.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data

        output_table = pd.concat([output_table, output_table6])

    if n == 7 :
        drop7 = dfCal[dropstart:dropend]  # start and end points of the drop in m/s^2
        drop7g = dfCalg[dropstart:dropend]  # start and end points of the drop in g
        drop7 = pd.DataFrame(drop7)  # makes dataframe including all data within the start and end points of the drop
        drop7g = pd.DataFrame(drop7g)
        dg = drop7g 
        d = drop7
        accPick(dg, d)  # chooses what accelerometer to use
        acc7 = acc
        acc7g = accg
        acc7Name = accName
        acc7Nameg = accNameg
        #findchangepts() #automatically identifies the moment of penetration
        peakpick() #interactive graph for user to choose start of penetration
        print("penstart: ", penstart)
        finddropend(n-1) #automatically find the end of the drop
        print("penend: ", penend)
        offset = acc7g[penstart] #sets the offset to be the acceleration at the moment of impact
        print("offset: ", offset)
        #print(drop)
        dg = dfCalg[penstart:penend] #chops the file to only include penetration 
        d = dfCal[penstart:penend]
        acc7g = accg[penstart:penend]
        acc7g = acc7g-offset
        acc7 = acc[penstart:penend]
        offsetms = offset * 9.80665
        acc7 = acc7-offsetms

        drop7 = d
        drop7g = dg
        integration(d, acc7) #double integration to solve for velocity and penetration depth
        drop7 = d
        areafind() #area calculations for bearing capacity calculations
        qdynfun(acc7) #calculation of dynamic bearing capacity
        qsbcfun(acc7) #calculation of quasi-static bearing capacities
        qdyn7 = qdyn
        integplot(drop7,acc7Name)
        print("Max penetration depth is: ", maxdep)
        drop7 = pd.DataFrame(drop7)
        drop7g = pd.DataFrame(drop7g)
        acc7 = pd.DataFrame(acc7)
        qdyn7 = pd.DataFrame(qdyn7)
        if droptype == "w":
            bctable = pd.DataFrame(bctable)

        drop7t = drop7 #.truncate(before=trunc_index) #this truncates each df such that only depths below cone height are considered
        drop7gt = drop7g #.truncate(before=trunc_index)
        acc7t = acc7 #.truncate(before=trunc_index)
        acc7t = np.array(acc7t)
        acc7gt = acc7g #.truncate(before=trunc_index)
        acc7gt = np.array(acc7gt)
        drop7t = drop7 #truncate(before=trunc_index)
        qdyn7t = qdyn7#.truncate(before=trunc_index)
        if droptype == "w":
            bctablet = bctable#.truncate(before=trunc_index)

        #bctruncate(drop1g, drop1, acc1, qdyn1)
        firmnessfactor(acc7t, drop7) #calculates the firmness factor
        if soiltype == "s":
            dr(acc7gt) #calculates relative density
            duncan_correlation() #calculates friction angle
        elif soiltype == "u":
            if maxdep <= .2:
                dr(acc7gt)
                duncan_correlation()
        pentime=(penend - penstart)/2 #ms #calculates the time of penetration

        dep_qsbc_comboplot(drop7g,drop7,acc7Nameg, qdyn7)
        plt.show()
        ignore_drop()
        if keep_drop == "y":
            duck_dep_qsbc_comboplot(drop7g, drop7, drop7t, acc7Nameg, qdyn7t)
            print("pentime= ", pentime)

            output_table7 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype, maxdep*100, maxacc, maxvel, maxbcdep*100, maxqdyn, maxHqsbc, maxLqsbc, maxAqsbc, ff, pentime, Dr, phi])
            output_table7 = output_table7.T
            output_table7.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data
        else:
            dep_qsbc_comboplot(drop7g,drop7,acc7Nameg, qdyn7)
            plotName = fileNum+" "+str(n)+" BC - REJECTED"
            plt.savefig(plotPath / plotName)
            output_table7 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype,'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'])
            output_table7 = output_table7.T
            output_table7.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data

        output_table = pd.concat([output_table, output_table7])

    if n == 8 :
        drop8 = dfCal[dropstart:dropend]  # start and end points of the drop in m/s^2
        drop8g = dfCalg[dropstart:dropend]  # start and end points of the drop in g
        drop8 = pd.DataFrame(drop8)  # makes dataframe including all data within the start and end points of the drop
        drop8g = pd.DataFrame(drop8g)
        dg = drop8g 
        d = drop8
        accPick(dg, d)  # chooses what accelerometer to use
        acc8 = acc
        acc8g = accg
        acc8Name = accName
        acc8Nameg = accNameg
        #findchangepts() #automatically identifies the moment of penetration
        peakpick() #interactive graph for user to choose start of penetration
        print("penstart: ", penstart)
        finddropend(n-1) #automatically find the end of the drop
        print("penend: ", penend)
        offset = acc8g[penstart] #sets the offset to be the acceleration at the moment of impact
        print("offset: ", offset)
        #print(drop)
        dg = dfCalg[penstart:penend] #chops the file to only include penetration 
        d = dfCal[penstart:penend]
        acc8g = accg[penstart:penend]
        acc8g = acc8g-offset
        acc8 = acc[penstart:penend]
        offsetms = offset * 9.80665
        acc8 = acc8-offsetms

        drop8 = d
        drop8g = dg
        integration(d, acc8) #double integration to solve for velocity and penetration depth
        drop8 = d
        areafind() #area calculations for bearing capacity calculations
        qdynfun(acc8) #calculation of dynamic bearing capacity
        qsbcfun(acc8) #calculation of quasi-static bearing capacities
        qdyn8 = qdyn
        integplot(drop8,acc8Name)
        print("Max penetration depth is: ", maxdep)
        drop8 = pd.DataFrame(drop8)
        drop8g = pd.DataFrame(drop8g)
        acc8 = pd.DataFrame(acc8)
        qdyn8 = pd.DataFrame(qdyn8)
        if droptype == "w":
            bctable = pd.DataFrame(bctable)

        drop8t = drop8 #.truncate(before=trunc_index) #this truncates each df such that only depths below cone height are considered
        drop8gt = drop8g #.truncate(before=trunc_index)
        acc8t = acc8 #.truncate(before=trunc_index)
        acc8t = np.array(acc8t)
        acc8gt = acc8g #.truncate(before=trunc_index)
        acc8gt = np.array(acc8gt)
        drop8t = drop8 #truncate(before=trunc_index)
        qdyn8t = qdyn8#.truncate(before=trunc_index)
        if droptype == "w":
            bctablet = bctable#.truncate(before=trunc_index)

        #bctruncate(drop1g, drop1, acc1, qdyn1)
        firmnessfactor(acc8t, drop8) #calculates the firmness factor
        if soiltype == "s":
            dr(acc8gt) #calculates relative density
            duncan_correlation() #calculates friction angle
        elif soiltype == "u":
            if maxdep <= .2:
                dr(acc8gt)
                duncan_correlation()
        pentime=(penend - penstart)/2 #ms #calculates the time of penetration

        dep_qsbc_comboplot(drop8g,drop8,acc8Nameg, qdyn8)
        plt.show()
        ignore_drop()
        if keep_drop == "y":
            duck_dep_qsbc_comboplot(drop8g, drop8, drop8t, acc8Nameg, qdyn8t)
            print("pentime= ", pentime)

            output_table8 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype, maxdep*100, maxacc, maxvel, maxbcdep*100, maxqdyn, maxHqsbc, maxLqsbc, maxAqsbc, ff, pentime, Dr, phi])
            output_table8 = output_table8.T
            output_table8.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data
        else:
            dep_qsbc_comboplot(drop8g,drop8,acc8Nameg, qdyn8)
            plotName = fileNum+" "+str(n)+" BC - REJECTED"
            plt.savefig(plotPath / plotName)
            output_table8 = pd.DataFrame([BD, fileNum, n, soiltype, atype, tiptype, droptype,'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN', 'NaN'])
            output_table8 = output_table8.T
            output_table8.columns = ["BlueDrop", "BD File", "Drop No.", "Soil Type", "Area Type", "Tip Type", "Drop Type", "Penetration Depth (cm)", "Maximum Deceleration (m/s^2)", "Impact Velocity (m/s)", "Depth of Max QSBC (cm)", "Max Qdyn (kPa)", "Max QSBC-Upper Bounds (kPa)", "Max QSBC-Lower Bounds (kPa)", "Max QSBC-Average (kPa)", "Firmness Factor (m-1)", "Penetration Time (ms)", "Relative Density (%)", "Friction Angle (deg.)"]  # Creates an empty Dataframe in pandas from the output data

        output_table = pd.concat([output_table, output_table8])
    
    n = n + 1
    a = a + 1
    print(output_table)

print(output_table)

exporttoexcel()

print("Run Complete")



'''

VARIABLE NAME INDEX

fileNum = bin file number you want to analyze (do not include 'bLog' or '.bin')
soiltype = type of soil drop is in; #s = sand, c=clay, m=mixed, u=unknown
atype = area type used for force calculation; # m = mantle area (best for sands), p = projected area (best for clays)
tiptype = tip used on penetrometer;  # c = cone, p = parabolic, b = blunt
offset = offset from deceleration readings in g such that the moment at impact is at 0g
droptype = whether the drop is in water or air; #w = water, #a = air
sign = effective unit weight value in kg/m^3 or "uk" if unknown '''

'''
HISTORICAL CALIBRATION FACTORS

blueDrop #1 calibration factors:
- for data collected from July 2020 onward, use the following calibration factors:
            g2g=(df['g2g']-42590.9)/1626361.1; %accelerometers are in g
            g18g=(df['g18g']-44492.9)/161125.5;
            g50g=(df['g50g']-171656.1)/64020.3;
            ppm=(df['ppm']+31776.1)/20679.7;  % this is psi
            g200g=(df['g200g']-723404.8)/32209.7);
            gX55g=(df['gX55g']-54881.1)/64858.6);
            gY55g=(df['gY55g']-28735.5)/63839.9;
            g250g=(df['g250g']+13299.7)/13697.1;
            ppm=ppm*6.89475729; % convert into kPa

blueDrop #2 calibration factors:

- for data collected between July 2019 and September 2021, use the following calibration factors:
   g2g=((double(F(:,3)))+31683)/1624691; %accelerometers are in g
            g18g=((double(F(:,4)))-27002.4)/159942.4;
            g50g=((double(F(:,5)))+91943.2)/63865;
            ppm=((double(F(:,6)))+37971.1)/12982;  % this is psi
            g200g=((double(F(:,7)))-17936.6)/26226.6;
            gX55g=((double(F(:,8)))+89249)/66161.6;
            gY55g=((double(F(:,9)))+11088.8)/64773.6;
            g250g=((double(F(:,10)))+18888.6)/13518.9;
            ppm=ppm*6.89475729; % convert into kPa

- for data collected following September 2021, use the following calibration factors:
   g2g=((double(F(:,3)))+37242.2)/1639250.2; %accelerometers are in g
            g18g=((double(F(:,4)))-26867.0)/160460.5;
            g50g=((double(F(:,5)))-213923.3)/64080.7;
            ppm=((double(F(:,6)))+55518.9)/18981.7;  % this is psi
            g200g=((double(F(:,7)))-171448.6)/30334.2;
            gX55g=((double(F(:,8)))-54242.6)/64767.7;
            gY55g=((double(F(:,9)))-40574.2)/66343.1;
            g250g=((double(F(:,10)))-40614.9)/13654.6;

blueDrop #3 calibration factors:

- for data collected before July 2019, use the following calibration factors:
   g2g=((double(F(:,3)))-38285.6)/1615800.9; %accelerometers are in g
   g18g=((double(F(:,4)))+13738)/163516.8;
     g50g=((double(F(:,5)))-238520.6)/63666.0;
   ppm=((double(F(:,6)))-139040.1)/20705.0;  % this is psi
   g200g=((double(F(:,7)))+12142.6)/27751.9;
   gX55g=((double(F(:,8)))-90237.0)/65351.5;
   gY55g=((double(F(:,9)))-57464.2)/65545.5;
   g250g=((double(F(:,10)))-40420.3)/13636.9;

- for data collected after July 2019, use the following calibration factors:
   g2g=((double(F(:,3)))-38285.6)/1615800.9; %accelerometers are in g
            g18g=((double(F(:,4)))+13738)/163516.8;
            g50g=((double(F(:,5)))-238520.6)/63666;
            ppm=((double(F(:,6)))-139040.1)/20705;  % this is psi
            g200g=(((double(F(:,7)))+12142.6)/27751.9);
            gX55g=((double(F(:,8)))-90237)/65351.5;
            gY55g=((double(F(:,9)))-57464.2)/65545.5;
            g250g=((double(F(:,10)))-40420.3)/13636.9;
            g2g(end) = []; g2g = [1;g2g];
            g200g(end) = []; g200g = [1;g200g];
            g18g(end) = []; g18g = [1;g18g];
            g200g(end) = []; g200g = [1;g200g];
            ppm=ppm*6.89475729; % convert into kPa
'''




