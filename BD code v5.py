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
from numpy import trapz
from scipy.integrate import cumtrapz

#  SETUP VARIABLES - USER INPUTS
project = "Insert Project Name"
station = "Insert Station Name"
date = "Insert date here"
BD = 3
atype = 'p'  # m = mantle area, p = projected area
tiptype = 'c'  # c = cone, p = parabolic, b = blunt
# paste the filepath to the folder where the bd data is stored
binFilepath = Path("H:\My Drive\CEE 5904 - Project & Report\FRF Data/test data")
# write the bin file number you want to analyze (do not include 'bLog' or '.bin')
fileNum = '02F4'
outputFile = 'data.xlsx'  #this currently doesn't do anything, but eventually all data will be printed out into an excel sheet
outputPath = Path("H:\My Drive\CEE 5904 - Project & Report\FRF Data/test data") #  Path for new files
offset = 1 # this value is subtracted from the accelerometer readings
droptype = 'w' #w = water, #a = air

def masslength(tiptype): #sets the mass and length of the pentrometer based on the tip
    global mass
    global length
    if tiptype == 'c':
        mass = 7.71 #kg
        length = 7.87
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

def finddropend(): #finds the location where the deceleration is 1-offset after the peak
    global num1
    global num2
    below0list = list()
    for i in range(dropstart+jindex, dropend, 1):
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
    vel = scipy.integrate.cumtrapz(accint, x=d["Time (s)"])  # solves for velocity
    vel = np.array(vel)
    vel = numpy.insert(vel, 0, 0) #not sure what to insert here, but it makes it the right size
    vel = np.flip(vel)
    maxvel = vel.max()
    dep = scipy.integrate.cumtrapz(vel, x=d["Time (s)"])  # solves for penetration depth
    dep = numpy.insert(dep, 0, 0)  # not sure what to insert here, but it makes it the right size
    maxdep = dep.max()
    d.insert(9, "Velocity (m/s)", vel)
    d.insert(10, "Penetration Depth (m)", dep)

def areafind(): #finds the embedded area based on the penetration depth, area type, and the tip
    global area
    a1 = list() #placeholder for the penetrated area at that specific depth
    d1 = dep*100 #penetration depth array, in cm
    print(len(d1))
    r = list() #placeholder for the radius at that specific depth
    if tiptype == 'c':
        if atype == 'm':
            for k in range(0,len(d1)):
                if d1[k]<length:
                    r.append(d1[k]*np.tan(30*pi/180))
                    a1.append(pi*r[k]*(sqrt((r[k]**2)+(d1[k]**2))))
                elif d1[k]>=length:
                    r.append(4.375)
                    a1.append(pi*r[k]*(sqrt((r[k]^2)+(length^2))))
                a1[k] = a1[k]/10000
            area = a1
        elif atype == 'p':
            for k in range(0,len(d1)):
                if d1[k]<length:
                    r.append(d1[k]*np.tan(30*pi/180))
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
                a1.append(pi*4.375^2)
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

def bc(acc): # calculates dynamic and quasi-staticbearing capacity
    global qdyn
    global bctable

    buoy = 1020*0.002473
    if droptype == "w": #water drops
        Fbe = (mass-buoy)*acc #drop force

    elif droptype =="a": #air drops
        Fbe = mass*acc

    qdyn = (Fbe/area)/1000 #Dynamic bearing capacity (kPa)

    srcv = np.log10(vel/0.02) #Velocity portion of the strain rate correction. 
    srfk = [0.2, 0.4, 1, 1.5] #list of strain rate factors to run

    bctable = pd.DataFrame()
    for i in range(0,len(srfk)):
        fsr = 1+srfk[i]*srcv
        qsbc = qdyn/fsr
        bctable.insert(i, "qsbc for srf = "+str(srfk[i]), qsbc)
    print(bctable)


#Set the penetrometer mass and length
masslength(tiptype)

# READ BD DATA IN
data_array = []  # creates an empty array for us to fill with bd data
fileName = 'bLog'+fileNum+".bin"
# print(fileName)
newPath = binFilepath / fileName
print(newPath)
file = open(newPath, 'rb')  # read file
element = file.read(3)  # create a byte list with each element having 3 bytes

while element:
    # Convert to signed integer before adding to data array
    iVAl = int.from_bytes(element, byteorder='big', signed=True)
    data_array.append(iVAl)  # adds the reshaped data from the bd file to the data frame
    element = file.read(3)

np_array = np.array(data_array)  # create numpy array from the list
np_array = np.reshape(np_array, (-1, 10))  # convert the 1d array to 2d array with 10 cols

print(np_array.shape)
# print(np_array)

df = pd.DataFrame(np_array) # Creates a Dataframe in pandas from the bd data
df.columns = ['Count', 'no clue', 'g2g', 'g18g', 'g50g', 'ppm', 'g200g', 'gX55g', 'gY55g', 'g250g'] # names columns
# print(dfCal)

# APPLY CALIBRATION FACTORS
if BD == 3:  # calibration factors from July 2019
    g2g = (df['g2g']-34426.5)/1615925.8 - offset# accelerometers are in g
    g18g = (df['g18g']+12322.1)/163530.7 - offset
    g50g = (df['g50g']-237384.9)/63651 - 0.1120 - offset
    ppm = ((df['ppm']+62496.7)/20583.0) * 6.89475729 # converts to kPa
    g200g = ((df['g200g'] -248943.7)/39009.4)+0.5518 - offset
    gX55g = (df['gX55g']-59093.7)/66674.3 - offset #check if lateral accelerometers also need to be offset
    gY55g = (df['gY55g']-140224.6)/66674.3- offset
    g250g = (df['g250g']-40536.1)/13631.6 - offset

if BD == 2: # calibration factors from Aug 26, 2021
    g2g = (df['g2g']+31384.7)/1624987.2-0.035 - offset# accelerometers are in g
    g18g = (df['g18g']-26631.0)/159945.4 - offset
    g50g = (df['g50g']+92987.0)/63783.5 - offset
    ppm = ((df['ppm']-35170.6)/12922.9) * 6.89475729 # converts to kPa
    g200g = (df['g200g']-16264.8)/26042.8 -0.277 - offset
    gX55g = (df['gX55g']+89890.3)/63897.1 - offset
    gY55g = (df['gY55g']+14993.0)/64118.0 - offset
    g250g = (df['g250g']+17362.1)/13533.5+0.0656 - offset

if BD == 1: # calibration factors from July 2020
    g2g = (df['g2g']+277743.2)/1637299.6 - offset # accelerometers are in g
    g18g = (df['g18g']-3755.9)/159932.2 - offset
    g50g = (df['g50g']+92817.6)/63237.1 - offset
    ppm = ((df['ppm']-33154.0)/14763.5) * 6.89475729 # this is kPa
    g200g = (df['g200g'] -1155309.9)/28368.5 - 1.464 - offset
    gX55g = (df['gX55g'] +97138.4)/62023.7 - offset
    gY55g = (df['gY55g']-9921.7)/62669.2 - offset
    g250g = (df['g250g']+59211.3)/13276.9 - offset

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
x = np.array(g250g)  # what accelerometer to get the peaks from
peaks, _ = find_peaks(x, height = 2, distance=10000, prominence=3)  # finds the largest peaks more than 2g spaced at least 10000 counts apart

peaksArray = np.array(peaks)  # prints a list of the count where the peaks occur
#print(peaksArray)
q = (peaksArray.shape) #gives number of peaks
nDrops = int(q[0]) #number of drops in the file
#print(nDrops)

# For each drop, find the start and end points and integrate to solve for velocity and acceleration

a = 0
n = 1

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
        findchangepts() #finds the start of the drop
        finddropend() #dinds the end of the drop
        #print(drop)
        d = d[jindex:num1] #shortens the dataframe to only include the data during penetration (jindex = start, num1 = end)
        dg = dg[jindex:num1]
        #print(d)
        #print(np.size(d))
        drop1 = d
        drop1g = dg
        integration(d) #solves for velocity and acceleration
        drop1 = d #this dataframe now includes velocity and acceleration data
        #print(drop1)
        areafind()
        acc1 = acc1[jindex:num1]
        bc(acc1)
        qdyn1 = qdyn


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
        finddropend()
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
        bc(acc2)
        qdyn2 = qdyn

    if n == 3 :
        drop3 = dfCal[dropstart:dropend]  # start and end points of the drop in m/s^2
        drop3g = dfCalg[dropstart:dropend]  # start and end points of the drop in g
        drop3 = pd.DataFrame(drop3)  # makes dataframe including all data within the start and end points of the drop
        drop3g = pd.DataFrame(drop3g)
        dg = drop3g 
        d = drop3
        accPick(dg, d)  # chooses what accelerometer to use
        acc3 = acc
        acc3Name = accName
        acc3Nameg = accNameg
        findchangepts() #finds the start of the drop
        finddropend() #dinds the end of the drop
        d = d[jindex:num1] #shortens the dataframe to only include the data during penetration (jindex = start, num1 = end)
        drop3 = d
        integration(d) #solves for velocity and acceleration
        drop3 = d #this dataframe now includes velocity and acceleration data

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
        #print(num1, num2)
        #print(acc)
        #print(acc.iloc[1])
        findchangepts()
        finddropend()
        #print(drop)
        d = d[jindex:num1]
        #print(np.size(d))
        drop4 = d
       # drop1plot = drop1.plot(y=accName, ylabel="Deceleration (g)", title="drop 1")
        #drop1plot = plt.plot(acc1Name, acc1Name[num1], "x")
        integration(d)
        drop4 = d

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
        finddropend()
        #print(num1, num2)
        #print(acc)
        #print(acc.iloc[1])
        findchangepts()
        #print(drop)
        d = d[jindex:num1]
        drop5 = d
       # drop1plot = drop1.plot(y=accName, ylabel="Deceleration (g)", title="drop 1")
        #drop1plot = plt.plot(acc1Name, acc1Name[num1], "x")
        integration(d)
        drop5 = d

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
        finddropend()
        #print(num1, num2)
        #print(acc)
        #print(acc.iloc[1])
        findchangepts()
        #print(drop)
        d = d[jindex:num1]
        drop6 = d
       # drop1plot = drop1.plot(y=accName, ylabel="Deceleration (g)", title="drop 1")
        #drop1plot = plt.plot(acc1Name, acc1Name[num1], "x")
        integration(d)
        drop6 = d

    n = n + 1
    a = a + 1

# GENERATE PLOTS

def overviewplot(): #Plot showing all accellerometers and pore pressure readings
    fig, (ax1, ax2) = plt.subplots(2)
    ax1.plot(time, g2g, label="2g")
    ax1.plot(time, g18g, label="18g")
    ax1.plot(time, g50g, label="50g")
    #plt.plot(time, ppm)
    #ax1.plot(time, g200g, label="200g")
    #plt.plot(time, gX55g, label="2g")
    #plt.plot(time, gY55g, label="2g")
    ax1.plot(time, g250g, label="250g")
    ax1.legend()
    ax1.set(ylabel="Deceleration (g)")
    ax1.set(xlabel="Time (s)")
    ax1.set_title("BD file "+fileNum)

    ax2.plot(time, ppm, label="Pore Pressure")
    ax2.set(ylabel="Pore Pressure (kPa)")
    ax2.set(xlabel="Time (s)")
    plt.show()

def qdynplot(drop, qdyn): #Plot showing dynamic bearing capacity vs time
    fig, (ax1) = plt.subplots(1)
    ax1.plot(qdyn, drop["Penetration Depth (m)"]*100, label="Qdyn") #color = "k", marker = 11
    ax1.plot(bctable.iloc[:,0], drop["Penetration Depth (m)"]*100, label=str(bctable.columns[0]))
    ax1.plot(bctable.iloc[:,len(bctable.columns)-1], drop["Penetration Depth (m)"]*100, label=str(bctable.columns[len(bctable.columns)-1]))
    for i in range(1,len(bctable.columns)-1):
        ax1.plot(bctable.iloc[:,i], drop["Penetration Depth (m)"]*100, label=str(bctable.columns[i]), color = "k")
    ax1.set(xlabel="Bearing Capacity (kPa)", ylabel="Penetration Depth (cm)")
    ax1.set_xlim(0,)
    ax1.invert_yaxis()
    ax1.legend(["Qdyn", str(bctable.columns[0]), str(bctable.columns[len(bctable.columns)-1])])
    ax1.set_title("Bearing Capacity- "+fileNum+ " "+str(n))
    plt.show()

def qsbcplot(drop):
    fig, (ax1) = plt.subplots(1)

    ax1.invert_yaxis()
    plt.show()


def peakplot(): # Plot showing peak deceleration
    peakplot = plt.plot(x)
    peakplot = plt.plot(peaks, x[peaks], "x")
    plt.show()

def integplot(drop, accName): #Deceleration,Velocity,and penetration depth vs time plots
    fig, (ax1, ax2, ax3) = plt.subplots(3)
    ax1.plot(drop["Time (s)"], drop[accName], color = "k", marker = 11)
    ax1.set(ylabel="Deceleration (m/s^2)", xlabel="Time(s)")
    ax2.plot(drop["Time (s)"], drop['Velocity (m/s)'] , color = "k", marker = 11)
    ax2.set(ylabel="Velocity (m/s)", xlabel="Time(s)")
    ax3.plot(drop["Time (s)"], drop["Penetration Depth (m)"] , color = "k", marker = 11)
    ax3.set(ylabel="Penetration Depth (m)", xlabel="Time(s)")
    plt.show()

def depthplot(dropg, drop, accNameg): #Velocity and develeration vs. penetration depth 
    fig, (ax1) = plt.subplots(1)
    ax1.plot(dropg[accNameg], drop["Penetration Depth (m)"]*100, color = "k", linestyle = "solid") #marker = 11
    ax1.plot(drop["Velocity (m/s)"], drop["Penetration Depth (m)"]*100, color = "k", linestyle = "dashed")
    ax1.set(xlabel="Deceleration (g) and Velocity (m/s)", ylabel="Penetration Depth (cm)")
    ax1.invert_yaxis()
    ax1.legend(["Deceleration (g)", "Velocity (m/s)"])
    plt.show()

'''overviewplot()
peakplot()
integplot(drop1,acc1Name)    
integplot(drop2,acc2Name)
integplot(drop3,acc3Name)
depthplot(drop2g,drop2,acc2Nameg)'''

qdynplot(drop2, qdyn2)









