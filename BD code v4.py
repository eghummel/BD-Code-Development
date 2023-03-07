import numpy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

import scipy.integrate
from scipy.signal import find_peaks
from numpy import trapz
from scipy.integrate import cumtrapz


#  SETUP VARIABLES - USER INPUTS
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

#def tellme(s):
 #   #print(s)
  #  ax1.set(title=s)
   # plt.draw()

def accPick(dg, d):
    # each of the following are the same; if an accelerometer breaks on a BD, can edit that section
    # the 200g accelerometer is ignored in all cases
    maxAcc = dg["250g (g)"].max()
    global acc
    global accName
    global accNameg
    if maxAcc < 5 - offset:
        if dg["2g (g)"].max() < 1.8 - offset:  # does an extra check for the 2g because of noise
            acc = d["2g (m/s^2)"]
            accName = "2g (m/s^2)"
            accNameg = "2g (g)"
        else:
            acc = d["18g (m/s^2)"]
            accName = "18g (m/s^2)"
            accNameg = "18g (g)"
    elif maxAcc < 18 - offset:
        acc = d["18g (m/s^2)"]
        accName = "18g (m/s^2)"
        accNameg = "18g (g)"
    elif maxAcc < 50 - offset:
        acc = d["50g (m/s^2)"]
        accName = "50g (m/s^2)"
        accNameg = "50g (g)"
    else:
        acc = d["250g (m/s^2)"]
        accName = "250g (m/s^2)"
        accNameg = "250g (g)"

'''def peakpicker(dnum, dse):
    global acc
    #global dnum
    #global dse
    pts = plt.ginput(2)
    (x0, y0), (x1, y1) = pts
    xmin, xmax = sorted([x0, x1])
    dropTimes[dnum, dse] = xmin
    xmin = round(xmin)
    print(xmin)
    dse = dse + 1
    dropTimes[dnum, dse] = xmax
    xmax = round(xmax)
    print(xmax)
    print(dropTimes)
    d = dfCal[xmin:xmax]
    acc = d[accName]
    print(d)'''

def integration(d):
    global vel
    global maxvel
    global dep
    global maxdep
    vel = scipy.integrate.cumtrapz(acc, x=d["Time (s)"])  # solves for velocity
    vel = np.array(vel)
    vel = numpy.insert(vel, 0, 0) #not sure what to insert here, but it makes it the right size
    vel = np.flip(vel)
    maxvel = vel.max()
    dep = scipy.integrate.cumtrapz(vel, x=d["Time (s)"])  # solves for penetration depth
    dep = numpy.insert(dep, 0, 0)  # not sure what to insert here, but it makes it the right size
    maxdep = dep.max()
    d.insert(9, "Velocity (m/s)", vel)
    d.insert(10, "Penetration Depth (m)", dep)



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
    gY55g = (df['gY55g']-140224.6)/66674.3 - offset
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
print(dfCalg)


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
print(dfCal)



#AUTOMATIC PEAK FINDING
x = np.array(g250g)  # what accelerometer to get the peaks from
peaks, _ = find_peaks(x, height = 2, distance=10000)  # finds the largest peaks more than 5g spaced at least 10000 counts apart

# CREATION OF INDIVIDUAL DROP FILES
peaksArray = np.array(peaks)  # prints a list of the count where the peaks occur
print(peaksArray)
q = (peaksArray.shape) #gives number of peaks
nDrops = int(q[0]) #number of drops in the file
print(nDrops)

a = 0
n = 1

dropTimes = np.array([[0]*2]*nDrops)
dnum = 0 #row index for drop; 0 = drop 1, etc
dse = 0 #column index for drop 0 = start, 1 = end

    # MAKE INDIVIDUAL DATAFRAMES FOR EACH DROP AND PLOT DECELERATION VS TIME
while n <= nDrops :
    b = int(peaksArray[a])  # count at the ath drop
    dropstart = b - 100  # offset in counts before impact to include in sliced df
    dropend = b + 100  # offset in counts after impact to include in sliced df

    if n == 1 :
        drop1 = dfCal[dropstart:dropend] #start and end points of the drop in m/s^2
        drop1g = dfCalg[dropstart:dropend] #start and end points of the drop in g
        drop1 = pd.DataFrame(drop1) # makes dataframe including all data within the start and end points of the drop
        drop1g = pd.DataFrame(drop1g)
        dg = drop1g # chooses what accelerometer to use based on the max g
        d = drop1
        accPick(dg, d)  # chooses what accelerometer to use
        acc1 = acc
        acc1Name = accName
        acc1Nameg = accNameg
        drop1.plot(y=accName, ylabel="Deceleration (g)", title='Select start and end of drop #' + str(n))
        #plt.show()
        #plt.set(title = 'Select start and end of drop #' + str(n))
        #peakpicker(dnum, dse)
        pts = plt.ginput(2)
        (x0, y0), (x1, y1) = pts
        xmin, xmax = sorted([x0, x1])
        dropTimes[dnum, dse] = xmin
        xmin = round(xmin)
        print(xmin)
        dse = dse + 1
        dropTimes[dnum, dse] = xmax
        xmax = round(xmax)
        print(xmax)
        print(dropTimes)
        d = dfCal[xmin:xmax]
        acc = d[accName]
        print(d)
        integration(d)
        drop1 = d
        dnum = dnum + 1
        dse = dse - 1

    if n == 2 :
        drop2 = dfCal[dropstart:dropend]
        drop2g = dfCalg[dropstart:dropend]
        drop2 = pd.DataFrame(drop2)
        drop2g = pd.DataFrame(drop2g)
        dg = drop2g
        d = drop2
        accPick(dg, d)  # chooses what accelerometer to use
        acc2 = acc
        acc2Name = accName
        acc2Nameg = accNameg
        drop2.plot(y=accName, ylabel="Deceleration (g)", title='Select start and end of drop #' + str(n))
        #plt.show()
        pts = plt.ginput(2)
        (x0, y0), (x1, y1) = pts
        xmin, xmax = sorted([x0, x1])
        dropTimes[dnum, dse] = xmin
        xmin = round(xmin)
        print(xmin)
        dse = dse + 1
        dropTimes[dnum, dse] = xmax
        xmax = round(xmax)
        print(xmax)
        print(dropTimes)
        d = dfCal[xmin:xmax]
        acc = d[accName]
        print(d)
        drop2 = d
        integration(d)
        dnum = dnum + 1
        dse = dse - 1

    if n == 3 :
        drop3 = dfCal[dropstart:dropend] #start and end points of the drop in m/s^2
        drop3g = dfCalg[dropstart:dropend] #start and end points of the drop in g
        drop3 = pd.DataFrame(drop3) # makes dataframe including all data within the start and end points of the drop
        drop3g = pd.DataFrame(drop3g)
        dg = drop3g # chooses what accelerometer to use based on the max g
        d = drop3
        accPick(dg, d)  # chooses what accelerometer to use
        acc3 = acc
        acc3Name = accName
        acc3Nameg = accNameg
        drop3.plot(y=accName, ylabel="Deceleration (g)", title='Select start and end of drop #' + str(n))
        pts = plt.ginput(2)
        (x0, y0), (x1, y1) = pts
        xmin, xmax = sorted([x0, x1])
        dropTimes[dnum, dse] = xmin
        xmin = round(xmin)
        print(xmin)
        dse = dse + 1
        dropTimes[dnum, dse] = xmax
        xmax = round(xmax)
        print(xmax)
        print(dropTimes)
        d = dfCal[xmin:xmax]
        acc = d[accName]
        print(d)
        #drop3.plot(x="Time (s)", y=accName, ylabel="Deceleration (g)")
        #plt.show()
        drop3 = d
        integration(d)
        dnum = dnum + 1
        dse = dse - 1

    if n == 4 :
        drop4 = dfCal[dropstart:dropend] #start and end points of the drop in m/s^2
        drop4g = dfCalg[dropstart:dropend] #start and end points of the drop in g
        drop4 = pd.DataFrame(drop4) # makes dataframe including all data within the start and end points of the drop
        drop4g = pd.DataFrame(drop4g)
        dg = drop4g # chooses what accelerometer to use based on the max g
        d = drop4
        accPick(dg, d)  # chooses what accelerometer to use
        acc4 = acc
        #drop4.plot(x="Time (s)", y=accName, ylabel="Deceleration (g)")
        #plt.show()
        integration(d)

    if n == 5 :
        drop5 = dfCal[dropstart:dropend] #start and end points of the drop in m/s^2
        drop5g = dfCalg[dropstart:dropend] #start and end points of the drop in g
        drop5 = pd.DataFrame(drop5) # makes dataframe including all data within the start and end points of the drop
        drop5g = pd.DataFrame(drop5g)
        dg = drop5g # chooses what accelerometer to use based on the max g
        d = drop5
        accPick(dg, d)  # chooses what accelerometer to use
        acc5 = acc
        #drop5.plot(x="Time (s)", y=accName, ylabel="Deceleration (g)")
        #plt.show()
        integration(d)

    if n == 6 :
        drop6 = dfCal[dropstart:dropend] #start and end points of the drop in m/s^2
        drop6g = dfCalg[dropstart:dropend] #start and end points of the drop in g
        drop6 = pd.DataFrame(drop6) # makes dataframe including all data within the start and end points of the drop
        drop6g = pd.DataFrame(drop6g)
        dg = drop6g # chooses what accelerometer to use based on the max g
        d = drop6
        accPick(dg, d)  # chooses what accelerometer to use
        acc6 = acc
        #drop6.plot(x="Time (s)", y=accName, ylabel="Deceleration (g)")
        #plt.show()
        integration(d)

    n = n + 1
    a = a + 1


#show plots

# GENERATE PLOTS

#PLot showing all accellerometers and pore pressure readings
fig, (ax1, ax2) = plt.subplots(2)

ax1.plot(time, g2g, label="2g" )
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

ax2.plot(time, ppm, label="Pore Pressure")
ax2.set(ylabel="Pore Pressure (kPa)")
ax2.set(xlabel="Time (s)")
plt.show()
#'''

#'''
# Plot showing peak deceleration
peakplot = plt.plot(x)
peakplot = plt.plot(peaks, x[peaks], "x")
plt.show()
#'''

#Deceleration,Velocity,and penetration depth vs time plots
#'''
fig, (ax1, ax2, ax3) = plt.subplots(3)
ax1.plot(drop1["Time (s)"], drop1[acc1Name], marker = 11, color = "k")
ax1.set(ylabel="Deceleration (m/s^2)", xlabel="Time(s)",title="Drop 1")
ax2.plot(drop1["Time (s)"], drop1['Velocity (m/s)'], marker = 11, color = "k")
ax2.set(ylabel="Velocity (m/s)", xlabel="Time(s)")
ax3.plot(drop1["Time (s)"], drop1["Penetration Depth (m)"], marker = 11, color = "k")
ax3.set(ylabel="Penetration Depth (m)", xlabel="Time(s)")
plt.show()

fig, (ax1, ax2, ax3) = plt.subplots(3)
ax1.plot(drop2["Time (s)"], drop2[acc2Name], marker = 11, color = "k")
ax1.set(ylabel="Deceleration (m/s^2)", xlabel="Time(s)")
ax2.plot(drop2["Time (s)"], drop2['Velocity (m/s)'], marker = 11, color = "k")
ax2.set(ylabel="Velocity (m/s)", xlabel="Time(s)")
ax3.plot(drop2["Time (s)"], drop2["Penetration Depth (m)"], marker = 11, color = "k")
ax3.set(ylabel="Penetration Depth (m)", xlabel="Time(s)")
plt.show()

fig, (ax1, ax2, ax3) = plt.subplots(3)
ax1.plot(drop3["Time (s)"], drop3[acc3Name], marker = 11, color = "k")
ax1.set(ylabel="Deceleration (m/s^2)", xlabel="Time(s)")
ax2.plot(drop3["Time (s)"], drop3['Velocity (m/s)'], marker = 11, color = "k")
ax2.set(ylabel="Velocity (m/s)", xlabel="Time(s)")
ax3.plot(drop3["Time (s)"], drop3["Penetration Depth (m)"], marker = 11, color = "k")
ax3.set(ylabel="Penetration Depth (m)", xlabel="Time(s)")
plt.show()

fig, (ax1) = plt.subplots(1)
ax1.plot(drop3g[acc3Nameg], drop3["Penetration Depth (m)"]*100, marker = 11, color = "k")
ax1.plot(drop3["Velocity (m/s)"], drop3["Penetration Depth (m)"]*100, marker = 11, color = "k")
ax1.set(xlabel="Deceleration (g) and Velocity (m/s)", ylabel="Penetration Depth (cm)")
ax1.invert_yaxis()
plt.show()

#'''