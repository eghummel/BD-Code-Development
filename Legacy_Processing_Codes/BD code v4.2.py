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
outputFile = 'data.xlsx' # this currently doesn't do anything, but eventually all data will be printed out into an excel sheet
outputPath = Path("H:\My Drive\CEE 5904 - Project & Report\FRF Data/test data") #  Path for new files

def accPick(dg):
    # each of the following are the same; if an accelerometer breaks on a BD, can edit that section
    # the 200g accelerometer is ignored in all cases
    if BD == 3:
        maxAcc = dg["250g (g)"].max()
        global acc
        global accName
        if maxAcc < 5:
            if dg["2g (g)"].max() < 1.8:  # does an extra check for the 2g because of noise
                acc = dg["2g (g)"]
                accName = "2g (g)"
            else:
                acc = dg["18g (g)"]
                accName = "18g (g)"
        elif maxAcc < 18:
            acc = dg["18g (g)"]
            accName = "18g (g)"
        elif maxAcc < 50:
            acc = dg["50g (g)"]
            accName = "50g (g)"
        else:
            acc = dg["250g (g)"]
            accName = "250g (g)"

    """def integration(d):
    v = scipy.integrate.cumtrapz(acc, x=d["Time (s)"])  # solves for velocity
    v = v
    maxv = v.max()
    depth = scipy.integrate.cumtrapz(v, x=d["Time (s)"])  # solves for penetration depth
    depth = depth
    maxd = depth.max()
    d.insert("Velocity", v)
    d.insert("Penetration Depth", depth)"""



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
    g2g = (df['g2g']-34426.5)/1615925.8 # accelerometers are in g
    g18g = (df['g18g']+12322.1)/163530.7
    g50g = (df['g50g']-237384.9)/63651 - 0.1120
    ppm = ((df['ppm']+62496.7)/20583.0) * 6.89475729 # converts to kPa
    g200g = ((df['g200g'] -248943.7)/39009.4)+0.5518
    gX55g = (df['gX55g']-59093.7)/66674.3
    gY55g = (df['gY55g']-140224.6)/66674.3
    g250g = (df['g250g']-40536.1)/13631.6

if BD == 2: # calibration factors from Aug 26, 2021
    g2g = (df['g2g']+31384.7)/1624987.2-0.035 # accelerometers are in g
    g18g = (df['g18g']-26631.0)/159945.4
    g50g = (df['g50g']+92987.0)/63783.5
    ppm = ((df['ppm']-35170.6)/12922.9) * 6.89475729 # converts to kPa
    g200g = (df['g200g']-16264.8)/26042.8 -0.277
    gX55g = (df['gX55g']+89890.3)/63897.1
    gY55g = (df['gY55g']+14993.0)/64118.0
    g250g = (df['g250g']+17362.1)/13533.5+0.0656

if BD == 1: # calibration factors from July 2020
    g2g = (df['g2g']+277743.2)/1637299.6 # accelerometers are in g
    g18g = (df['g18g']-3755.9)/159932.2
    g50g = (df['g50g']+92817.6)/63237.1
    ppm = ((df['ppm']-33154.0)/14763.5) * 6.89475729 # this is kPa
    g200g = (df['g200g'] -1155309.9)/28368.5 - 1.464
    gX55g = (df['gX55g'] +97138.4)/62023.7
    gY55g = (df['gY55g']-9921.7)/62669.2
    g250g = (df['g250g']+59211.3)/13276.9

time = (df['Count']-df['Count'].iloc[0]+1)/2000 # gives time in s

dfCalg = pd.DataFrame([time, g2g, g18g, g50g, g200g, g250g, gX55g, gY55g, ppm]) # copies the data frame; this version will be edited with calibration factors
dfCalg = dfCalg.T
dfCalg.columns = ['Time (s)', '2g (g)', '18g (g)', '50g (g)', '200g (g)', '250g (g)', 'X55g (g)', 'Y55g (g)', 'Pore Pressure (kPa)']  # names columns
#print(dfCalg)

dfCal = dfCalg
dfCal.columns = ['Time (s)', '2g (m/s^2)', '18g (m/s^2)', '50g (m/s^2)', '200g (m/s^2)', '250g (m/s^2)', 'X55g (m/s^2)', 'Y55g (m/s^2)', 'Pore Pressure (kPa)']  # names columns
dfCal['2g (m/s^2)'] = dfCal['2g (m/s^2)'] * 9.80665
dfCal['18g (m/s^2)'] = dfCal['18g (m/s^2)'] * 9.80665
dfCal['50g (m/s^2)'] = dfCal['50g (m/s^2)'] * 9.80665
dfCal['200g (m/s^2)'] = dfCal['200g (m/s^2)'] * 9.80665
dfCal['250g (m/s^2)'] = dfCal['250g (m/s^2)'] * 9.80665
dfCal['X55g (m/s^2)'] = dfCal['X55g (m/s^2)'] * 9.80665
dfCal['Y55g (m/s^2)'] = dfCal['Y55g (m/s^2)'] * 9.80665
print(dfCal)


'''
# GENERATE PLOTS

fig, (ax1, ax2) = plt.subplots(2)

# plot of all deceleration data
ax1.plot(time, g2g, label="2g")
ax1.plot(time, g18g, label="18g")
ax1.plot(time, g50g, label="50g")
#plt.plot(time, ppm)
ax1.plot(time, g200g, label="200g")
#plt.plot(time, gX55g, label="2g")
#plt.plot(time, gY55g, label="2g")
ax1.plot(time, g250g, label="250g")
ax1.legend()
ax1.set(ylabel="Deceleration (g)")
ax1.set(xlabel="Time (s)")

# Plot pore pressure data
ax2.plot(time, ppm, label="Pore Pressure")
ax2.set(ylabel="Pore Pressure (kPa)")
ax2.set(xlabel="Time (s)")
plt.show()

'''

#AUTOMATIC PEAK FINDING
x = np.array(g250g)  # what accelerometer to get the peaks from
peaks, _ = find_peaks(x, height = 5, distance=10000)  # finds the largest peaks more than 5g spaced at least 10000 counts apart
plt.plot(x)
plt.plot(peaks, x[peaks], "x")
#plt.show()

# CREATION OF INDIVIDUAL DROP FILES
peaksArray = np.array(peaks)  # prints a list of the count where the peaks occur
print(peaksArray)
q = (peaksArray.shape) #gives number of peaks
nDrops = int(q[0]) #number of drops in the file
print(nDrops)

a = 0
n = 1

    # MAKE INDIVIDUAL DATAFRAMES FOR EACH DROP AND PLOT DECELERATION VS TIME
while n <= nDrops :
    b = int(peaksArray[a])  # count at the ath drop
    dropstart = b - 100  # offset in counts before impact to include in sliced df
    dropend = b + 100  # offset in counts after impact to include in sliced df

    if n == 1 :
        #drop1 = dfCal[dropstart:dropend] #start and end points of the drop in m/s^2
        drop1g = dfCalg[dropstart:dropend] #start and end points of the drop in g
        #drop1 = pd.DataFrame(drop1) # makes dataframe including all data within the start and end points of the drop
        dg = drop1g # chooses what accelerometer to use based on the max g
        #d = drop1
        accPick(dg)  # chooses what accelerometer to use
        v = scipy.integrate.cumtrapz(acc, x=dg["Time (s)"])  # solves for velocity
        v = np.array(v)
        #v = v.T
        print(v.shape)
        print(v)
        #integration(d)
        drop1g.plot(x="Time (s)", y=accName, ylabel="Deceleration (m/s^2)")
        plt.show()
        print(drop1g)

    """ if n == 2 :
        drop2 = dfCalMs[dropstart:dropend]
        drop2 = pd.DataFrame(drop2)  # makes dataframe including all data within the start and end points of the drop
        d = drop2
        accPick(d)
        drop2.plot(x="Time (s)", y=accName, ylabel="Deceleration (g)")
        plt.show()
        #print(drop2)

    if n == 3 :
        drop3 = dfCalMs[dropstart:dropend]
        drop3 = pd.DataFrame(drop3)  # makes dataframe including all data within the start and end points of the drop
        d = drop3
        accPick(d)
        drop3.plot(x="Time (s)", y=accName, ylabel="Deceleration (g)")
        plt.show()
        #print(drop3)

    if n == 4 :
        drop4 = dfCalMs[dropstart:dropend]
        drop4 = pd.DataFrame(drop4)  # makes dataframe including all data within the start and end points of the drop
        d = drop4
        accPick(d)
        drop4.plot(x="Time (s)", y=accName, ylabel="Deceleration (g)")
        plt.show()
        #print(drop4)

    if n == 5 :
        drop5 = dfCalMs[dropstart:dropend]
        drop5 = pd.DataFrame(drop5)  # makes dataframe including all data within the start and end points of the drop
        d = drop5
        accPick(d)
        drop5.plot(x="Time (s)", y=accName, ylabel="Deceleration (g)")
        plt.show()
        #print(drop5)

    if n == 6 :
        drop6 = dfCalMs[dropstart:dropend]
        drop6 = pd.DataFrame(drop6)  # makes dataframe including all data within the start and end points of the drop
        d = drop6
        accPick(d)
        drop6.plot(x="Time (s)", y=accName, ylabel="Deceleration (g)")
        plt.show()
        #print(drop6)
    """
    n = n + 1
    a = a + 1
    #outputName = fileNum + " Drop " + str(n)
    #print(outputName)










