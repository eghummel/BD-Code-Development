import numpy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.signal import find_peaks

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

def accPick(d):
    # each of the following are the same; if an accelerometer breaks on a BD, can edit that section
    # the 200g accelerometer is ignored in all cases
    if BD == 3:
        maxAcc = d["250g (g)"].max()
        global acc
        if maxAcc < 5:
            if d["2g (g)"].max() < 1.8:  # does an extra check for to 200g because of noise
                acc = g2g  #d["2g (g)"]
            else:
                acc = d["18g (g)"]
        elif maxAcc < 18:
            acc = g18g  #d["18g (g)"]
        elif maxAcc < 50:
            acc = g50g  #d["50g (g)"]
        else:
            acc = d["250g (g)"]

    if BD == 2:
        maxAcc = d["250g (g)"].max()
        if maxAcc < 5:
            if d["2g (g)"].max() < 1.8:  # does an extra check for to 200g because of noise
                acc = d["2g (g)"]
            else:
                acc = d["18g (g)"]
        elif maxAcc < 18:
            acc = d["18g (g)"]
        elif maxAcc < 50:
            acc = d["50g (g)"]
        else:
            acc = d["250g (g)"]

    if BD == 1:
        maxAcc = d["250g (g)"].max()
        if maxAcc < 5:
            if d["2g (g)"].max() < 1.8:  # does an extra check for to 200g because of noise
                acc = d["2g (g)"]
            else:
                acc = d["18g (g)"]
        elif maxAcc < 18:
            acc = d["18g (g)"]
        elif maxAcc < 50:
            acc = d["50g (g)"]
        else:
            acc = d["250g (g)"]


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

#print(np_array.shape)
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

dfCal = pd.DataFrame([time, g2g, g18g, g50g, g200g, g250g, gX55g, gY55g, ppm]) # creates a calibrated data frame
dfCal = dfCal.T
dfCal.columns = ['Time (s)', '2g (g)', '18g (g)', '50g (g)', '200g (g)', '250g (g)', 'X55g (g)', 'Y55g (g)', 'Pore Pressure (kPa)']  # names columns

#print(dfCal)

#AUTOMATIC PEAK FINDING
x = np.array(g200g)  # what accelerometer to get the peaks from
peaks, _ = find_peaks(x, height = 5, distance=10000)  # finds the largest peaks more than 5g spaced at least 10000 counts apart
plt.plot(x)
plt.plot(peaks, x[peaks], "x")
plt.show()

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
        drop1 = dfCal[dropstart:dropend]
        drop1 = pd.DataFrame(drop1) # makes dataframe including all data within the start and end points of the drop
        d = drop1
        accPick(d)
        print(acc)
        drop1.plot(x="Time (s)", y=acc, ylabel="Deceleration (g)")
        plt.show()
        #print(drop1)

    if n == 2 :
        drop2 = dfCal[dropstart:dropend]
        drop2 = pd.DataFrame(drop2)  # makes dataframe including all data within the start and end points of the drop
        d = drop2
        accPick(d)
        print(acc)
        drop2.plot(x="Time (s)", y="50g (g)", ylabel="Deceleration (g)")
        plt.show()
        #print(drop2)

    if n == 3 :
        drop3 = dfCal[dropstart:dropend]
        drop3 = pd.DataFrame(drop3)  # makes dataframe including all data within the start and end points of the drop
        d = drop3
        accPick(d)
        print(acc)
        drop3.plot(x="Time (s)", y="50g (g)", ylabel="Deceleration (g)")
        plt.show()
        #print(drop3)

    if n == 4 :
        drop4 = dfCal[dropstart:dropend]
        drop4 = pd.DataFrame(drop4)  # makes dataframe including all data within the start and end points of the drop
        d = drop4
        drop4.plot(x="Time (s)", y="50g (g)", ylabel="Deceleration (g)")
        plt.show()
        #print(drop4)

    if n == 5 :
        drop5 = dfCal[dropstart:dropend]
        drop5 = pd.DataFrame(drop5)  # makes dataframe including all data within the start and end points of the drop
        d = drop5
        drop5.plot(x="Time (s)", y="50g (g)", ylabel="Deceleration (g)")
        plt.show()
        #print(drop5)

    if n == 6 :
        drop6 = dfCal[dropstart:dropend]
        drop6 = pd.DataFrame(drop6)  # makes dataframe including all data within the start and end points of the drop
        d = drop6
        drop6.plot(x="Time (s)", y="50g (g)", ylabel="Deceleration (g)")
        plt.show()
        #print(drop6)

    #outputName = fileNum + " Drop " + str(n)
    #print(outputName)
    n = n + 1
    a = a + 1

