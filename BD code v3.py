import numpy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def tellme(s):
    #print(s)
    ax1.set(title=s)
    plt.draw()

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

dfCal = pd.DataFrame([time, g2g, g18g, g50g, g200g, g250g, gX55g, gY55g, ppm]) # copies the data frame; this version will be edited with calibration factors
dfCal= dfCal.T
dfCal.columns = ['Time (s)', '2g (g)', '18g (g)', '50g (g)', '200g (g)', '250g (g)', 'X55g (g)', 'Y55g (g)', 'Pore Pressure (kPa)']  # names columns
#print(dfCal)

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
#plt.show()

#  SELECT START AND END POINTS
tellme("How many drops?")
#plt.show()

nDrops = int(input("Input the number of Drops in file:"))
#print(nDrops)

n = 1
b = 0
a = 0

dropTimes = np.array([[0]*2]*nDrops)

while n <= nDrops :
    tellme('Select start and end of drop #' + str(n))
    pts = plt.ginput(2)
    (x0, y0), (x1, y1) = pts
    xmin, xmax = sorted([x0, x1])
    dropTimes[b,a] = xmin
    a = a + 1
    dropTimes[b,a] = xmax
    b = b + 1
    a = a - 1
    n = n + 1

print(dropTimes)
'''
#Make a new dataframe for each drop
n = 1
x1r = 0
x1c = 0
x2r = 0
x2c = 1

def slicedDf(dStart,dEnd):
    dfCal.iloc[(dfCal["Time (s)"] == dStart : dfCal["Time (s)"] == dEnd),:]

while n<nDrops:
    dStart = dropTimes[x1r,x1c]
    dEnd = dropTimes[x2r, x2c]
    drop1 = slicedDf(dStart,dEnd)
    x1r+1
    x2r+1
    n=n+1

print(drop1)

'''

'''This file provides all calibration constants that were used in this study. 
Depending on the date of the survey, the calibration constants may need to be adjusted.
This refers to lines 233-277 in the code callled "BDMasterBeach.m" within the 
subfunction called gdata. The three cases refer to the three blueDrops.

blueDrop #1 calibration factors:
	- for data collected from July 2020 onward, use the following calibration factors:
	    g2g=((double(F(:,3)))-42590.9)/1626361.1; %accelerometers are in g
            g18g=((double(F(:,4)))-44492.9)/161125.5;
            g50g=((double(F(:,5)))-171656.1)/64020.3;
            ppm=((double(F(:,6)))+31776.1)/20679.7;  % this is psi
            g200g=(((double(F(:,7)))-723404.8)/32209.7);
            gX55g=((double(F(:,8))-54881.1)/64858.6);
            gY55g=((double(F(:,9)))-28735.5)/63839.9;
            g250g=((double(F(:,10)))+13299.7)/13697.1;
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
            ppm=ppm*6.89475729; % convert into kPa'''