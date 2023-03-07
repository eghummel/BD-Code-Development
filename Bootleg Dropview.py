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
atype = 'p'  # m = mantle area (best for sands), p = projected area (best for clays)
tiptype = 'c'  # c = cone, p = parabolic, b = blunt
offset = 1 # this value is subtracted from the accelerometer readings
droptype = 'w' #w = water, #a = air
sign = "uk" #enter an effective unit weight value in kg/m^3 or "uk" if unknown 
# paste the filepath to the folder where the BD data is stored
binFilepath = Path("H:\My Drive\CEE 5904 - Project & Report/2021 FRF Data\Bluedrop\October/14 October 2021 AM\Pier - BD2")
#paste the filepath that you want the files containing drops copied to 
outputPath = Path("H:\My Drive\CEE 5904 - Project & Report/2021 FRF Data\Bluedrop\October/14 October 2021 AM\Analysis Results-14 October 2021 AM.xlsx") #  Path to pre-existing Excel File
plotPath = Path("H:\My Drive\CEE 5904 - Project & Report/2021 FRF Data\Bluedrop\October/14 October 2021 AM\Analysis Figures- 14 October 2021 AM")
#if applicable, paste the filepath to an excel file that troubleshooting data will be printed in
troubleshootingPath = Path("H:\My Drive\CEE 5904 - Project & Report/2021 FRF Data\Bluedrop\October/14 October 2021 AM\Troubleshooting.xlsx")

def overviewplot(): #Plot showing all accellerometers and pore pressure readings
    fig, (ax1, ax2, ax3) = plt.subplots(3)
    fig.set_size_inches(14, 7)
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
    #plotName = fileNum+" Overview.png"
    #plt.savefig(plotPath / plotName)
    plt.show()

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
'''
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

overviewplot()

print("Does the file contain a drop? (y/n)")
drop_input=input()
    if drop_input=='y':   
    files = Path(directory).glob('*')
    for file in files:
    print(file)
'''       
name = [*fileNum]
print(name)
print(name[3])

def newName():
    a = int(name[3])
    print(a)
    if a <=8:
        b = a + 1
        b = str(b)
    elif a == 9:
        b = "A"
    elif a == "A":
        b = "B"
    elif a == "B":
        b = "C"
    elif a == "C":
        b = "D"
    elif a == "D":
        b = "E"

        newName = name
        newName[3] = b
        print(newName)

newName()
