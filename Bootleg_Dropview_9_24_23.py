import numpy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import math
from math import pi
from math import sqrt
import shutil

# SETUP VARIABLES - USER INPUTS
BD = 3 #Bluedrop file is from 
fileNum = '00A0'
# write the bin file number for the earliest drop in the folder (do not include 'bLog' or '.bin')
numFiles = 20 #write the number of bin files in the folder of interest (if only looking at 1 drop, write 2)
offset = 0 # this value is subtracted from the accelerometer readings
# paste the filepath to the folder where the BD data is stored
binFilepath = Path("C:/Users/elise/Desktop/Sequim_2021/")
#paste the filepath that you want the files containing drops copied to 
outputPath = Path("C:/Users/elise/Desktop/Sequim_2021/Drops Only")
plotPath = Path("C:/Users/elise/Desktop/Sequim_2021/Drops Only/Dec-Time Plots")

def overviewplot(): #Plot showing all accelerometers and pore pressure readings
    fig, (ax1) = plt.subplots(1)
    fig.set_size_inches(4, 4)
    ax1.plot(time, g250g, label="250g")
    ax1.plot(time, g50g, label="50g")
    ax1.plot(time, g18g, label="18g")
    ax1.plot(time, g2g, label="2g")
    ax1.legend(loc = "upper right")
    ax1.set(ylabel="Deceleration (g)")
    ax1.set_title("BD file "+ newName)
    ax1.grid()

    plt.draw()
    print("Is there a drop? y/n")
    plt.waitforbuttonpress(2)
    drop_input=input()

    if drop_input=='y':   
        origFile = 'bLog'+newName+".bin"
        origFilePath = binFilepath / origFile
        newFilePath = outputPath / origFile
        plotFilePath = plotPath / newName
        shutil.copyfile(origFilePath, newFilePath)
        plt.savefig(plotPath / newName)

    plt.close(fig)

def newNamefun(): #find the next sequential bin name a=existing last variable, b=new last variable
    global newName
    a = name[3]
    c = name[2]
    e = name[1]
    g = name[0]
    if a == "A":
        b = "B"
    elif a == "B":
        b = "C"
    elif a == "C":
        b = "D"
    elif a == "D":
        b = "E"
    elif a == "E":
            b = "F"
    elif a == "F":
        b = "0"
        if c == "A":
            d = "B"
        elif c == "B":
            d = "C"
        elif c == "C":
            d = "D"
        elif c == "D":
            d = "E"
        elif c == "E":
            d = "F"
        elif c == "F":
            d = "0"
            if e == "A":
                f = "B"
            elif e == "B":
                f = "C"
            elif e == "C":
                f = "D"
            elif e == "D":
                f = "E"
            elif e == "E":
                f = "F"
            elif e == "F":
                f = "0"
                if g == "A":
                    h = "B"
                elif g == "B":
                    h = "C"
                elif g == "C":
                    h = "D"
                elif g == "D":
                    h = "E"
                elif g == "E":
                    h = "F"
                elif g == "F":
                    h = "0"
                    f = "0"
                    d = "0"
                    b = "0"
                elif g == "9":
                    h = "A"
                else: #a <=8:
                    g = int(name[0]) #the only other options are numbers
                    h = g + 1
                    h = str(b)    

            elif e == "9":
                f = "A"
            else: #a <=8:
                e = int(name[1]) #the only other options are numbers
                f = e + 1
                f = str(f)    


        elif c == "9":
            d = "A"
        else: #a <=8:
            c = int(name[2]) #the only other options are numbers
            d = c + 1
            d = str(d)    


    elif a == "9":
        b = "A"
    else: #a <=8:
        a = int(name[3]) #the only other options are numbers
        b = a + 1
        b = str(b)    


    newNameList = name
    newNameList[3] = b #b always will change
    if a == "F":
        newNameList[2] = d
        if c == "F":
            newNameList[1] = f
            if e == "F":
                newNameList[0] = h
                    
    else:
        newNameList[2] = c
        newNameList[1] = e
        newNameList[0] = g

    print(newNameList)
    newName = ""
    newName = newName.join(newNameList)
    print(newName)

newName = fileNum #for the first drop, select the earliest drop defined earlier

for i in range(1, numFiles):

    # READ BD DATA IN
    data_array = []  # creates an empty array for us to fill with bd data
    fileName = 'bLog'+newName+".bin"
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


    df = pd.DataFrame(np_array) # Creates a Dataframe in pandas from the bd data
    df.columns = ['Count', 'no clue', 'g2g', 'g18g', 'g50g', 'ppm', 'g200g', 'gX55g', 'gY55g', 'g250g'] # names columns

    # APPLY CALIBRATION FACTORS
    if BD == 3:  # calibration factors from March 2023
        g2g = (df['g2g']-33570.1)/1614577.9 - offset# accelerometers are in g
        g18g = (df['g18g']+13495)/163387.2 - offset
        g50g = (df['g50g']-238817.4)/63779.5 - offset
        ppm = ((df['ppm']-139040.1)/20705) * 6.89475729 # converts to kPa
        g200g = ((df['g200g'] -262332.4/38888.7)) - offset
        gX55g = (df['gX55g']-70406)/59754.3  
        gY55g = (df['gY55g']-69421.1)/141871.5
        g250g = (df['g250g']-39077.4)/13746.9 - offset

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

    name = [*newName]

    newNamefun()

print("all files checked")