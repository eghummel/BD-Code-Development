#change python environment to BlueDrop

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ruptures as rpt
from pathlib import Path 
import math
import statistics
import os

import scipy.integrate
from scipy.signal import find_peaks
from numpy import trapz
from scipy.integrate import cumtrapz

#Sets the plots to higher default DPI
plt.rcParams['figure.dpi'] = 600

#%%
#  SETUP VARIABLES - USER INPUTS
BD = 1
atype = 'p'  # m = mantle area, p = projected area
tiptype = 'c'  # c = cone, p = parabolic, b = blunt
# paste the filepath to the folder where the bd data is stored
#Change to external hard drive depending on which computer is used
binFilepath = Path("H:\My Drive\CEE 5904 - Project & Report/2021 FRF Data\Bluedrop\October/14 October 2021 AM\BD2 - Pier\Pier")
# write the bin file number you want to analyze (do not include 'bLog' or '.bin')

fileNum = '02B7'

# input based file entering system for ease of running the program for numerour files. Can be commented out and the original method used 
#of desired.
#print('Type the file number')
#fileNum=input()

#I added a section to combine the dataframes for the individual drops and export them as an excel sheet_FTF
# outputFile = 'data.xlsx'  #this currently doesn't do anything, but eventually all data will be printed out into an excel sheet
# outputPath = Path("H:\My Drive\CEE 5904 - Project & Report\FRF Data/test data") #  Path for new files
offset = 1 # this value is subtracted from the accelerometer readings
droptype = 'a' #w = water, #a = air

def masslength(tiptype): #sets the mass and length of the pentrometer based on the tip type 
    global mass
    global length
    if tiptype == 'c':
        mass = 7.71 #units?
        length = 0.0787 #m
    elif tiptype == 'p':
        mass = 9.15
        length = 8.26
    elif tiptype == 'b':
        mass = 10.30
        length = 8.57
    return mass, length

#def areafind(atype, tiptype, length):
    

def areafind():#this currenly only works for conical tips, can be updapdated for other tips if they are needed_FTF
    r=np.zeros(len(dep))
    A=np.zeros(len(dep))
    for i in range(1,len(dep)):
        if dep[i]<length:#tip length
            r[i]=dep[i]*np.tan(np.radians(30))
            A[i]=np.pi*r[i]**2
        elif dep[i]>=length:
            r[i]=0.04375
            A[i]=np.pi*r[i]**2 
        #A[i]=A[i]/10000
    return A
        

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
        return dropstart,dropend


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
        accg = dg["250g (g)"] #this was set to 50g (g), I assume this should be 250g to match the rest?_FTF
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
        jlist.append(j)#adds the value of j to the end of the list

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
    #jindex = int(jindex)
    #print("jindex is")
    #print(jindex)

def finddropend(): #finds the location where the deceleration is 1-offset after the peak
    global num1
    global num2
    below0list = list()
    for i in range(dropstart+jindex, dropend):
        if accg[i] < 1 - offset:
            num1 = i - dropstart
            #num2 = i-jindex-11
            below0list = np.append(below0list, num1)
    num1=int(min(below0list))
    #print("num 1 is ")
    #print(num1)



def integration(d): #integrates the deceleration data to solve for velocity and penetration depth
    global vel
    global maxvel
    global dep
    global maxdep
    accint = acc[jindex:num1]
    vel = scipy.integrate.cumtrapz(accint, x=d["Time (s)"])  # solves for velocity
    vel = np.array(vel)
    vel = np.insert(vel, 0, 0) #not sure what to insert here, but it makes it the right size
    vel = np.flip(vel)
    maxvel = vel.max()
    dep = scipy.integrate.cumtrapz(vel, x=d["Time (s)"])  # solves for penetration depth
    dep = np.insert(dep, 0, 0)  # not sure what to insert here, but it makes it the right size
    maxdep = dep.max()
    d.insert(9, "Velocity (m/s)", vel)
    d.insert(10, "Penetration Depth (m)", dep)

#This function is added to check the values of the start and end of the drop for softer drops where the autofinding
#fucntion is not currently working FTF 12/19/2022
#This function can be commented out so that it does not do anything if the autofinding function is working for the series of drops being worked on
#This function is placed within the while loop for each drop
def CheckingFunction():
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
#these values should be slightly before the start and slightly after the end             
        print('What is the approximate drop start and drop end?')
        approx_start=input()
        approx_end=input()
        approx_end=int(approx_end)
        approx_start=int(approx_start)
        test=pd.DataFrame.copy(temp_dg[approx_start:approx_end])
        #this is a little weird, im sure there is a better way to do it but this will pick the start and end of the penetration
        #after the correct jindex and num1 values are input_FTF.
        for i in range(0,len(test)):
            if  test[i+test.index[0]]<0:
                test[i+test.index[0]]=np.nan
            else: test[i+test.index[0]]=test[i+test.index[0]]        
        jindex_add=pd.Series.first_valid_index(test)
        num1_add=pd.Series.last_valid_index(test)
        jindex=jindex_add-temp_dg.index[0]
        num1=num1_add-temp_dg.index[0]
      
        plt.plot(xvalues[jindex:num1],temp_dg[jindex:num1])
        plt.grid(visible=True,which='both')
        plt.show()
        #This is what I originally had for the corrections, this can be removed, it is not accurate enough. I am going to keep it
        #in my file for now_FTF
        # plt.figure(num=2)
        # plt.plot(xvalues[approx_start:approx_end],temp_dg[approx_start:approx_end])
        # plt.grid(visible=True)
        # plt.show()
        # print('Select the correct drop start')
        # jindex=input()
        # jindex=int(jindex)
        # plt.figure(num=1)
        # print ('Select the correct drop end')
        # num1=input()
        # num1=int(num1)
        return jindex,num1
#The following variable and function are to calculate the QSBC for the drops. The QSBC values do not always look correct
#I need to look into this more_FTF  

#define the cone factors to use
srfk =[0.2, 0.4, 1.0, 1.5]
#computed the quasistartic bearing capacity for the drops for each cone factor selected
def QSBC():
    global qsbc_1
    global qsbc_2
    global qsbc_3
    global qsbc_4
    
    force=(mass*d[accName])
    qdyn=force/A #dynamic bearing capacity [Pa]
    srcv=np.array(np.log10(d['Velocity (m/s)']/0.02))#velocity portion of the strain correction, 0.02 is the 2cm/s push speed of a CPT
    # srfK = [0.2, 0.4 ,1 ,1.5]
    frs_1=1+srfk[0]*srcv
    qsbc_1=(qdyn/frs_1)/1000
    frs_2=1+srfk[1]*srcv
    qsbc_2=(qdyn/frs_2)/1000    
    frs_3=1+srfk[2]*srcv
    qsbc_3=(qdyn/frs_3)/1000
    frs_4=1+srfk[3]*srcv
    qsbc_4=(qdyn/frs_4)/1000
    #This section removes negative values of bearing capacity. This is probably not the best way to do it, I need to check the 
    #matlab code a littler closer to figure out how this function should be changed_FTF
    for i in range (qsbc_1.index[0],qsbc_1.index[0]+len(qsbc_1)):
        if qsbc_1[i]<0 or qsbc_1[i]==-np.inf:
           qsbc_1[i]=0
    for i in range (qsbc_2.index[0],qsbc_2.index[0]+len(qsbc_2)):
        if qsbc_2[i]<0 or qsbc_2[i]==-np.inf:
            qsbc_2[i]=0
    for i in range (qsbc_3.index[0],qsbc_3.index[0]+len(qsbc_3)):
        if qsbc_3[i]<0 or qsbc_3[i]==-np.inf:
            qsbc_3[i]=0
    for i in range (qsbc_4.index[0],qsbc_4.index[0]+len(qsbc_4)):
        if qsbc_4[i]<0 or qsbc_4[i]==-np.inf:
            qsbc_4[i]=0
                  
    return qsbc_1, qsbc_2, qsbc_3, qsbc_4
   
# READ BD DATA IN
data_array = []  # creates an empty array for us to fill with bd data
fileName = 'Copy of bLog'+fileNum+".bin"
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
#These calibration factors do not match what julie has in her bluedrop code
if BD == 1: # calibration factors from July 2020
    # trying with values from Julies File
    g2g = (df['g2g']-42590.9)/1626361.1 - offset # accelerometers are in g
    g18g = (df['g18g']-44492.9)/161125.5 - offset
    g50g = (df['g50g']-171656.1)/64020.3 - offset
    ppm = ((df['ppm']-33154.0)/14763.5) * 6.89475729 # this is kPa
    g200g = (df['g200g'] -723404.8)/32209.7  - offset# - 1.464
    gX55g = (df['gX55g'] -54881.1)/64858.6 - offset
    gY55g = (df['gY55g']-28735.5)/63839.9 - offset
    g250g = (df['g250g']+13299.7)/13697.1 - offset
    
    # g2g = (df['g2g']+277743.2)/1637299.6 - offset # accelerometers are in g
    # g18g = (df['g18g']-3755.9)/159932.2 - offset
    # g50g = (df['g50g']+92817.6)/63237.1 - offset
    # ppm = ((df['ppm']-33154.0)/14763.5) * 6.89475729 # this is kPa
    # g200g = (df['g200g'] -1155309.9)/28368.5 - 1.464 - offset
    # gX55g = (df['gX55g'] +97138.4)/62023.7 - offset
    # gY55g = (df['gY55g']-9921.7)/62669.2 - offset
    # g250g = (df['g250g']+59211.3)/13276.9 - offset

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

peaks, _ = find_peaks(x, height = 5, distance=10000)  # finds the largest peaks more than 5g spaced at least 10000 counts apart

peaksArray = np.array(peaks)  # prints a list of the count where the peaks occur
#print(peaksArray)
q = (peaksArray.shape) #gives number of peaks
nDrops = int(q[0]) #number of drops in the file
#print(nDrops)

# For each drop, find the start and end points and integrate to solve for velocity and acceleration


# #%%

# a = 0
# n = 1
# for i in range(nDrops):
#     peak=peaksArray[i]
#     start=dropstartend(peak)
#     print('dropstart=',dropstart,'dropend=',dropend)
#     print('peak=',peak)
#%%
masslength(tiptype)
a=0
n=1
while n <= nDrops :
    peak = int(peaksArray[a])  # count at the ath drop
    dropstartend(peak) #zooms in the drop file to only consider 500 counts before and 1500 counts after the peak deceleration
    print('drop start =',dropstart,'drop end=', dropend)

    if n == 1 :
        print('a=',a)
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
        finddropend() #finds the end of the drop
        #print(drop)

        CheckingFunction()
        
        d = d[jindex:num1] #shortens the dataframe to only include the data during penetration (jindex = start, num1 = end)
        dg = dg[jindex:num1]
        #print(d)
        #print(np.size(d))
        drop1 = d
        drop1g = dg
        integration(d) #solves for velocity and displacement
        drop1 = d #this dataframe now includes velocity and acceleration data
        A=areafind()
        drop1['Cone Area (m^2)']=A #add cone area to the dataframe
      
        QSBC()#run the QSBC function
        #add the QSBC values to the dataframe
        drop1['QSBC (kPa) K=0.2']=qsbc_1
        drop1['QSBC (kPa) K=0.4']=qsbc_2
        drop1['QSBC (kPa) K=1.0']=qsbc_3
        drop1['QSBC (kPa) K=1.5']=qsbc_4
        #print(drop1)
        


    if n == 2 :
        print('a=',a)
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

        CheckingFunction()
        
        d = d[jindex:num1]
        dg = dg[jindex:num1]
        drop2 = d
        drop2g = dg
        # drop1plot = drop1.plot(y=accName, ylabel="Deceleration (g)", title="drop 1")
        # drop1plot = plt.plot(acc1Name, acc1Name[num1], "x")
        integration(d)
        drop2 = d
        A=areafind()
        drop2['Cone Area (m^2)']=A #add cone area to the dataframe
        QSBC()#run the QSBC function
        #add the QSBC values to the dataframe
        drop2['QSBC (kPa) K=0.2']=qsbc_1
        drop2['QSBC (kPa) K=0.4']=qsbc_2
        drop2['QSBC (kPa) K=1.0']=qsbc_3
        drop2['QSBC (kPa) K=1.5']=qsbc_4
        
       
    if n == 3 :
        print ('a=', a)
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
        finddropend() #finds the end of the drop
        #this section of code goes with the ChangeFunction to confirm that the points are selected correctly. 

        CheckingFunction()
        
        d = d[jindex:num1] #shortens the dataframe to only include the data during penetration (jindex = start, num1 = end)
        dg = dg[jindex:num1]
        drop3 = d
        drop3g = dg
        integration(d) #solves for velocity and acceleration
        drop3 = d #this dataframe now includes velocity and acceleration data
        A=areafind()
        drop3['Cone Area (m^2)']=A #add cone area to the dataframe
        QSBC()#run the QSBC function
        #add the QSBC values to the dataframe
        drop3['QSBC (kPa) K=0.2']=qsbc_1
        drop3['QSBC (kPa) K=0.4']=qsbc_2
        drop3['QSBC (kPa) K=1.0']=qsbc_3
        drop3['QSBC (kPa) K=1.5']=qsbc_4

    if n == 4 :
        print('a=',a)
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
        #this section of code goes with the ChangeFunction to confirm that the points are selected correctly. 

        CheckingFunction()
        
        d = d[jindex:num1]
        dg = dg[jindex:num1]
        #print(np.size(d))
        drop4 = d
        drop4g = dg
        # drop1plot = drop1.plot(y=accName, ylabel="Deceleration (g)", title="drop 1")
        #drop1plot = plt.plot(acc1Name, acc1Name[num1], "x")
        integration(d)
        drop4 = d
        A=areafind()
        drop4['Cone Area (m^2)']=A #add cone area to the dataframe
        QSBC()#run the QSBC function
        #add the QSBC values to the dataframe
        drop4['QSBC (kPa) K=0.2']=qsbc_1
        drop4['QSBC (kPa) K=0.4']=qsbc_2
        drop4['QSBC (kPa) K=1.0']=qsbc_3
        drop4['QSBC (kPa) K=1.5']=qsbc_4

    if n == 5 :
        print('a=',a)
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

        CheckingFunction()
        
        #print(drop)
        d = d[jindex:num1]
        dg = dg[jindex:num1]
        drop5 = d
        drop5g = dg
        # drop1plot = drop1.plot(y=accName, ylabel="Deceleration (g)", title="drop 1")
        #drop1plot = plt.plot(acc1Name, acc1Name[num1], "x")
        integration(d)
        drop5 = d
        A=areafind()
        drop5['Cone Area (m^2)']=A #add cone area to the dataframe
        QSBC()#run the QSBC function
        #add the QSBC values to the dataframe
        drop5['QSBC (kPa) K=0.2']=qsbc_1
        drop5['QSBC (kPa) K=0.4']=qsbc_2
        drop5['QSBC (kPa) K=1.0']=qsbc_3
        drop5['QSBC (kPa) K=1.5']=qsbc_4
        

    if n == 6 :
        print('a=',a)
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

        CheckingFunction()
        
        d = d[jindex:num1]
        dg = dg[jindex:num1]
        drop6 = d
        drop6g = dg
        # drop1plot = drop1.plot(y=accName, ylabel="Deceleration (g)", title="drop 1")
        #drop1plot = plt.plot(acc1Name, acc1Name[num1], "x")
        integration(d)
        drop6 = d
        A=areafind()
        drop6['Cone Area (m^2)']=A #add cone area to the dataframe
        QSBC()#run the QSBC function
        #add the QSBC values to the dataframe
        drop6['QSBC (kPa) K=0.2']=qsbc_1
        drop6['QSBC (kPa) K=0.4']=qsbc_2
        drop6['QSBC (kPa) K=1.0']=qsbc_3
        drop6['QSBC (kPa) K=1.5']=qsbc_4
        
        
        if n == 7 :
            print('a=',a)
            drop7 = dfCal[dropstart:dropend]  # start and end points of the drop in m/s^2
            drop7g = dfCalg[dropstart:dropend]  # start and end points of the drop in g
            drop7 = pd.DataFrame(drop7)  # makes dataframe including all data within the start and end points of the drop
            drop7g = pd.DataFrame(drop7g)
            dg =drop7g # chooses what accelerometer to use based on the max g
            d = drop7
            accPick(dg, d)  # chooses what accelerometer to use
            acc7 = acc
            acc7Name = accName
            acc7Nameg = accNameg
            finddropend()
            #print(num1, num2)
            #print(acc)
            #print(acc.iloc[1])
            findchangepts()
            #print(drop)

            CheckingFunction()
            
            d = d[jindex:num1]
            drop7 = d
            drop7g = dg
            # drop1plot = drop1.plot(y=accName, ylabel="Deceleration (g)", title="drop 1")
            #drop1plot = plt.plot(acc1Name, acc1Name[num1], "x")
            integration(d)
            drop7 = d
            A=areafind()
            drop7['Cone Area (m^2)']=A #add cone area to the dataframe
            QSBC()#run the QSBC function
            #add the QSBC values to the dataframe
            drop7['QSBC (kPa) K=0.2']=qsbc_1
            drop7['QSBC (kPa) K=0.4']=qsbc_2
            drop7['QSBC (kPa) K=1.0']=qsbc_3
            drop7['QSBC (kPa) K=1.5']=qsbc_4
            

    n = n + 1
    a = a + 1
#%%
#combine the dataframes and export to excel csv
if nDrops==1:
    frames=[drop1]
elif nDrops==2:
    frames=[drop1,drop2]
elif nDrops==3:
    frames=[drop1, drop2, drop3]
elif nDrops==4:
    frames=[drop1, drop2, drop3, drop4]
elif nDrops==5:
    frames=[drop1, drop2, drop3, drop4, drop5]
elif nDrops==6:
    frames=[drop1, drop2, drop3, drop4, drop5, drop6]
elif nDrops==7:
    frames=[drop1, drop2, drop3, drop4, drop5, drop6, drop7]
data_file=pd.concat(frames)
data_file.to_csv(fileNum+'.csv')
#show plots
#%%
# GENERATE PLOTS
# If there are less than 6 drops higher drop plots will plot blank
# Add more plot sections if there are more than 6 drops


#PLot showing all accellerometers and pore pressure readings
fig, (ax1, ax2) = plt.subplots(2, sharex=True)

ax1.plot(time, g2g, label="2g", )
ax1.plot(time, g18g, label="18g")
ax1.plot(time, g50g, label="50g")
plt.plot(time, ppm)
ax1.plot(time, g200g, label="200g")
#plt.plot(time, gX55g, label="2g")
#plt.plot(time, gY55g, label="2g")
#ax1.plot(time, g250g, label="250g")
ax1.legend()
ax1.set(ylabel="Deceleration (g)")
# ax1.set(xlabel="Time (s)")

ax2.plot(time, ppm, label="Pore Pressure")
ax2.set(ylabel="Pore Pressure (kPa)")
ax2.set(xlabel="Time (s)")
plt.show()
#%%
# Plot showing peak deceleration
peakplot = plt.plot(x)
peakplot = plt.plot(peaks, x[peaks], "x")
plt.show()

#I edited the following plot outputs from what was originally in the file for what I wanted the plots to output as_FTF
#%%
#Deceleration,Velocity,and penetration depth vs time plots for Drop 1

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
ax1.plot(drop1["Time (s)"], drop1[acc1Name], color = "k", marker = 11)
ax1.set(ylabel="Deceleration (m/s^2)")
ax2.plot(drop1["Time (s)"], drop1['Velocity (m/s)'] , color = "k", marker = 11)
ax2.set(ylabel="Velocity (m/s)")
ax3.plot(drop1["Time (s)"], drop1["Penetration Depth (m)"] , color = "k", marker = 11)
ax3.set(ylabel="Penetration Depth (m)", xlabel="Time(s)")
ax1.plot(drop1["Time (s)"], drop1["250g (m/s^2)"], color = "k", marker = 11)
plt.suptitle(fileNum+' Drop 1')

plt.show()
#%%
#Deceleration,Velocity,and penetration depth vs time plots for Drop 2
fig, (ax1, ax2, ax3) = plt.subplots(3,sharex=True)
ax1.set(ylabel="Deceleration (m/s^2)")
ax2.plot(drop2["Time (s)"], drop2['Velocity (m/s)'], color = "k", marker = 11)
ax2.set(ylabel="Velocity (m/s)",)
ax3.plot(drop2["Time (s)"], drop2["Penetration Depth (m)"], color = "k", marker = 11)
ax3.set(ylabel="Penetration Depth (m)", xlabel="Time(s)")
ax1.plot(drop2["Time (s)"], drop2["250g (m/s^2)"], color = "k", marker = 11)
plt.suptitle(fileNum+' Drop 2')

plt.show()
#%%
#Deceleration,Velocity,and penetration depth vs time plots for Drop 3
fig, (ax1, ax2, ax3) = plt.subplots(3,sharex=True)
ax1.plot(drop3["Time (s)"], drop3["250g (m/s^2)"], marker = 11, color = "k")
ax1.set(ylabel="Deceleration (m/s^2)")

ax2.plot(drop3["Time (s)"], drop3['Velocity (m/s)'], marker = 11, color = "k")
ax2.set(ylabel="Velocity (m/s)")
ax3.plot(drop3["Time (s)"], drop3["Penetration Depth (m)"], marker = 11, color = "k")
ax3.set(ylabel="Penetration Depth (m)", xlabel="Time(s)")
ax1.plot(drop3["Time (s)"], drop3["250g (m/s^2)"], color = "k", marker = 11)
plt.suptitle(fileNum+' Drop 3')
plt.show()

#%%
#Deceleration,Velocity,and penetration depth vs time plots for Drop 4
fig, (ax1, ax2, ax3) = plt.subplots(3,sharex=True)
ax1.plot(drop4["Time (s)"], drop4["250g (m/s^2)"], marker = 11, color = "k")
ax1.set(ylabel="Deceleration (m/s^2)")

ax2.plot(drop4["Time (s)"], drop4['Velocity (m/s)'], marker = 11, color = "k")
ax2.set(ylabel="Velocity (m/s)", xlabel="Time(s)")
ax3.plot(drop4["Time (s)"], drop4["Penetration Depth (m)"], marker = 11, color = "k")
ax3.set(ylabel="Penetration Depth (m)", xlabel="Time(s)")
ax1.plot(drop4["Time (s)"], drop4["250g (m/s^2)"], color = "k", marker = 11)
plt.suptitle(fileNum+' Drop 4')
plt.show()

#%%
#Deceleration,Velocity,and penetration depth vs time plots for Drop 5
fig, (ax1, ax2, ax3) = plt.subplots(3,sharex=True)
ax1.plot(drop5["Time (s)"], drop5["250g (m/s^2)"], marker = 11, color = "k")
ax1.set(ylabel="Deceleration (m/s^2)")

ax2.plot(drop5["Time (s)"], drop5['Velocity (m/s)'], marker = 11, color = "k")
ax2.set(ylabel="Velocity (m/s)")
ax3.plot(drop5["Time (s)"], drop5["Penetration Depth (m)"], marker = 11, color = "k")
ax3.set(ylabel="Penetration Depth (m)", xlabel="Time(s)")
ax1.plot(drop5["Time (s)"], drop5["250g (m/s^2)"], color = "k", marker = 11)
plt.suptitle(fileNum+' Drop 5')
plt.show()

#%%
#Deceleration,Velocity,and penetration depth vs time plots for Drop 6
fig, (ax1, ax2, ax3) = plt.subplots(3,sharex=True)
ax1.plot(drop6["Time (s)"], drop6["250g (m/s^2)"], marker = 11, color = "k")
ax1.set(ylabel="Deceleration (m/s^2)")

ax2.plot(drop6["Time (s)"], drop6['Velocity (m/s)'], marker = 11, color = "k")
ax2.set(ylabel="Velocity (m/s)")
ax3.plot(drop6["Time (s)"], drop6["Penetration Depth (m)"], marker = 11, color = "k")
ax3.set(ylabel="Penetration Depth (m)", xlabel="Time(s)")
ax1.plot(drop6["Time (s)"], drop6["250g (m/s^2)"], color = "k", marker = 11)
plt.suptitle(fileNum+' Drop 6')
plt.show()


#%%
#Deceleration and Velocity profile vs penetration depth for Drop 1


fig, (ay1, ay2) = plt.subplots(1, 2, sharey=True)
ay1.plot(drop1g[acc1Nameg], drop1["Penetration Depth (m)"]*100, color = "k", linestyle = "solid") #marker = 11
ay1.plot(drop1["Velocity (m/s)"], drop1["Penetration Depth (m)"]*100, color = "k", linestyle = "dashed")
ay1.set(xlabel="Deceleration (g) and Velocity (m/s)", ylabel="Penetration Depth (cm)")
ay1.invert_yaxis()
ay1.legend(["Deceleration (g)", "Velocity (m/s)"],prop={'size': 8})
ay1.title.set_text('Deceleration/Velocity')
ay1.grid(visible=True)

ay2.plot(drop1['QSBC (kPa) K=0.2'],drop1["Penetration Depth (m)"]*100)
ay2.plot(drop1['QSBC (kPa) K=0.4'],drop1["Penetration Depth (m)"]*100)
ay2.plot(drop1['QSBC (kPa) K=1.0'],drop1["Penetration Depth (m)"]*100)
ay2.plot(drop1['QSBC (kPa) K=1.5'],drop1["Penetration Depth (m)"]*100)
ay2.set(xlabel="QSBC (kPa)")
ay2.grid(visible=True)
ay2.title.set_text('QSBC')
ay2.legend(["K=0.2", "K=0.4", "K=1.0", "K=1.5"],prop={'size': 8})

plt.suptitle(fileNum+' Drop 1', fontsize=10)
plt.show()

title=fileNum+' Drop 1 Deceleration+Velocity Profile'
fig.savefig(title)

#%%
#Deceleration and Velocity profile vs penetration depth for Drop 2
fig, (ay1, ay2) = plt.subplots(1, 2, sharey=True)
ay1.plot(drop2g[acc1Nameg], drop2["Penetration Depth (m)"]*100, color = "k", linestyle = "solid") #marker = 11
ay1.plot(drop2["Velocity (m/s)"], drop2["Penetration Depth (m)"]*100, color = "k", linestyle = "dashed")
ay1.set(xlabel="Deceleration (g) and Velocity (m/s)", ylabel="Penetration Depth (cm)")
ay1.invert_yaxis()
ay1.legend(["Deceleration (g)", "Velocity (m/s)"],prop={'size': 8})
ay1.title.set_text('Deceleration/Velocity')
ay1.grid(visible=True)

ay2.plot(drop2['QSBC (kPa) K=0.2'],drop2["Penetration Depth (m)"]*100)
ay2.plot(drop2['QSBC (kPa) K=0.4'],drop2["Penetration Depth (m)"]*100)
ay2.plot(drop2['QSBC (kPa) K=1.0'],drop2["Penetration Depth (m)"]*100)
ay2.plot(drop2['QSBC (kPa) K=1.5'],drop2["Penetration Depth (m)"]*100)
ay2.set(xlabel="QSBC (kPa)")
ay2.grid(visible=True)
ay2.title.set_text('QSBC')
ay2.legend(["K=0.2", "K=0.4", "K=1.0", "K=1.5"],prop={'size': 8})

plt.suptitle(fileNum+' Drop 2', fontsize=10)
plt.show()

title=fileNum+' Drop 2 Deceleration+Velocity Profile'
fig.savefig(title)

#%%
#Deceleration and Velocity profile vs penetration depth for Drop 3
fig, (ay1, ay2) = plt.subplots(1, 2, sharey=True)
ay1.plot(drop3g[acc1Nameg], drop3["Penetration Depth (m)"]*100, color = "k", linestyle = "solid") #marker = 11
ay1.plot(drop3["Velocity (m/s)"], drop3["Penetration Depth (m)"]*100, color = "k", linestyle = "dashed")
ay1.set(xlabel="Deceleration (g) and Velocity (m/s)", ylabel="Penetration Depth (cm)")
ay1.invert_yaxis()
ay1.legend(["Deceleration (g)", "Velocity (m/s)"],prop={'size': 8})
ay1.title.set_text('Deceleration/Velocity')
ay1.grid(visible=True)

ay2.plot(drop3['QSBC (kPa) K=0.2'],drop3["Penetration Depth (m)"]*100)
ay2.plot(drop3['QSBC (kPa) K=0.4'],drop3["Penetration Depth (m)"]*100)
ay2.plot(drop3['QSBC (kPa) K=1.0'],drop3["Penetration Depth (m)"]*100)
ay2.plot(drop3['QSBC (kPa) K=1.5'],drop3["Penetration Depth (m)"]*100)
ay2.set(xlabel="QSBC (kPa)")
ay2.grid(visible=True)
ay2.title.set_text('QSBC')
ay2.legend(["K=0.2", "K=0.4", "K=1.0", "K=1.5"],prop={'size': 8})

plt.suptitle(fileNum+' Drop 3', fontsize=10)
plt.show()

title=fileNum+' Drop 3 Deceleration+Velocity Profile'
fig.savefig(title)


#%%
#Deceleration and Velocity profile vs penetration depth for Drop 4
fig, (ay1, ay2) = plt.subplots(1, 2, sharey=True)
ay1.plot(drop4g[acc1Nameg], drop4["Penetration Depth (m)"]*100, color = "k", linestyle = "solid") #marker = 11
ay1.plot(drop4["Velocity (m/s)"], drop4["Penetration Depth (m)"]*100, color = "k", linestyle = "dashed")
ay1.set(xlabel="Deceleration (g) and Velocity (m/s)", ylabel="Penetration Depth (cm)")
ay1.invert_yaxis()
ay1.legend(["Deceleration (g)", "Velocity (m/s)"],prop={'size': 8})
ay1.title.set_text('Deceleration/Velocity')
ay1.grid(visible=True)

ay2.plot(drop4['QSBC (kPa) K=0.2'],drop4["Penetration Depth (m)"]*100)
ay2.plot(drop4['QSBC (kPa) K=0.4'],drop4["Penetration Depth (m)"]*100)
ay2.plot(drop4['QSBC (kPa) K=1.0'],drop4["Penetration Depth (m)"]*100)
ay2.plot(drop4['QSBC (kPa) K=1.5'],drop4["Penetration Depth (m)"]*100)
ay2.set(xlabel="QSBC (kPa)")
ay2.grid(visible=True)
ay2.title.set_text('QSBC')
ay2.legend(["K=0.2", "K=0.4", "K=1.0", "K=1.5"],prop={'size': 8})

plt.suptitle(fileNum+' Drop 4', fontsize=10)
plt.show()

title=fileNum+' Drop 4 Deceleration+Velocity Profile'
fig.savefig(title)

#%%
#Deceleration and Velocity profile vs penetration depth for Drop 5
fig, (ay1, ay2) = plt.subplots(1, 2, sharey=True)
ay1.plot(drop5g[acc1Nameg], drop5["Penetration Depth (m)"]*100, color = "k", linestyle = "solid") #marker = 11
ay1.plot(drop5["Velocity (m/s)"], drop5["Penetration Depth (m)"]*100, color = "k", linestyle = "dashed")
ay1.set(xlabel="Deceleration (g) and Velocity (m/s)", ylabel="Penetration Depth (cm)")
ay1.invert_yaxis()
ay1.legend(["Deceleration (g)", "Velocity (m/s)"],prop={'size': 8})
ay1.title.set_text('Deceleration/Velocity')
ay1.grid(visible=True)

ay2.plot(drop5['QSBC (kPa) K=0.2'],drop5["Penetration Depth (m)"]*100)
ay2.plot(drop5['QSBC (kPa) K=0.4'],drop5["Penetration Depth (m)"]*100)
ay2.plot(drop5['QSBC (kPa) K=1.0'],drop5["Penetration Depth (m)"]*100)
ay2.plot(drop5['QSBC (kPa) K=1.5'],drop5["Penetration Depth (m)"]*100)
ay2.set(xlabel="QSBC (kPa)")
ay2.grid(visible=True)
ay2.title.set_text('QSBC')
ay2.legend(["K=0.2", "K=0.4", "K=1.0", "K=1.5"],prop={'size': 8})

plt.suptitle(fileNum+' Drop 5', fontsize=10)
plt.show()

title=fileNum+' Drop 5 Deceleration+Velocity Profile'
fig.savefig(title)

#%%
#Deceleration and Velocity profile vs penetration depth for Drop 6
fig, (ay1, ay2) = plt.subplots(1, 2, sharey=True)
ay1.plot(drop6g[acc1Nameg], drop6["Penetration Depth (m)"]*100, color = "k", linestyle = "solid") #marker = 11
ay1.plot(drop6["Velocity (m/s)"], drop6["Penetration Depth (m)"]*100, color = "k", linestyle = "dashed")
ay1.set(xlabel="Deceleration (g) and Velocity (m/s)", ylabel="Penetration Depth (cm)")
ay1.invert_yaxis()
ay1.legend(["Deceleration (g)", "Velocity (m/s)"],prop={'size': 8})
ay1.title.set_text('Deceleration/Velocity')
ay1.grid(visible=True)

ay2.plot(drop6['QSBC (kPa) K=0.2'],drop6["Penetration Depth (m)"]*100)
ay2.plot(drop6['QSBC (kPa) K=0.4'],drop6["Penetration Depth (m)"]*100)
ay2.plot(drop6['QSBC (kPa) K=1.0'],drop6["Penetration Depth (m)"]*100)
ay2.plot(drop6['QSBC (kPa) K=1.5'],drop6["Penetration Depth (m)"]*100)
ay2.set(xlabel="QSBC (kPa)")
ay2.grid(visible=True)
ay2.title.set_text('QSBC')
ay2.legend(["K=0.2", "K=0.4", "K=1.0", "K=1.5"],prop={'size': 8})

plt.suptitle(fileNum+' Drop 6', fontsize=10)
plt.show()

title=fileNum+' Drop 6 Deceleration+Velocity Profile'
fig.savefig(title)


#%%
#Deceleration and Velocity profile vs penetration depth for Drop 7
fig, (ay1, ay2) = plt.subplots(1, 2, sharey=True)
ay1.plot(drop7g[acc1Nameg], drop7["Penetration Depth (m)"]*100, color = "k", linestyle = "solid") #marker = 11
ay1.plot(drop7["Velocity (m/s)"], drop7["Penetration Depth (m)"]*100, color = "k", linestyle = "dashed")
ay1.set(xlabel="Deceleration (g) and Velocity (m/s)", ylabel="Penetration Depth (cm)")
ay1.invert_yaxis()
ay1.legend(["Deceleration (g)", "Velocity (m/s)"],prop={'size': 8})
ay1.title.set_text('Deceleration/Velocity')
ay1.grid(visible=True)

ay2.plot(drop7['QSBC (kPa) K=0.2'],drop7["Penetration Depth (m)"]*100)
ay2.plot(drop7['QSBC (kPa) K=0.4'],drop7["Penetration Depth (m)"]*100)
ay2.plot(drop7['QSBC (kPa) K=1.0'],drop7["Penetration Depth (m)"]*100)
ay2.plot(drop7['QSBC (kPa) K=1.5'],drop7["Penetration Depth (m)"]*100)
ay2.set(xlabel="QSBC (kPa)")
ay2.grid(visible=True)
ay2.title.set_text('QSBC')
ay2.legend(["K=0.2", "K=0.4", "K=1.0", "K=1.5"],prop={'size': 8})

plt.suptitle(fileNum+' Drop 7', fontsize=10)
plt.show()

title=fileNum+' Drop 7 Deceleration+Velocity Profile'
fig.savefig(title)


xvalues=np.linspace(0,len(dg)-1,len(dg))

temp_dg=dg['250g (g)']
test=np.zeros(len(temp_dg))
