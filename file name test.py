import numpy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

fileNum = "0FF1"

name = [*fileNum]
print(name)

def newName(): #find the next sequential bin name a=existing last variable, b=new last variable
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

newName()