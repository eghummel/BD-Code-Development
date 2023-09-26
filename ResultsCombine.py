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


outputPath = Path("H:/.shortcut-targets-by-id/1aF9t2aiRGWTftJMZFAOBixqvQniFBjnb/Duck  2023/Data/Intertidal/BlueDrop, Samples & Moisture Gage/08March23/BlueDrop Processing 3.8.23 - July 2023.xlsx") #  Path to pre-existing Excel File

df = pd.concat(pd.read_excel(outputPath, sheet_name=None), ignore_index=True)

print(df)

with pd.ExcelWriter(outputPath, mode="a") as writer:
        df.to_excel(writer, sheet_name = 'All Data', index=False)
