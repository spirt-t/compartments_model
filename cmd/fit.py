from scipy.optimize import curve_fit
import numpy as np
from numpy import array
from tkinter import filedialog as fd

def peppas(x, k, n):
    return 3.8*k*pow(x,n)

def read_data():
    filename = fd.askopenfilename()
    x=[]
    y=[]
    with open(filename, 'r') as f:
        for line in f:
            x.append(float(line.split('\t')[0]))
            y.append(float(line.split('\t')[1]))
    x=array(x) # probabilities from data
    y=array(y)
    f.close()
    return x, y

x, y = read_data()
pars, cov = curve_fit(f=peppas, xdata=x, ydata=y, p0=[1, 0.5], bounds=(-np.inf, np.inf))
print(pars)