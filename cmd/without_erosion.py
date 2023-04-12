# -*- coding: utf-8 -*-

### INSTRUCTIONS ###

# 1. Download the "Enthought Canopy" IDE from the 
#    https://www.enthought.com/product/canopy/ and install it.
# 2. Run the Enthought Canopy and click to the "Package Manager" icon.
# 3. Go to the tab "Installed" at the left and check for the following packages: 
#   matplotlib, scipy, numpy, tk_headers.
# 4. If you can not find them, go to the tab "Available" at the left.
#    If you find them, go to item 6.
# 5. Find and install the following packages in the list: 
#    matplotlib, scipy, numpy, tk_headers.
# 6. Close the Package Manager and open the "Editor" icon.
# 7. Select the program file from the computer or create new file and copy 
#    text of the program. 
# 8. Run the project.
# 9. Select the geometry of the system in the window that opens.
# 10. Enter the other initial paremeters (additional information 
#    about each parameter when you click on the button with a question mark "?").
#    Test distribution and release data for each geometry you can also find 
#    in the supplementary materials.
# 11. Press the "ENTER" button.
# 12. Wait for the result.
# 13. Save the results, if necessary.


# -*- coding: utf-8 -*-
import matplotlib
import matplotlib.pyplot as plt
from scipy.special import iv
from numpy import tanh
from scipy.optimize import leastsq
from numpy import array, any, sqrt, pi, sin, tan, exp, linspace, cosh, sinh
from tkinter import *
from tkinter import ttk
from tkinter.ttk import *
from tkinter import messagebox as mb
from tkinter import filedialog as fd
import tkinter.font as tkFont

### CHOOSE EXPERIMENTAL DATA 
### AND DEFINE STARTING VALUES 




# Laplace function: Qdi(s) (amount desorbed for a single radius fibre) to be defined below

def save_file(rm, rad, sn, cv, td, Dd, st):
    f = fd.asksaveasfile(mode='w', defaultextension=".txt")
    if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
        return
    f.write("Rmean = %.2f %s\n" % (rm, rad))
    f.write("Standard deviation of a random variable = %.2f \n" % (sn))
    f.write("Coefficient of variation (CV) = %.2f \n" % (cv))
    f.write("Characteristic time of diffusion (td) = %.2f %s\n" % (td, st))
    f.write("Diffusion coefficient = %.8f %s^2/%s\n" % (Dd, rad, st))
    f.close()


def saveSimulationPoints(xlist, ylist):
    fl = fd.asksaveasfile(mode='w', defaultextension=".txt")
    if fl is None:  # asksaveasfile return `None` if dialog closed with "cancel".
        return

    for x, y in zip(xlist, ylist):
        fl.write("%s\t%s\n" % (x, y))

    fl.close()

def next_step():
    #print("linspace start")
    tt=linspace(t0,tmax,100)
    #print("linspace end")
    
    # Ri from data 
    
    
    cv=[] 	# coefficient of variation cv=sn/rm
    sn=[] 	# the standard deviation of a random variable
    
    sn=sqrt(sum(p*(ri-rm)**2))
    cv=sn/rm
    pname = (['td','m0','q0'])
    par0 = array([td_0 , m0_0, q0_0]) #values for original guess
    plsq = leastsq(residuals, par0, args=(y, x), maxfev=2000) # Perform fitting 
    
    Dd = (((rm**2))/(plsq[0][0]))
    plt.rc('font',family='Times New Roman')
    plt.plot(x,y,'ko', label='%s' % (leg_exp)) # here you can change the parameter 'k'
                    # which which identifies the color of the experimental points. 
                    # For example, 'r' is red, 'g' is green, 'b' is blue,
                    # 'c' is azure, 'k' is black.
                    
                    # 'o' near the 'k' indicates the shape of the experimental point.
                    # You can change it, for example:
                    # '>', '<', 'v' and '^' are the various triangles, 'D' is the rhombus,
                    # 's' is the squares, '*' is the stars, 'x' is the X
    yy = peval(tt,plsq[0],p)

    plt.plot(tt,yy,'k-', label='%s' % (leg_mod)) # Here you can also
                    # change the black color 'k' like in the previous case,
                    # and also you can use different tipes of lines (here symbol next to 'k'):
                    # '-' is the solid line, '--' is the dashed line, ':' is the dotted line,
                    # '-.' is the dashed-dotted line

    plt.xlabel('Time, %s' %st)
    plt.ylabel('Drug release')
    plt.show()
    plt.legend()
    
    #print("new TK window start")
    root1 = Tk()
    if var.get() == 1:
        name = "membrane"
    elif var.get() == 2:
        name = "cylinder"
    elif var.get() == 3:
        name = "sphere"
    name = "cylinder"
    root1.title("Results for %s" %name)
    xw = (root1.winfo_screenwidth()) / 6
    yw = (root1.winfo_screenheight()) / 4
    root1.geometry("%dx%d" % (xw, yw))
    x22 = (root1.winfo_screenwidth() - root1.winfo_reqwidth()) / 1.5
    y22 = (root1.winfo_screenheight() - root1.winfo_reqheight()) / 4
    root1.wm_geometry("+%d+%d" % (x22, y22))
    
    label1 = Label(root1, text="Rmean = %.2f %s" % (rm, rad))
    label1.grid(row=0, column=0, pady=6)
    label2 = Label(root1, text="Standard deviation of a random variable = %.2f " % (sn))
    label2.grid(row=1, column=0, pady=6)
    label3 = Label(root1, text="Coefficient of variation (CV) = %.2f " % (cv))
    label3.grid(row=2, column=0, pady=6)
    label5 = Label(root1, text="Characteristic time of diﬀusion (td) = %.2f %s" % (plsq[0][0], st))
    label5.grid(row=3, column=0, pady=6) 
    label4 = Label(root1, text="Diffusion coefficient = %.8f %s^2/%s" % (Dd, rad, st))
    label4.grid(row=4, column=0, pady=6)
    but = ttk.Button(root1, text='Save', style = 'exit.TButton', 
                    command=lambda:[save_file(rm, rad, sn, cv, plsq[0][0], Dd, st)])
    but.grid(row=5, column=0, pady=6)

    saveSimulationPoints(tt, yy)
    #print("new TK window end")
    root1.mainloop()
    

def quit():
    #print("quit start")
    root.destroy()
    #print("quit end")
    next_step()

def Qdi(s,tdi):
   # print("Qdi start")
    sr=sqrt(s*tdi)
    if var.get() == 1:
        q = (2/(s*sr))*((cosh(sr)-1)/sinh(sr))
        name = "membrane"
    elif var.get() == 2:
        q = 2*iv(1,sr)/(iv(0,sr)*s*sr)
        name = "cylinder"
    elif var.get() == 3:
        q = (3/(s*sr))*((1/tanh(sr))-(1/sr))
        name = "sphere"
    #print("Qdi end")
    return (q)
# taking r distribution into account (summing all spheres with different radiuses):
def Qdl(s,td,m0,q0,p): 
    #print("Qdl start")
    if var.get() == 1:
        tm2=ri_rm # = (r[i]/rm)
    elif var.get() == 2:
        tm2=ri_rm**2 # = (r[i]/rm)^2
    elif var.get() == 3:
        tm2=ri_rm**3 	# = (r[i]/rm)^3
    tdi=td*tm2	 # produce n tdi values based on ri ( tdi[i]=(r[i]/rm)^3*rm^2/D 
    # = td*(r[i]/rm)^3)
    pr2w=tm2*p
    pr2w=pr2w/sum(pr2w) # pr2w=tm2*p/sum(tm2*p)
    tm=Qdi(s,tdi)*pr2w # weighted laplace function based on p(ri)*(ri/rw)^3=pr2w(ri); rw^2=sum(ri^2*p(ri)) 
    #print("Qdl end")
    return sum(tm)*m0
    
    # function in time domain Q(t); Numerical Laplace transform 
def Qd(t,td,m0,q0,p):
    #print("Qd start")
    #These constants are dependent on the function and might be different if func changes
    shift=0.01
    N = 24
    if any(t == 0.):
        mb.showerror("Error", "Inverse transform can not be calculated for t=0!")
    # Initiate the stepsize
    h = 2*pi/N
    ans = 0.0
    # parameters from
    # T. Schmelzer, L.N. Trefethen, SIAM J. Numer. Anal. 45 (2007) 558-571
    c1 = 0.5017
    c2 = 0.6407
    c3 = 0.6122
    c4 = 0+0.2645j
    # The for loop is evaluating the Laplace inversion at each point theta i
    # which is based on the trapezoidal rule
    for k in range(N):
        theta = -pi + (k+0.5)*h
        z = shift + N/t*(c1*theta/tan(c2*theta) - c3 + c4*theta)
        dz = N/t * (-c1*c2*theta/sin(c2*theta)**2 + c1/tan(c2*theta)+c4)
        ans += exp(z*t)*Qdl(z,td,m0,q0,p)*dz
    #print("Qd end")
    return ((h/(2j*pi))*ans).real+q0

def residuals(par, y, x):
    #print("res start") 
    err = y-peval(x,par,p);
    #print("res end") 
    return array(err)

def peval(x, par,p):
    #print("pev start")
    tm=[]
    for i in range(len(x)): tm.append(Qd(x[i],par[0],par[1],par[2],p))
    #print("pev end")
    return array(tm) 
        
def check(entry_td, entry_m0, entry_q0, entry_st, 
    entry_exp, entry_mod, entry_rad,rm):
    #print("check start")
    try:
        global td_0
        global D_0
        D_0 = entry_td.get()
        D_0 = float(D_0)

        td_0 = rm**2/D_0
        global m0_0
        m0_0 = entry_m0.get()

        global q0_0
        q0_0 = entry_q0.get()

        td_0 = float(td_0)
        m0_0 = float(m0_0)
        q0_0 = float(q0_0)
        global st
        st = entry_st.get()

        global leg_exp
        leg_exp = entry_exp.get()

        global leg_mod
        leg_mod = entry_mod.get()

        global rad
        rad = entry_rad.get()
    except ValueError:
        mb.showerror("Error", "A number must be entered!") 
    #print("check end")

def insertText1():
    #print("text1 start")
    file_name = fd.askopenfilename()
    global ri
    global p
    ri=[]
    p=[]
    with open(file_name, 'r') as f:
        for line in f:
            print(line)
            print(line.split('\t')[0])
            ri.append(float(line.split('\t')[0]))
            p.append(float(line.split('\t')[1]))
    p=array(p) # probabilities from data
    p=p/sum(p)
    ri=array(ri)
    f.close()
    global rm
    rm=sum(p*ri)
    global ri_rm
    ri_rm=ri/rm
    #print("text1 end")
    
def insertText2():
    #print("text2 start")
    file_name = fd.askopenfilename()
    global x
    global y
    x1=[]
    y1=[]
    with open(file_name, 'r') as f:
        for line in f:
            x1.append(float(line.split('\t')[0]))
            y1.append(float(line.split('\t')[1]))
    x=array(x1) # probabilities from data
    y=array(y1)
    global t0
    global tmax
    t0 = x[0]
    if t0 == 0: t0=0.01
    tmax = x[-1]
    f.close()
    #print("text2 end")

def Help_info(text):
    note = Toplevel(root) 
    xx = (note.winfo_screenwidth() - note.winfo_reqwidth()) / 2
    yy = (note.winfo_screenheight() - note.winfo_reqheight()) / 2
    note.wm_geometry("+%d+%d" % (xx, yy))   
    info = Label(note, text=text)
    info.pack()
            
root = Tk()
root.style = ttk.Style()
root.title("Size Distribution Model")
x22 = (root.winfo_screenwidth() - root.winfo_reqwidth()) / 2
y22 = (root.winfo_screenheight() - root.winfo_reqheight()) / 4
root.wm_geometry("+%d+%d" % (x22, y22))

arial = tkFont.Font(family='Times New Roman', size=11, weight=tkFont.BOLD)
infofont = tkFont.Font(family='Times New Roman', size=9, weight=tkFont.BOLD)

note = Notebook(root)
tab1 = Frame(note)
tab2 = Frame(note)
tab3 = Frame(note)
note.add(tab1, text = "Fitting", compound=LEFT)
note.add(tab2, text = "Information")
note.add(tab3, text = "Help")
note.grid(row=0, column=0, columnspan=2)

######### TAB 1 Fitting ##########


label2 = Label(tab1, text='Please, select the system geometry:', font=arial)
label2.grid(row=1, column=0, columnspan=3, pady=5)
root2 = Label(tab1)
root2.grid(row=2, column=0, columnspan=3, pady=5)
var=IntVar()
var.set(0)
rbutton1=Radiobutton(root2,text='Membrane', variable=var, value=1, command=lambda:[radioClick(1)])
rbutton2=Radiobutton(root2,text='Cylinder', variable=var, value=2, command=lambda:[radioClick(2)])
rbutton3=Radiobutton(root2,text='Sphere', variable=var, value=3, command=lambda:[radioClick(3)])
rbutton1.grid(row=0, column=0)
rbutton2.grid(row=0, column=1)
rbutton3.grid(row=0, column=2)
def radioClick(val):
    var.set(val)
    #print(var.get())
    

main_label = Label(tab1, text='Please, enter the estimated initial values.', font=arial)
main_label.grid(row=3, column=0, columnspan=3, pady=5)


td_label = Label(tab1, text='Enter estimated diffusion\n coefficient'+ 
                            '(in dimension\ncorresponding to time\nand radius dimension):', 
                            anchor='center')
td_label.grid(row=7, column=0)
td_info = ttk.Button(tab1, text="?", width = 3, command=lambda:[Help_info(
                    '''Estimated value,
if you do not know this value,
enter the 1.''')])
td_info.grid(row=7, column=1)

m0_label = Label(tab1, text='Enter loaded \n amount of drug (M0):')
m0_label.grid(row=8, column=0)
m0_info = ttk.Button(tab1, text="?", width = 3, command=lambda:[Help_info(
                    '''Estimated amount of loaded drug. 
You can use the relative value 
in the persentage (%), then the 
data will be presented on the
graph also in a percentage.''')])
m0_info.grid(row=8, column=1)

q0_label = Label(tab1, text='Enter amount of drug \n on the surface (Q0):')
q0_label.grid(row=9, column=0)
q0_info = ttk.Button(tab1, text="?", width = 3, command=lambda:[Help_info(
                    '''Estimated amount of drug, located
on the vehicle surface.
If you do not know this value,
you can set it to zero (0).
If the M0 was in the %, this value
also must be in the %.''')])
q0_info.grid(row=9, column=1)
entry_td = Entry(tab1)
entry_td.grid(row=7, column=2)
entry_m0 = Entry(tab1)
entry_m0.grid(row=8, column=2)
entry_q0 = Entry(tab1)
entry_q0.grid(row=9, column=2)

b1 = ttk.Button(tab1, text="Open size distribution \n data file (.txt)", 
                command=lambda:[insertText1()])
b1.grid(row=6, column = 0, pady=12)
b2 = ttk.Button(tab1, text="Open drug release \n data file (.txt)", 
command=lambda:[insertText2()]).grid(row=6, column = 2,pady=12)

tu_label = Label(tab1, text='Time dimension:')
tu_label.grid(row=4, column=0, pady=4)
tu_info = ttk.Button(tab1, text="?", width = 3, command=lambda:[Help_info(
'Example: days, h, min, s...')])
tu_info.grid(row=4, column=1)

entry_st = Entry(tab1)
entry_st.grid(row=4, column=2,pady=4)
exp_label = Label(tab1, text='Experimental label:')
exp_label.grid(row=11, column=0, pady=4)
exp_info = ttk.Button(tab1, text="?", width = 3, command=lambda:[Help_info(
'''The signature of the experimental data,
that will be displayed in the legend 
on the graph. 
Example: Experimental data''')])
exp_info.grid(row=11, column=1)

mod_label = Label(tab1, text='Simulation label:')
mod_label.grid(row=12, column=0, pady=4)
mod_info = ttk.Button(tab1, text="?", width = 3, command=lambda:[Help_info(
'''The signature of the simulation data,
that will be displayed in the legend 
on the graph.
Example: Simulation data''')])
mod_info.grid(row=12, column=1)

entry_exp = Entry(tab1)
entry_exp.grid(row=11, column=2, pady=4)
entry_mod = Entry(tab1)
entry_mod.grid(row=12, column=2, pady=4)
rad_label = Label(tab1, text='Radius (thickness) dimension :')
rad_info = ttk.Button(tab1, text="?", width = 3, command=lambda:[Help_info(
'''Units of your size distribution. 
Example: nm, mkm...''')])
rad_info.grid(row=5, column=1)

rad_label.grid(row=5, column=0, pady=4)
entry_rad = Entry(tab1)
entry_rad.grid(row=5, column=2, pady=4)


ttk.Style().configure('exit.TButton', background='green')

######### TAB 2 INFORMATION ##########

info_label = Label(tab2, text="This program based on the size distribution model "+ 
                    "and was developed by the Griffith University and "+
                    "Tomsk Polytechnic University research team \n "+
                    "(Anissimov Y.G., Spiridonova T.I., Petlin D.G.)", 
                    wraplength = root.winfo_reqwidth()*1.3).pack(pady = 10)
info1_label = Label(tab2, text='The contact e-mails:\ny.anissimov@griffith.edu.au'+ 
                                                       '\nspiridonovatis2@gmail.com', 
                    wraplength = root.winfo_reqwidth()*1.3, anchor='w').pack(pady = 10)


######### TAB 3 HELP ##########

info_label = Label(tab3, text='Instruction', font = arial).pack(pady=10)
info2_label = Label(tab3, text='1. Select the geometry of your model: '+
                    'membrane, cylinder (in case of fibrous scaffolds)'+
                    'or sphere (for micro-and nanoparticles).', 
                    wraplength = root.winfo_reqwidth()*1.4, anchor='w').pack(pady=2)
info8_label = Label(tab3, text='2. Enter the time dimension (for example, h or min), '+
                    'as well as the radius or thickness dimension (for example, nm or mkm).', 
                    wraplength = root.winfo_reqwidth()*1.4, anchor='w').pack(pady=2)                   
info3_label = Label(tab3, text='3. Enter the initial values in the appropriate fields: '+
                    'supposed coefficient of diffusion (at least approximately)'+
                    'in dimension  corresponding to the previously entered'+
                    'time and radius (or thickness) dimensions. '+
                    'Then enter the supposed amount of drug into your delivery vehicle and '+
                    'on the vehicle surface.\n If you will enter the last two values as ' + 
                    'a percentage, the graph will also show the percentage of the drug released.', 
                    wraplength = root.winfo_reqwidth()*1.4, anchor='w').pack(pady=2)
info4_label = Label(tab3, text='3. Press the "Open size distribution data file" button '+
                    'and select the .txt file with the size distribution data on your computer.',
                    wraplength = root.winfo_reqwidth()*1.4, anchor='w').pack(pady=2)
info5_label = Label(tab3, text='4. Press the "Open drug release data file" button '+
                    'and select the .txt file with the experimental data of '+
                    'time-dependance drud release on your computer.',
                    wraplength = root.winfo_reqwidth()*1.4, anchor='w').pack(pady=2)
info7_label = Label(tab3, text='5. Enter the units of time, '+
                    'texts of legends for experimental and simulation data, '+
                    'as well as the units of size.',
                    wraplength = root.winfo_reqwidth()*1.4, anchor='w').pack(pady=2)
info7_label = Label(tab3, text='6. Press the "ENTER" button to perform the fitting.',
                    wraplength = root.winfo_reqwidth()*1.4, anchor='w').pack(pady=2)
info7_label = Label(tab3, text='If you have some quastions, please, '+
                    'write to the e-meil: spiridonovatis2@gmail.com.',
                    wraplength = root.winfo_reqwidth()*1.4, anchor='w').pack(pady=2)
b5 = ttk.Button(tab1, text='ENTER', style = 'exit.TButton', command=lambda:
[check(entry_td, entry_m0, entry_q0, entry_st, 
entry_exp, entry_mod, entry_rad,rm),quit()])
b5.grid(row=14, column=2, pady=12)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
root.mainloop()

