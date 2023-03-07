import math

#USER INPUTS

moist = 180.4 #moisture in mV
gammad = 15.3 #dry unit weight in kN/m^3
d60 = .413 # D60 in mm
cu = 1.69 # coefficient of uniformity (unitless)

emin = 0.468 #minimum void ratio
emax = 0.753 #maximum void ratio
gammadmin = 15.196 #min dry unit weight (kN/m^3)
gammadmax = 18.143 #max dry unit weight (kN/m^3)
gammasea = 10.05 #unit weight of seawater (kN/m^3)
Gs = 2.65 #specific gravity

#weight volume calculations
def wtvol():
    global S
    global gammabulk
    global n
    global Dr
    V = moist/1000 #moisture in volts
    theta = 2.832*V**4 - 3.6426*V**3 + 1.3985*V**2 + 0.4112*V - 0.0149 #Calibrated volumetric water content
    e = Gs*gammasea/gammad - 1 #void ratio
    n = e/(1+e) #porosity
    Dr = (emax-e)/(emax-emin) #relative density (%)
    S = theta/n*100 #degree of saturation (%)
    w = (S*e)/Gs #gravitational water content (%)
    gammabulk = gammad*(1+w/100) #bulk unit weight (kN/m^3)

#matric suction calc
def matric_suction():
    global ms
    thetar = 2.1 #residual volumetric water content (%)
    c1 = 1.07
    c2 = 12.07
    sal = 30 #salinity (g/kg)
    t = 20 #temperature (deg. C)
    st1 = 0.073 #surface tension (N/m)
    st2 = st1*(1+(3.766*10**-4)*sal+(2.347*10**-6)*sal*t) #surface tension (N/m)
    Sr = thetar/n #residual degress of saturation (%)
    coefn = (c1/math.log10(cu))+1
    coefa = (c2*st2/(d60/1000))/1000 #kPa
    Se = (S-Sr)/(1-Sr/100) #effective degree of saturation (%)
    ms = coefa*((Se/100)**(coefn/(1-coefn)) - 1)**(1/coefn) #matric suction (kPa)

def Duncan():#Duncan Correlation (phi')
    global phi
    A = 34
    B = 10
    C = 3
    D = 2
    sigN= .05 * gammabulk #normal stress for duncan correlation, assuming D = .05m (kPa)
    pa = 1.03 #atmospheric pressure, kPa
    phi = A + B* Dr - (C + D*Dr)*math.log10(sigN/pa)

def Vanapalli(): #Bearing Capacity (Vanapalli & Mohammed)
    c = 0 #effective cohesion, kPa
    psi = 1
    D = 0.05
    B = 0.0875
    Kp = math.tan(math.radians(45+phi/2))**2
    Nq = math.e**(math.pi*math.tan(math.radians(phi)))*Kp
    Nc = (Nq-1)*(1/math.tan(math.radians(phi)))
    Ngam = 2*(Nq+1)*math.tan(math.radians(phi))
    Sc = 1+Nq/Nc
    Sgam = 0.6
    qu = (c + ms*(S/100**psi)*math.tan(math.radians(phi)))*Nc*Sc+0.5*gammabulk*B*Ngam*Sgam #partially saturated bearing capacity, kPa
    print(qu)


wtvol()
matric_suction()
Duncan()
Vanapalli()



