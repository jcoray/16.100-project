# Aircraft.py
# Created by Team 2 - Jakob, Jose, and Austin on 11/27/18

import math
import warnings
from skaero.atmosphere import coesa
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np

class Aircraft(object):
    '''An aircraft container'''
    def __init__(self, components, ap=None, fp=None):
        ''' 
        fp (dict): SI units. Includes: 
        woe (operating empty MASS), wfuelland (fuel to land MASS), 
        wpay (Payload MASS), wfuel (Fuel MASS for cruise),
        R (range), Minf (Mach), ainf (speed o' sound), 
        g (like in F = m*g), rho (density inf), nu (kinematic visc)
        TSFC (Thrust Spec Fuel Consumption; Here because the notes gave us a value to use)

        '''

        # _, T, _, rho = coesa.table(alt)
        # alt -> rho, nu = <mu> / rho, ainf = sqrt(<1.4*287.5>*T)
        
        # By default, assume a 737 max

        # Aircraft parameters (all assumed constant)
        if ap == None:
            self.ap = {'woe':45e3, 'wfuelland':2300, 'wpay':20e3, \
            'R':6500e3, 'g':9.81, 'Sref':127}
        else: self.ap = ap

        # Flight parameters (all assumed variable)
        if fp == None:
            self.fp = {'wfuel':15800, 'Minf':0.78, 'alt': 10000, \
            'TSFC': 1.42e-5}
        else: 
            self.fp = fp

        self.update_fp(self.fp)
            
        self.components = {c.name: c for c in components}
        return

    def update_fp(self, fp):
        self.fp = fp

        #print("New altitude:", self.fp['alt'])
        _, Temp, _, rho = coesa.table(self.fp['alt'])
         
        if True: # Dumb if statement, but keeps indentation   
        #if 'rho'  not in self.fp:
            self.fp['rho'] = rho
        #if 'nu'   not in self.fp:
            mu = 1.983e-5 # Dynamic Viscosity
            self.fp['nu']  = mu / self.fp['rho']
        #if 'ainf' not in self.fp:
            self.fp['ainf']= math.sqrt(1.4 * 287.058 * Temp)
        #if 'Vinf' not in self.fp:
            self.fp['Vinf'] = self.fp['Minf']*self.fp['ainf']
        #if 'qinf' not in self.fp:
            self.fp['qinf'] = 1/2*self.fp['rho']*self.fp['Vinf']**2
        #if 'T' not in self.fp:
            self.fp['T'] = self.ap['R'] / self.fp['Vinf'] # Total cruise time
        #if 'wfinal' not in self.fp:
            self.fp['wfinal'] = self.ap['woe'] + self.ap['wfuelland'] + self.ap['wpay']
        #~ if 'winit'  not in self.fp:
            #~ self.fp['winit'] = self.fp['wfinal'] + self.fp['wfuel']
        # TODO replace all instances of winit with wfinal + wfuel
        #print("New fp:", self.fp)

    
    def lift_at_time(self, t, wfuel = None):
        ''' Find lift at time t after start of cruise flight '''
        if wfuel == None: 
            wfuel = self.fp['wfuel']
        if (t < 0) or (t > self.fp['T']):
            raise ValueError("t must be positive and less than total cruise time.")
        lift = self.ap['g'] * (wfuel+self.fp['wfinal']  - (wfuel/self.fp['T'])*t)
        return lift

    def averageLift(self,wfuel = None):
        ''' Calculates average lift using initial and final weights '''
        if wfuel == None: 
            wfuel = self.fp['wfuel']
        Lavg = self.ap['g']/2*(2*self.fp['wfinal'] + wfuel) 
        return Lavg
    
    def profileDrag(self):
        ''' Calculates profile drag using Shevell's method. Dp ~ Dform + Dfriction '''
        Dp = 0
        for component in self.components.values():
            Rei = self.fp['Vinf'] * component.length2 /self.fp['nu'] # length2 = li or ci

            
            Cfi = 0
            if Rei > 1e6: # if Turbulent 
                Cfi = 0.455 * math.log10(Rei)**(-2.58)
            else: # if laminar
                warnings.warn('Flow assumed to be laminar. Component Re = %g' % Rei)
                Cfi = 1.328 * Rei**(-.5)
                 # FIX ME assumes all components transition at the same Re
                 # FIXED: Arthur said Re > 1e6 is very likely to transition

            Kfi = 0
            if component.isairfoil:
                tici = component.t/component.c
                Kfi = 1 + 2.0*tici + 60*(tici)**4
            else:
                dili = component.d/component.l
                Kfi = 1 + 1.5*(dili)**1.5 + 7*(dili)**3
            Dpi = self.fp['qinf'] * Kfi * Cfi * component.awet
            Dp += Dpi 
        
        return Dp
        
    def inducedDrag(self, lift=None):
        '''
        Calculates induced drag during cruise
        Assumes: 
            1) wings generate all lift. 
            2) e = 1
        '''
        if lift == None:
            lift = self.averageLift();
        e = 1
        Di = (lift/self.components['wing'].span)**2 / (math.pi*e*self.fp['qinf'])
        return Di
        
    def waveDrag(self):
        '''
        Calculates wave drag. Assumes 0
        '''
        Dw = 0
        return Dw
        
    def drag_at_time(self, t, wfuel = None):
        ''' Find drag at time t after start of cruise flight '''
        if (t < 0) or (t > self.fp['T']):
            raise ValueError("t must be positive and less than total cruise time.")
        lift = self.lift_at_time(t, wfuel)
        drag = self.profileDrag() + self.waveDrag() + self.inducedDrag(lift=lift)
        return drag

    def averageDrag(self, wfuel= None, nsamples=50):
        # Values of drag average when using averaged lift vs iterating through time
        # are not that different.
        #return self.profileDrag() + self.inducedDrag() + self.waveDrag()
        if wfuel == None: 
            wfuel = self.fp['wfuel']
        drag_samples = []
        for i in range(nsamples):
            t = self.fp['T'] * (i / (nsamples-1))
            drag = self.drag_at_time(t, wfuel)
            drag_samples.append(drag)
        mean_drag = sum(drag_samples) / nsamples
        return mean_drag
    
    def averageLD(self, wfuel = None, nsamples=50):
        # Similar to averageDrag, iterating through L/D values
        # vs using L and D's averages hadly changes anything
        #return self.averageLift()/self.averageDrag()
        ld_samples = []
        for i in range(nsamples):
            t = self.fp['T'] * (i / (nsamples-1))
            ld = self.lift_at_time(t, wfuel) / self.drag_at_time(t, wfuel)
            ld_samples.append(ld)
        mean_ld = sum(ld_samples) / nsamples
        return mean_ld
        
    # TSFC is assumed using the same logic with which we assume wfuel
    # In other words, the notes told us to assume TSFC = 1.42e-5
    # NOTE: winit and wfinal are now in the self.fp dictionary
    def bregue_wfuel(self, ld=None):
        ''' Estimates required wfuel from Bregue Range Equation model. '''
        if ld == None: 
            ld = self.averageLD()
        return self.fp['wfinal'] * (math.exp(self.fp['TSFC']*self.ap['g']*self.ap['R']
                                    / (self.fp['Vinf']*ld)   ) - 1)
    
    def findWfuel(self):
        wld = self.fp['wfuel']
        wbr = self.bregue_wfuel()
        itercount = 0
        while abs(wbr - wld) > .01:
            # Hope convergance 
            wld = wbr
            ld = self.averageLD(wld)
            wbr = self.bregue_wfuel(ld)
            itercount += 1
            #~ print('ld',ld)
            #~ print('wbr',wbr)
        #print('itercount', itercount)
        self.fp['wfuel'] = wbr
        return wbr

class Component(object):
    def __init__(self, name, isairfoil, length1, length2, awet, span = None):
        ''' name (str): e.g. Main wing, horizontal tail
            isairfoil (bool):   True if Airfoil, False if axisymmetric body
            length1 (float):    t (max thickness) if airfoil, d (max diameter) if axisym body
            length2 (float):    c (chord) if airfoil, l (length) if axisym body
            awet (float):       Wetted Area
            span (float):       Wing span. Used only for wings.
        '''
        self.name = name
        self.isairfoil = isairfoil
        self.length1 = length1; self.length2 = length2;
        if self.isairfoil:
            self.t = length1; self.c = length2;
        else:
            self.d = length1; self.l = length2;
        self.awet = awet
        if span != None:
            self.span = span
        return
       
       
def alt_fuel_plot(aircraft, max_alt, max_wfuel, max_thrust):
    alt = np.linspace(0, max_alt, 100)
    
    wfuel_list = list()
    for alti in alt:
        aircraft.fp['alt'] = alti 
        aircraft.fp['wfuel'] = max_wfuel/2
        aircraft.update_fp(aircraft.fp)
        wfuel_list.append(aircraft.findWfuel())
    
    plt.plot(wfuel_list, alt)
    plt.ylabel("Altitude (m)")
    plt.xlabel("Wfuel (kg)")
    plt.title("Altitude-Fuel-Burn")
    plt.grid(True)
    plt.show()
    
    pass
    
def alt_drag_plot(aircraft, max_alt, max_wfuel, max_thrust):
    alt = np.linspace(0, max_alt, 100)
    
    drag_list = list()
    for alti in alt:
        aircraft.fp['alt'] = alti 
        aircraft.update_fp(aircraft.fp)
        drag_list.append(aircraft.drag_at_time(0))
    
    plt.plot(drag_list, alt)
    plt.ylabel("Altitude (m)")
    plt.xlabel("Drag (N)")
    plt.title("Altitude-Drag")
    plt.grid(True)
    plt.show()
    
    pass    
    
def alt_cl_plot(aircraft, max_alt, max_wfuel, max_thrust):
    alt = np.linspace(0, max_alt, 100)
    
    cl_list = list()
    for alti in alt:
        aircraft.fp['alt'] = alti 
        aircraft.update_fp(aircraft.fp)
        cl_list.append(aircraft.lift_at_time(0)/(0.5*aircraft.fp['rho']*aircraft.ap['Sref']*aircraft.fp['Vinf']**2))

    
    plt.plot(cl_list, alt)
    plt.ylabel("Altitude (m)")
    plt.xlabel("$C_L$")
    plt.title("Altitude-$C_L$")
    plt.grid(True)
    plt.show()
    
    pass      
    
    

def optimize_fuel_burn(aircraft, max_alt, max_wfuel, max_thrust):
    # State: (altitude, fuel)
    x0 = (max_alt/2, max_wfuel/2)
    bounds = [(0, max_alt), (0, max_wfuel)]

    def set_fp_from_x(x):
        aircraft.fp['alt'] = x[0] # Alt
        aircraft.fp['wfuel']=x[1] # wfuel
        aircraft.update_fp(aircraft.fp)

    def thrust_obj(x):
        set_fp_from_x(x)
        T = aircraft.averageDrag()
        return max_thrust - T

    def wfuel_obj(x):
        set_fp_from_x(x)
        wfuel = aircraft.findWfuel()
        return max_wfuel - wfuel
    
    thrust_constraint = {'type':'ineq', 'fun': thrust_obj}
    wfuel_constraint  = {'type':'ineq', 'fun': wfuel_obj}
    constraints = (thrust_constraint, wfuel_constraint)

    def fuel_burn_obj(x):
        set_fp_from_x(x)
        wfuel = aircraft.findWfuel()
        return wfuel

    res = minimize(fuel_burn_obj, x0, bounds=bounds,\
        constraints=constraints)
    print(res)
    return res

def altitude_study_data(aircraft, altitudes):
    ''' 
    Collects data for altitude study, including:
    altitude, wfuel, qinf, total drag, induced drag, profile drag, CL, L/D
    '''
    dash = '-' * 40

    data = []
    for a in altitudes:
        aircraft.fp['alt'] = a
        aircraft.update_fp(aircraft.fp)
        wfuel = aircraft.findWfuel()
        qinf = aircraft.fp['qinf']

        lift = aircraft.lift_at_time(0)
        drag = aircraft.drag_at_time(0)
        indDrag = aircraft.inducedDrag(lift=lift)
        profDrag = aircraft.profileDrag()
        CL = lift / (0.5*aircraft.fp['rho']*aircraft.ap['Sref']*aircraft.fp['Vinf']**2)
        LD = lift / drag

        data.append([a, wfuel, qinf, drag, indDrag, profDrag, CL, LD])

    headers = ['Altitude (m)', 'wfuel (kg)', 'qinf (Pa)', 'drag (N)', 'indDrag (N)', 'profDrag (N)', 'CL (1)', 'LD (1)']

    print("ALTITUDE STUDY DATA")
    print(*headers, sep='\t')
    for d in data:
        print(*d, sep='\t')
    return data

# NOTE:    Assuming wfuel is 15800, as said in lecture, to calculate L/D from Drag

# FIX ME - Need to measure airfoil chords from diagram to calculate accurate Re numbers.
# FIXED  - Jakob measured the chords and found thicknesses using the t/c ratios.

s3smaxFP   = None
s3sWing    = Component('wing',    True,  .13*3, 3,  210, span=35.79) # RE should be the cord not the span
s3sFuse    = Component('fuse',    False, 3.76, 38,  415)
s3sNacellL = Component('nacellL', True, 0.69,  3.38, 59.0/2) # awet is assumed to be for both
s3sNacellR = Component('nacellR', True, 0.69,  3.38, 59.0/2) # The nacell length was eyeballed
s3sTails   = Component('tails',   True,  .1*3.25, 3.25, 115) # Assume all the tails are one

s3smax = Aircraft( [s3sWing, s3sFuse, s3sNacellL, s3sNacellR, s3sTails])
print('Dp', round(s3smax.profileDrag()))
print('Di', round(s3smax.inducedDrag()))
print('D', round(s3smax.averageDrag()))
print('L', round(s3smax.averageLift()))
print('L/D', s3smax.averageLD())
print("Guessed wfuel", s3smax.fp['wfuel'])
print('Bregue wfuel', s3smax.bregue_wfuel())
wfuel = s3smax.findWfuel()
print("wfuel updated", wfuel)
print('L/D updated', s3smax.averageLD(wfuel))
print()

print("Fuel burn optimization:")
res = optimize_fuel_burn(s3smax, 15e3, 26e3, 2*128e3)
print("Found optimal. Setting aircraft fp to optimal:")
optimal_alt = res.x[0]
s3smax.fp['alt'] = optimal_alt
s3smax.update_fp(s3smax.fp)
optimal_wfuel = s3smax.findWfuel()
s3smax.fp['wfuel'] = optimal_wfuel
print("Optimal altitude (m), optimal wfuel (kg):", optimal_alt, optimal_wfuel)
print("Optimal fp:")
print(s3smax.fp)
# One Nacell:   L/D 18.22712568809953
# Two Nacells:  L/D 16.73512500441754

#alt_fuel_plot(s3smax, 15e3, 26e3, 2*128e3)
# alt_drag_plot(s3smax, 15e3, 26e3, 2*128e3)
# alt_cl_plot(s3smax, 15e3, 26e3, 2*128e3)

print()
altitudes = list(np.linspace(3e3, 15e3, 5)) + [round(optimal_alt)]
altitude_study_data(s3smax, altitudes)


