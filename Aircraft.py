# Aircraft.py
# Created by Team 2 - Jakob, Jose, and Austin on 11/27/18

import math
import warnings

class Aircraft(object):
    '''An aircraft container'''
    def __init__(self, components, fp=None):
        ''' 
        fp (dict): SI units. Includes: 
        woe (operating empty MASS), wfuelland (fuel to land MASS), 
        wpay (Payload MASS), wfuel (Fuel MASS for cruise),
        R (range), Minf (Mach), ainf (speed o' sound), 
        g (like in F = m*g), rho (density inf), nu (kinematic visc)
        TSFC (Thrust Spec Fuel Consumption; Here because the notes gave us a value to use)
        '''
        
        # By default, assume a 737 max
        if fp == None:
            self.fp = {'woe':45e3, 'wfuelland':2300, 'wpay':20e3, 'wfuel':15800, \
            'R':6500e3, 'Minf':0.78, 'ainf':300, 'g':9.8, 'rho':0.41351, \
            'nu':0.000035251, 'TSFC': 1.42e-5}
        else: 
            self.fp = fp
            
        if 'Vinf' not in self.fp:
            self.fp['Vinf'] = self.fp['Minf']*self.fp['ainf']
        if 'qinf' not in self.fp:
            self.fp['qinf'] = 1/2*self.fp['rho']*self.fp['Vinf']**2
        if 'T' not in self.fp:
            self.fp['T'] = self.fp['R'] / self.fp['Vinf'] # Total cruise time
        if 'wfinal' not in self.fp:
            self.fp['wfinal'] = self.fp['woe'] + self.fp['wfuelland'] + self.fp['wpay']
        if 'winit'  not in self.fp:
            self.fp['winit'] = self.fp['wfinal'] + self.fp['wfuel']
            
        self.components = {c.name: c for c in components}
        return
    
    def lift_at_time(self, t):
        ''' Find lift at time t after start of cruise flight '''
        if (t < 0) or (t > self.fp['T']):
            raise ValueError("t must be positive and less than total cruise time.")
        lift = self.fp['g'] * (self.fp['winit']  - (self.fp['wfuel']/self.fp['T'])*t)
        return lift

    def averageLift(self):
        ''' Calculates average lift using initial and final weights '''
        Lavg = self.fp['g']/2*(2*self.fp['wfinal'] + self.fp['wfuel']) 
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
        
    def drag_at_time(self, t):
        ''' Find drag at time t after start of cruise flight '''
        if (t < 0) or (t > self.fp['T']):
            raise ValueError("t must be positive and less than total cruise time.")
        lift = self.lift_at_time(t)
        drag = self.profileDrag() + self.waveDrag() + self.inducedDrag(lift=lift)
        return drag

    def averageDrag(self, nsamples=50):
        # Values of drag average when using averaged lift vs iterating through time
        # are not that different.
        #return self.profileDrag() + self.inducedDrag() + self.waveDrag()
        drag_samples = []
        for i in range(nsamples):
            t = self.fp['T'] * (i / (nsamples-1))
            drag = self.drag_at_time(t)
            drag_samples.append(drag)
        mean_drag = sum(drag_samples) / nsamples
        return mean_drag
    
    def averageLD(self, nsamples=50):
        # Similar to averageDrag, iterating through L/D values
        # vs using L and D's averages hadly changes anything
        #return self.averageLift()/self.averageDrag()
        ld_samples = []
        for i in range(nsamples):
            t = self.fp['T'] * (i / (nsamples-1))
            ld = self.lift_at_time(t) / self.drag_at_time(t)
            ld_samples.append(ld)
        mean_ld = sum(ld_samples) / nsamples
        return mean_ld
        
    # TSFC is assumed using the same logic with which we assume wfuel
    # In other words, the notes told us to assume TSFC = 1.42e-5
    # NOTE: winit and wfinal are now in the self.fp dictionary
    def bregue_wfuel(self):
        ''' Estimates required wfuel from Bregue Range Equation model. '''
        return self.fp['wfinal'] * (math.exp(self.fp['TSFC']*self.fp['g']*self.fp['R']
                                    / (self.fp['Vinf']*self.averageLD())   ) - 1)
        

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
        
# NOTE:    Assuming wfuel is 15800, as said in lecture, to calculate L/D from Drag

# FIX ME - Need to measure airfoil chords from diagram to calculate accurate Re numbers.
# FIXED  - Jakob measured the chords and found thicknesses using the t/c ratios.

s3smaxFP   = None
s3sWing    = Component('wing',    True,  .13*3, 3,  210, span=35.79) # RE should be the cord not the span
s3sFuse    = Component('fuse',    False, 3.76, 38,  415)
s3sNacellL = Component('nacellL', False, 2.44,  4.0, 59.0/2) # awet is assumed to be for both
s3sNacellR = Component('nacellR', False, 2.44,  4.0, 59.0/2) # The nacell length was eyeballed
s3sTails   = Component('tails',   True,  .1*3.25, 3.25, 115) # Assume all the tails are one

s3smax = Aircraft( [s3sWing, s3sFuse, s3sNacellL, s3sNacellR, s3sTails])
print('Dp', round(s3smax.profileDrag()))
print('Di', round(s3smax.inducedDrag()))
print('D', round(s3smax.averageDrag()))
print('L', round(s3smax.averageLift()))
print('L/D', s3smax.averageLD())
print("Guessed wfuel", s3smax.fp['wfuel'])
print('Bregue wfuel', s3smax.bregue_wfuel())

# One Nacell:   L/D 18.22712568809953
# Two Nacells:  L/D 16.73512500441754



