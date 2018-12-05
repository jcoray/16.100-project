from Aircraft import *
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
import math
from skaero.atmosphere import coesa
import copy

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
        cl_list.append(aircraft.lift_at_time(0)/(0.5*aircraft.fp['rho']*aircraft.fp['Sref']*aircraft.fp['Vinf']**2))

    
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

def optimize_fuel_burn_again(aircraft, max_alt):
    # State: (altitude, Sref, Aeng)
    x0 = (max_alt/2, 127, math.pi*(2.44**2)/4)
    bounds = [(0, max_alt)] # TODO: Figure out how to plug-in partially defined bounds

    # TODO: Create function to set state in aircraft
    def set_state(x):
        # Set state in aircraft
        # TODO: Figure out how to input wfuel here and in subsequent funs.
        pass

    # TODO: Create constraint objective functions:
    # Drag divergence, Fuel capacity

    def mdd_obj(x):
        set_state(x)
        mdd = aircraft.mdd()
        m_perp = aircraft.m_perp()
        return mdd - m_perp

    def fuel_cap_obj(x):
        set_state(x)
        

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
        CL = lift / (0.5*aircraft.fp['rho']*aircraft.fp['Sref']*aircraft.fp['Vinf']**2)
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

#~ s3smaxFP   = None
#~ s3sWing    = Component('wing',    True,  .13*3, 3,  210, span=35.79) # RE should be the cord not the span
#~ s3sFuse    = Component('fuse',    False, 3.76, 38,  415)
#~ s3sNacellL = Component('nacellL', True, 0.69,  3.38, 59.0/2) # awet is assumed to be for both
#~ s3sNacellR = Component('nacellR', True, 0.69,  3.38, 59.0/2) # The nacell length was eyeballed
#~ s3sTails   = Component('tails',   True,  .1*3.25, 3.25, 115) # Assume all the tails are one

# base_ap = {'wfuelland':2300, 'wpay':20e3, 'R':6500e3, 'g':9.81}
# base_fp = {'wfuel':15800, 'Minf':0.78, 'alt': 10000, \
#             'TSFC': 1.42e-5, 'Sref':127}

ap = {'wfuelland':2300, 'wpay':20e3, 'R':6500e3, 'g':9.81}
fp = {'wfuel':15800, 'Minf':0.78, 'alt': 10000, \
            'TSFC': 1.42e-5, 'Sref':127}

s3sWing    = Wing(     .13*3,    3, span=35.79, K_wing=0.71, Lambda=25*math.pi/180.0, rho_box=2700, omega_box=2.1e8, taper=0.3) 
s3sFuse    = Fuselage(  3.76,   38,  415, 19200 )
s3sNacell  = Engine(    0.69, 3.38, 59.0, 11000, math.pi*2.44**2 / 4) # awet is assumed to be for both # The nacell length was eyeballed
s3sTails   = Tail(   .1*3.25, 3.25,  115, .2 * s3sWing.mass(s3sFuse, ap['wpay'])) # Assume all the tails are one


s3smax = Aircraft( [s3sWing, s3sFuse, s3sNacell, s3sTails], ap=copy.deepcopy(ap), fp=copy.deepcopy(fp))
baseline = Aircraft( [s3sWing, s3sFuse, s3sNacell, s3sTails], ap=copy.deepcopy(ap), fp=copy.deepcopy(fp))

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

# alt_fuel_plot(s3smax, 15e3, 26e3, 2*128e3)
# alt_drag_plot(s3smax, 15e3, 26e3, 2*128e3)
# alt_cl_plot(s3smax, 15e3, 26e3, 2*128e3)

print()
altitudes = list(np.linspace(3e3, 15e3, 5)) + [round(optimal_alt)]
altitude_study_data(s3smax, altitudes)


