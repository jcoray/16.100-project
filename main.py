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

def alt_fuel_plot_again(aircraft, baseline, max_alt):
    alt = np.linspace(0, max_alt, 100)

    wfuel_list = list()
    for alti in alt:
        res = optimize_fuel_burn_again(aircraft, baseline, alti)
        wfuel_list.append(res.fun)

    plt.plot(wfuel_list, alt)
    plt.ylabel("Altitude (m)")
    plt.xlabel("Wfuel (kg)")
    plt.title("Altitude-Fuel-Burn (varying Sref)")
    plt.grid(True)
    plt.show()

def all_plots_sref_part(aircraft, baseline, max_alt):
    alt = np.linspace(0, max_alt, 100)

    sref_list = list()
    wfuel_list = list()
    mdd_m_perp_diff_list = list()
    wmax_wfuel_diff_list = list()
    for alti in alt:
        res = optimize_fuel_burn_again_plotter(aircraft, baseline, alti)
        wfuel = res.fun
        sref_list.append(res.x[0])
        wfuel_list.append(wfuel)
        mdd_m_perp_diff_list.append(aircraft.mdd() - aircraft.m_perp())
        wmax = baseline.ap['rho_fuel'] * aircraft.components['wing'].vtank()
        wmax_wfuel_diff_list.append(wmax - wfuel)


    plt.plot(sref_list, alt)
    plt.ylabel("Altitude (m)")
    plt.xlabel("Sref (m^2)")
    plt.title("Altitude-Sref (varying Sref)")
    plt.grid(True)
    plt.show()

    plt.plot(wfuel_list, alt)
    plt.ylabel("Altitude (m)")
    plt.xlabel("Wfuel (kg)")
    plt.title("Altitude-Fuel-Burn (varying Sref)")
    plt.grid(True)
    plt.show()

    plt.plot(mdd_m_perp_diff_list, alt)
    plt.ylabel("Altitude (m)")
    plt.xlabel("Mdd - M_perp")
    plt.title("Mach Number Constraint (Positive means within constraint)")
    plt.grid(True)
    plt.show()

    plt.plot(wmax_wfuel_diff_list, alt)
    plt.ylabel("Altitude (m)")
    plt.xlabel("Wmax - Wfuel")
    plt.title("Fuel Capacity Constraint (Positive means within constraint)")
    plt.grid(True)
    plt.show()
    
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

def optimize_fuel_burn_again_plotter(aircraft, baseline, altitude):
    # State: (Sref)
    x0 = (baseline.fp['Sref'],)
    bounds = [(0.1,1e3)] # TODO: Figure out how to plug-in partially defined bounds

    # TODO: Create function to set state in aircraft
    def set_state(x):
        # Set state in aircraft
        # TODO: Figure out how to input wfuel here and in subsequent funs.
        # fp: altitude, Sref, Aeng
        aircraft.fp['alt'] = altitude
        aircraft.fp['Sref']= x[0]
        # aircraft.fp['Aeng']= x[2]
        aircraft.components = update_components(baseline, aircraft.fp)
        aircraft.update_fp(aircraft.fp)
        aircraft.findAeng(baseline)

    # TODO: Create constraint objective functions:
    # Drag divergence, Fuel capacity

    def mdd_obj(x):
        set_state(x)
        mdd = aircraft.mdd()
        m_perp = aircraft.m_perp()
        return mdd - m_perp

    def fuel_cap_obj(x):
        set_state(x)
        wfuel_max = baseline.ap['rho_fuel'] * aircraft.components['wing'].vtank()
        wfuel = aircraft.findWfuel()
        return wfuel_max - wfuel

    mdd_constraint = {'type':'ineq', 'fun': mdd_obj}
    fuel_cap_constraint  = {'type':'ineq', 'fun': fuel_cap_obj}
    constraints = (mdd_constraint, fuel_cap_constraint)

    def fuel_burn_obj(x):
        set_state(x)
        wfuel = aircraft.findWfuel()
        return wfuel

    res = minimize(fuel_burn_obj, x0, bounds=bounds,\
        constraints=constraints)
    #print(res)
    return res
    
def optimize_fuel_burn_sref_alt(aircraft, baseline, maxalt, initalt=None):
    if initalt == None: initalt = baseline.fp['alt']
    # State: (alt, Sref)
    x0 = (initalt,baseline.fp['Sref'],)
    bounds = [(10000,maxalt),(100,200)] # TODO: Figure out how to plug-in partially defined bounds

    # TODO: Create function to set state in aircraft
    def set_state(x):
        # Set state in aircraft
        # TODO: Figure out how to input wfuel here and in subsequent funs.
        # fp: altitude, Sref, Aeng
        aircraft.fp['alt'] = x[0]
        aircraft.fp['Sref']= x[1]
        #~ aircraft.fp['Aeng']= x[2]
        aircraft.components = update_components(baseline, aircraft.fp)
        aircraft.update_fp(aircraft.fp)
        aircraft.findAeng(baseline)

    # TODO: Create constraint objective functions:
    # Drag divergence, Fuel capacity

    def mdd_obj(x):
        set_state(x)
        mdd = aircraft.mdd()
        m_perp = aircraft.m_perp()
        return mdd - m_perp

    def fuel_cap_obj(x):
        set_state(x)
        wfuel_max = baseline.ap['rho_fuel'] * aircraft.components['wing'].vtank()
        wfuel = aircraft.findWfuel()
        return wfuel_max - wfuel

    mdd_constraint = {'type':'ineq', 'fun': mdd_obj}
    fuel_cap_constraint  = {'type':'ineq', 'fun': fuel_cap_obj}
    constraints = (mdd_constraint, fuel_cap_constraint)

    def fuel_burn_obj(x):
        set_state(x)
        wfuel = aircraft.findWfuel()
        return wfuel

    res = minimize(fuel_burn_obj, x0, bounds=bounds,\
        constraints=constraints)
    #print(res)
    return res

        

#~ def altitude_study_data(aircraft, altitudes):
    #~ ''' 
    #~ Collects data for altitude study, including:
    #~ altitude, wfuel, qinf, total drag, induced drag, profile drag, CL, L/D
    #~ '''
    #~ dash = '-' * 40
#~ 
    #~ data1 = []
    #~ data2 = []
    #~ for a in altitudes:
        #~ aircraft.fp['alt'] = a
        #~ aircraft.update_fp(aircraft.fp)
        #~ wfuel = aircraft.findWfuel()
        #~ qinf = aircraft.fp['qinf']
#~ 
        #~ lift = aircraft.lift_at_time(0)
        #~ drag = aircraft.drag_at_time(0)
        #~ indDrag = aircraft.inducedDrag(lift=lift)
        #~ profDrag = aircraft.profileDrag()
        #~ CL = lift / (0.5*aircraft.fp['rho']*aircraft.fp['Sref']*aircraft.fp['Vinf']**2)
        #~ LD = lift / drag
        #~ CLperp = CL/math.cos(aircraft.components['wing'].Lambda())**2
#~ 
        #~ data1.append([a, res.x[0], aircraft.components['wing'].span(), aircraft.components['engine'].Aeng(), \
            #~ aircraft.components['wing'].mass(), aircraft.components['tail'].mass(), aircraft.woe(), \
            #~ aircraft.components['wing'].vtank(), res.fun])
        #~ data2.append([a, qinf, drag, Dind, Dprof, LD, CL, CLperp, aircraft.mdd()])
#~ 
    #~ headers1 = ['Altitude',   'S_{ref}' , 'b' , 'A_{eng}' , 'W_{wing}', 'W_{tail}', 'W_{OE}', 'W_{fuel}', 'W_{fuel}^{max}']
    #~ headers2 = ['Altitude', 'q', 'D', 'Dind', 'Dprof', 'LD', 'CL', 'CLperp', 'Mdd']
   #~ 
    #~ print("ALTITUDE STUDY DATA")
    #~ print(*headers1, sep='\t')
    #~ for d in data1:
        #~ print(*d, sep=' & \t')
    #~ 
    #~ print(*headers2, sep='\t')
    #~ for d in data2:
        #~ print(*d, sep=' & \t')
    #~ 
    #~ return data1, data2
    #~ 
    
def altitude_study_data_part4(aircraft, baseline, altitudes):
    ''' 
    Collects data for altitude study, including:
    '''
    dash = '-' * 40

  
    data1 = []
    data2 = []
    for a in altitudes:
        aircraft.fp['alt'] = a
        aircraft.update_fp(aircraft.fp)
        res = optimize_fuel_burn_again_plotter(aircraft, baseline, a)
        wfuel = aircraft.findWfuel()
        qinf = aircraft.fp['qinf']

        lift = aircraft.lift_at_time(0)
        drag = aircraft.drag_at_time(0)
        indDrag = aircraft.inducedDrag(lift=lift)
        profDrag = aircraft.profileDrag()
        CL = lift / (0.5*aircraft.fp['rho']*aircraft.fp['Sref']*aircraft.fp['Vinf']**2)
        LD = lift / drag
        CLperp = CL/math.cos(aircraft.components['wing'].Lambda())**2

        data1.append([a, \
            round(res.x[0],2), \
            round(aircraft.components['wing'].span(),2), \
            round(aircraft.components['engine'].Aeng(),2), \
            int(aircraft.components['wing'].mass(aircraft.components['fuse'],aircraft.ap['wpay'])), \
            int(aircraft.components['tail'].mass()), \
            int(aircraft.woe()), \
            int(res.fun), \
            int(aircraft.components['wing'].vtank()*baseline.ap['rho_fuel']), \
            '''\\\ \hline'''])
        data2.append([a, round(qinf/1000,2), round(drag/1000,1), round(indDrag/1000,1), round(profDrag/1000,1), round(LD,2), round(CL,3), round(CLperp,3), round(aircraft.mdd(),3), \
            '''\\\ \hline'''])

    headers1 = ['Altitude',   'S_{ref}' , 'b' , 'A_{eng}' , 'W_{wing}', 'W_{tail}', 'W_{OE}', 'W_{fuel}', 'W_{fuel}^{max}']
    headers2 = ['Altitude', 'q', 'D', 'Dind', 'Dprof', 'LD', 'CL', 'CLperp', 'Mdd']
   
    print("ALTITUDE STUDY DATA")
    print(*headers1, sep='\t')
    for d in data1:
        print(*d , sep=' & \t')
    
    print(*headers2, sep='\t')
    for d in data2:
        print(*d, sep=' & \t')
    
    return data1, data2




def part1():
    
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

    Sref = 127
    b = 35.79
    AR = b**2 / Sref
    c = math.sqrt(Sref / AR)
    print("c = ", c)
    Aeng = math.pi*2.44**2 / 4

    ap = {'wfuelland':2300, 'wpay':20e3, 'R':6500e3, 'g':9.81, 'rho_fuel': 800}
    fp = {'wfuel':15800, 'Minf':0.78, 'alt': 10000, \
                'TSFC': 1.42e-5, 'Sref':Sref, 'AR':AR, 'Aeng': Aeng}

    s3sWing    = Wing(     .13*c,    c, span=b, K_wing=0.71, Lambda=25*math.pi/180.0, rho_box=2700, omega_box=2.1e8, taper=0.3) 
    s3sFuse    = Fuselage(  3.76,   38,  415, 19200 )
    s3sNacell  = Engine(    0.69, 3.38, 59.0, 11000, Aeng) # awet is assumed to be for both # The nacell length was eyeballed
    s3sTails   = Tail(   .1*3.25, 3.25,  115, .2 * s3sWing.mass(s3sFuse, ap['wpay'])) # Assume all the tails are one


    s3smax = Aircraft( [s3sWing, s3sFuse, s3sNacell, s3sTails], ap=copy.deepcopy(ap), fp=copy.deepcopy(fp))

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
    
    return 0

def part3():
    
    Sref = 127
    b = 35.79
    AR = b**2 / Sref
    c = math.sqrt(Sref / AR)
    print("c = ", c)
    Aeng = math.pi*2.44**2 / 4

    ap = {'wfuelland':2300, 'wpay':20e3, 'R':6500e3, 'g':9.81, 'rho_fuel': 800}
    fp = {'wfuel':15800, 'Minf':0.78, 'alt': 10000, \
                'TSFC': 1.42e-5, 'Sref':Sref, 'AR':AR, 'Aeng': Aeng}

    s3sWing    = Wing(     .13*c,    c, span=b, K_wing=0.71, Lambda=25*math.pi/180.0, rho_box=2700, omega_box=2.1e8, taper=0.3) 
    s3sFuse    = Fuselage(  3.76,   38,  415, 19200 )
    s3sNacell  = Engine(    0.69, 3.38, 59.0, 11000, Aeng) # awet is assumed to be for both # The nacell length was eyeballed
    s3sTails   = Tail(   .1*3.25, 3.25,  115, .2 * s3sWing.mass(s3sFuse, ap['wpay'])) # Assume all the tails are one


    s3smax = Aircraft( [s3sWing, s3sFuse, s3sNacell, s3sTails], ap=copy.deepcopy(ap), fp=copy.deepcopy(fp))
    baseline = Aircraft( [s3sWing, s3sFuse, s3sNacell, s3sTails], ap=copy.deepcopy(ap), fp=copy.deepcopy(fp))

    # print("Fuel burn optimization:")
    # res = optimize_fuel_burn(s3smax, 15e3, 26e3, 2*128e3)
    # print("Found optimal. Setting aircraft fp to optimal:")
    # optimal_alt = res.x[0]
    # s3smax.fp['alt'] = optimal_alt
    # s3smax.update_fp(s3smax.fp)
    # optimal_wfuel = s3smax.findWfuel()
    # s3smax.fp['wfuel'] = optimal_wfuel
    # print("Optimal altitude (m), optimal wfuel (kg):", optimal_alt, optimal_wfuel)
    # print("Optimal fp:")
    # print(s3smax.fp)
    # One Nacell:   L/D 18.22712568809953
    # Two Nacells:  L/D 16.73512500441754

    # alt_fuel_plot(s3smax, 15e3, 26e3, 2*128e3)
    # alt_drag_plot(s3smax, 15e3, 26e3, 2*128e3)
    # alt_cl_plot(s3smax, 15e3, 26e3, 2*128e3)
    

    print()
    altitudes = list(np.linspace(3e3, 15e3, 5)) #+ [round(optimal_alt)]
    altitude_study_data(s3smax, altitudes)
    
def part4():
    
    Sref = 127
    b = 35.79
    AR = b**2 / Sref
    print('AR', AR)
    print('Lambda', 25*math.pi/180.0)

    c = math.sqrt(Sref / AR)
    #~ print("c = ", c)
    Aeng = math.pi*2.44**2 / 4

    ap = {'wfuelland':2300, 'wpay':20e3, 'R':6500e3, 'g':9.81, 'rho_fuel': 800}
    fp = {'wfuel':15800, 'Minf':0.78, 'alt': 10000, \
                'TSFC': 1.42e-5, 'Sref':Sref, 'AR':AR, 'Aeng': Aeng}

    s3sWing    = Wing(     .13*c,    c, span=b, K_wing=0.71, Lambda=25*math.pi/180.0, rho_box=2700, omega_box=2.1e8, taper=0.3) 
    s3sFuse    = Fuselage(  3.76,   38,  415, 19200 )
    s3sNacell  = Engine(    0.69, 3.38, 59.0, 11000, Aeng) # awet is assumed to be for both # The nacell length was eyeballed
    s3sTails   = Tail(   .1*3.25, 3.25,  115, .2 * s3sWing.mass(s3sFuse, ap['wpay'])) # Assume all the tails are one


    s3smax = Aircraft( [s3sWing, s3sFuse, s3sNacell, s3sTails], ap=copy.deepcopy(ap), fp=copy.deepcopy(fp))
    baseline = Aircraft( [s3sWing, s3sFuse, s3sNacell, s3sTails], ap=copy.deepcopy(ap), fp=copy.deepcopy(fp))

    maxalt = 15e3
    #~ print("Test alt:", test_alt)
    #~ res = optimize_fuel_burn_sref_alt(s3smax, baseline, maxalt)
    #~ print('res',res)
    #~ print(s3smax.fp)
    alts = [a*1000+10000 for a in range(6)]


    def initalConditionSensitivityStudy():
        for a in alts:
            res = optimize_fuel_burn_sref_alt(s3smax, baseline, maxalt)
            print('Alt', a, '   X',res.x, '   Wfuel', res.fun)
            #Alt 10000    X [ 13422.52597943    178.42461114]    Wfuel 13925.468604663049
            #Alt 11000    X [ 13289.9721146     177.50699479]    Wfuel 13924.927986742347
            #Alt 12000    X [ 13629.80167927    184.62443866]    Wfuel 13926.53501829638
            #Alt 13000    X [ 13677.19563179    180.01826931]    Wfuel 13930.24982839971
            #Alt 14000    X [ 13660.4809263     180.20236997]    Wfuel 13929.423210925443
    
    #~ initalConditionSensitivityStudy()
    
    #~ altitude_study_data_part4(s3smax, baseline, alts)
    #~ print("Plotting Sref Altitude-Fuel-Burn (this will take a while...)")
    all_plots_sref_part(s3smax, baseline, maxalt)
    pass

def main():
    part4()
    return 0
    
if __name__ == '__main__':
     main()
