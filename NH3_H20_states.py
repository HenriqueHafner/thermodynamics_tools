import matplotlib.pyplot as plt
import numpy
from scipy.optimize import minimize

# plt.ion()
# fig, ax = plt.subplots()
# ax.set_xlim(0.1,0.9)
# ax.set_ylim(-90+273,190+273)
# fig.set_dpi(200)
# fig.set_size_inches(4,6)

print_output = True

def printf(data):
    if print_output:
        print(data)
    else:
        return data

def ln(number):
    number = float(number)
    number = max(number,0.0000001)
    number = numpy.log(number)
    printf(number)
    return number

def exp(number):
    number = float(number)
    number = numpy.exp(number)
    printf(number)
    return number

def NH3_molar_to_massic_frac(NH3_molar_frac):
  NH3_molar_wheight = 17.031 # kg/kmol
  H2O_molar_wheight = 18.015 # kg/kmol
  H2O_molar_frac = 1-NH3_molar_frac
  NH3_mass = NH3_molar_frac*NH3_molar_wheight
  H2O_mass = H2O_molar_frac*H2O_molar_wheight
  NH3_massic_frac = NH3_mass/(NH3_mass+H2O_mass)
  printf(NH3_massic_frac)
  return NH3_massic_frac

def NH3_massic_to_molar_frac(NH3_mass_frac):
  NH3_molar_wheight = 17.031 # kg/kmol
  H2O_molar_wheight = 18.015 # kg/kmol
  H2O_mass_frac = 1-NH3_mass_frac
  NH3_mols = NH3_mass_frac/NH3_molar_wheight
  H2O_mols = H2O_mass_frac/H2O_molar_wheight
  NH3_molar_frac = NH3_mols/(NH3_mols+H2O_mols)
  printf(NH3_molar_frac)
  return NH3_molar_frac

def T_px(p,x):
    p0 = 2   # MPa
    t0 = 100 # K
    table_1 = [
    [0,  0, +0.322302e+1],
    [0,  1, -0.384206e+0],
    [0,  2, +0.460965e-1],
    [0,  3, -0.378945e-2],
    [0,  4, +0.135610e-3],
    [1,  0, +0,487755e+0],
    [1,  1, -0.120108e+0],
    [1,  2, +0.106154e-1],
    [2,  3, -0.533589e-3],
    [4,  0, +0.785041e+1],
    [5,  0, -0.115941e+2],
    [5,  1, +0.523150e-1],
    [6,  0, +0.489596e+1],
    [13, 1, +0.421059e-1]]
    ser_i = 0
    for i in table_1:
        mi = i[0]
        ni = i[1]
        ai = i[2]
        value_i = ai*((1-x)**mi)*(ln(p0/p)**ni)
        ser_i += value_i 
    T = t0*ser_i
    printf(T)
    return T

def T_py(p,y):
    p0 =  2  # MPa
    t0 =  100# K
    table_7 = [
    [0, 0, +0.324004e+1],
    [0, 1, -0.395920e+0],
    [0, 2, +0.435624e-1],
    [0, 3, -0.218943e-2],
    [1, 0, -0.143526e+1],
    [1, 1, +0.105256e+1],
    [1, 2, -0.719281e-1],
    [2, 0, +0.122362e+2],
    [2, 1, -0.224368e+1],
    [3, 0, -0.201780e+2],
    [3, 1, +0.110834e+1],
    [4, 0, +0.145399e+2],
    [4, 2, +0.644312e+0],
    [5, 0, -0.221246e+1],
    [5, 2, -0.756266e+0],
    [6, 0, -0.135529e+1],
    [7, 2, +0.183541e+0]]
    ser_i = 0
    for i in table_7:
        mi = i[0]
        ni = i[1]
        ai = i[2]
        value_i = ai*((1-y)**(mi/4))*(ln(p0/p)**ni)
        ser_i += value_i 
    T = t0*ser_i
    printf(T)
    return T

def y_px(p,x):
    p0 =  2  # MPa
    table_3 = [
    [0, 0, +1.98022017e+1],
    [0, 1, -1.18092669e+1],
    [0, 6, +2.77479980e+1],
    [0, 7, -2.88634277e+1],
    [1, 0, -5.91616608e+1],
    [2, 1, +5.78091305e+2],
    [2, 2, -6.21736743e+0],
    [3, 2, -3.42198402e+3],
    [4, 3, +1.19403127e+4],
    [5, 4, -2.45413777e+4],
    [6, 5, +2.91591865e+4],
    [7, 6, -1.84782290e+4],
    [7, 7, +2.34819434e+1],
    [8, 7, +4.80310617e+3]]
    ser_i = 0
    for i in table_3:
        mi = i[0]
        ni = i[1]
        ai = i[2]
        value_i = ai*((p/p0)**mi)*(x**(ni/3))
        ser_i += value_i 

    y = 1 - exp(ln(1-x)*ser_i)
    printf(y)
    return y

def hl(T,x):
    h0 = 100.00*0.6842 # kJ/kg
    T0 = 273.16 # K
    table_4 = [
    [0,  1, -0.761080e+1],
    [0,  4, +0.256905e+2],
    [0,  8, -0.247092e+3],
    [0,  9, +0.325952e+3],
    [0, 12, -0.158854e+3],
    [0, 14, +0.619084e+2],
    [1,  0, +0.114314e+2],
    [1,  1, +0.118157e+1],
    [2,  1, +0.284179e+1],
    [3,  3, +0.741609e+1],
    [5,  3, +0.891844e+3],
    [5,  4, -0.161309e+4],
    [5,  5, +0.622106e+3],
    [6,  2, -0.207588e+3],
    [6,  4, -0.687393e+1],
    [8,  0, +0.350716e+1]]
    ser_i = 0
    for i in table_4:
        mi = i[0]
        ni = i[1]
        ai = i[2]
        value_i = ai*((T/T0-1)**mi)*(x**ni)
        ser_i += value_i 
    hl = h0*ser_i
    printf(hl)
    return hl

def hg(T,y):
    h0 = 1000.00 # kJ/kg
    T0 =  324.00 # K
    table_5 = [
    [0, 0,  +0.128827e+1],
    [1, 0,  +0.125247e+0],
    [2, 0,  -0.208748e+1],
    [3, 0,  +0.217696e+1],
    [0, 2,  +0.235687e+1],
    [1, 2,  -0.886987e+1],
    [2, 2,  +0.102635e+2],
    [3, 2,  -0.237440e+1],
    [0, 3,  -0.670515e+1],
    [1, 3,  +0.164508e+2],
    [2, 3,  -0.936849e+1],
    [0, 4,  +0.842254e+1],
    [1, 4,  -0.858807e+1],
    [0, 5,  -0.277049e+1],
    [4, 6,  -0.961248e+0],
    [2, 7,  +0.988009e+0],
    [1, 10, +0.308482e+0]]
    ser_i = 0
    for i in table_5:
        mi = i[0]
        ni = i[1]
        ai = i[2]
        value_i = ai*((1-T/T0)**mi)*((1-y)**(ni/4))
        ser_i += value_i 
    hg = h0*ser_i
    printf(hg)
    return hg

def p_Tx(T,x,p0=0.1):
    def T_residue(p,T,x):
        return (T-T_px(p,x))**2
    pf = minimize(T_residue,p0,tol=0.001,args=(T,x))
    print()
    pf = pf.x[0]
    pf = round(pf,5)
    printf(p_Tx)
    return pf

def draw_diagram():
    plt.figure(figsize=(36,36*1.62))
    colors_table = ['b','g','r','c','m','y','k','purple','b','g','r','c','b']
    massic_fraction_horizontal_values = numpy.arange(0,1.005,step=0.005)

    molar_fraction_horizontal_values =numpy.zeros_like(massic_fraction_horizontal_values)
    for i in range(len(molar_fraction_horizontal_values)):
        molar_fraction_horizontal_values[i] = NH3_massic_to_molar_frac(massic_fraction_horizontal_values[i])

    isobars = [2.000,1.000,0.750,0.500,0.300,0.200,0.100,0.075,0.050,0.030,0.020,0.010]
    color_index = 0
    for curr_p_isobar in isobars:
        p = curr_p_isobar

        hl_isobaric_vertical_values =numpy.zeros_like(molar_fraction_horizontal_values)
        for i in range(len(hl_isobaric_vertical_values)):
            Ti = T_px(p,molar_fraction_horizontal_values[i])
            hl_isobaric_vertical_values[i] = hl(Ti,molar_fraction_horizontal_values[i])

        hg_isobaric_vertical_values = numpy.zeros_like(molar_fraction_horizontal_values)
        for i in range(len(hg_isobaric_vertical_values)):
            xi = molar_fraction_horizontal_values[i]
            Ti = T_py(p, xi)
            hg_isobaric_vertical_values[i] = hg(Ti, xi)
        
        hgy_isobaric_vertical_values = numpy.zeros_like(molar_fraction_horizontal_values)
        for i in range(len(hgy_isobaric_vertical_values)):
            xi = molar_fraction_horizontal_values[i]
            yi = y_px(p, xi)
            Ti = T_py(p, yi)
            hgy_isobaric_vertical_values[i] = hg(Ti, yi)

        curr_color = colors_table[color_index]
        color_index += 1
        label_iteration = 'at '+str(p*1000)+' kPa'
        plt.plot(massic_fraction_horizontal_values, hl_isobaric_vertical_values, color=curr_color, label=label_iteration)
        plt.plot(massic_fraction_horizontal_values, hg_isobaric_vertical_values, color=curr_color)
        plt.plot(massic_fraction_horizontal_values, hgy_isobaric_vertical_values, color=curr_color)

    plt.legend(fontsize=50)
    plt.grid()
    plt.xticks(numpy.arange(0,1.02,step=0.02))
    plt.yticks(numpy.arange(-400,2840,step=40))
    plt.xlabel("Mass Fraction")
    plt.ylabel("Enthalpy [kJ/kg]")
    plt.title('h [kJ/kg] (liquid and vapor) X mass fraction of NH3.')
    #plt.show()
    return True
