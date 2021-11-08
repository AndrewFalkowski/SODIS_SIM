import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
import numba
import time
from scipy.integrate import odeint



# a sample differential equation dy/dx = (x-y)/2

def dydx(x,y):
    return ((x-y)/2)

# find the value of y for a given x using step size h
# and an initial value y0 at x0

def rungeKutta(x0, y0, x, h):
    #count num iteratings using step size or step height h
    n = int(((x - x0)/h))
    # iterate for number of iterations
    y = y0
    for i in range(1, n + 1):
        # apply runge kutta formulas to find the next value of y
        k1 = h * dydx(x0, y)
        k2 = h * dydx(x0 + 0.5 * h, y + 0.5 * k1)
        k3 = h * dydx(x0 + 0.5 * h, y + 0.5 * k2)
        k4 = h * dydx(x0 + h, y + k3)

        # update the next value of y
        y = y + (1.0 / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

        # update the next value of x
        x0 = x0 + h

    return y


# driver method
x0 = 0
y = 1
x = 2
h = 0.2
print('The value of y at x is:', rungeKutta(x0, y, x, h))

#%%

def box_dim(A_c, h, prct_f):
    # all dimensions in meters
    box_vol = A_c * h
    vol_f = box_vol * prct_f # L
    m_a = box_vol * (1-prct_f) * 1.225
    m_f = vol_f * 997 # kg
    print('Contained Water: ', m_f, 'Liters')
    A_s = 4 * h * np.sqrt(A_c)
    return m_f, m_a, A_s

m_f, m_a, A_s = box_dim(0.25, 0.15, 0.9)


def boxODE(x, t, m_f, m_a, A_s):

    # constants
    A_c = 0.25 # square meters
    A_s = A_s
    A_f = A_c # square meters
    T_amb = 298 # kelvin
    T_sky = T_amb - 6 # kelvin
    alpha_g = 0.02 # %
    alpha_p = 0.98
    t_g = 0.85 # %
    t_f = 0.85 # %
    Irr = 0.0426*(t+27000) + 1.38E-6*(t+27000)**2 - 7.94E-11*(t+27000)**3 + 7.3E-16*(t+27000)**4
    # Irr = 600
    x_b = 0.065 # insulation thickness meters
    x_s = 0.065 # insulation thickness meters

    k_i = 1.0 # thermal conductivity of side materials, foamed glass # W/mK
    h_rad_g2_g1 = 10
    h_cov_g2_g1 = 5
    h_rad_g1_sky = 10
    h_rad_g1_amb = 10
    h_rad_p_g2 = 25
    h_cov_a_g2 = 10
    h_cov_f_a = 10
    h_cov_p_f = 20
    h_cov_g1_amb = 100

    M_f = m_f * 4.187
    M_g1 = 1150 * (A_c * 0.001) * 1.67 # assuming acrylic
    M_g2 = M_g1
    M_p = 8960 * (A_c * 0.065) * 0.5 # assuming coper
    M_a = 0.718 * m_a

    # assign each ODE to a vector element
    T_g1 = x[0]
    T_g2 = x[1]
    T_a = x[2]
    T_p = x[3]
    T_f = x[4]

    Q_rad_g2_g1 = h_rad_g2_g1 * A_c * (T_g2 - T_g1)
    Q_cov_g2_g1 = h_cov_g2_g1 * A_c * (T_g2 - T_g1)
    Q_rad_g1_sky = h_rad_g1_sky * A_c * (T_g1 - T_sky)
    Q_cov_g1_amb = h_rad_g1_amb * A_c * (T_g1 - T_amb)
    Q_rad_p_g2 = h_rad_p_g2 * A_c * (T_p - T_g2)
    Q_cov_a_g2 = h_cov_a_g2 * A_c * (T_a - T_g2)
    Q_cov_f_a = h_cov_f_a * (A_c) * (T_f - T_a)
    Q_cov_p_f = h_cov_p_f * A_c * (T_p - T_f)
    U_base = ((x_b/k_i) + 1/(h_cov_g1_amb))**(-1)
    U_side = ((x_s/k_i) + 1/(h_cov_g1_amb))**(-1)
    Q_amb_loss = (U_base*A_c + U_side*A_s)*(T_p - T_amb)



    # define each ODE
    dT_g1dt = (Irr * alpha_g * A_c + Q_rad_g2_g1 + Q_cov_g2_g1 - Q_rad_g1_sky - Q_cov_g1_amb) / M_g1
    dT_g2dt = (Irr * alpha_g * t_g * A_c + Q_rad_p_g2 + Q_cov_a_g2 - Q_rad_g2_g1) / M_g2
    dT_adt = (Q_cov_f_a - Q_cov_a_g2)/M_a
    dT_pdt = (Irr * alpha_p * t_g**2 * t_f * A_c - Q_rad_p_g2 - Q_amb_loss - Q_cov_p_f) / M_p
    dT_fdt = (Q_cov_p_f - Q_cov_f_a) / M_f

    return [dT_g1dt, dT_g2dt, dT_adt, dT_pdt, dT_fdt]

x0 = [298, 298, 298, 298, 280]


# test the defined ODES
# print(boxODE(x=x0, t=0, m_f=m_f, m_a=m_a, A_s=A_s))


# declare a time vector (time window)
t = np.linspace(0,27000,1000)
x = odeint(boxODE,x0,t, args=(m_f, m_a, A_s))

Tf_2 = x[:,4]
Tp = x[:,3]

# plot the results
plt.plot((t/3600)+5.8,Tf_2, label='fluid')
# plt.plot(t/3600,Tp, label='plate')
plt.legend()
plt.ylim(298, 338)
plt.show()

#%%

xs = np.arange(27000,28201,1)
ys = 0.0226*xs - 295

#%%

fig = plt.figure(figsize=(5,5))
fig, ax1 = plt.subplots()

plt.plot((t/3600)+5.8,Tf, color='r')
plt.plot(xs/3600 + 5.8, ys, color='r')
plt.plot(np.arange(27000,27601,1)/3600+5.8, )
plt.hlines(338, -100, 100, linestyle=':', color='k')
plt.text(6.5, 339, 'Pasteurization Temperature')

ax1.tick_params(direction='in', length=7,top=True, right=True, left=True)
minor_locator_x = AutoMinorLocator(2)
minor_locator_y = AutoMinorLocator(2)
ax1.get_xaxis().set_minor_locator(minor_locator_x)
ax1.get_yaxis().set_minor_locator(minor_locator_y)
# rotate and align the tick labels so they look better
plt.tick_params(which='minor',
                direction='in',
                length=4,
                right=True,
                left=True,
                top=True)
plt.xlim(6,21)
plt.xlabel('Hour of Day')
plt.ylim(298, 350)
plt.ylabel('Water Temperature (K)')

plt.savefig('Figures/comb_img.png', dpi=300)