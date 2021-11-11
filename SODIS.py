import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import seaborn as sns
from Models.TubeThermal import TubeThermal
from Models.BoxThermal import BoxThermal
from Models.SolarIrradiance import ClearSkyIrradiance
from Models.FluidFlow import FluidFlow
from Utils.utils import build_mirror
from BoxThermal import *
import pandas as pd
from tqdm import tqdm
import numpy as np
import itertools
from scipy.integrate import odeint

#%%

# single model param prediction

m_f, m_a, A_s = box_dim(0.25, 0.15, 0.9) # calculate relavent box info
x0 = [298, 298, 298, 298, 298] # initial conditions for box receiver

# calculate pretreatment module




D_in = 0.001
v_len = 0.14
T_in = 298

model_params = {
    'D_ri': D_in, # meters
    'D_ro': D_in+(D_in*0.1), # meters
    'D_ci': D_in+(D_in*0.2), # meters
    'D_co': D_in+(D_in*0.3), # meters
    'mod_width': 1.0, # meters
    'mod_height': 0.75, # meters
    'tilt' : 0,
    'v_len': v_len,
    'date': '2021-07-15',
    'time': '13:05:00',
    'Tam':298,
}

mp = pd.Series(model_params)

CSI = ClearSkyIrradiance()
irrad_total = CSI.get_irradiance(date=mp.date, tilt=mp.tilt, surface_azimuth=180)
mp['Irr'] = CSI.irrad_at_time(mp.date, mp.time)['GHI']
trange = np.arange(0,1440,1)
irrad_total.index = trange

Touts = np.zeros_like(Tf, dtype='float64')
for i,t in enumerate(trange):
    mp['Irr'] = irrad_total['GHI'].loc[t]

# calculate fluid flow

    FF = FluidFlow(D_co=mp.D_co, D_ri=mp.D_ri, v_len=mp.v_len, mod_width=mp.mod_width, mod_height=mp.mod_height)
    mp['tube_vol'] = FF.volume()
    mp['Re'] = FF.reynolds_num()/2
    mp['Pr'] = FF.prandtl_num()
    mp['n_tube_units'] = FF.n_units
    mp['total_tube_len'] = FF.TP_len
    mp['V'] = FF.velocity()

    mp = build_mirror(mp) # update model parameters with mirror info


    TT = TubeThermal(TP_len=mp.total_tube_len, v_len = mp.v_len, n_tube_units=mp.n_tube_units,
                     mod_width=mp.mod_width, D_ri=mp.D_ri, D_ro=mp.D_ro, D_ci=mp.D_ci,
                     D_co=mp.D_co, T_am=mp.Tam, Irr=mp.Irr, V=mp.V, T_in=T_in, Re=mp.Re, Pr=mp.Pr,
                     Wa=mp.Wa, C=mp.C)

    Touts[i] = TT.thermal_prop.T_out - 298
    # prop = TT.thermal_prop
#%%

t = np.linspace(0,63000,900)
Tf = odeint(boxODE,x0, t, args=(m_f, m_a, A_s))[:,4]
before = 298 * np.ones((330,))
Tf = np.insert(Tf, 0, before)
after = 293.251 * np.ones(((1440 - len(Tf)),))
Tf = np.append(Tf, after)

#heat decay control
td = np.arange(720,1440,1)
y = 0.02*np.cos(0.0043633*td)+1.02
decay = np.insert(y, 0, np.ones(720,))
Tf = Tf * decay
Tf[1190:1440] = 302.4

plt.plot((trange / 60),((Tf+Touts)))
plt.plot((trange / 60),Tf)
plt.plot((trange / 60),Touts+298)
plt.hlines(338, -10,1500, color='k')
plt.xlim(0,24)
#%%

# method for iterating over many possible combinations of parameters

D_in = np.arange(0.001,0.051, 0.01)
v_lens = np.arange(0.01, 0.401, 0.01)
combs = itertools.product(D_in, v_lens)
combs_len = len(list(combs))
values = np.zeros((combs_len, 3))

for i, params in tqdm(enumerate(itertools.product(D_in, v_lens)), total=combs_len):
    D_in , v_len = params
    model_params = {
        'D_ri': D_in, # meters
        'D_ro': D_in+(D_in*0.1), # meters
        'D_ci': D_in+(D_in*0.2), # meters
        'D_co': D_in+(D_in*0.3), # meters
        'mod_width': 1.0, # meters
        'mod_height': 1.0, # meters
        'tilt' : 0,
        'v_len': v_len,
        'date': '2021-07-15',
        'time': '12:05:00',
        'Tam':298,
    }


    mp = pd.Series(model_params)

    CSI = ClearSkyIrradiance()
    irrad_total = CSI.get_irradiance(date=mp.date, tilt=mp.tilt, surface_azimuth=180)
    mp['Irr'] = CSI.irrad_at_time(mp.date, mp.time)['GHI']

    # calculate pretreatment module

    FF = FluidFlow(D_co=mp.D_co, D_ri=mp.D_ri, v_len=mp.v_len, mod_width=mp.mod_width, mod_height=mp.mod_height)
    mp['tube_vol'] = FF.volume()
    mp['Re'] = FF.reynolds_num()/2
    mp['Pr'] = FF.prandtl_num()
    mp['n_tube_units'] = FF.n_units
    mp['total_tube_len'] = FF.TP_len
    mp['V'] = FF.velocity()

    mp = build_mirror(mp) # update model parameters with mirror info


    TT = TubeThermal(TP_len=mp.total_tube_len, v_len = mp.v_len, n_tube_units=mp.n_tube_units,
                     mod_width=mp.mod_width, D_ri=mp.D_ri, D_ro=mp.D_ro, D_ci=mp.D_ci,
                     D_co=mp.D_co, T_am=mp.Tam, Irr=mp.Irr, V=mp.V, T_in=298, Re=mp.Re, Pr=mp.Pr,
                     Wa=mp.Wa, C=mp.C)

    prop = TT.thermal_prop
    values[i,:] = [D_in, v_len, prop.T_out]

#%%
V = FF.velocity()
vol = FF.volume()
Re = FF.reynolds_num()
Pr = FF.prandtl_num()
# determine time it will take to move through pipes


# find way to calculate time in residence.

#%%

d_range = np.arange(0.005, 0.0501, 0.005)
heat_df = pd.DataFrame()
heat_df.index = np.arange(600, 1050, 50)
for d in tqdm(d_range):
    rates = []
    L=1.0
    D_ri = d
    D_ro = D_ri + 0.002
    D_ci = D_ro + 0.005
    D_co = D_ci + 0.002
    G_b_range = np.arange(600, 1050, 50)
    for Gb in G_b_range:
        TT = TubeThermal(W=0.05, L=L, f=0.1, A_a=0.5, C=2.0, D_ri=D_ri,
                         D_ro=D_ro, D_ci=D_ci, D_co=D_co, T_am=298, G_b=Gb,
                         V=1, T_in=298)
        prop = TT.thermal_prop
        del_T = 338 - TT.mp.T_in
        ts = (4200 * del_T) / prop.Q_u
        # print(f'iter: {del_T}  {prop.Q_u}')
        rates.append(ts)
    heat_df[f'{d}'] = rates

#%%
d_range = np.arange(0.005, 0.0501, 0.005)
heat_df = pd.DataFrame()
heat_df.index = np.arange(0.001, 1.01, 0.001)
for d in tqdm(d_range):
    rates = []
    D_ri = d
    D_ro = D_ri + 0.002
    D_ci = D_ro + 0.005
    D_co = D_ci + 0.002
    v_range = np.arange(0.001, 1.01, 0.001)
    for V in v_range:
        TT = TubeThermal(W=0.05, L=0.3, f=0.1, A_a=0.5, C=2.0, D_ri=D_ri,
                         D_ro=D_ro, D_ci=D_ci, D_co=D_co, T_am=298, G_b=950,
                         V=V, T_in=298)
        prop = TT.thermal_prop
        del_T = prop.T_out - TT.mp.T_in
        rate = del_T/0.3
        rates.append(rate)
    heat_df[f'{d}'] = rates

#%%
colors = sns.color_palette("rocket", 10)
d_range = np.arange(0.005, 0.0501, 0.005)
fig = plt.figure(figsize=(5,5))
fig, ax1 = plt.subplots()
for i, d in enumerate(d_range):
    plt.semilogx(heat_df.index, heat_df[f'{d}'], linestyle='-', marker=None,
             color=colors[i], alpha=0.9, label = f'{round(d, 3)} (m)')
plt.legend(loc='upper right', framealpha=0.95)
ax1.yaxis.set_label_position("left")
# ax1.yaxis.tick_right()
ax1.tick_params(direction='in', length=7,top=True, right=True, left=True)
# minor_locator_x = AutoMinorLocator(2)
# minor_locator_y = AutoMinorLocator(2)
# ax1.get_xaxis().set_minor_locator(minor_locator_x)
# ax1.get_yaxis().set_minor_locator(minor_locator_y)
# import matplotlib.dates as mdates
# ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
# # rotate and align the tick labels so they look better
# fig.autofmt_xdate()
plt.tick_params(which='minor',
                direction='in',
                length=4,
                right=True,
                left=True,
                top=True)
plt.xlabel('Velocity (m/s)')
plt.ylabel('Heating Rate ($^o$C/m)')

plt.show()




#%%

BT = BoxThermal(vol=0.1, A_g=0.38, A_a=0.38, T_am=298, G_b=950, V=0.001, T_in=298)
prop = BT.thermal_prop
print(prop.T_fm)


#%%

CSI = ClearSkyIrradiance()
july_irrad = CSI.get_irradiance(date='07-15-2021', tilt=25, surface_azimuth=180)
aug_irrad = CSI.get_irradiance(date='08-15-2021', tilt=25, surface_azimuth=180)
sept_irrad = CSI.get_irradiance(date='09-15-2021', tilt=25, surface_azimuth=180)

irrad_list = [july_irrad, aug_irrad, sept_irrad]
months =['July', 'August', 'September']

colors = sns.color_palette("rocket_r", 3)

fig = plt.figure(figsize=(5,5))
fig, ax1 = plt.subplots()
for i, month in enumerate(irrad_list):
    plt.plot(july_irrad.index, month['GHI'], linestyle='-', marker=None,
             color=colors[i], alpha=1.0, label = f'{months[i]}')
plt.legend(loc='upper left', framealpha=0.95)
ax1.yaxis.set_label_position("right")
ax1.yaxis.tick_right()
ax1.tick_params(direction='in', length=7,top=True, right=True, left=True)
minor_locator_x = AutoMinorLocator(2)
minor_locator_y = AutoMinorLocator(2)
ax1.get_xaxis().set_minor_locator(minor_locator_x)
ax1.get_yaxis().set_minor_locator(minor_locator_y)
import matplotlib.dates as mdates
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
# rotate and align the tick labels so they look better
fig.autofmt_xdate()
plt.tick_params(which='minor',
                direction='in',
                length=4,
                right=True,
                left=True,
                top=True)
plt.ylim(-30,1000)
plt.xlabel('Time')
plt.ylabel('GHI Irradiance (W/m$^2$)')

plt.draw()

#%%

CSI = ClearSkyIrradiance()


irrads = pd.DataFrame()
irrads['0'] = CSI.get_irradiance(date='09-15-2021', tilt=0, surface_azimuth=180)['GHI']
tilts = np.arange(3,40,3)

for tilt in tilts:
    july_irrad = CSI.get_irradiance(date='09-15-2021', tilt=tilt, surface_azimuth=180)
    irrads[f'{tilt}'] = july_irrad['POA']

tilts = np.arange(0,40,3)
colors = sns.color_palette("icefire", 14)

fig = plt.figure(figsize=(5,5))
fig, ax1 = plt.subplots()
for i, tilt in enumerate(tilts):
    plt.plot(irrads[f'{tilt}'], linestyle='-', marker=None,
             color=colors[i], alpha=0.6, label = f'{tilts[i]}$^o$')
plt.legend(loc='upper left', framealpha=0.95)
ax1.yaxis.set_label_position("right")
ax1.yaxis.tick_right()
ax1.tick_params(direction='in', length=7,top=True, right=True, left=True)
minor_locator_x = AutoMinorLocator(2)
minor_locator_y = AutoMinorLocator(2)
ax1.get_xaxis().set_minor_locator(minor_locator_x)
ax1.get_yaxis().set_minor_locator(minor_locator_y)
import matplotlib.dates as mdates
ax1.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
# rotate and align the tick labels so they look better
fig.autofmt_xdate()
plt.tick_params(which='minor',
                direction='in',
                length=4,
                right=True,
                left=True,
                top=True)
plt.ylim(750,1000)
plt.xlabel('Time')
plt.ylabel('POA Irradiance (W/m$^2$)')

plt.draw()


# fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
# summer_irrad['GHI'].plot(ax=ax1, label='GHI')
# summer_irrad['POA'].plot(ax=ax1, label='POA')
# winter_irrad['GHI'].plot(ax=ax2, label='GHI')
# winter_irrad['POA'].plot(ax=ax2, label='POA')
# ax1.set_xlabel('Time of day (Summer)')
# ax2.set_xlabel('Time of day (Winter)')
# ax1.set_ylabel('Irradiance ($W/m^2$)')
# ax1.legend()
# ax2.legend()
# plt.show()

#%%

# df = pd.read_clipboard()
# df2 = pd.read_clipboard()

colors = sns.color_palette('coolwarm', 2)
fig, axs = plt.subplots(2)
axs[0].plot(df['hour'], df['temp'], color=colors[1], mec='k', linestyle='None', alpha=0.9, marker='s', label='Internal Temperature')
axs[0].plot(df['hour'], df['trend'], color='red', label='Temperature Trend')
axs[1].plot(df2['hour'], df2['irrad'], color=colors[0], mec='k', linestyle='None', alpha=0.9, marker='s')
axs[0].legend(loc='best')


axs[0].set_ylabel('Internal Temperature ($^o$C)')
# axs[0].tick_params(axis='y')
axs[0].tick_params(direction='in',
                    length=7,top=True, right=True)
axs[0].set_ylim(0,60)
axs[0].set_xlim(0,80)
minor_locator_x = AutoMinorLocator(2)
minor_locator_y = AutoMinorLocator(2)
axs[0].get_xaxis().set_minor_locator(minor_locator_x)
axs[0].get_yaxis().set_minor_locator(minor_locator_y)
axs[0].tick_params(which='minor',
                direction='in',
                length=4,
                right=True,
                top=True)


axs[1].set_xlabel('Time (hours)')
axs[1].set_ylabel('Irradiance (W/m$^2$)')
# axs[0].tick_params(axis='y')
axs[1].tick_params(direction='in',
                    length=7,top=True, right=True)
axs[1].set_ylim(0,750)
axs[1].set_xlim(0,80)
minor_locator_x = AutoMinorLocator(2)
minor_locator_y = AutoMinorLocator(2)
axs[1].get_xaxis().set_minor_locator(minor_locator_x)
axs[1].get_yaxis().set_minor_locator(minor_locator_y)
plt.tick_params(which='minor',
                direction='in',
                length=4,
                right=True,
                top=True)


#%%


color = colors[1]
ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('Irradiance', color=color)
ax1.plot(df2['hour'], df2['irrad'], color=color, mec='k', alpha=0.5, marker='s')
ax1.tick_params(axis='y', labelcolor=color)
ax1.tick_params(direction='in',
                    length=7,top=True, right=True)

minor_locator_x = AutoMinorLocator(2)
minor_locator_y = AutoMinorLocator(2)
ax1.get_xaxis().set_minor_locator(minor_locator_x)
ax1.get_yaxis().set_minor_locator(minor_locator_y)
plt.tick_params(which='minor',
                direction='in', labelcolor=color,
                length=4,
                right=True,
                top=True)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = colors[0]
ax2.set_ylabel('Internal Temperature ($^o$C)', color=color)
ax2.plot(df['hour'], df['temp'], color=color, mec='k', alpha=0.5, marker='s')
ax2.tick_params(axis='y', labelcolor=color, direction='in',
                    length=7)
# ax2.set_ylim(120,170)
minor_locator_x = AutoMinorLocator(2)
minor_locator_y = AutoMinorLocator(2)
ax2.get_xaxis().set_minor_locator(minor_locator_x)
ax2.get_yaxis().set_minor_locator(minor_locator_y)
plt.tick_params(which='minor',
                direction='in', labelcolor=color,
                length=4,
                right=True,
                top=True)

plt.draw()


#%%

# df = pd.read_clipboard()

colors = sns.color_palette("rocket_r", 3)

fig = plt.figure(figsize=(5,5))
fig, ax1 = plt.subplots()

plt.scatter(df['Irrad'], df['Temp'], color=colors[0], label='Experimental Data')
plt.plot(df['Irrad'], df['Trend'], label='Linear Fit', color=colors[2])
plt.text(200, 46, 'R$^2$ = 0.321')


plt.legend(loc='upper left', framealpha=0.95)
ax1.yaxis.set_label_position("left")
# ax1.yaxis.tick_right()
ax1.tick_params(direction='in', length=7,top=True, right=True, left=True)
minor_locator_x = AutoMinorLocator(2)
minor_locator_y = AutoMinorLocator(2)
ax1.get_xaxis().set_minor_locator(minor_locator_x)
ax1.get_yaxis().set_minor_locator(minor_locator_y)
plt.tick_params(which='minor',
                direction='in',
                length=4,
                right=True,
                left=True,
                top=True)
plt.xlabel('Irradiance (W/m$^2$)')
plt.ylabel('Internal Temperature ($^o$C)')
plt.ylim(30,50)
# plt.draw()
