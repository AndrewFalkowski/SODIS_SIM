import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import seaborn as sns
from Models.TubeThermal import TubeThermal
from Models.BoxThermal import BoxThermal
from Models.SolarIrradiance import ClearSkyIrradiance
from Models.FluidFlow import FluidFlow
import pandas as pd
from tqdm import tqdm
import numpy as np

#%%
'''
model inputs
D_ri
D_ro
D_ci
D_co
mod_width
mod_height
v_len
Gloabl Irradiance
Temp Ambient
'''

D = 0.00025
D_ro = D + 0.05
D_ci = D_ro + 0.005
D_co = D_ci+0.005

# update the module config to account for inner and outer diameter
FF = FluidFlow(D=D, v_len=0.03, mod_width=1.0, mod_height=1.0)
V = FF.velocity()
vol = FF.volume()
Re = FF.reynolds_num()
Pr = FF.prandtl_num()
# determine time it will take to move through pipes

TT = TubeThermal(W=0.05, L=FF.TP_len, f=0.1, A_a=0.5, C=2.0, D_ri=D,
                 D_ro=D_ro, D_ci=D_ci, D_co=D_co, T_am=298, G_b=1000,
                 V=V, T_in=298, Re=Re, Pr=Pr)
prop = TT.thermal_prop

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
