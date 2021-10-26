#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 10:57:22 2021

@author: andrewf
"""

import numpy as np
import pandas as pd
from Utils.utils import *
#%%

class BoxThermal:
    def __init__(self, vol=0.015, A_g=0.5, A_a=0.5, T_am=298, G_b=1000, V=0.001, T_in=298):

        self.model_params = {
            "vol": vol, # cubic meters
            'A_ri': A_a, # reciever inner A (m^2)
            'A_ro': A_a, # reciever outer A (m^2)
            'A_ci': A_g, # cover inner A (m^2)
            'A_co': A_g, # cover outer A (m^2)
            'E_r': 0.2, # reciever emittance
            'E_c': 0.9, # cover emittance
            'a_a': 0.96, # absorber absorbance
            't_c': 0.95, # cover transmittance
            'sig': 5.67E-8, # W / m^2 K
            'T_am': T_am, # ambient temperature (K)
            'G_b': G_b, # solar direct beam irradiation (W/m^2)
            'Q_s': 0.0, # Solar heat flux (W)
            'h_out': 10, # transfer coefficent for cover (W/m^2 K)
            'h_in': 0.0, # heat transfer coefficient for absorber-fluid (W/m^2 K)
            'V': V , # volumetric flow rate (m^3/s)
            'm': 0.0, # mass flow rate (kg/s) # calc
            'C_p': 4.169, # specific heat J/kg K
            'T_in': T_in, # inlet temperature (K)
            'n_opt': 0.0 # max optical efficiency # calc
            }
        self.mp = pd.Series(self.model_params)
        self.mp = calc_box_params(self.mp)

        self.calc_properties()

    def calc_properties(self):
        thermal_prop = {
            'Q_loss' : calc_Q_loss(self.mp),
            'T_r' : calc_T_r(self.mp),
            'T_c' : calc_T_c(self.mp),
            'T_out': calc_T_out(self.mp),
            'T_fm': calc_T_fm(self.mp)
            }

        self.thermal_prop = pd.Series(thermal_prop)

        return self.thermal_prop

#%%

