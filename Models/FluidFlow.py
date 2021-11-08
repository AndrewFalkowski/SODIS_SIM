import numpy as np
import fluids # Good practice
from fluids.fittings import *

#%%

# def module_config(mod_width, mod_height, D, v_len):
#     '''
#     mod_width = module width in meters
#     mod_height = module height in meters
#     D = pipe diameter in meters
#     v_len = length of vertical section of pipe
#     '''
#     try:
#         assert v_len > D
#     except:
#         print('\nv_len is too short, redefining v_len as D+0.01 and continuing calculation...')
#         v_len = D + 0.01
#     h_len = mod_width - D
#     unit_len = h_len + v_len
#     num_units = np.floor(mod_height / v_len)
#     total_pipe_len = (num_units * unit_len)
#     vert_pipe_len = num_units * v_len
#     return num_units, h_len, v_len, total_pipe_len

# def fluid_velocity(D, L=None, C=140):
#     # calculate fluid velocity with the HazenWilliams Equation
#     A = (D/2)**2 * np.pi
#     P = 2 * np.pi * (D/2)
#     k = 0.849
#     hydR = A / P
#     n_units, h_len, v_len, total = module_config(1.0, 1.0, D, 0.06)
#     if L is not None:
#         total = L
#         print(L)
#     S = (n_units * v_len) / total
#     v = k * C * hydR**0.63 * S**0.54 # m/s
#     return v

# def fittings_loss(num_units, D, total_len):
#     Re = fluids.Reynolds(V=fluid_velocity(D), D=D, rho=997, mu=0.38)
#     fd = fluids.friction.friction_factor(Re)
#     K = fluids.K_from_f(fd=fd, L=total_len, D=D)
#     K = entrance_sharp()

#     for n in np.arange(num_units):
#         K += 2 * bend_rounded(D, 90, method='Swamee')
#     return K

# def friction_loss(D):
#     Re = fluids.Reynolds(V=fluid_velocity(D), D=D, rho=997, mu=0.00038)
#     fd = fluids.friction.friction_factor(Re)
#     return fd

# D = 0.05
# nunits, h_len, v_len, total_pipe_len = module_config(1.0,1.0,0.05,0.0)
# K = fittings_loss(nunits, D, total_pipe_len)
# fd = friction_loss(D)
# Leq = fluids.L_from_K(K, D, fd)
# v = fluid_velocity(D, L=Leq)
# v = fluid_velocity(D )


#%%


class FluidFlow:
    def __init__(self, D_co, D_ri, v_len, mod_width, mod_height, C=150):
        self.D_co = D_co
        self.D_ri = D_ri
        self.v_len = v_len
        self.mod_width = mod_width
        self.mod_height = mod_height
        self.C = C
        self.n_units, self. h_len, self.v_len, self.TP_len = self.module_config()
        self.v = self.velocity()

    def module_config(self):
        '''
        mod_width = module width in meters
        mod_height = module height in meters
        D = pipe diameter in meters
        v_len = length of vertical section of pipe
        '''
        try:
            assert self.v_len > self.D_co
        except:
            print('\nv_len is too short, redefining v_len as D+0.01 and continuing calculation...')
            self.v_len = self.D_co + 0.01
        h_len = self.mod_width - self.D_co
        unit_len = h_len + self.v_len
        num_units = np.floor(self.mod_height / self.v_len)
        total_pipe_len = (num_units * unit_len)
        vert_pipe_len = num_units * self.v_len
        return num_units, h_len, self.v_len, total_pipe_len

    def velocity(self):
        # calculate fluid velocity with the HazenWilliams Equation
        A = (self.D_ri/2)**2 * np.pi
        P = 2 * np.pi * (self.D_ri/2)
        k = 0.849
        hydR = A / P
        S = (self.n_units * self.v_len) / self.TP_len
        self.v = k * self.C * hydR**0.63 * S**0.54 # m/s
        return self.v

    def volume(self):
        self.vol = self.D_ri * self.TP_len
        return self.vol

    def reynolds_num(self):
        Re = fluids.Reynolds(self.v, D=self.D_ri, rho=997, mu=0.00038)
        return Re

    def prandtl_num(self):
        Pr = fluids.Prandtl(Cp=4190, k=0.5899, mu=0.00038)
        return Pr

