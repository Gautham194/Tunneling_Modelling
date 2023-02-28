import numpy as np
from matplotlib import pyplot as plt

# Constants

q = 1.6e-19
pi = np.pi
h = 6.6261e-34
hbar = 1.0546e-34
m_0 = 9.109e-31
m_hbn = m_0 * 2.21

hbn_bg = 5.955  # Bandgap in eV
hbn_ea = 2.0  # Electron Affinity in eV
hbn_dc = 3.76  # Dielectric Constant

grp_dp = 4.5  # Dirac Point
grp_dc = 3  # epsilon out-of-plane in 2.5 for SL, 2.6 BL, estimate 3 for multilayer

inse_bg = 1.4
inse_ea = 4.55
inse_dc = 8  # see Politano et al, 2016

barrier_height = (inse_ea - hbn_ea) * q

# Variables

s1 = 2e-10  # Approx contact area between InSe + hBN for 20x10um flake
s2 = 1e-10  # Approx contact area between Gate + hBN for 10x10um gate

d_hbn_cap = 20e-9  # Thickness of encapsulating layer - Assume 0.4nm for hBN
d_grp = 3e-9  # Assume 0.35nm for SL
d_inse = 5e-9  # Assume 0.8nm for SL

V = 5


# Use Capacitor Approximation
# Area = 4micron * 10micron area of top gate

def get_E_tunneling(d_hbn_tunnel=5e-9, volt=V):
    top = hbn_dc * grp_dc * inse_dc
    b_1 = d_hbn_cap * hbn_dc * inse_dc * grp_dc
    b_2 = d_grp * hbn_dc * hbn_dc * inse_dc
    b_3 = d_hbn_tunnel * hbn_dc * inse_dc * grp_dc
    b_4 = d_inse * hbn_dc * hbn_dc * grp_dc
    e_tunnel = top * volt / (b_1 + b_2 + b_3 + b_4)
    return e_tunnel


def get_fn_current(e_tunnel):
    c_1 = (q ** 3) / (8 * pi * h * barrier_height)
    c_2 = (-4 * (2 * m_hbn) ** 2 * barrier_height ** (3 / 2)) / (3 * hbar * q)
    j = c_1 * (e_tunnel ** 2) * np.exp(c_2 * (1 / e_tunnel))
    return j


def loop_over_t(num_layers, volt = 5):
    arr_1 = np.linspace(1, num_layers, num_layers)
    d_tunnels = arr_1 * 0.4e-9
    j_arr = []

    for thickness in d_tunnels:
        var1 = get_E_tunneling(d_hbn_tunnel=thickness, volt=volt)
        var2 = get_fn_current(var1)
        j_arr.append(var2)

    plt.plot(d_tunnels, j_arr)
    plt.show()


def loop_over_v(peak_v, thickness=5e-9):
    arr_1 = np.linspace(1, 1000, 1000)
    arr_2 = arr_1 * (peak_v / 1000)
    j_arr = []

    for voltage in arr_2:
        var1 = get_E_tunneling(d_hbn_tunnel=thickness, volt=voltage)
        var2 = get_fn_current(var1)
        j_arr.append(var2)

    plt.plot(arr_2, j_arr)
    plt.show()

loop_over_t(1000)
