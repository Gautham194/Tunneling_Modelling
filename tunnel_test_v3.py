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
A = (25e-6 ** 2)

d_hbn_cap = 20e-9  # Thickness of encapsulating layer - Assume 0.4nm for hBN
d_grp = 3e-9  # Assume 0.35nm for SL
d_inse = 5e-9  # Assume 0.8nm for SL

V = 5


# Use Capacitor Approximation
# Area = 10micron * 10micron area of top gate

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


def loop_over_t(num_layers, volt=5):
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


def get_charge_stored(input_v, input_t):
    e_tunnel = get_E_tunneling(volt=input_v)
    current = get_fn_current(e_tunnel)
    charge = current * input_t
    return charge


def get_vs(v_in, d_hbn=5e-9):
    base = (d_inse * hbn_dc * grp_dc * hbn_dc) + (d_hbn * inse_dc * hbn_dc * grp_dc) + (
                d_grp * inse_dc * hbn_dc * hbn_dc) + (d_hbn_cap * inse_dc * hbn_dc * grp_dc)
    v3 = (((hbn_dc * hbn_dc * grp_dc) * d_inse) / base) * v_in
    v2 = ((inse_dc * hbn_dc * grp_dc) * d_hbn / base) * (v_in)
    v1 = ((inse_dc * hbn_dc * hbn_dc) * d_grp / base) * (v_in)
    v0 = ((inse_dc * hbn_dc * grp_dc) * d_hbn_cap / base) * (v_in)
    return [v3, v2, v1, v0]

def get_q_top(v_in, d_hbn=5e-9):
    c_top = (inse_dc * hbn_dc * grp_dc * hbn_dc) * A
    c_bot = (d_inse * hbn_dc * grp_dc * hbn_dc) + (d_hbn * inse_dc * hbn_dc * grp_dc) + (
                d_grp * inse_dc * hbn_dc * hbn_dc) + (d_hbn_cap * inse_dc * hbn_dc * grp_dc)
    c = c_top/c_bot
    q = c * v_in
    return q

def get_v_gate(q_gate, d_hbn = 5e-9):
    c_top = (hbn_dc * inse_dc) * A
    c_bot = (d_hbn * inse_dc) + (d_inse * hbn_dc)
    c_gate = c_top/c_bot
    v_gate = q_gate / c_gate
    return v_gate

def get_v_gate_alt(q_gate, d_hbn = 5e-9):
    c = (hbn_dc * A) / d_hbn
    vg = q_gate / c
    return vg

def get_e_gate(q_gate, d_hbn=5e-9):
    e_gate = -1 * q_gate/(hbn_dc * A)
    return e_gate


a = get_E_tunneling(5e-9, 5e-3)
b = get_fn_current(a)
c = get_v_gate_alt(b*5e-9)
d = get_v_gate(b*5e-9)
e = get_e_gate(b*5e-9)
f = get_fn_current(a+e)
print(a, b, c, d, e, f)

# Set Loop

q_gate = 0
pulse_time = 5e-7
pulse_int = 5e-10
t = 0
app_v = 5e-3
d_t = 5e-9

times = np.linspace(1, 1000, 1000) * pulse_int
q_gates = []

for time in times:
    e1 = get_E_tunneling(d_t, app_v)
    e2 = get_e_gate(q_gate, d_t)
    fn = get_fn_current(e1+e2)
    q_int = fn * pulse_int
    q_gate += q_int
    q_gates.append(q_gate)

print(q_gates)
plt.plot(q_gates)
plt.show()

dum1 = get_v_gate(q_gates[-1])
dum2 = get_v_gate_alt(q_gates[-1])
print(dum1, dum2)
