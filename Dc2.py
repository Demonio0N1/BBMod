import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Constantes
c = 3e10  # Velocidad de la luz en cm/s
h = 6.626e-27  # Constante de Planck en erg*s
k = 1.38e-16  # Constante de Boltzmann en erg/K
beta = 1.5  # Índice de emisividad
nu_ref = 250e9  # Frecuencia de referencia en Hz
kappa_ref = 0.4  # Coeficiente de absorción de referencia en cm^2/g
mu_H = 1.67e-24  # Masa del hidrógeno en gramos

# Función Planck
def planck_function(nu, T):
    return (2 * h * nu**3 / c**2) / (np.exp(h * nu / (k * T)) - 1)

# Ecuación del modelo de emisión de polvo
def dust_emission(nu, T_d, d, R_s, NH):
    Omega_s = np.pi * R_s**2 / d**2  # Ángulo sólido
    tau_nu = kappa_ref * (nu / nu_ref)**beta * NH * mu_H / 100  # Profundidad óptica
    B_nu = planck_function(nu, T_d)  # Función de Planck
    return B_nu * (1 - np.exp(-tau_nu)) * Omega_s

# Configuración inicial
nu = np.logspace(10, 14, 500)  # Frecuencia en Hz
T_d_init = 30  # Temperatura inicial en K
d_init = 1e27  # Distancia inicial en cm
R_s_init = 1e19  # Radio inicial en cm
NH_init = 1e21  # Densidad de columna inicial en cm^-2

# Crear figura y gráfico inicial
fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.35)
l, = plt.plot(nu, dust_emission(nu, T_d_init, d_init, R_s_init, NH_init), label="Emisión de polvo")
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Frecuencia (Hz)')
ax.set_ylabel('Intensidad ($S_\\nu$)')
ax.legend()

# Crear sliders
ax_Td = plt.axes([0.25, 0.25, 0.65, 0.03])
ax_d = plt.axes([0.25, 0.2, 0.65, 0.03])
ax_Rs = plt.axes([0.25, 0.15, 0.65, 0.03])
ax_NH = plt.axes([0.25, 0.10, 0.65, 0.03])

slider_Td = Slider(ax_Td, 'T_d (K)', 10, 100, valinit=T_d_init)
slider_d = Slider(ax_d, 'd (cm)', 1e26, 1e28, valinit=d_init, valstep=1e26)
slider_Rs = Slider(ax_Rs, 'R_s (cm)', 1e18, 1e20, valinit=R_s_init, valstep=1e18)
slider_NH = Slider(ax_NH, 'NH (cm⁻²)', 1e21, 9e22, valinit=NH_init, valstep=1e21)

# Actualizar función
def update(val):
    T_d = slider_Td.val
    d = slider_d.val
    R_s = slider_Rs.val
    NH = slider_NH.val
    l.set_ydata(dust_emission(nu, T_d, d, R_s, NH))
    fig.canvas.draw_idle()

# Conectar sliders con la función de actualización
slider_Td.on_changed(update)
slider_d.on_changed(update)
slider_Rs.on_changed(update)
slider_NH.on_changed(update)

plt.show()