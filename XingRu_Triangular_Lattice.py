##Python Script for Triangular Lattice
from scipy import * 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm  
import pylab as pl
import scipy.sparse as sp
from numpy import linalg as LA

def tri_lattice_plot(Gamma, a):
    sq_3 = sqrt(3)
    
    d_x = d_y =1000
    
    ##For plotting of E(k) verus k along the trajectory Gamma to M to K to Gamma 
    ##Gamma to M
    k_x_1 = zeros(d_x)
    k_y_1 = linspace(0, ((2*pi)/(sqrt(3)*a)), d_y)
    
    E_k_11 = 8*Gamma - 4*Gamma*(cos((1/2)*k_x_1*a))*(cos((1/2)*k_x_1*a) + cos((sq_3/2)*k_y_1*a))
    E_k_21 = 8*Gamma - 4*Gamma*(cos((1/2)*k_x_1*a))*(cos((1/2)*k_x_1*a) + cos((sq_3/2)*k_y_1*a))
    
    ##M to K
    k_x_2 = linspace(0, ((2*pi)/(3*a)), d_x)
    k_y_2 = ((2*pi)/(sqrt(3)*a))*ones(d_y)
    
    E_k_12 = 8*Gamma - 4*Gamma*(cos((1/2)*k_x_2*a))*(cos((1/2)*k_x_2*a) + cos((sq_3/2)*k_y_2*a))
    E_k_22 = 8*Gamma - 4*Gamma*(cos((1/2)*k_x_2*a))*(cos((1/2)*k_x_2*a) + cos((sq_3/2)*k_y_2*a))
    
    
    ##K to Gamma
    k_x_3 = linspace(((2*pi)/(3*a)), 0, d_x)
    k_y_3 = linspace(((2*pi)/(sqrt(3)*a)), 0, d_y)
    
    E_k_13 = 8*Gamma - 4*Gamma*(cos((1/2)*k_x_3*a))*(cos((1/2)*k_x_3*a) + cos((sq_3/2)*k_y_3*a))
    E_k_23 = 8*Gamma - 4*Gamma*(cos((1/2)*k_x_3*a))*(cos((1/2)*k_x_3*a) + cos((sq_3/2)*k_y_3*a))
    
    E_k1 = concatenate((E_k_11, E_k_12, E_k_13))
    E_k2 = concatenate((E_k_21, E_k_22, E_k_23))
    
    plt.figure()
    plt.title(r'$E(\vec{k})$ vs $\vec{k}$ graph')
    plt.plot(E_k2,'b', label = 'Energy band')
    plt.xlabel(r'$\vec{k}$')
    plt.ylabel(r'$E(\vec{k})$ in eV')
    plt.legend(loc='upper right')
    ax = plt.gca() # grabbing the current axis
    plt.xticks([])
    ax.set_xticks([0,(1000-1), (2000-1), (3000-1)]) 
    ax.set_xticklabels([r"$\Gamma$", r"$M$", r"$K$", r"$\Gamma$"]) 
    
Gamma=1.0
a = 1.0
tri_lattice_plot(Gamma, a)
sq_3 = sqrt(3)

# Setting up meshgrid for plotting of 3d graph for energy surface of triangular lattice
kvec = linspace(-pi,pi,(100+1))
kx,ky = meshgrid(kvec,kvec)
E = 8*Gamma - 4*Gamma*(cos((1/2)*kx*a))*(cos((1/2)*kx*a) + cos((sq_3/2)*ky*a))

# Plot 3d graph for energy surface of triangular lattice
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(kx,ky,E,cmap='viridis',vmin=-E.max(),vmax=E.max(),rstride=1,cstride=1)
plt.title(r"Energy dispersion graph of $E(\vec{k})$ for triangular lattice")
ax.set_xlabel(r"$k_{x}$")
ax.set_ylabel(r"$k_{y}$")
ax.set_zlabel(r"$E(\vec{k})$")