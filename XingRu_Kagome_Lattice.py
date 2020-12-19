##Python Script for Kagome Lattice
from scipy import * 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm  
import pylab as pl
import scipy.sparse as sp
from numpy import linalg as LA

def Kagome_Lattice_plot(Gamma, a):
    sq_3 = sqrt(3)
    
    d_x = d_y =1000
    
    ##For plotting of E(k) verus k along the trajectory Gamma to M to K to Gamma 
    ##Gamma to M
    k_x_1 = zeros(d_x)
    k_y_1 = linspace(0, ((2*pi)/(sqrt(3)*a)), d_y)
    
    E_k_11 = -Gamma + Gamma*sqrt( 4*( cos((1/2)*k_x_1*a)*cos((1/2)*k_x_1*a) + cos((1/4)*k_x_1*a + (sq_3/4)*k_y_1*a)*cos((1/4)*k_x_1*a + (sq_3/4)*k_y_1*a) + cos((1/4)*k_x_1*a - (sq_3/4)*k_y_1*a)*cos((1/4)*k_x_1*a - (sq_3/4)*k_y_1*a)   )  -3    )
    E_k_21 = -Gamma - Gamma*sqrt( 4*( cos((1/2)*k_x_1*a)*cos((1/2)*k_x_1*a) + cos((1/4)*k_x_1*a + (sq_3/4)*k_y_1*a)*cos((1/4)*k_x_1*a + (sq_3/4)*k_y_1*a) + cos((1/4)*k_x_1*a - (sq_3/4)*k_y_1*a)*cos((1/4)*k_x_1*a - (sq_3/4)*k_y_1*a)   )  -3    )
    E_k_31 = 2*Gamma*ones(d_x)
    
    ##M to K
    k_x_2 = linspace(0, ((2*pi)/(3*a)), d_x)
    k_y_2 = ((2*pi)/(sqrt(3)*a))*ones(d_y)
    
    E_k_12 = -Gamma + Gamma*sqrt( 4*( cos((1/2)*k_x_2*a)*cos((1/2)*k_x_2*a) + cos((1/4)*k_x_2*a + (sq_3/4)*k_y_2*a)*cos((1/4)*k_x_2*a + (sq_3/4)*k_y_2*a) + cos((1/4)*k_x_2*a - (sq_3/4)*k_y_2*a)*cos((1/4)*k_x_2*a - (sq_3/4)*k_y_2*a)   )  -3    )
    E_k_22 = -Gamma - Gamma*sqrt( 4*( cos((1/2)*k_x_2*a)*cos((1/2)*k_x_2*a) + cos((1/4)*k_x_2*a + (sq_3/4)*k_y_2*a)*cos((1/4)*k_x_2*a + (sq_3/4)*k_y_2*a) + cos((1/4)*k_x_2*a - (sq_3/4)*k_y_2*a)*cos((1/4)*k_x_2*a - (sq_3/4)*k_y_2*a)   )  -3    )
    E_k_32 = 2*Gamma*ones(d_x)
    
    
    ##K to Gamma
    k_x_3 = linspace(((2*pi)/(3*a)), 0, d_x)
    k_y_3 = linspace(((2*pi)/(sqrt(3)*a)), 0, d_y)
    
    E_k_13 = -Gamma + Gamma*sqrt( 4*( cos((1/2)*k_x_3*a)*cos((1/2)*k_x_3*a) + cos((1/4)*k_x_3*a + (sq_3/4)*k_y_3*a)*cos((1/4)*k_x_3*a + (sq_3/4)*k_y_3*a) + cos((1/4)*k_x_3*a - (sq_3/4)*k_y_3*a)*cos((1/4)*k_x_3*a - (sq_3/4)*k_y_3*a)   )  -3    )
    E_k_23 = -Gamma - Gamma*sqrt( 4*( cos((1/2)*k_x_3*a)*cos((1/2)*k_x_3*a) + cos((1/4)*k_x_3*a + (sq_3/4)*k_y_3*a)*cos((1/4)*k_x_3*a + (sq_3/4)*k_y_3*a) + cos((1/4)*k_x_3*a - (sq_3/4)*k_y_3*a)*cos((1/4)*k_x_3*a - (sq_3/4)*k_y_3*a)   )  -3    )
    E_k_33 = 2*Gamma*ones(d_x)
    
    E_k1 = concatenate((E_k_11, E_k_12, E_k_13))
    E_k2 = concatenate((E_k_21, E_k_22, E_k_23))
    E_k3 = concatenate((E_k_31, E_k_32, E_k_33))
    
    plt.figure()
    plt.title(r'$E(\vec{k})$ vs $\vec{k}$ graph')
    plt.plot(E_k3,'k', label = 'Second (flat) conduction band')
    plt.plot(E_k1,'r', label = 'First conduction band')
    plt.plot(E_k2,'b', label = 'Valence band')
    plt.xlabel(r'$\vec{k}$')
    plt.ylabel(r'$E(\vec{k})$ in eV')
    plt.legend(loc='upper right')
    ax = plt.gca() # grabbing the current axis
    plt.xticks([])
    ax.set_xticks([0,(1000-1), (2000-1), (3000-1)]) 
    ax.set_xticklabels([r"$\Gamma$", r"$M$", r"$K$", r"$\Gamma$"]) 
    
Gamma=1.0
a = 1.0
Kagome_Lattice_plot(Gamma, a)
sq_3 = sqrt(3)

# Setting up meshgrid for plotting of 3d graph for energy surface of Kagome Lattice
kvec = linspace(-2*pi,2*pi,(100+1))
kx,ky = meshgrid(kvec,kvec)
E1 = -Gamma + Gamma*sqrt( 4*( cos((1/2)*kx*a)*cos((1/2)*kx*a) + cos((1/4)*kx*a + (sq_3/4)*ky*a)*cos((1/4)*kx*a + (sq_3/4)*ky*a) + cos((1/4)*kx*a - (sq_3/4)*ky*a)*cos((1/4)*kx*a - (sq_3/4)*ky*a)   )  -3    )
E2 = -Gamma - Gamma*sqrt( 4*( cos((1/2)*kx*a)*cos((1/2)*kx*a) + cos((1/4)*kx*a + (sq_3/4)*ky*a)*cos((1/4)*kx*a + (sq_3/4)*ky*a) + cos((1/4)*kx*a - (sq_3/4)*ky*a)*cos((1/4)*kx*a - (sq_3/4)*ky*a)   )  -3    )
# Plot 3d graph for energy surface of Kagome Lattice
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(kx,ky,E1,cmap='viridis',vmin=-E1.max(),vmax=E1.max(),rstride=1,cstride=1)
ax.plot_surface(kx,ky,E2,cmap='viridis',vmin=-E1.max(),vmax=E1.max(),rstride=1,cstride=1)
ax.set_title(r"Energy dispersion graph for Kagome Lattice excluding flat band")
ax.set_xlabel(r"$k_{x}$")
ax.set_ylabel(r"$k_{y}$")
ax.set_zlabel(r"$E(\vec{k})$")