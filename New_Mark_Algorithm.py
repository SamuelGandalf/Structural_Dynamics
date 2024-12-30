import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy import integrate
#from scipy.integrate import cumtrapz
import warnings
warnings.filterwarnings(action="ignore", category=RuntimeWarning)

data_frame = pd.read_excel("./RSN12_KERN.PEL_PEL180.xlsx")

def New_Mark_Algorithm(mass, damping, stiffness, initial_displacement, initial_velocity, B, Y,):
    g = 9.81
    delta_t = 0.01
    data_frame["Force"] = -mass * data_frame["Ground Acceleration (in G)"]*g
    Load_P = data_frame["Force"]
    
    #initialise arrays
    n = data_frame.shape[0]
    Displacement, Displacement[0] = [0]*data_frame.shape[0], initial_displacement
    Velocity, Velocity[0] = [0]*data_frame.shape[0], initial_velocity
    Acceleration = [0]*data_frame.shape[0]
    P_hat = [0]*data_frame.shape[0]

    #initial calculations
    Acceleration[0] = (Load_P[0] - damping*Velocity[0] - stiffness*Velocity[0])/mass
    a1 = (mass/(B*delta_t**2)) + ((Y*damping)/(B*delta_t))
    a2 = (mass/(B*delta_t)) + (((Y/B )-1) * damping)
    a3 = (((1/(2*B)) -1 )*mass) + (((Y/2*B)-1)*damping*delta_t)
    K_hat = stiffness + a1

    #the time step loop// the iterations
    for i in range (n-1):
        P_hat[i+1] = Load_P[i+1] + a1*Displacement[i] + a2*Velocity[i] + a3*Acceleration[i]
        Displacement[i+1] = P_hat[i+1]/K_hat
        Velocity[i+1] = (Y / (B * delta_t)) * (Displacement[i+1] - Displacement[i]) + (1 - (Y / B)) * Velocity[i] + delta_t * (1 - Y / (2 * B)) * Acceleration[i]
        Acceleration[i+1] = ((Displacement[i+1] - Displacement[i]) / (B * delta_t**2)) - (Velocity[i] / (B * delta_t)) - ((1 / (2 * B)) - 1) * Acceleration[i]
    
    return Displacement, Velocity, Acceleration

data_frame['Displacement']= New_Mark_Algorithm(2,10,1000,0,0,0.25,0.5)[0]
data_frame['Velocity']= New_Mark_Algorithm(2,10,1000,0,0,0.25,0.5)[1]
data_frame['Acceleration']= New_Mark_Algorithm(2,10,1000,0,0,0.25,0.5)[2]

fig, axes = plt.subplots(figsize=(10,10),nrows=3,ncols=1) 
axes[0].set_xlabel('time')
axes[0].set_ylim([-3,3])
axes[1].set_xlabel('time')
axes[1].set_ylim([-10,10])
axes[2].set_xlabel('time')
axes[2].set_ylim([-1,1])

axes[0].plot(data_frame['delta t (sec)'],5*data_frame['Acceleration'],label='a(t)',linewidth = 0.9,color='red',linestyle ='-')
axes[1].plot(data_frame['delta t (sec)'],400*data_frame['Velocity'],label='v(t)',linewidth = 0.9,color='blue',linestyle ='-')
axes[2].plot(data_frame['delta t (sec)'],400*data_frame['Displacement'],label='u(t)',linewidth = 0.9, color='black',linestyle ='-')

axes[0].set_title('The structure acceleration-time history')
axes[0].legend(loc='upper left')
axes[1].set_title('The structure velocity-time history')
axes[1].legend(loc='upper left')

axes[2].set_title('The structure displacement-time history')
axes[2].legend(loc='upper left')

plt.show()