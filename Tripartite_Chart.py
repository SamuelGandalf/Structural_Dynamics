import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy import integrate
#from scipy.integrate import cumtrapz
import warnings
warnings.filterwarnings(action="ignore", category=RuntimeWarning)

def Tripartite(mass,damping,stiffness):
    g = 9.81
    data_frame = pd.read_excel("./RSN12_KERN.PEL_PEL180.xlsx")
    data_frame["Force"] = -mass*data_frame["Ground Acceleration (in G)"]*g
    
    # Initializing vectors of known quantities
    Displacement = [0]*data_frame.shape[0]
    Velocity = [0]*data_frame.shape[0]
    Acceleration = [0]*data_frame.shape[0]
    Load_P = data_frame["Force"]
    delta_t = 0.01
    P_hat =[0]*data_frame.shape[0]
   
    # Implementing the Central Difference Algorithm
    Displacement[-1]= Displacement[0] - delta_t*Velocity[0]+ 0.5*Acceleration[0]*delta_t**2
    K_hat = mass/delta_t**2 + damping/(2*delta_t)
    a = mass/delta_t**2 - damping/(2*delta_t)
    b = stiffness - 2*mass/delta_t**2

    #the iterations
    for i in range(0, data_frame.shape[0]-1):
        P_hat[i]= data_frame["Force"][i]-a*Displacement[i-1]-b*Displacement[i]
        Displacement[i+1]= P_hat[i]/K_hat
        Velocity[i] = (Displacement[i+1] - Displacement[i-1])/(2*delta_t)
        Acceleration[i] = (Displacement[i+1]-2*Displacement[i]+Displacement[i-1])/delta_t**2

    return max(Displacement),max(Velocity),max(Acceleration)


Damping_Ratios= [0.05,0.10,0.20,0.50]
Period = np.arange(0,0.4,0.001)
stiffness=1000

Displacement = []
Velocity = []
Acceleration = []
Pseudo_Spectra_Velocity=[]
Pseudo_Spectra_Acceleration = []

for zeta in Damping_Ratios:
    for Time_n in Period:
        mass = stiffness*Time_n/4*np.pi**2
        damping = 2*zeta*np.sqrt(stiffness*mass)
        #print(mass)
        Peak_Displacement = Tripartite(mass,damping,stiffness)[0]*10000
        Peak_Velocity =  Tripartite(mass,damping,stiffness)[1]*10000
        Peak_Acceleration =  Tripartite(mass,damping,stiffness)[2]*10000

        Pseudo_Spectra_Acceleration.append(Peak_Displacement*(2*np.pi)/Time_n)
        Pseudo_Spectra_Velocity.append(Peak_Displacement*((2*np.pi)/Time_n)**2)

        Displacement.append(Peak_Displacement)
        Velocity.append(Peak_Velocity)
        Acceleration.append(Peak_Acceleration)

# fig, axes = plt.subplots(figsize=(15,10),nrows=2,ncols=2) 

# axes[0][0].plot(period, Acceleration[0:400],label='5%',linewidth = 1,color='green',linestyle ='-')
# axes[0][0].plot(period, Acceleraion[400:800],label='10%',linewidth = 1,color='blue',linestyle ='-')
# axes[0][0].plot(period, Acceleration[800:1200],label='20%',linewidth = 2,color='red',linestyle ='-')
# axes[0][0].plot(period, Acceleration[1200:1600],label='50%',linewidth = 2,color='yellow',linestyle ='-')
# axes[0][0].legend(loc='lower right')
# axes[0][0].set_title('The structure acceleration response spectra')

# axes[0][1].plot(Period, Pseudo_Spectra_Acceleration[0:400],label='5%',linewidth = 1,color='green',linestyle ='-')
# axes[0][1].plot(Period, Pseudo_Spectra_Acceleration[400:800],label='10%',linewidth = 1,color='blue',linestyle ='-')
# axes[0][1].plot(Period, Pseudo_Spectra_Acceleration[800:1200],label='20%',linewidth = 1,color='red',linestyle ='-')
# axes[0][1].plot(Period, Pseudo_Spectra_Acceleration[1200:1600],label='50%',linewidth = 1,color='yellow',linestyle ='-')
# axes[0][1].legend(loc='upper right')
# axes[0][1].set_title('The Pseudo-Acceleration Response Spectrum')

# axes[1][0].plot(Period, Velocity[0:400],label='5%',linewidth = 1,color='green',linestyle ='-')
# axes[1][0].plot(Period, Velocity[400:800],label='10%',linewidth =1,color='blue',linestyle ='-')
# axes[1][0].plot(Period, Velocity[800:1200],label='20%',linewidth =1,color='red',linestyle ='-')
# axes[1][0].plot(Period, Velocity[1200:1600],label='50%',linewidth =1,color='yellow',linestyle ='-')
# axes[1][0].legend(loc='lower right')
# axes[1][0].set_title('The structure velocity response spectra')

# axes[1][1].plot(Period, Displacement[0:400],label='5%',linewidth= 1, color='green',linestyle ='-')
# axes[1][1].plot(Period, Displacement[400:800],label='10%',linewidth =1, color='blue',linestyle ='-')
# axes[1][1].plot(Period, Displacement[800:1200],label='20%',linewidth = 1, color='red',linestyle ='-')
# axes[1][1].plot(Period, Displacement[1200:1600],label='50%',linewidth = 1, color='yellow',linestyle ='-')
# axes[1][1].legend(loc='lower right')
# axes[1][1].set_title('The structure displacment response spectra')
# Tripartite(2, 10, 100)

fig, axes = plt.subplots(figsize=(25,18), nrows=2, ncols=2)

def init():
    for ax in axes.flatten():
        ax.clear()
    return axes.flatten()

def update(frame):
    # Clear previous data
    for ax in axes.flatten():
        ax.clear()

    # Update plots
    axes[0][0].plot(Period[:frame], Acceleration[:frame], label='20%', linewidth=2, color='green', linestyle='-')
    axes[0][0].plot(Period[:frame], Acceleration[400:frame+400], label='30%', linewidth=2, color='blue', linestyle='-')
    axes[0][0].plot(Period[:frame], AccelerationA[800:frame+800], label='40%', linewidth=2, color='red', linestyle='-')
    axes[0][0].plot(Period[:frame], Acceleration[1200:frame+1200], label='50%', linewidth=2, color='yellow', linestyle='-')
    axes[0][0].legend(loc='lower right')
    axes[0][0].set_title('System Acceleration Response Spectra')

    axes[1][0].plot(Period[:frame], Velocity[:frame], label='20%', linewidth=2, color='green', linestyle='-')
    axes[1][0].plot(Period[:frame], Velocity[400:frame+400], label='30%', linewidth=2, color='blue', linestyle='-')
    axes[1][0].plot(Period[:frame], Velocity[800:frame+800], label='40%', linewidth=2, color='red', linestyle='-')
    axes[1][0].plot(Period[:frame], Velocity[1200:frame+1200], label='50%', linewidth=2, color='yellow', linestyle='-')
    axes[1][0].legend(loc='lower right')
    axes[1][0].set_title('System Velocity Response Spectra')

    axes[1][1].plot(Period[:frame], Displacement[:frame], label='20%', linewidth=2, color='green', linestyle='-')
    axes[1][1].plot(Period[:frame], Displacement[400:frame+400], label='30%', linewidth=2, color='blue', linestyle='-')
    axes[1][1].plot(Period[:frame], Displacement[800:frame+800], label='40%', linewidth=2, color='red', linestyle='-')
    axes[1][1].plot(Period[:frame], Displacement[1200:frame+1200], label='50%', linewidth=2, color='yellow', linestyle='-')
    axes[1][1].legend(loc='lower right')
    axes[1][1].set_title('System Displacement Response Spectra')

    axes[0][1].plot(Period[:frame], Pseudo_Spectra_Acceleration[:frame], label='20%', linewidth=2, color='green', linestyle='-')
    axes[0][1].plot(Period[:frame], Pseudo_Spectra_Acceleration[400:frame+400], label='30%', linewidth=2, color='blue', linestyle='-')
    axes[0][1].plot(Period[:frame], Pseudo_Spectra_Acceleration[800:frame+800], label='40%', linewidth=2, color='red', linestyle='-')
    axes[0][1].plot(Period[:frame], Pseudo_Spectra_Acceleration[1200:frame+1200], label='50%', linewidth=2, color='yellow', linestyle='-')
    axes[0][1].legend(loc='upper right')
    axes[0][1].set_title('Pseudo-Acceleration Response Spectrum')

    return axes.flatten()

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=len(period), init_func=init, blit=False, interval=100)
ani.save('ResponseSpectr1.gif')

plt.tight_layout()
plt.show()
      
