import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy import integrate
#from scipy.integrate import cumtrapz
import warnings
warnings.filterwarnings(action="ignore", category=RuntimeWarning)

data_frame = pd.read_excel("./RSN12_KERN.PEL_PEL180.xlsx")

#visualise the earthquake data
# data_frame["Ground_Velocity"] = integrate.cumulative_trapezoid(9.81*data_frame["Ground Acceleration (in G)"], data_frame["delta t (sec)"], initial=0)
# data_frame["Ground_Displacement"] = integrate.cumulative_trapezoid(data_frame["Ground_Velocity"], data_frame["delta t (sec)"], initial=0)

# fig, axes = plt.subplots(figsize=(25,18), nrows=3, ncols=1)
# axes[0].set_xlabel("time(s)")
# axes[0].set_ylim([-1,1])

# axes[1].set_xlabel("time(s)")
# axes[1].set_ylim([-80,80])

# axes[2].set_xlabel("time(s)")
# axes[2].set_ylim([-200,200])

# axes[0].plot(data_frame["delta t (sec)"], 9.81*data_frame["Ground Acceleration (in G)"], label=r'$\dot{\dot{u}}_g(t)$', linewidth=2, color="green", linestyle="-")
# axes[1].plot(data_frame["delta t (sec)"], 1000*data_frame["Ground_Velocity"], label=r'$\dot{u}_g(t)$', linewidth=2, color="blue", linestyle="-")
# axes[2].plot(data_frame["delta t (sec)"], 4000*data_frame["Ground_Displacement"], label=r'u$_g$(t)', linewidth=2, color="green", linestyle="-")

# axes[0].set_title("Ground Acceleration time history")
# axes[0].legend(loc="upper left")
# axes[1].set_title("Ground Velcotiy time history")
# axes[1].legend(loc="upper left")
# axes[2].set_title("Ground Displacement time history")
# axes[2].legend(loc="upper left")

# def Annotate_Peak(ax, x, y, factor, label):
#     peak_index=np.argmax(y)
#     peak_value=y[peak_index]
#     peak_time=x[peak_index]
#     ax.scatter(peak_time, peak_value, color="red", s=100, edgecolor="black", zorder=5)
#     ax.annotate(f"{peak_value: .2f}", xy=(peak_time, peak_value), xytext=(peak_time, peak_value), fontsize=12, color="red")

# #apply the annotations here
# Annotate_Peak(axes[0], data_frame["delta t (sec)"], 9.81*data_frame["Ground Acceleration (in G)"], 1, r'$\dot{\dot{u}}_g(t)$')
# Annotate_Peak(axes[0], data_frame["delta t (sec)"], 1000*data_frame["Ground_Velocity"], 10, r'$\dot{u}_g(t)$')
# Annotate_Peak(axes[0], data_frame["delta t (sec)"], 4000*data_frame["Ground_Displacement"], 100, r'u$_g$(t)')

# plt.show()

#THE CENTRAL DIFFERENCE ALGROITHM
def Central_Difference_Algorithm(mass, damping, stiffness, Initial_Displacement, Initial_Velocity):
    Displacement,  Displacement[0] = [0]*data_frame.shape[0], Initial_Displacement
    Velocity, Velocity[0] = [0]*data_frame.shape[0], Initial_Velocity
    Acceleration = [0]*data_frame.shape[0]
    g=9.81
    data_frame["Force"] = -mass*data_frame["Ground Acceleration (in G)"]*g
    Lood_P = data_frame["Force"]
    delta_T = 0.01
    P_hat = [0]*data_frame.shape[0]

    #Implementing the central difference algorithm
    Displacement[-1] = Displacement[0]-delta_T*Velocity[0] + 0.5*Acceleration[0]*delta_T**2
    K_hat = mass/delta_T**2 + damping/(2*delta_T)
    a = mass/ delta_T**2 - damping/(2*delta_T)
    b = stiffness - 2*mass/delta_T**2

    #the iterations
    for i in range(0, data_frame.shape[0]-1):
        P_hat[i] = Lood_P[i] - a*Displacement[i-1] - b*Displacement[i]
        Displacement[i+1] = P_hat[i]/K_hat
        Velocity[i] = (Displacement[i+1]-Displacement[i-1])/(2*delta_T)
        Acceleration[i] = (Displacement[i+1]- 2*Displacement[i] + Displacement[i-1])/delta_T**2

    return Displacement, Velocity, Acceleration

data_frame["CD_Displacement"] = Central_Difference_Algorithm(2, 10, 1000, 0, 0)[0]
data_frame["CD_Velocity"] = Central_Difference_Algorithm(2, 10, 1000, 0, 0)[1]
data_frame["CD_Acceleration"] = Central_Difference_Algorithm(2, 10, 1000, 0, 0)[2]

fig, axes = plt.subplots(figsize=(15,15), nrows=3, ncols=1)
axes[0].set_xlabel("Time (s)")
axes[0].set_ylim([-10, 10])
axes[1].set_xlabel("Time (s)")
axes[1].set_ylim([-5, 5])
axes[2].set_xlabel("Time (s)")
axes[2].set_ylim([-1, 1])

axes[0].plot(data_frame["delta t (sec)"], 15*data_frame["CD_Acceleration"], label = "a(t)", linewidth = 0.9, color="red", linestyle="-")
axes[1].plot(data_frame["delta t (sec)"], 100*data_frame["CD_Velocity"], label = "v(t)", linewidth = 0.9, color="blue", linestyle="-")
axes[2].plot(data_frame["delta t (sec)"], 400*data_frame["CD_Displacement"], label = "u(t)", linewidth = 0.9, color="black", linestyle="-")

axes[0].set_title("The structure acceleration-time history")
axes[0].legend(loc="upper left")
axes[1].set_title("The structure velocity-time history")
axes[1].legend(loc="upper left")
axes[2].set_title("The structure displacemen-time history")
axes[2].legend(loc="upper left")

plt.show()