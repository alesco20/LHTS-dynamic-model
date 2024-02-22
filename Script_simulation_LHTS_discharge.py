# -*- coding: utf-8 -*-
"""
@author: alessandro.colangelo

Purpose: Dynamic simulation of Longitudinal LHTS discharging process

The LHTS unit is made of 96 pipes immersed in a PCM. Each pipe is made of Copper and surrounded by Aluminum fins 
                                
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.sparse
import scipy.sparse.linalg
import time


#Default plots parameters setting
plt.style.use('seaborn')


#SIMULATION DATA
#Pipe geometry
di = 2*7.45e-3         #[m] pipe inner diameter
L_tube = 1.29          #[m] pipe length
A_tube = np.pi*di**2/4  #[m2] pipe cross section
n_pipes = 96;           #number of pipes contained in the LHTS used for the discharging phase

#Water properties @65°C
rho = 980.6            #[kg/m3]
cp = 4184              #[J/(kg*K)]
k = 0.6455             #[W/(m*K)]
mu = 433.4e-6          #[Pa/s] 
alpha = k/(rho*cp)     #[m2/s]

#Operational (initial) values
T_in = 50             #[°C] HTF inlet temperature
T_initial = 80        #[°C] Average PCM Initial temperature 
G = 0.006             #[kg/s] HTF mass flow rate in a single pipe
u = G/(rho*A_tube)    #[m/s] HTF velocity in each pipe

#Adimensional numbers
Pr = mu*cp/k
Re = G/A_tube*di/mu
Nu_t = 0.023*Re**(4/5)*Pr**(0.4)           #valid for Re>10000


#LHTS features:
#PCM:    
rho_PCM = 858          #[kg/m3] PCM density
cp_s = 1800            #[J/(kg*K)] 
cp_l = 2000            #[J/(kg*K)]
T_pc = 71              #[°C] Mean phase change temperature
dT_pc = 3
T_s = T_pc - dT_pc/2
T_l = T_pc + dT_pc/2
Dh_pc = 224e3            #[J/kg] Latent heat
k_PCM = 0.28             #[W/mK] PCM thermal conductivity
#HCM:
rho_Cu = 8978            #[kg/m3]
rho_Al = 2701            #[kg/m3]
cp_Cu = 381              #[J/(kg*K)]
cp_Al = 871              #[J/(kg*K)]   
#Geometry:
V_pcm_tot = 0.573448     #[m3]
V_Cu_tot = 4.56648e-3    #[m3]
V_Al_tot = 0.0443597     #[m3]
V_pcm = V_pcm_tot/L_tube/n_pipes         #[m3/m] Vol. of PCM around each meter of pipe
V_Cu = V_Cu_tot/L_tube/n_pipes           #[m3/m] Vol. of Copper (Cu) around each meter of pipe
V_Al = V_Al_tot/L_tube/n_pipes           #[m3/m] Vol. of Aluminum (Al) around each meter of pipe

A_exchange = 82.536                       #[m2] Total exchange surface (i.e. overall contact surface between PCM and fins)
l_c = V_pcm_tot/A_exchange                #[m] Characteristic length of the specific LHTS design 

M_PCM = rho_PCM*V_pcm               #[kg/m] Mass of PCM around each meter of pipe
M_HCM = rho_Cu*V_Cu + rho_Al*V_Al   #[kg/m] Mass of high conducting material around each meter of pipe

T_ref = 20                          #[°C] Chosen reference temperature for calculation of LHTS energy content


#Energy content of PCM-HCM assembly "around" one tube (as a function of temperature)  #[J/m]
f_En = lambda T: (rho_Cu*V_Cu*cp_Cu*(T-T_ref) + rho_Al*V_Al*cp_Al*(T-T_ref) + rho_PCM*V_pcm*(cp_s*(T-T_ref)*(T<=T_s) + 
    (cp_s*(T-T_ref)+Dh_pc*(T-T_s)/dT_pc)*((T>T_s) & (T<T_pc)) + 
    (cp_s*(T_pc-T_ref)+cp_l*(T-T_pc)+Dh_pc*(T-T_s)/dT_pc)*((T>=T_pc) & (T<T_l)) + 
    (cp_s*(T_pc-T_ref)+cp_l*(T-T_pc)+Dh_pc)*(T>=T_l)))  

#Characteristic curve of thermal power (physical dimensions depend on how delta energy is expressed)
beta = M_PCM/(M_PCM + M_HCM)   #[-]
qq = lambda tau,soc,de: (soc*beta/tau*(np.log(1/soc + 1e-10))**((beta-1)/beta))*de

#Function for calculating the time scale of the heat transfer (tau)
f_tau = lambda soc0,T_w: 0.632*(rho_PCM*(l_c*soc0)**2*
    ((T_w<=T_s)*(cp_s*(T_s-T_w) + cp_l*(T_initial-T_l) + Dh_pc) +
     ((T_w>T_s) & (T_w<T_l))*(cp_l*(T_initial-T_l) + Dh_pc*(T_l-T_w)/dT_pc) +
     (T_w>=T_l)*(cp_l*(T_initial-T_w))))/(k_PCM*(T_initial-T_w))


#Discretization and simulation parameters
t_end = 3600                                    #[s]
dt = 30                                         #[s]
Nt = int(t_end/dt)                              #number of timesteps
t = np.linspace(0,t_end,Nt)                     #[s] time vector
dz = 0.01                                       #[m]
n_centroids = int(L_tube/dz)                    #number of centroids with finite volume discretization
z = np.linspace(dz/2,L_tube-dz/2,n_centroids)   #[m] space vector
z = np.append(np.insert(z,0,0),[L_tube])        #[m] If the space vector is built in this way it is consistent with a finite volume discretization (centroids + 2 BCs)
Nz = np.size(z)                                 #Number of centroids + 2 boundaries (inlet and outlet)
A_lat = np.pi*di*dz                             #[m2] Lateral area of each small cylinder in which the whole pipe is divided

Nu_l = np.mean(3.66 + (0.0668*di/z[1:-1]*Re*Pr)/(1 + 0.04*(di/z[1:-1]*Re*Pr)**(2/3)))  #This formula considers the thermal entry length in laminar flow (see pag. 513 Incropera - "Fundamentals of heat and mass transfer")
Nu = (Re<=2300)*Nu_l + (Re>=2300)*Nu_t   #Nusselt number is evaluated choosing between the laminar and the turbulent value depending on Re number (the "region" between laminar and turbulent is approximated with the turbulent value, which in principle should be applicable only fo Re>10000)
h_conv = k*Nu/di                         #[W/(m2*K)]


#Initial values
T_old = T_initial*np.ones(Nz)                   #[°C] Initialization of axial HTF temperature
SOC_0 = np.ones(Nz-2)                           #[-] Initial LHTS state of charge 
SOC = SOC_0                                     #[-] Initialization of LHTS axial states of charge
E_0 = f_En(T_initial)*dz*np.ones(Nz-2)          #[J] Initial energy content in each small volume of PCM-HCM assembly
E_t = E_0                                       #[J] Initialization of the energy content in each small volume of PCM-HCM assembly

#Construction of problem matrix (A) with BCs
a_low = -u*dt/dz*np.ones(Nz)
a_diag = (1 + u*dt/dz)*np.ones(Nz)
a_diag[0] = 1
a_diag[-1] = 1     #[-1] indicates last element of a_diag
a_low[-2] = -1     #The last element in the matrix is [-2], because the last element of the array a_low will be cut in the matrix construction 
A = sp.sparse.spdiags([a_low, a_diag], [-1,0], Nz, Nz).tocsc()    #sparse matrix allows to reduce computational time 

#Right-hand side vector with BCs
b = np.zeros(Nz)
b[0] = T_in
b[-1] = 0

#Inizialization of output variables
T_out = T_in*np.ones(Nt)             #[°C] HTF outlet temperature
T_axial = T_in*np.ones((Nz,Nt))      #[°C] HTF axial temperature
q_axial = np.zeros((Nz-2,Nt))        #[W/m2] axial heat flux
SOC_axial = np.ones((Nz-2,Nt))       #[-] axial state of charge of each pipe section
T_wall = T_in*np.ones((Nz-2,Nt))     #[°C] Wall temperature at the contact between HTF and PCM-HCM assembly
G_t = G*np.ones(Nt)                  #[kg/s] mass flow rate in a single tube in time G(t)
Tin_t = T_in*np.ones(Nt)             #[°C] inlet temperature in time T_in(t)
q_tot = np.zeros(Nt)                 #[W/tube] Total heat transfer rate from each LHTS tube in time 
state_charge = np.ones(Nt)           #[-] Average storage state of charge

tic = time.time()

for j in range(1,Nt):
    tau_ref = f_tau(SOC_0,T_ref)                              #reference time constant for linearization of the problem (could also be calculated outside the "for loop")
    DE_ref = (f_En(T_initial)-f_En(T_ref))*dz*np.ones(Nz-2)   #reference delta energy for linearization of the problem (could also be calculated outside the "for loop")
    q_ref = qq(tau_ref,SOC,DE_ref)                            #reference thermal power value for linearization of the problem (must be calculated inside the "for loop" because it depends on the state of charge)
    
    Tw = (q_ref - T_ref*q_ref/(T_ref-T_initial) + h_conv*A_lat*T_old[1:-1])/(h_conv*A_lat - q_ref/(T_ref-T_initial))  #[°C] Temperature at the contact wall between the PCM and the pipe
    q_wall = h_conv*A_lat*(Tw - T_old[1:-1])                 #[W] Thermal power exchanged at the contact wall between the PCM and the pipe
    
    b[1:-1] = T_old[1:-1] + dt/(rho*cp*A_tube*dz)*q_wall     #update of right-end side vector inner elements
    
    T_z = sp.sparse.linalg.spsolve(A,b)                      #[°C] calculation of axial HTF temperature T(z)
    
    #Update variables
    T_old = T_z
    T_in = 50                           #[°C] here the value of the HTF inlet temperature can be updated at each time-step if desired
    E_t = E_t - q_wall*dt               #[J] new energy content of each small volume around a pipe
    E_end = f_En(T_in)*dz               #[J] update of the "final" LHTS energy content (i.e. the ideal LHTS energy content at the end of the discharging process) depending on the new value of the HTF inlet temperature
    SOC = (E_t - E_end)/(E_0 - E_end)   #Update of the state of charge of each small volume around the pipe
    G = 0.006                           #[kg/s] ere the value of the HTF mass flow rate in a single tube can be updated at each time-step if desired
    u = G/(rho*A_tube)                  #[m/s] update of the HTF velocity
    
    #Adimensional numbers and h_conv update
    Re = G/A_tube*di/mu
    Nu_t = 0.023*Re**(4/5)*Pr**(0.4)           
    Nu_l = np.mean(3.66 + (0.0668*di/z[1:-1]*Re*Pr)/(1 + 0.04*(di/z[1:-1]*Re*Pr)**(2/3)))  
    Nu = (Re<=2300)*Nu_l + (Re>=2300)*Nu_t
    h_conv = k*Nu/di

    #Update problem matrix and right-end side vector
    b[0] = T_in
    a_low = -u*dt/dz*np.ones(Nz)
    a_diag = (1 + u*dt/dz)*np.ones(Nz)
    a_diag[0] = 1
    a_diag[-1] = 1     #[-1] indicates last element of a_diag
    a_low[-2] = -1     #The last element in the matrix is [-2], because the last element of the array a_low will be cut in the matrix construction 
    A = sp.sparse.spdiags([a_low, a_diag], [-1,0], Nz, Nz).tocsc()    #sparse matrix allows to reduce computational time 

    
    #Output variables
    T_out[j] = T_z[-1]                  
    T_axial[:,j] = T_z
    q_axial[:,j] = q_wall/(np.pi*di*dz)  #[W/m2] axial heat flux exchanged at the contact wall between the pipe and the PCM
    SOC_axial[:,j] = SOC
    T_wall[:,j] = Tw
    G_t[j] = G        
    Tin_t[j] = T_in
    q_tot[j] = sum(q_wall)              
    state_charge[j] = np.mean(SOC)
    
    #Print number of timestep
    print(j)
    

toc = time.time()
print(toc-tic, 'sec elapsed')




#Plots: T_in, T_out
fig,ax = plt.subplots()
ax.plot(t/60, T_out, linewidth=4, zorder = 2, label='$T_{out, sim}$')
ax.plot(t/60, Tin_t, marker = 'o', linestyle ='-', label='$T_{in}$')
fig.legend(loc='upper right',fontsize = 16, bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
ax.set_xlabel('Time [min]')
ax.set_ylabel('Temperature [°C]')
# fig.savefig('Fig_Long_Tin_Tout.png', dpi=fig.dpi)
plt.show()

#Plots: LHTS Thermal power
fig,ax = plt.subplots()
ax.plot(t/60, q_tot*n_pipes/1000, linewidth=3, zorder = 2, label='$q_{LHTS, sim}$')
fig.legend(loc='upper right',fontsize = 16, bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
ax.set_xlabel('Time [min]')
ax.set_ylabel('Thermal power [kW]')
# fig.savefig('Fig_Long_qLHTS.png', dpi=fig.dpi)
plt.show()

#Plots: State of charge
fig,ax = plt.subplots()
ax.plot(t/60, state_charge, linewidth=3, zorder = 2, label='$SOC$')
fig.legend(loc='upper right',fontsize = 16, bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
ax.set_xlabel('Time [min]')
ax.set_ylabel('SOC [-]')
# fig.savefig('Fig_Long_Edischarged.png', dpi=fig.dpi)
plt.show()

