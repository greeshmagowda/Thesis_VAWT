import numpy as np
import math as mt
import scipy.io
import scipy.integrate as integrate
import matplotlib.pyplot as plt
N = 36
f = 0.98
R = 1
solidity = 0.1
B = 2
c = solidity * 2 * R / B
rps = 3.0 / 60
TSR = 4
Vinf = [5]  # m/s
rho = 1.225
Pitch = 0
nu = 1.460e-5
alpha,Cl,Cd = np.loadtxt('airfoildatanew.txt', skiprows=1,usecols=(0, 1, 2), 
                         unpack=True)
corr1 = np.polyfit(alpha,Cl,5)
corr2 = np.polyfit(alpha,Cd,5)                         
#%%
# Rotor Sections
theta1 = []
deltheta = 2 * mt.pi / N
theta = np.arange(deltheta/2, 2*mt.pi, deltheta)
for val in theta: 
    theta1.append(mt.degrees(val))
x = -f * np.sin(theta)
y = f * np.cos(theta)  
#Influence Coefficients
Rwx = np.zeros((N,N))
Rwy = np.zeros((N,N))
rwx = np.zeros((N,N))
rwy = np.zeros((N,N))
# Functions to calculate the value of integral
def intX(deltheta1,x,y):
    funX = (-((x + mt.sin(deltheta1))*mt.sin(deltheta1)) + ((y - mt.cos(deltheta1)) * mt.cos(deltheta1))) / ((x+mt.sin(deltheta1))**2 + (y-mt.cos(deltheta1))**2)
    return funX
def intY(deltheta1,x,y):
    funY = (-((x + mt.sin(deltheta1))*mt.cos(deltheta1)) - ((y - mt.cos(deltheta1)) * mt.sin(deltheta1))) / ((x+mt.sin(deltheta1))**2 + (y-mt.cos(deltheta1))**2)
    return funY

for j in range(N):
    deltheta1 = np.arange((theta[j]- (0.5*deltheta)),(theta[j]+ (0.5*deltheta)),(0.1*deltheta))
    for i in range(N): theta = np.arange(1,5,0.5)
        rwx[j][i] = intX (deltheta1,x[i],y[i])
        rwy[j][i] = intY (deltheta1,x[i],y[i])
    #scipy.io.savemat('Rwx.mat', mdict={'Rwx': Rwx})
    #scipy.io.savemat('Rwy.mat', mdict={'Rwy': Rwy})
#%%
#Initial guess for induced velocity
Vinf = [5]
Wx = np.zeros(N)  
Wy = np.zeros(N) 
Wxn = np.zeros(N)
Wyn = np.zeros(N)
Qn = np.zeros(N)
Qt = np.zeros(N)
Fn = np.zeros(N)
Ft = np.zeros(N)
ct = np.zeros(N)
cp = np.zeros(N)
alfa1 = np.zeros(N)

for i in range(len(Vinf)):
    omega = TSR * Vinf[j] / R
    finished = False
    count = 0 
    while not (finished):
        print(count)
        print(Wx) 
        print(Wy)
        for j in range(N):
          Vx = Vinf[i] * (1 + Wx[j]) - (omega * R * mt.sin(theta[j]))
          #print(Vx)
          Vy = (Vinf[i] * Wy[j])  + (omega * R* mt.cos(theta[j]))
          #print(Vy)
          Vn = Vx * mt.cos(theta[j]) + Vy * mt.sin(theta[j])
          Vt = - (Vx * mt.sin(theta[j])) + (Vy * mt.cos(theta[j]))
          Vr = mt.sqrt((Vn ** 2) + (Vt ** 2))
          alfa = mt.atan(Vn/Vt)
          alfa1[j] = mt.degrees(alfa)
          #print("Alfa", mt.degrees(alfa))
          CL = np.polyval (corr1,mt.degrees(alfa))
          CD = np.polyval (corr2,mt.degrees(alfa))
          Cn = CD * mt.sin(alfa) + CL * mt.cos(alfa)
          Ct = CL * mt.sin(alfa) - CD * mt.cos(alfa)
          Fn[j] = 0.5 * rho * c * (Vr ** 2)  * Cn
          Ft[j] = 0.5 * rho * c * (Vt ** 2)  * Ct
          Qn[j] = B * ((Fn[j] * mt.cos(Pitch)) - (Ft[j] * mt.sin(Pitch))) / (2* mt.pi * R)
          Qt[j] = B * ((Ft[j] * mt.cos(Pitch)) + (Fn[j] * mt.sin(Pitch))) / (2* mt.pi * R)
        for k in range(N):
            ct[k] = Qn[k] * mt.sin(theta[k]) - Qt[k] * mt.cos(theta[k])
            cp[k] = ((B * Ft[k] * mt.cos(Pitch) + Fn[k] * mt.sin(Pitch)) * omega)
            sumx = 0
            sumy = 0
            #Y = y
            for l in range(N):
                sumx = sumx + (Rwx[l][k] * Qn[l])
                sumy = sumy + (Rwy[l][k] * Qn[l])
            Wxn[k] = (1 / (2 * mt.pi)) * sumx - Qn[k] 
            Wyn[k] = (1 / (2 * mt.pi)) * sumy
        Wx = Wxn
        Wy = Wyn 
        CT = np.trapz(ct,theta)
        CP = (1/( 2 * mt.pi)) * np.trapz(cp,theta) /(rho * Vinf[i] ** 3)
          
        count = count + 1
        if count == 1:
            finished = True

plt.figure(1)         
plt.plot(theta1,alfa1,'bo-')
plt.xlabel(r'$\Theta$')
plt.ylabel(r'$\alpha$')
plt.title('Azimuthal variation of Angle of attack')
plt.figure(2)
plt.plot(theta1,Wx,'bo')
plt.xlabel(r'$\theta$')
plt.ylabel(r'$W_x$')
plt.title('Variation of Induced Velocity ('r'$W_x$)')
plt.figure(3)
plt.plot(theta1,Wy,'bo')
plt.xlabel(r'$\theta$')
plt.ylabel(r'$W_y$')
plt.title('Variation of Induced Velocity ('r'$W_y$)')           