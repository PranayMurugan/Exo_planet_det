""" This program takes the data in the form of radial velocity and the time (in Julian days) at which the data was recorded,
and applies Lomb-Scargle Periodogram, and finding the orbital time period,
after finding the time period the radial velocity is folded to get the radial_velocity curve
"""
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
import astropy.units as u
import math
from scipy.signal import find_peaks
from scipy import optimize
import numpy as np
import csv
import pandas as pd

t = []
r = []
dr = []
new_t = []
m_star = 1.11
G = 6.67 * 10**-11

a = []
b = []

K = semi_amplitude
P = period
e = 0.
w = 0.
tau = t[np.argmax(rv)]
vr0 = voffset
guess = (K,P,e,w,tau,vr0)

# Read the data file as a csv

with open('/home/safesoccer/Projects/Planetary/RV_51Pegasi_Data.txt','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=' ')
    for row in plots:
        t.append(float(row[0]))
        r.append(float(row[1]))
        dr.append(float(row[2]))

# simple scatter plot of the raw data

plt.scatter(t,r, label='Radial velocity')
plt.xlabel('Time period (in Julian days)')
plt.ylabel('Radial Velocity')
plt.title('Scatter plot')
plt.show()

# Apply LombScargle Periodogram to find orbital time period

frequency, power = LombScargle(t, r, dr).autopower()
plt.plot(frequency, power)
plt.xlabel('frequency')
plt.ylabel('power/probablity')
plt.title('LombScargle Periodogram')
plt.show()

# Find the peak in the LombScargle probablity distribution

peaks, height = find_peaks(power, height = 0.9)  # height = 0.9 refers the 90% probablity
tp = 1/(frequency[peaks])
tp2 = np.round(tp,2)

for q in t:
    new_t.append(q % tp)

new_t2 = np.array(new_t).ravel().tolist()

# Radial-velocity curve plotted

plt.errorbar(new_t, r, yerr = dr, fmt="o")
plt.xlabel('Time period (in Julian days)')
plt.ylabel('Radial velocity')
plt.title('R_V curve for %s days' % tp2)
plt.show()

list_c = list(zip(new_t2, r))
list_c.sort()
a_sorted = [a for b, a in list_c]

#print(val_Y)
df = pd.DataFrame(list_c)
df.to_csv("/home/safesoccer/Projects/my_array_pandas.csv",header=None, index=None)

with open('/home/safesoccer/Projects/my_array_pandas.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        a.append(float(row[0]))
        b.append(float(row[1]))

def keplerian_fit(t,K,P,e,w,tau,vr0):
  e_anomaly = solve_kepler((t-tau)*2*np.pi/P,e)
  theta = 2*np.arctan2(np.sqrt(1.+e)*np.sin(0.5*e_anomaly),
                       np.sqrt(1.-e)*np.cos(0.5*e_anomaly))
  return K*(np.cos(theta+w)+e*np.cos(w))+vr0

def solve_kepler(M,e):
  eanom = np.zeros(M.shape)
  for i,mi in enumerate(M):
    # do iterative root solve with e=0 giving E=M as guess
    tmp,=fsolve(lambda E: E-e*np.sin(E)-mi,mi)
    eanom[i] = tmp
  return eanom

rvfit = keplerian_fit(t,K,P,e,w,tau,vr0)
chisq = np.sum(((rv-rvfit)/er)**2)
print("Chi-squared of initial guess is %10.5f" % chisq)

popt,pcov = curve_fit(keplerian_fit,t,rv,sigma=er,absolute_sigma=True,p0=guess)

(K,P,e,w,tau,vr0) = popt
print(popt)
rvfit = keplerian_fit(t,K,P,e,w,tau,vr0)
chisq = np.sum(((rv-rvfit)/er)**2)
print("Chi-squared of least-squares fit is %10.5f" % chisq)


plt.errorbar(a, b, yerr = dr, fmt="o", label = 'Scatter plot')
plt.plot(a, rv_func(a, params[0], params[1], params[2]),label='Fitted function')
plt.xlabel('Time period (in Julian days)')
plt.ylabel('Radial Velocity')
plt.title('Scatter plot')
plt.show()



