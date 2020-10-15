# Advanced Material Balance Calculation

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimization
import math
from waterinflux import water_influx
from zfactor import z_factor

# Import properties into the model 
prop = input("Enter the file name for properties: ")
inp_prop = open(prop, "r")
line = inp_prop.readline()
line = line.rstrip()
tokens = line.split()
cf = float(tokens[0])
cw = float(tokens[1])
sgi = float(tokens[2])
swi = float(tokens[3])
phi = float(tokens[4])
ct = float(tokens[5])
h = float(tokens[6])
ra = float(tokens[7])
rr = float(tokens[8])
theta = float(tokens[9])
k = float(tokens[10])
vis = float(tokens[11])
aq = tokens[12]
bnd = tokens[13]
Pinit = float(tokens[14])
Temp = float(tokens[15])
Bw = float(tokens[16])
haq = float(tokens[17]) 
inp_prop.close()

# Import pressure vs time data into the model 
pres = input("Enter file name for pressure: ")
inp = open(pres, "r")
line1 = inp.readline()
line1 = line1.rstrip()
tokens1 = line1.split()
m = int(tokens1[0])
t = np.zeros((m), int)
P = np.zeros((m), float)
for i in range (m):
	line2 = inp.readline()
	line2 = line2.rstrip()
	tokens2 = line2.split()
	t[i] = int(tokens2[0])
	P[i] = float(tokens2[1])
inp.close()

T = np.ones((m))*Temp
We1, We2 = water_influx(m, phi, ct, h, ra, rr, theta, k, vis, aq, bnd, Pinit, P, t)
plt.figure()
plt.plot(t, We1, "ro", label ="$VanEverdingen-Hurst$")
plt.plot(t, We2, "g--", label= "$Carter-Tracy$")
#plt.plot(t, We3, "b--", label = "$Fetkovich$")
plt.xlabel("$Time,$ $Days$")
plt.ylabel("$Cummulative$ $Water$ $Influx$, $MMBBL$")
plt.title("$Water$ $Influx$")
plt.legend()
plt.grid()

plt.show()
choice = input("Choose which method you want to retain. For VEH type V and for Cartar-Tracy type C: ")
if choice == "V":
	We = We1
#elif choice == "F":
#	We = We3
else:
	We = We2
    
z, Bg = z_factor(P) 
Bgi = Bg[0]
zi = z[0]
# Import production data into the model
prod = input("Enter file name for production: ")
intake = open(prod,"r")
line3 = intake.readline()
line3 = line3.rstrip()
tokens3 = line3.split()
n = int(tokens3[0])
days = np.zeros((n), int)
cum_gas = np.zeros((n), float)
cum_water = np.zeros((n), float)
cum_cond = np.zeros((n), float)
for i in range (n):
	line4 = intake.readline()
	line4 = line4.rstrip()
	tokens4 = line4.split()
	days[i] = int(tokens4[0])
	cum_gas[i] = float(tokens4[1])
	cum_water[i] = float(tokens4[2])
	cum_cond[i] = float(tokens4[3])
intake.close()

G_p = np.zeros((m), float)
Gp = np.zeros((m), float)
Wp = np.zeros((m), float)
Cp = np.zeros((m), float)

for i in range (m):
	temp = np.zeros((n))
	for j in range (n):
		temp[j] = abs(t[i] - days[j])
	x = np.where(temp == temp.min())
	G_p[i] = cum_gas[x[0][0]]
	Wp[i] = cum_water[x[0][0]]
	Cp[i] = cum_cond[x[0][0]]
	Gp[i] = G_p[i] + (Cp[i]*0.0011)

plt.figure()
plt.plot(P,z, "b.")
plt.title("$z$ $Factor$ $(Peng-Robinson$ $EOS)$")
plt.xlabel("$Pressure,$ $Psia$")
plt.ylabel("$z$")
plt.grid()

plt.show()
We = We*5.615

def best_fit(m, x, y):
	theta = np.array([0,1])
	sigma = np.ones((m), float)        
	def linear(x, a, b):
		return (a+b*x)
	a, cov = optimization.curve_fit(linear, x, y, theta, sigma)
	GIIP = (-(a[0]/a[1]))
	return (GIIP, a[0], a[1])   

def best_fit_mod(m, x, y):
	a = np.zeros((m))
	for i in range(m):
		if y[i] == 0:
			y[i] = 0.0001
		elif y[i] == 1:
			y[i] = 0.99999
		else:
			y[i] = y[i]
	a[:] = x[:]*((1+np.log(y[:]))/np.log(y[:]))
	A = min(a)
	return A 
    
def best_fit_expo(m, x, y, G):
	theta = np.array([G])
	sigma = np.ones((m), float)        
	def expo(x, b):
		return (-x/(b-x))
	a, cov = optimization.curve_fit(expo, x, y, theta, sigma)
	return (a[0]) 
    
def PZ(P, z, Gp):    
	pz = np.zeros((m), float)
	pz[:] = P[:]/z[:]
	plt.figure()
	plt.plot(Gp, pz, "bo", label = "$P/z$")
	G, intercept, slope = best_fit(m, Gp, pz)
	plt.plot(Gp,(intercept +(Gp*slope)), "r--", label = "$Best$ $Fit$")
	plt.xlabel("$Cummulative$ $Production,$ $MMSCf$")
	plt.ylabel("$P/z,$ $Psia$")
	plt.title("$Classical$ $Gas$ $Material$ $Balance$ $(P/z$ $Method)$")
	plt.legend()
	plt.grid()
    
	plt.show()
	return G

def PZ_advanced(We, Wp, Bg, Bgi, Bw, Sgi, swi, Pinit, P, z, G):
	cwip = np.zeros((m), float)
	cwip[:] = (We[:] - (Wp[:]*Bw))*sgi/(G*Bgi) 
	cep = np.zeros((m), float)
	cep[:] = (cf + (cw*swi))*(Pinit - P[:])
	exp_init = np.zeros((m), float)
	exp = np.zeros((m), float)
	exp_init[:] = (sgi - cep[:])
	exp[:] = (sgi - cwip[:] - cep[:])
	zdstar_init = np.zeros((m))
	zdstar = np.zeros((m))
	zdstar_init[:] = P[:]/((1/sgi)*(P[:]/z[:])*(exp_init[:]))
	zdstar[:] = P[:]/((1/sgi)*(P[:]/z[:])*(exp[:]))
	plt.figure()
	plt.plot(Gp, (P/zdstar_init), "bo", label = "$P/z**$")
	G_mod1, intercept_mod1, slope_mod1 = best_fit(m, Gp, (P/zdstar_init))
	G_mod2, intercept_mod2, slope_mod2 = best_fit(m, Gp,(P/zdstar))
	plt.plot(Gp,(intercept_mod1 +(Gp*slope_mod1)), "r--", label = "$Best$ $Fit$")
	plt.xlabel("$Cummulative$ $Production,$ $MMSCf$")
	plt.ylabel("$P/z**,$ $Psia$")
	plt.title("$Advanced$ $Gas$ $Material$ $Balance$ $(Without$ $considering$ $Water-Influx)$")
	plt.legend()
	plt.grid()

	plt.show()
	plt.figure()
	plt.plot(Gp, (P/zdstar), "bo")
	plt.plot(Gp,(intercept_mod2 +(Gp*slope_mod2)), "r--", label = "$Best$ $Fit$")
	plt.xlabel("$Cummulative$ $Production,$ $MMSCf$")
	plt.ylabel("$P/z**,$ $Psia$")
	plt.title("$Advanced$ $Gas$ $Material$ $Balance$ $(Considering$ $Water-Influx)$")
	plt.legend()
	plt.grid()

	plt.show()
	temp = np.zeros((m))
	temp[:] = ((intercept_mod1 +(Gp[:]*slope_mod1)))
	We_mod1= np.zeros((m))
	cwip_mod1 = np.zeros((m), float)
	cwip_mod1[:] = -temp[:]*sgi*(z[:]/P[:]) - cep[:] + sgi 
	#cwip_mod1[:] = (We_mod1[:] - (Wp[:]*Bw))*sgi/(G_mod1*Bgi) 
	exp_mod1 = np.zeros((m), float)
	exp_mod1[:] = (sgi - cwip_mod1[:] - cep[:])
	zdstar_mod1 = np.zeros((m))
	zdstar_mod1[:] = P[:]/((1/sgi)*(P[:]/z[:])*(exp_mod1[:]))
	G_final, intercept_final, slope_final = best_fit(m, Gp, (P/zdstar_mod1))
	We_mod1 = (((G_final*Bgi)/sgi)*cwip_mod1) +(Wp*Bw)
	We_mod1[0] = 0
	plt.figure()
	plt.plot(Gp, (P/zdstar_mod1), "bo", label = "$P/z**$")
	plt.plot(Gp,(intercept_final +(Gp*slope_final)), "r--", label = "$Best$ $Fit$")
	plt.xlabel("$Cummulative$ $Production,$ $MMSCf$")
	plt.ylabel("$P/z**,$ $Psia$")
	plt.title("$Advanced$ $Gas$ $Material$ $Balance$ $(Balancing$ $Water-Influx)$")
	plt.legend()
	plt.grid()

	plt.show()
	return (G_final, We_mod1)

def HavlenaOdeh(m, Wp, Bw, We, cw, swi, P, Pinit, Gp, cf, z, T):
	Bg =  0.02827*(z*T/P)
	Bgi = Bg[0]
	Eg =  np.zeros((m)) 
	Efw =  np.zeros((m))
	F =  np.zeros((m))
	Eg[:] = Bg[:] - Bgi
	Efw[:] = Bgi*(((Pinit - P[:])*(cf + (cw*swi)))/(1 - swi))
	for i in range(m):
		if (We[i]*Bw) > ((Gp[i]*Bg[i]) + (Wp[i]*Bw)):
			We[i] = 0
		else:
			We[i] = We[i]
	F[:] = (Gp[:]*Bg[:]) + (Wp[:]*Bw) - (We[:]*Bw)
	junk, intercept, slope = best_fit(m, (Eg+ Efw), F)
	plt.figure()
	plt.plot(Gp, (F/(Eg + Efw)), "bo", label= "$Havlena-Odeh$")
	plt.xlabel("$Cummulative$ $Production,$ $MMSCf$")
	plt.ylabel("$F/(Eg+Efw)$")
	plt.title("$Havlena-Odeh$ $Gas$ $Material$ $Balance$")
	plt.ylim([0, max(F/(Eg + Efw))])
	plt.xlim([0, max(Gp)])
	plt.legend()
	plt.grid()

	plt.show()
	plt.figure()
	plt.plot((Eg + Efw), F, "bo", label= "$Havlena-Odeh$")
	plt.plot((Eg+Efw), (intercept + (slope*(Eg+Efw))),"r--", label="$interpolated$")
	plt.xlabel("$Eg+Efw,$ $MMSCf/MMSCf$")
	plt.ylabel("$F,$ $MMSCf$")
	plt.ylim([0, max(F)])
	plt.xlim([0, max(Eg + Efw)])
	plt.title("$Havlena-Odeh$ $Gas$ $Material$ $Balance$")
	plt.legend()
	plt.grid()
	plt.show()
	G = slope
	return G
    
def Wed_inf(Td):
	if (Td < 1):
		Wd = (2*(np.sqrt(Td/np.pi))) + (Td/2) - ((Td/6)*(np.sqrt(Td/np.pi))) + ((Td**2)/16)
	elif (1 < Td < 100):
		a = np.array([0.81638, 0.85373, (-2.7455 *(10**(-2))), (1.0284*(10**(-3))), (-2.2740*(10**(-5))), (2.8354*(10**(-7))), (-1.8436*(10**(-9))), (4.8534*(10**(-12)))])
		Wd = a[0] + (a[1]*Td) + (a[2]*(Td**2)) + (a[3]*(Td**3)) + (a[4]*(Td**4)) + (a[5]*(Td**5)) + (a[6]*(Td**6)) + (a[7]*(Td**7))
	else:
		Wd = (2*Td)/(np.log(Td))
	return Wd

def WI_Func(aq, td, delP, red):
	We = np.zeros((m))
	if (aq == "I"):
		WD = np.zeros((m))
		for l in range (1, m):
			temp = 0
			sum = 0
			for j in range(l):
				temp = td[l] - td[j]
				Wd = Wed_inf(temp)
				sum += delP[j+1]*Wd
			WD[l] = Wd
			We[l] = Wd*sum
	elif (aq == "F"):
		WD = np.zeros((m))
		tdstar = 0.4*((red**2) - 1)
		for i in range (1, m):
			if tdstar > td[i]:
				temp = 0
				sum = 0
				for j in range(i):
					temp = td[i] - td[j]
					Wd = Wed_inf(temp)
					sum += delP[j+1]*Wd
				WD[i] = Wd
				We[l] = Wd*sum
			else:
				Jstar = ((red**4)*(np.log(red))/((red**2) - 1)) + (0.25 * (1 - (3*(red**2))))
				Wd = 0.5*((red**2) - 1)*(1-np.exp(-(2*td[i]/Jstar)))
				WD[i] = Wd
				We[l] = Wd*sum
	else:
		WD = np.zeros((m))
		for l in range (1, m):
			temp = 0
			sum = 0
			for j in range(l):
				temp = td[l] - td[j]
				Wd = Wed_inf(temp)
				sum += delP[j+1]*Wd
			WD[l] = Wd
			We[l] = Wd*sum
	return WD, We
    
def CARET(m, G, We, Wp, Bw, Gp, P, Pinit, cf, cw, ha, hr, swi, z, T, phi, vis, rr, ra, aq, k, t):
	B = 1.119*phi*ct*h*(rr**2)*(theta/360)
	Pbar = np.ones((m))
	delP = np.zeros((m))
	Pbar[1:] = (P[1:] + P[:-1])/2
	Pbar[0] = Pinit
	delP[1:] = Pbar[:-1] - Pbar[1:]
	delP[0] = 0
	td = np.ones((m))
	td[:] = t[:]*(6.33*k/(phi*vis*cf*(rr**2)))
	red = ra/rr
	WD, WI = WI_Func(aq, td, delP, red)
	Bg =  .005035*(z*T/P)
	Bgi = Bg[0]
	Eg =  np.zeros((m)) 
	Efw =  np.zeros((m))
	F =  np.zeros((m))
	Ecaret = np.zeros((m))
	Vpr = G*Bgi/(1-swi)
	We = We/5.615
	Eg[:] = (Bg[:] - Bgi)
	Efw[:] = Bgi*(((Pinit - P[:])*(cf + (cw*swi)))/(1 - swi))
	F[:] = (Gp[:]*Bg[:]) + (Wp[:]*Bw)
	Ecaret[:] = (2*cf*WD[:]*Bgi*ha/((1-swi)*hr)) + Eg[:] + Efw[:]
	junk, intercept, slope = best_fit(m, Ecaret, F)
	Ffit = np.zeros((m))
	Ffit[:] = Ecaret[:]*slope
	sFE = np.sqrt(sum(((F[:] - Ffit[:])**2))/(m-1))
	vcoeff = (sFE/np.mean(F))*100
	plt.figure()
	plt.plot(Ecaret, F, "bo", label= "$CARET$")
	plt.plot(Ecaret,((Ecaret*slope)), "r--", label = "$Interpolated$")
	plt.xlabel("$Ecaret,$ $MMRB/MMSCf$")
	plt.ylabel("$F,$ $MMRB$")
	plt.ylim([0,max(F)])
	plt.xlim([0, max(Ecaret)])
	plt.title("$CARET$")
	plt.legend()
	plt.grid()

	plt.show()
	#print ("Error :",vcoeff, "%")
	#print ("UnCorrected Gas Inplace by CARET method (SPE-28630) : ", slope/35.31464, "MMm3")
	if vcoeff > 3:
		WD[:] = (((1-swi)*hr)/(2*cf*G*Bgi*ha))*(F[:] - (G*(Eg[:] + Efw[:])))
		Ecaret[:] = (2*cf*WD[:]*Bgi*ha/((1-swi)*hr)) + Eg[:] + Efw[:]
		junk, intercept, slope = best_fit(m, Ecaret, F)
		Ffit = np.zeros((m))
		Ffit[:] = Ecaret[:]*slope
		sFE = np.sqrt(sum(((F[:] - Ffit[:])**2))/(m-1))
		vcoeff = (sFE/np.mean(F))*100
	plt.figure()
	plt.plot(Ecaret, F, "bo", label= "$CARET$")
	plt.plot(Ecaret,((Ecaret*slope)), "r--", label = "$Interpolated$")
	plt.xlabel("$Ecaret,$ $MMRB/MMSCf$")
	plt.ylabel("$F,$ $MMRB$")
	plt.title("$CARET$ $(reduced error)$")
	plt.ylim([0,max(F)])
	plt.xlim([0, max(Ecaret)])
	plt.legend()
	plt.grid()
	plt.show()
	GIIP = slope
	print ("Error :",vcoeff, "%")
	return GIIP
    
def modified(m, P, z, Gp, G):
	pz = np.zeros((m), float)
	pz[:] = P[:]/z[:]
	pzi = P[0]/z[0]
	pzd = np.zeros((m))
	pzd[:] = pz[:]/pzi 
	plt.figure()
	plt.plot(Gp, pzd, "bo", label = "$P/z$")
	b = best_fit_expo(m, Gp, np.log(pzd), G)
	plt.plot(Gp,(np.exp(-(Gp/(b-Gp)))), "r--", label = "$Best$ $Fit$")
	plt.xlabel("$Cummulative$ $Production,$ $MMSCf$")
	plt.ylabel("$Dimensionless P/z,$ $Psia$")
	plt.title("$Modified$ $Gas$ $Material$ $Balance$ $(Exponential$ $Method)$")
	plt.legend()
	plt.grid()

	plt.show()
	return b
def driveindex(m, G, Bg, Bgi, We, Gp, cf, cw, P, Pinit):
	DDI = np.zeros((m))
	WDI = np.zeros((m))
	CDI = np.zeros((m))
	cep = np.zeros((m), float)
	cep[:] = (cf + (cw*swi))*(Pinit - P[:])
	for i in range(m):
		if Gp[i]*Bg[i] == 0. :
			DDI[i] = 0
			WDI[i] = 0
			CDI[i] = 0
		else:
			DDI[:] = (G*(Bg[:] - Bgi))/(Gp[:]*Bg[:])
			WDI[:] = (We[:] - (Wp[:]*Bw))/(Gp[:]*Bg[:])
			CDI[:] = cep[:]/(Gp[:]*Bg[:])
	Energy = DDI + WDI + CDI
	for j in range (m):
		if DDI[j] < 0. :
			DDI[j] = 0
			CDI[j] = CDI[j]/Energy[j]
			WDI[j] = WDI[j]/Energy[j]
			Energy[j] = CDI[j] + WDI[j]
		elif WDI[j] < 0.:
			WDI[j] = 0
			CDI[j] = CDI[j]/Energy[j]
			DDI[j] = DDI[j]/Energy[j]
			Energy[j] = CDI[j] + DDI[j]
		elif CDI[j] < 0:
			CDI[j] = 0
			DDI[j] = DDI[j]/Energy[j]
			WDI[j] = WDI[j]/Energy[j]
			Energy[j] = DDI[j] + WDI[j]
	for j in range(m):
		if Energy[j] > 1:
			DDI[j] = DDI[j]/Energy[j]
			WDI[j] = WDI[j]/Energy[j]
			CDI[j] = CDI[j]/Energy[j]
			Energy[j] = DDI[j] + WDI[j] + CDI[j]
		else:
			continue

	plt.figure()
	plt.plot(t, DDI, "r--", label = "$Depletion$ $Drive$ $Index$")
	plt.plot(t, WDI, "b--", label = "$Water$ $Drive$ $Index$")
	plt.plot(t, CDI, "c--", label = "$Compaction$ $Drive$ $Index$")
	plt.ylabel("$Drive$ $Indices$")
	plt.ylim([0,1])
	plt.xlim([0, max(t)])
	plt.xlabel("$Days$")
	plt.legend()
	plt.title("$Drive$ $Indices$ $Plot$")
	plt.grid()
	plt.show()
	return (DDI, WDI, CDI)

def write_data(G,Gadv,Gho,Gcaret,Gmod, file):
	out = open(file, "w")
	out.write("Calculated GIIP by different methods (MMm3)" + "\n")
	out.write("Classical" + "                    " + "Advanced Gas Material Balance (SPE-139428)" + "                    " + "Havlena-Odeh Method" + "                    " + "CARET Method (SPE-28630)" + "                    " + "Exponential P/z Method" + "\n")
	out.write(str(G/35.31464) + "                    " + str(Gadv/35.31464) + "                    " + str(Gho/35.31464) + "                    " + str(Gcaret/35.31464) + "                    " + str(Gmod/35.31464))
	out.close()
G = PZ(P, z, Gp)
print ("Inital Free Gas by mormal P/z method :", G/35.31464, "MMm3"  )

Gadv, Wef = PZ_advanced(We, Wp, Bg, Bgi, Bw, sgi, swi, Pinit, P, z, G)
print ("Corrected Gas Inplace by Advanced Material Balance (SPE-139428) : ", Gadv/35.31464, "MMm3")

Gho = HavlenaOdeh(m, Wp, Bw, Wef, cw, swi, P, Pinit, Gp, cf, z, T)
print ("Corrected Gas Inplace by Havlena-Odeh method : ", Gho/35.31464, "MMm3"  ) 

Gcaret = CARET(m, Gadv, Wef, Wp, Bw, Gp, P, Pinit, cf, cw, haq, h, swi, z, T, phi, vis, rr, ra, aq, k, t)
print ("Corrected Gas Inplace by CARET method (SPE-28630) : ", Gcaret/35.31464, "MMm3" )

Gmod = ((modified(m, P, z, Gp, Gcaret)))
print ("Corrected Gas Inplace by modified method : ", (Gmod/35.31464), "MMm3" )

DDI, WDI, CDI = driveindex(m, Gcaret, Bg, Bgi, We, Gp, cf, cw, P, Pinit)

GIIP = (G/35.31464, Gadv/35.31464, Gho/35.31464, Gcaret/35.31464, Gmod/35.31464)
ind = np.arange(5)
width = 0.4

fig, ax = plt.subplots()
rect =ax.bar(ind, GIIP, width, color = "r")
ax.set_ylabel("$GIIP,$ $BCM$")
ax.set_xticks(ind+width/2)
ax.set_xticklabels(("$Classical$", "$Advanced$", "$Havlena-Odeh$", "$CARET$", "$Exponential$"))
plt.show()


"""
plt.figure()
plt.plot(Wef, We, "bo", Wef, We1*5.615, "b--", Wef, We2*5.615, "r--", Wef, We3*5.615, "k--")
plt.show()
fig, ax1 = plt.subplots()
ax1.plot(t, P, "bd", label = "$Pressure,$ $Psia$")
ax1.set_ylim([0,max(t)])
ax2 = ax1.twinx()
ax2.plot(t, Wef, "b-", label = "$Cummulative$ $Water$ $Influx,$ $MMBBL$")
ax2.plot(t, Gp, "r-", label = "$Cummulative$ $Gas$ $MMBBL$")
ax2.set_ylim([0,50])
ax1.set_xlim([0,3000])
plt.show()
"""
write_data(G,Gadv,Gho,Gcaret,Gmod, "GIIP.txt")# Advanced Material Balance Calculation
