# Water Influx Calculation
def water_influx(m, phi, ct, h, ra, rr, theta, k, vis, aq, bnd, Pinit, P, days):
	import numpy as np
	import matplotlib.pyplot as plt

	def Wed_inf(Td):
		if (Td < 1):
			Wd = (2*(np.sqrt(Td/np.pi))) + (Td/2) - ((Td/6)*(np.sqrt(Td/np.pi))) + ((Td**2)/16)
		elif (1 < Td < 100):
			a = np.array([0.81638, 0.85373, (-2.7455 *(10**(-2))), (1.0284*(10**(-3))), (-2.2740*(10**(-5))), (2.8354*(10**(-7))), (-1.8436*(10**(-9))), (4.8534*(10**(-12)))])
			Wd = a[0] + (a[1]*Td) + (a[2]*(Td**2)) + (a[3]*(Td**3)) + (a[4]*(Td**4)) + (a[5]*(Td**5)) + (a[6]*(Td**6)) + (a[7]*(Td**7))
		else:
			Wd = (2*Td)/(np.log(Td))
		return Wd

	def VEH(B, delP, td, red, aq, m):
		We = np.zeros((m))
		if (aq == "I"):
			for k in range (1, m):
				temp = 0
				sum = 0
				for j in range(k):
					temp = td[k] - td[j]
					Wd = Wed_inf(temp)
					sum += delP[j+1]*Wd
				We[k] = B*sum
		elif (aq == "F"):
			tdstar = 0.4*((red**2) - 1)
			for i in range (1, m):
				if tdstar > td[i]:
					temp = 0
					sum = 0
					for j in range(i):
						temp = td[i] - td[j]
						Wd = Wed_inf(temp)
						sum += delP[j+1]*Wd
					We[i] = B*sum
				else:
					Jstar = ((red**4)*(np.log(red))/((red**2) - 1)) + (0.25 * (1 - (3*(red**2))))
					Wd = 0.5*((red**2) - 1)*(1-np.exp(-(2*td[i]/Jstar)))
					We[i] = B*sum
		else:
			for k in range (1, m):
				temp = 0
				sum = 0
				for j in range(k):
					temp = td[k] - td[j]
					Wd = Wed_inf(temp)
					sum += delP[j+1]*Wd
				We[k] = B*sum
		return We

	def CT(B, td, delP, m):
		We = np.zeros((m))
		Pd = np.zeros((m))
		dP = np.zeros((m))
		for i in range (1, m):
			if td[i] < 100:
				Pd[i] = ((380.529*np.sqrt(td[i])) + (137.528*td[i]) + (5.60549*(td[i]**1.5)))/(328.834 + (265.488*np.sqrt(td[i])) + (45.2157*td[i]) + (td[i]**1.5))
				E = 716.441 + (46.7984*(td[i]**0.5)) + (270.038*td[i]) + (71.0098*(td[i]**1.5))
				F = (1296.86*(td[i]**0.5)) + (1204.73*td[i]) + (618.618*(td[i]**1.5)) + (538.072*(td[i]**2)) + (142.41*(td[i]**2.5))
				dP[i] = E/F
				We[i] = We[i-1] + (td[i] - td[i-1])*(((B*delP[i]) - We[i-1]*dP[i])/(Pd[i] - td[i]*dP[i]))
			else:
				Pd[i] = 0.5*(np.log(td[i]) + 0.80907)
				dP[i] = 1/(2*td[i])
				We[i] = We[i-1] + (td[i] - td[i-1])*(((B*delP[i]) - We[i-1]*dP[i])/(Pd[i] - td[i]*dP[i]))
		return We

	def FETKOVICH(P, Pinit, days, k, h, phi, theta, vis, ct, ra, rr, red, aq, bnd):
		We = np.zeros((m))
		dWe = np.zeros((m))
		J = np.zeros((m))
		A = np.zeros((m))
		Pr = np.zeros((m))
		Pa = np.zeros((m))
		Pr[0] = Pinit
		sum = 0
		Wi = np.pi*((ra**2) - (rr**2))*h*phi/5.615
		Wei = ct*Wi*Pinit*(theta/360)
		for i in range(1, m):
			if aq == "F" and bnd == "NF":
				J[i] = (0.00708*k*1000*h*(theta/360))/(vis*(np.log(red) - 0.75))
			elif aq == "F" and bnd == "CP":
				J[i] = (0.00708*k*1000*h*(theta/360))/(vis*np.log(red))
			else:
				A[i] = np.sqrt(0.0142*k*1000*days[i]/(vis*ct*(theta/360)))
				J[i] = (0.00708*k*1000*h*(theta/360))/(vis/np.log(A[i]/rr))
			Pr[i] = (P[i] + P[i-1])/2
			Pa[i-1] = Pinit*(1-(We[i-1]/Wei))
			dWe[i] = (Wei/Pinit)*(Pa[i-1] - Pr[i])*(1- np.exp(-J[i]*Pinit*(days[i] - days[i-1])/Wei))
			sum += dWe[i]
			We[i] = sum
		return We
        
	def write_data(m, t, P, We1, We2, file):
		out = open(file, "w")
		out.write(str(m) + "\n")
		out.write("time(days)" + "                " + "Pressure(Psia)" + "                " + "Van Everdingen-Hurst" + "                " + "Carter-Tracy" + "\n")
		for i in range(0, m):
			out.write(str(t[i]) + "               " + str(P[i]) + "               " + str(We1[i]) + "               " + str(We2[i]) + "\n")
		out.close()

	B = 1.119*phi*ct*h*(rr**2)*(theta/360)
	Pbar = np.ones((m))
	delP = np.zeros((m))
	Pbar[1:] = (P[1:] + P[:-1])/2
	Pbar[0] = Pinit
	delP[1:] = Pbar[:-1] - Pbar[1:]
	delP[0] = 0
	td = np.ones((m))
	td[:] = days[:]*(6.33*k/(phi*vis*ct*(rr**2)))
	red = ra/rr
    
	Influx_VEH = VEH(B, delP, td, red, aq, m)/1000000
	Influx_CT = CT(B, td, delP, m)/1000000
	#Influx_FETKOVICH = FETKOVICH(P, Pinit, days, k, h, phi, theta, vis, ct, ra, rr, red, aq, bnd)/1000000
   
	write_data(m, days, P, Influx_VEH, Influx_CT, "waterinflux.txt")
	return (Influx_VEH, Influx_CT)