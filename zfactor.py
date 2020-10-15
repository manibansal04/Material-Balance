# Program to calculate z factor using Hall Yaarborough EOS

def z_factor(P):

	import matplotlib.pyplot as plt
	import numpy as np

	#Initialisation Section 
	Tc = {"CH4": 343.00, "C2H6": 549.59, "C3H8": 665.73, "i-C4H10": 734.13, "n-C4H10": 765.29, "i-C5H12": 828.77, "n-C5H12": 845.47, "n-C6H14": 913.27, "n-C7H16": 972.37, "n-C8H18": 1023.89, "n-C9H20":1070.35, "n-C10H22": 1111.67,"H2": 59.36, "Water": 1164.85, "N2": 227.16, "O2": 278.24, "H2S": 672.35, "CO2": 547.48}
	Tb = {"CH4": 200.94, "C2H6": 332.18, "C3H8": 415.92, "i-C4H10": 470.45, "n-C4H10": 490.75, "i-C5H12": 541.79, "n-C5H12": 556.59, "n-C6H14": 615.39, "n-C7H16": 668.83, "n-C8H18": 717.88, "n-C9H20":763.14, "n-C10H22": 805.15,"H2": 36.715, "Water": 617.67, "N2": 139.219, "O2": 162.338, "H2S": 383.173, "CO2": 350.413}
	Pc = {"CH4": 666.40, "C2H6": 706.50, "C3H8": 616.00, "i-C4H10": 527.90, "n-C4H10": 550.60, "i-C5H12": 490.40, "n-C5H12": 488.60, "n-C6H14": 436.90, "n-C7H16": 396.80, "n-C8H18": 360.70, "n-C9H20":331.80, "n-C10H22": 305.20, "H2": 187.5, "Water": 3200.10, "N2": 493.10, "O2": 731.40, "H2S": 1306.0, "CO2": 1071.0}
	Mg = {"CH4": 16.043, "C2H6": 30.070, "C3H8": 44.097, "i-C4H10": 58.123, "n-C4H10": 58.123, "i-C5H12": 72.150, "n-C5H12": 72.150, "n-C6H14": 86.177, "n-C7H16": 100.204, "n-C8H18": 114.231, "n-C9H20":128.258, "n-C10H22": 142.285, "H2": 2.1059, "Water": 18.0153, "N2": 28.0134, "O2": 31.9988, "H2S": 34.08, "CO2": 44.010}
	omega = {"CH4": 0.0115, "C2H6": 0.0995, "C3H8": 0.1523, "i-C4H10": 0.1770, "n-C4H10": 0.2002, "i-C5H12": 0.2275, "n-C5H12": 0.2515, "n-C6H14": 0.3013, "n-C7H16": 0.3495, "n-C8H18": 0.3996, "n-C9H20":0.4435, "n-C10H22": 0.4923, "H2": 0.2202, "Water": 0.3443, "N2": 0.0372, "O2": 0.0216, "H2S": 0.0948, "CO2": 0.2667}

	# Functions
	def normalized(composition):
		sum = 0
		for j in composition:
			sum += composition[j]
        
		if sum == 1.:
			correction_factor = 0
		else:
			correction_factor = 1.0/sum
        
		for l in composition:
			if composition[l] > 0.:
				composition[l] = composition[l] + correction_factor
			else:
				composition[l] = composition[l]
		return composition
	
	def pseudocritical(composition, P_c, T_c, T_b, M_g, Om):
		p_temp = 0.
		t_temp = 0.
		m_temp = 0.
		tb_temp = 0.
		Om_temp = 0.
		for k in P_c:
			p_temp = p_temp + (P_c[k]*composition[k])
			t_temp = t_temp + (T_c[k]*composition[k])
			print(k, P_c[k], composition[k], p_temp, T_c[k], composition[k],  t_temp)
			m_temp = m_temp + (M_g[k]*composition[k])
			tb_temp = tb_temp + (T_b[k]*composition[k])
			Om_temp = Om_temp + (Om[k]*composition[k])
        
		return (p_temp, t_temp, tb_temp, m_temp, Om_temp)
        
	def accentric(P_c, T_c, T_b):
		om = ((3.0*(np.log(P_c/14.70)))/(7.0*((T_c/T_b)-1))) - 1
		return om
        
	def pseudoreduced(Pr, temp, P_pc, T_pc):
		Ppr = np.zeros((len(Pr)))
		Tpr = np.zeros((len(Pr)))
		for i in range (len(Pr)):
			#print Pr[i]/P_pc, "                      |                          ", temp[i]/T_pc
			Ppr[i] =(Pr[i]/P_pc)
			Tpr[i] = (temp[i]/T_pc)
		return Ppr, Tpr

	def comp():
		gas_composition = {}
		file = input("Please gas composition enter file name: ")
		inp = open(file, "r")
		line = inp.readline()
		line = line.rstrip()
		tokens = line.split()
		m = int(tokens[0])
		for i in range(m):
			line = inp.readline()
			line = line.rstrip()
			tokens = line.split()
			k = tokens[0]
			gas_composition[k] = gas_composition.get(k, float(tokens[1]))
		inp.close()
		print(gas_composition)
		#gas_composition = normalized(gas_composition)
		return gas_composition
    
	def Pressure():
		file = input("Enter pressure file name: ")
		inp = open(file, "r")
		line = inp.readline()
		line = line.rstrip()
		tokens = line.split()
		m = int(tokens[0])
		P = np.zeros((m), float)
		for i in range (m):
			line = inp.readline()
			line = line.rstrip()
			tokens = line.split()
			P[i] = float(tokens[0])
		inp.close()
		return P
    
	def Temperature(Pr):
		Ti = float(input("Enter Temperature:"))
		T = np.ones((len(Pr)))*Ti
		return T
       
	def boiling(grav, mw):
		var = np.array((6.77857, 0.401673, -1.58262, (3.77409*10**(-3)), 2.984036, (-4.25288*10**(-3))))
		Tb = var[0]*(mw**var[1])*(grav**var[2])*np.exp((var[3]*mw)+(var[4]*grav)+(var[5]*mw*grav))
		return Tb
        
	def sutton_pseudocritical(grav, P_crit, T_crit, M_weight):
		if grav > 0.57 and grav <1.68:
			P_c = 756.8 - (131.0*grav) - (3.6*(grav)**2)
			T_c = 169.2 + (349.5*grav) - (74.0*(grav)**2)
			return (P_c, T_c)
		elif grav < 0.57:
			print ("Sutton correlation is not applicable as gs gravity is not in range 0.56<gamma<1.68. Therefore using Standing Correlation")
			P_c = 667.0 + (15.0*grav) - (37.5*(grav)**2)
			T_c = 168.0 + (325.0*grav) - (12.5*(grav)**2)
			return (P_c, T_c)
		else:
			print ("Sutton correlation is not applicable as gs gravity is not in range 0.56<gamma<1.68. Therefore gas composition is must. Exit the program start with gas composition")
			pass
			#P_c, T_c =  pseudocritical(comp, P_crit, T_crit, M_weight)
            
	def zcalc(P_r, T_r, T_pc, P_pc, Pr, Temp, acc):
		z = np.ones((len(P_r)))
		for i in range(len(P_r)):
			z[i] = zfactor_PR(P_r[i], T_r[i], T_pc, P_pc, Pr[i], Temp[i], acc)
		#z[:] = zfactor_PR(P_r[:], T_r[:], T_pc, P_pc, Pr[:], Temp[:], acc)
		return z

	def zfactor_PR(P_r, T_r, T_pc, P_pc, Pr, Temp, omega):
		z = 1
		R = 10.731
		error = 10
		while error >= 0.0001 :
			if omega <= 0.49:
				k = 0.37464 + (1.54226*omega) - (0.26992*(omega**2))
			else:
				k = 0.379642 + (1.48503*omega) - (0.164423*(omega**2)) + (0.016666*(omega**3))
			ac = 0.457236 *((R**2)*(T_pc**2)/ P_pc)
			a = ((1+(k*(1-(T_r**0.5))))**2)*ac
			b = 0.077796 *(R*T_pc/P_pc)
			A = (a*Pr)/((R**2)*(Temp**2))
			B = (b*Pr)/(R*Temp)
			F = (z**3) + ((B-1)*z) + ((A-(2*B)-(3*(B**2)))*z) + ((B**3) + (B**2) - (A*B))
			dF = (3*(z**2)) + (2*(B-1)*z) + (A - (2*B) - (3*(B**2)))
			z1 = (z - (F/dF))
			error = abs(z1 - z)
			z = z1
            
		return z
		
	def zcalc_HY(P_r, T_r):
		if isinstance(P_r, np.ndarray) == False:
			z = zfactor_HY(P_r,T_r)
		else:
			z = np.ones((len(P_r)))
			for i in range(len(P_r)):
				z[i] = zfactor_HY(P_r[i], T_r[i])
		return z
		
	def zfactor_HY(p_pr,t_pr):
		delta = 10.0
		t = (1/t_pr)
		X1 = -0.06125*p_pr*t*np.exp(-1.2*(1-t)**2)
		X2 = 14.76*t - 9.76*t**2 + 4.58*t**3 
		X3 = 90.7*t - 242.2*t**2 + 42.4*t**3
		X4 = 2.18 + 2.82*t
		Y = 0.0125*p_pr*t*np.exp(-1.2*(1-t)**2)
		
		while delta >= 0.0001:
			f = X1 + ((Y + Y**2 + Y**3 + Y**4)/((1-Y)**3))-(X2*(Y**2)) + (X3*(Y**(X4)))
			df = (1+(4*Y) + (4*(Y**2)) - (4*(Y**3)) + (Y**4))/((1-Y)**4)
			Y1 = (Y - (f/df))
			delta = abs(Y1 - Y)
			Y = Y1
		
		z = 0.06125*p_pr*t*np.exp(-1.2*(1-t)**2)/Y
		return z
		
	def visc(P, T, z, M):
		A = np.array((0.00149406, 0.00094, 0.000002, 1.5, 209, 19, 3.5, 986, 0.01, 2.4, 0.2))
		rho = A[0]*(P*M)/(z*T)
		K1 = (A[1] + A[2]*M)*(T**A[3])/(A[4]+(A[5]*M)+T)
		X = (A[6] + (A[7]/T) + (A[8]*M))
		Y = (A[9] - (A[10]*X))
		mug = K1*np.exp(X*rho**Y)
		return mug
		
	def write_data(m, P, z, Bg, mug, outfile):
		out = open(outfile, "w")
		out.write("{0} \n".format(m))
		for i in range(0, m):
			#out.write("{4:10.5f} {4:10.5f} {4:10.5f} \n".format(P[i], z[i], Bg[i]))
			out.write(str(P[i])+ "     "+ str(z[i])+ "     "+ str(Bg[i])+ "     "+ str(mug[i])+ "\n")
		out.close()

	# Input Section 
	resp = input("Is gas composition available. If yes press y else n:")
	if resp == "y":
		gas_comp = comp()
		#gas_comp = {"H2": 0.0, "Water": 0.0, "N2": 0.0138, "O2": 0.0, "H2S": 0.0, "CO2": 0.0, "CH4": 0.9302, "C2H6": 0.0329, "C3H8": 0.0136, "i-C4H10": 0.0023, "n-C4H10": 0.0037, "i-C5H12": 0.0012, "n-C5H12": 0.0010, "n-C6H14": 0.0008, "n-C7H16": 0.0005, "n-C8H18": 0.0, "n-C9H20":0.0, "n-C10H22": 0.0}
	else:
		gamma = float(input("Please enter Gas Gravity:"))
        
	#P = Pressure()
	T = Temperature(P)
	gas = input("Whether gas is sweet(does not contain H2S or CO2 or both)? press y sweet and n for sour: ")
    
	# Evaluation section   
	if resp == "y":
		Ppc, Tpc, Tpb, mol, O = pseudocritical(gas_comp, Pc, Tc, Tb, Mg, omega)
		if (gas == "n"):
			A = gas_comp["H2S"]
			B = gas_comp["CO2"]
			e = (120*((A**0.9) - (A**1.6))) + (15*((B**0.5) - (B**4)))
			Tpc_corrected = Tpc - e
			Ppc_corrected = (Ppc*Tpc_corrected)/(Tpc +(B*(1 - B)*e))
			Tpc = Tpc_corrected 
			Ppc = Ppc_corrected
		else:
			pass
		print ("Critical Pressure, Pc(psia):", Ppc)
		print ("Critical temperature, Tc(deg Rankin):", Tpc)
		print ("Boiling temperature, Tb(deg Rankin):", Tpb)
		print ("Accentric Factor :", O)
		print ("Average Molecular Weight:", mol)
	else:
		print ("YOU CAN GO AHEAD BUT STRONGLY IT IS NOT SUGGESTED!!!!!")
		Ppc, Tpc = sutton_pseudocritical(gamma, Pc, Tc, Mg)
		mol = gamma*28.97
		Tpb = boiling(gamma, mol)
		if (gas == "n"):
			A = float(input("Enter H2S fraction: "))
			B = float(input("Enter CO2 fraction: "))
			e = (120*((A**0.9) - (A**1.6))) + (15*((B**0.5) - (B**4)))
			Tpc_corrected = Tpc - e
			Ppc_corrected = (Ppc*Tpc_corrected)/(Tpc +(B*(1 - B)*e))
			Tpc = Tpc_corrected 
			Ppc = Ppc_corrected
		else:
			pass
		O = accentric(Ppc, Tpc, Tpb)
		print ("Critical Pressure(Sutton's Correlation), Pc(psia):", Ppc)
		print ("Critical temperature (Sutton's Correlation), Tc(deg Rankin):", Tpc)
		print ("Boiling temperature, Tb(deg Rankin):", Tpb)
		print ("Accentric Factor :", O)

	#print "Pr", "                      |                          ", "Tr"
	#print "---------------------------------------------------------------"
	Pr ,Tr = pseudoreduced(P, T, Ppc, Tpc)
	z = zcalc_HY(Pr,Tr)
	#z = zcalc(Pr, Tr, Tpc, Ppc, P, T, O)
	Bg = (0.02827*z*T)/P
	mug = visc(P, T, z, mol)
    
	write_data(len(P), P, z, Bg, mug, "eos.txt")
    
	#z, Bg, P = main()

	plt.plot(P,z)
	plt.title("$z$ $Factor$ $(Hall Yaarborough$ $EOS)$")
	plt.xlabel("$Pressure,$ $Psia$")
	plt.ylabel("$z$")
	plt.show()
	
	return (z, Bg)
	