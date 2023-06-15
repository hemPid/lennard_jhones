import matplotlib.pyplot as plt
import numpy as np
import mnk

def read_dir(d):
	ret_data = []
	with open(d+'/pyout.txt', 'r') as f:
		for line in f:
			data = line.rstrip().split("=")
			ret_data.append(float(data[1]))
	return ret_data


dirs = ["t1", "t1.5", "t1.8", "t2", "t2.8", "t3", "t3.5", "t4"] #t1.8, t2.5, t2.8?

T = []
T_err = []
D = []
D_err = []
Sigma = []
Sigma_err = []

for d in dirs:
	ret = read_dir(d)
	T.append(ret[0])
	T_err.append(ret[1])
	D.append(ret[2])
	D_err.append(ret[3])
	Sigma.append(ret[4])
	Sigma_err.append(ret[5])

print(T)
print(T_err)
print(D)
print(D_err)
print(Sigma)
print(Sigma_err)

T = np.array(T)
T_err = np.array(T_err)
D = np.array(D)
D_err = np.array(D_err)
Sigma = np.array(Sigma)
Sigma_err = np.array(Sigma_err)

Tobr = 1 / T
Tobr_err = T_err / (T*T)

m = mnk.solve_mnk(Tobr, Sigma)
print(m)

saz = m[0]/m[2]
s_saz = saz * ((((m[1]/m[0])**2)+((m[3]/m[2])**2))**0.5)
print("S=",saz)
print("s_S=", s_saz)

epsk = [120, 95.1]

for x in epsk:
	print(x*saz, "+-", x*s_saz)

fig, ax = plt.subplots()

ax.errorbar(Tobr, Sigma, xerr=Tobr_err, yerr=Sigma_err, fmt='.', ecolor='red')
ax.plot(Tobr, m[2] + m[0]*Tobr)
ax.grid()
ax.set_xlabel('$1/T^*$')
ax.set_ylabel(r'$\sigma$')
plt.show()


k = 1.38e-23
c_sig = 3.7e-10
c_eps = 95.1*k
c_m = 4.65e-26
alph1 = c_sig* ((c_eps/c_m) ** 0.5)
Patm = 101300
P = (2000/((27*c_sig) ** 3)*T*c_eps)
print(P)
alph2 = P/Patm
Tex = T*c_eps/k
print(Tex)
T0 = 273
alph3 = (T0/Tex) ** 1.5
Dex = D*alph1*alph2*alph3*1e4
s_Dex = D_err*alph1*alph2*alph3*1e4
print(Dex)
print(s_Dex)

mD = Dex.mean()

s = 0

for x in Dex:
	s += (mD-x)**2

s_mD = (s/56) ** 0.5

print("mD=", mD)
print("s_mD", s_mD)