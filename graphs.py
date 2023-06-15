import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import numpy as np
import mnk
import sys

d = sys.argv[1]
parts = 2000
dt = 0.001
outstring = ""

def axisdata(axis, vels, show_graphs):
	s = 10
	dv = (max(vels) - min(vels))/s
	xaxis = np.linspace(min(vels), max(vels), 50)
	phi = []
	for x in xaxis:
		phi.append(0)
		for y in vels:
			if x-dv <= y <= x+dv:
				phi[-1] += 1
			if y > x+dv:
				break
	phi = np.array(phi)/n
	vels_x = np.array(vels)
	phi_nz = []
	xaxis_nz = []

	for i in range(len(xaxis)):
		if phi[i]:
			phi_nz.append(phi[i])
			xaxis_nz.append(xaxis[i])

	phi_nz = np.array(phi_nz)
	xaxis_nz = np.array(xaxis_nz)


	factor = 0.3
	begin = int(len(xaxis_nz)*factor)
	end = int(len(xaxis_nz)*(1-factor))

	xaxis_nz_mnk = xaxis_nz[begin:end]
	phi_nz_mnk = phi_nz[begin:end]
	ln_fi_mnk = -1*np.log(phi_nz_mnk)
	sqr_x_mnk = xaxis_nz_mnk*xaxis_nz_mnk
	k, s_k, b, s_b = mnk.solve_mnk(sqr_x_mnk, ln_fi_mnk)

	print(axis, ": ")
	print("k = ", k)
	print("s_k = ", s_k)
	print("b = ", b)
	print("s_b = ", s_b)
	print("T* = ", 1/(2*k))
	print("sT* = ", s_k/(2*k*k))

	if show_graphs:
		fig, ax = plt.subplots()
		ax.scatter(xaxis_nz, phi_nz, s=2)
		ax.plot(xaxis_nz, np.exp(-1*(b+k*(xaxis_nz*xaxis_nz))), color="red")
		ax.grid()
		ax.set_xlabel('$v_'+axis+'$')
		ax.set_ylabel(r'$\varphi(v_'+axis+')$')
		plt.show()

		fig, ax = plt.subplots()
		ax.scatter(sqr_x_mnk, ln_fi_mnk, s=2)
		ax.plot(sqr_x_mnk, k*sqr_x_mnk + b, color="red")
		ax.grid()
		ax.set_xlabel('$v_'+axis+'^2$')
		ax.set_ylabel(r'$-\ln\varphi$')
		plt.show()
	return (1/(2*k), s_k/(2*k*k))


tick = []
energy = []
penergy = []
kenergy = []

with open(d+'/energy.txt', 'r') as f:
	for line in f:
		l = line.rstrip().split(", ")
		t_str = l[0].split(": ")
		e_str = l[1].split(": ")
		p_str = l[2].split(": ")
		k_str = l[3].split(": ")
		tick.append(float(t_str[1]))
		energy.append(float(e_str[1]))
		penergy.append(float(p_str[1]))
		kenergy.append(float(k_str[1]))
tick = np.array(tick)
energy = np.array(energy)/parts
penergy = np.array(penergy)/parts
kenergy = np.array(kenergy)/parts
mpenergy = (penergy[8000:]).mean()
mkenergy = (kenergy[8000:]).mean()

print((2*mkenergy)/3)

fig, ax = plt.subplots()
ax.grid()
ax.plot(tick, energy)
ax.set_xlabel('steps')
ax.set_ylabel(r'Full energy, $E/\varepsilon$')
plt.show()

fig, ax = plt.subplots()
ax.grid()
ax.plot(tick, energy, tick, penergy, tick, kenergy)#, tick[2:], mpenergy*(tick[2:]/tick[2:]), tick[2:], mkenergy*(tick[2:]/tick[2:])
ax.legend(("Full Energy", "Potential energy", "Kinetic Energy"),
          shadow=True, loc=(0.01, 0.48), handlelength=1.5, fontsize=16)
ax.set_xlabel('steps')
ax.set_ylabel('Energy, $E/\\varepsilon$')
#plt.ylim(-50,300)
plt.show()

n= 0;
vels = []
vels_x = []
vels_y = []
vels_z = []

with open(d+'/vels.txt', 'r') as f:
	for line in f:
		if (line[0] == "n"):
			n = int(line.rstrip()[2:])
			print(n)
			continue;
		l = line.rstrip().split(" ")
		vels.append(float(l[1]))
		vels_x.append(float(l[2]))
		vels_y.append(float(l[3]))
		vels_z.append(float(l[4]))

vels_x.sort()
vels_y.sort()
vels_z.sort()
vels = np.array(vels)

flag = False

xT, s_xT = axisdata("x", vels_x, flag)
yT, s_yT = axisdata("y", vels_y, flag)
zT, s_zT = axisdata("z", vels_z, flag)
mT = (xT + yT + zT)/3
print("mean T*=", mT)
outstring += "mean T*=" + str(mT) + '\n'
smT = ((((mT-xT)**2) + ((mT-yT)**2) + ((mT-zT)**2))/6)**0.5
print("sT*=", smT)
outstring += "sT*=" + str(smT) + '\n'

ticks = []
drs = []

with open(d+'/diffusion.txt', 'r') as f:
	for line in f:
		l = line.rstrip().split(", ")
		t_str = l[0].split(": ")
		e_str = l[1].split(": ")
		ticks.append(float(t_str[1]))
		drs.append(float(e_str[1]))

ticks = np.array(ticks)
drs = np.array(drs)
l = len(ticks)
ticks = ticks*dt
t0 = 3000
tend = 4500

m = mnk.solve_mnk(ticks[t0:tend], drs[t0:tend])


fig, ax = plt.subplots()
ax.plot(ticks, drs)
ax.plot(ticks, m[2] + m[0]*ticks)
ax.grid()
ax.set_xlabel('$t$')
ax.set_ylabel(r'$\overline{r^2}$')
plt.show()

D = m[0]/6
sD = m[1]/6
print(D)
outstring += "D=" + str(D) + '\n'
print(sD)
outstring += "sD=" + str(sD) + '\n'
mv = vels.mean()

dmv2_sum = 0
count = 0
for x in vels:
	dmv2_sum += (mv - x)**2
	count += 1
s_mv = (dmv2_sum/(count*(count-1)))**0.5
print("mv=",mv)
print("s_mv=",s_mv)
lam = (3*D)/(2*mv)
s_lam = lam*((((sD/D)**2)+((s_mv/mv)**2))**0.5)
print("lam=", lam)
print("s_lam=", s_lam)
n = parts / (27 ** 3)
sigma = 1/(n*lam*(2**0.5))
s_sigma = sigma*(s_lam/lam)
print("sigma=", sigma)
outstring += "sigma=" + str(sigma) + '\n'
print("s_sigma=", s_sigma)
outstring += "s_sigma=" + str(s_sigma) + '\n'



k = 1.38e-23
c_sig = 3.4e-10
c_eps = 0.0104*1.6e-19
c_m = 66.3e-27
alph1 = c_sig* ((c_eps/c_m) ** 0.5)
Patm = 101300
P = (parts/((27*c_sig) ** 3)*mT*c_eps)
print(P)
alph2 = P/Patm
Tex = mT*c_eps/k
print(Tex)
T0 = 273
alph3 = (T0/Tex) ** 1.5
Dex = D*alph1*alph2*alph3
s_Dex = sD*alph1*alph2*alph3
print(Dex*1e4)
print(s_Dex*1e4)
"""
lticks = np.log(ticks[2000:8000])
ldrs = np.log(drs[2000:8000])

plt.plot(lticks, ldrs)
plt.grid()
plt.show()
"""
"""
with open(d + '/pyout.txt', 'w') as f:
	f.write(outstring)"""