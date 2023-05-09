#!/usr/bin/env python3
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

rc0 = 0.5
mus0 = 0.5
mur0 = 0.1
bo0 = 0.0
rc = ["0.3", "0.5", "0.8"]
mus = ["0", "0.5", "1"]
mur = ["0", "0.1", "0.5"]
bo = ["0.0", "1", "5", "25"]
fig = plt.figure(figsize=(26 / 2.54, 14 / 2.54))
fig.suptitle('LinearViscoelasticFrictionReversibleAdhesiveSpecies', fontsize=16)
sMax = 500
alim=[20, 80]
dlim=[10, 90]

# Restitution coefficient
splot = 1
x = [float(b) for b in rc]
aor = np.zeros(len(rc))
daor = np.zeros(len(rc))
shear = np.zeros([len(rc), 4])
for i in range(len(rc)):
	name = "CalibrationHeap_%s_%s_%s_%s.txt" % (rc[i], mus0, mur0, bo0)
	aor[i] = float(open(name, 'r').readline().split()[0])
	name = "CalibrationDrum_%s_%s_%s_%s.txt" % (rc[i], mus0, mur0, bo0)
	daor[i] = float(open(name, 'r').readline().split()[0])
	name = "CalibrationShearCell_%s_%s_%s_%s.txt" % (rc[i], mus0, mur0, bo0)
	shear[i,:] = np.array([float(v) for v in open(name, 'r').readline().split()])
ax = fig.add_subplot(3, 4, splot)
ax.plot(x, aor,'.-')
ax.set_ylim(alim)
ax.set_ylabel("AoR")
ax = fig.add_subplot(3, 4, splot+4)
ax.plot(x, daor,'.-')
ax.set_ylim(dlim)
ax.set_ylabel("dAoR")
ax = fig.add_subplot(3, 4, splot+8)
ax.plot(x, shear,'.-')
ax.set_ylim([0, sMax])
ax.set_xlabel("rc\nmu_s=%s mu_r=%s Bo=%s" % (mus0, mur0, bo0))
ax.set_ylabel("shear")

# Sliding friction
splot += 1
aor = np.zeros(len(mus))
daor = np.zeros(len(mus))
shear = np.zeros([len(mus), 4])
for i in range(len(mus)):
	name = "CalibrationHeap_%s_%s_%s_%s.txt" % (rc0, mus[i], mur0, bo0)
	aor[i] = float(open(name, 'r').readline().split()[0])
	name = "CalibrationDrum_%s_%s_%s_%s.txt" % (rc0, mus[i], mur0, bo0)
	daor[i] = float(open(name, 'r').readline().split()[0])
	name = "CalibrationShearCell_%s_%s_%s_%s.txt" % (rc0, mus[i], mur0, bo0)
	shear[i, :] = np.array([float(v) for v in open(name, 'r').readline().split()])
ax = fig.add_subplot(3, 4, splot)
ax.plot(x, aor,'.-')
ax.set_ylim(alim)
ax = fig.add_subplot(3, 4, splot+4)
ax.plot(x, daor,'.-')
ax.set_ylim(dlim)
ax = fig.add_subplot(3, 4, splot+8)
ax.plot(x, shear,'.-')
ax.set_ylim([0, sMax])
ax.set_xlabel("mu_s\nrc=%s mu_r=%s Bo=%s" % (rc0, mur0, bo0))

# Rolling friction
splot += 1
x = [float(b) for b in mur]
aor = np.zeros(len(mur))
daor = np.zeros(len(mur))
shear = np.zeros([len(mur), 4])
for i in range(len(mur)):
	name = "CalibrationHeap_%s_%s_%s_%s.txt" % (rc0, mus0, mur[i], bo0)
	aor[i] = float(open(name, 'r').readline().split()[0])
	name = "CalibrationDrum_%s_%s_%s_%s.txt" % (rc0, mus0, mur[i], bo0)
	daor[i] = float(open(name, 'r').readline().split()[0])
	name = "CalibrationShearCell_%s_%s_%s_%s.txt" % (rc0, mus0, mur[i], bo0)
	shear[i, :] = np.array([float(v) for v in open(name, 'r').readline().split()])
ax = fig.add_subplot(3, 4, splot)
ax.plot(x, aor,'.-')
ax.set_ylim(alim)
ax = fig.add_subplot(3, 4, splot+4)
ax.plot(x, daor,'.-')
ax.set_ylim(dlim)
ax = fig.add_subplot(3, 4, splot+8)
ax.plot(x, shear,'.-')
ax.set_ylim([0, sMax])
ax.set_xlabel("mu_r\nrc=%s mu_s=%s Bo=%s" % (rc0, mus0, bo0))

# Bond number
splot += 1
x = [float(b) for b in bo]
aor = np.zeros(len(bo))
daor = np.zeros(len(bo))
shear = np.zeros([len(bo), 4])
for i in range(len(bo)):
	name = "CalibrationHeap_%s_%s_%s_%s.txt" % (rc0, mus0, mur0, bo[i])
	aor[i] = float(open(name, 'r').readline().split()[0])
	name = "CalibrationDrum_%s_%s_%s_%s.txt" % (rc0, mus0, mur0, bo[i])
	daor[i] = float(open(name, 'r').readline().split()[0])
	name = "CalibrationShearCell_%s_%s_%s_%s.txt" % (rc0, mus0, mur0, bo[i])
	shear[i, :] = np.array([float(v) for v in open(name, 'r').readline().split()])
ax = fig.add_subplot(3, 4, splot)
ax.plot(x, aor,'.-')
ax.set_ylim(alim)
ax = fig.add_subplot(3, 4, splot+4)
ax.plot(x, daor,'.-')
ax.set_ylim(dlim)
ax = fig.add_subplot(3, 4, splot+8)
ax.plot(x, shear,'.-')
ax.set_ylim([0, sMax])
ax.set_xlabel("bo\nrc=%s mu_r=%s mu_s=%s" % (rc0, mur0, mus0))
plt.show()
