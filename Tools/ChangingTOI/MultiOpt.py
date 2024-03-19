#Copyright (c) 2013-2022, The MercuryDPM Developers Team. All rights reserved.
#For the list of developers, see <http://www.MercuryDPM.org/Team>.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name MercuryDPM nor the
#    names of its contributors may be used to endorse or promote products
#    derived from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE MERCURYDPM DEVELOPERS TEAM BE LIABLE FOR ANY
#DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


# This expansion of SimpleOpt tool uses a custom optimization instruments from scipy.optimize
# to achieve the desirable system parameters in a simulation-guided optimization

import sys
import os
import numpy as np
from src.color_lib import colorClass
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline
from noisyopt import minimizeSPSA
from noisyopt import minimizeCompass


# Global parameters
SourceDir = "/home/iostanin/Desktop/MercuryGit/MercurySource"				# Mercury source dir
BuildDir = "/home/iostanin/Desktop/MercuryGit/MercuryBuild"					# Mercury Build dir

N = 5						# Number of maneuver reference points
TEST_MODE = False			# True for self-test on toy functional, False for real simulation-guided optimization
TEST_OPTIONS = {'xtol':1e-3, 'ftol':1e-4, 'maxiter':3} 	# Optimizer options (test mode)
REAL_OPTIONS = {'xtol':1e-3, 'ftol':1e-3, 'maxiter':100}	# Optimizer options (real mode)
QMin = 0.5					# Minimum allowed value of principal moment of inertia
QMax = 1.5					# Maximum allowed value of principal moment of inertia
I0 = 10 					# Initial (spherical) moment of inertia

ProgDuration = 100 					# Maneuver duration (absolute time)
SymDuration = 105 					# Maneuver duration (absolute time)
BaseAngvel = 5					# Initial angular velocity of the rotation

Timestep = 5.56833e-05 * 50

NTimesteps = ProgDuration / Timestep				# Number of timesteps in <duration> - should be sync W number of program timesteps!

QInitValues = np.array([1, 1])  # Initial/final conditions for control parameters
QFinalValues = np.array([1, 1])

x0 = np.ones(2*N) #Initial guess for control vector sequence
print(x0)
xc = np.random.rand(3*N) # Known solution (for toy functional only)


"""
# Transition 1-2
initial_theta = np.pi/2 	  # Initial orientation is along z axis
initial_phi = np.pi/4

final_theta = np.pi/4 	  # Final orientation is along x axis
final_phi = np.pi/2

# Transition 2-3
initial_theta = np.pi/4 	  # Initial orientation is along z axis
initial_phi = np.pi/2

final_theta = np.pi/4 	  # Final orientation is along x axis
final_phi = 0

# Transition 3-1
initial_theta = np.pi/4 	  # Initial orientation is along z axis
initial_phi = 0

final_theta = np.pi/2 	  # Final orientation is along x axis
final_phi = np.pi/4

# Transition 1-4
initial_theta = np.pi/2 	  # Initial orientation is along z axis
initial_phi = np.pi/4

final_theta = np.arccos(1/3**0.5) 	  # Final orientation is along y axis
final_phi = np.pi/4

# Transition 5-8
delta = np.pi/1000.
initial_theta = delta 	  # Initial orientation is along z axis
initial_phi = 0

final_theta = np.pi 	  # Final orientation is along x axis
final_phi = 0

# Transition 5-7
delta = np.pi/1000.
initial_theta = delta 	  # Initial orientation is along z axis
initial_phi = 0

final_theta = np.pi/2 	  # Final orientation is along x axis
final_phi = np.pi/2
"""

# Transition 4-5
Delta = np.pi / 1000.
InitialTheta = np.arccos(1 / 3 ** 0.5)	  # Initial orientation is along x=y=z
InitialPhi = np.pi/4

FinalTheta = Delta 	  # Final orientation is along z axis
FinalPhi = 0

def ComputeInertiaProfiles(nodal_values, duration, N_t):
	# pre-computes time evolution of tensor of inertia for MercuryDPM driver file
	# nodal_values is an array N X 3  equispaced values of three principal moments of inertia
	N = len(nodal_values)//3
	n_I1 = nodal_values[:N]
	n_I2 = nodal_values[N:2*N]
	n_I3 = nodal_values[2*N:]
	x = np.arange(0, duration * (N/(N-1)), duration/(N-1))
	cs1 = CubicSpline(x, n_I1, bc_type = 'clamped')
	cs2 = CubicSpline(x, n_I2, bc_type = 'clamped')
	cs3 = CubicSpline(x, n_I3, bc_type = 'clamped')
	xc = np.arange(0, duration * (N_t/(N_t-1)), duration/(N_t-1))
	I1 = cs1(xc); dI1 = cs1(xc, 1)
	I2 = cs2(xc); dI2 = cs2(xc, 1)
	I3 = cs3(xc); dI3 = cs3(xc, 1)
	return x, n_I1, n_I2, n_I3, xc, I1, dI1, I2, dI2, I3, dI3


FunCount = 0
FunBatch = False
# This function computes the value of the functional being optimized based on the simulation run in MercuryDPM.
def RealFun(x):
	import os
	import subprocess
	global SourceDir, BuildDir, N
	global ProgDuration, SymDuration, NTimesteps, FunCount, FunBatch
	global QMax, QMin, I0, QInitValues, QFinalValues

	ret = 0 # Return zero by default

	# Project the unconstrained vector of optimization parameters x into (QMin,QMax)
	q = (QMax + QMin) / 2. + (QMax - QMin) / 2. * np.cos(x)

	# Expand the set of control parameters with initial and final values
	q_ext_1 = np.hstack((QInitValues[0], q[:N], QFinalValues[0]))
	q_ext_2 = np.hstack((QInitValues[1], q[N:], QFinalValues[1]))

	# Compute principal components of inertia
	I_1 = q_ext_2 * I0
	I_2 = q_ext_1 * I0
	I_3 = q_ext_1 * q_ext_2 * I0
	nodal_values = np.hstack((I_1, I_2, I_3))
	print("Nodal values of principal moments of inertia:", nodal_values)

	# Compute cubic splines approximating the evolution of TOI
	xx, n_I1, n_I2, n_I3, xc, I1, dI1, I2, dI2, I3, dI3 = ComputeInertiaProfiles(nodal_values, ProgDuration, NTimesteps)

	arr = np.transpose(np.vstack((xc, I1, I2, I3, dI1, dI2, dI3)))

	# Expanded evolution incluidng the final piece with constant TOI
	I1_ext = np.hstack((I1, I1[len(I1)-1] * np.ones((int) (NTimesteps * (SymDuration - ProgDuration) / SymDuration))))
	I2_ext = np.hstack((I2, I2[len(I2)-1] * np.ones((int) (NTimesteps * (SymDuration - ProgDuration) / SymDuration))))
	I3_ext = np.hstack((I3, I3[len(I3)-1] * np.ones((int) (NTimesteps * (SymDuration - ProgDuration) / SymDuration))))
	T_I = np.arange(0, SymDuration, SymDuration / len(I1_ext)) # Timestamp

	# Save inertia control parameters of the simulation
	np.savetxt(BuildDir + '/Drivers/Clump/ChangingTOI/opt/inertia_profiles.txt', arr)
	np.savetxt(BuildDir + '/Drivers/Clump/ChangingTOI/opt/I_1.txt', I1_ext)
	np.savetxt(BuildDir + '/Drivers/Clump/ChangingTOI/opt/I_2.txt', I2_ext)
	np.savetxt(BuildDir + '/Drivers/Clump/ChangingTOI/opt/I_3.txt', I3_ext)
	np.savetxt(BuildDir + '/Drivers/Clump/ChangingTOI/opt/T_I.txt', T_I)
	np.savetxt(BuildDir + '/Drivers/Clump/ChangingTOI/opt/nodal_values.txt', nodal_values)

	# Run Mercury Driver file
	# (feed program duration, simulation duration, and initial angular velocity
	# as command line arguments)
	os.chdir(BuildDir + '/Drivers/Clump/ChangingTOI')
	command  = [BuildDir + '/Drivers/Clump/ChangingTOI/ChangingTOI']
	command.append('-p1')
	command.append(str(ProgDuration))
	command.append('-p2')
	command.append(str(SymDuration))
	command.append('-p3')
	command.append(str(BaseAngvel))
	print(command)
	subprocess.run(command)

	# Get the functional value
	f = open(BuildDir + '/Drivers/Clump/ChangingTOI/opt/functional.txt', "r")
	ret = float(f.read())

	# Save the functional value in the in the hystory batch
	if (FunCount == 0):
		FunBatch = ret
	else:
		FunBatch = np.append(FunBatch, ret)
		np.savetxt(BuildDir + '/Drivers/Clump/ChangingTOI/opt/fun_batch.txt', FunBatch)

	FunCount+=1 # Number of functional evaluations
	return ret

# Toy functional necessary for testing purposes
def ToyFun(x):
	global xc
	ret = 0
	for m in range(len(x)):
		ret = ret + (x[m] - xc[m])**2
	return ret


def main():
	global SourceDir, BuildDir, initial_phi, InitialTheta, FinalPhi, FinalTheta

	# Set up terminal output
	clr = colorClass()
	clr.set_use_colors(True)
	clr.def_colors()

	# Hello tag
	print(clr.END + clr.BOLD + clr.VIOLET + "-------------------------------------------")
	print(clr.END + clr.BOLD + clr.VIOLET + "     MercuryDPM MultiOpt experiment        ")
	print(clr.END + clr.BOLD + clr.VIOLET + "-------------------------------------------")


	# If the values for sourceDir and Build dir are provided in the command, then overwrite the default values
	if (len(sys.argv) == 3):  # Full format: "run source_dir build_dir"
		SourceDir = sys.argv[1]
		BuildDir = sys.argv[2]
	print("Source dir: ", SourceDir)
	print("Build dir: ", BuildDir)

	# Path to optimization directory
	OptPath = BuildDir + "/Drivers/Clump/ChangingTOI/opt"
	if ( not os.path.exists( OptPath )):
		os.mkdir( OptPath )
	
	# Let Mercury code know what is the goal orientation
	np.savetxt(BuildDir + '/Drivers/Clump/ChangingTOI/opt/init_orientation.txt', (InitialTheta, InitialPhi), delimiter =' ')
	np.savetxt(BuildDir + '/Drivers/Clump/ChangingTOI/opt/final_orientation.txt', (FinalTheta, FinalPhi), delimiter =' ')

	print(clr.END)

	if (TEST_MODE):
		res = minimize(ToyFun, x0, method='Powell', options = TEST_OPTIONS)
	else:
		res = minimize(RealFun, x0, method='Powell', options = REAL_OPTIONS)
		# Alternative optimization tools
		#res = minimizeSPSA(real_fun, x0=x0, niter=1000, paired=False)
		#res = minimizeCompass(real_fun, x0=x0, deltatol=0.01, paired=False)

	print(res)

	if (TEST_MODE): print(clr.END + "Converged to the known solution: ", np.allclose(res["x"], xc))

# Main script that runs the maneuver optimization procedure
if __name__ == "__main__":
    sys.exit(main())