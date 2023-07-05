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


# This tool uses a custom optimization instruments from scipy.optimize 
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
source_dir = "/home/iostanin/Desktop/MercuryGit/MercurySource"				# Mercury source dir
build_dir = "/home/iostanin/Desktop/MercuryGit/MercuryBuild"					# Mercury Build dir
N = 100						# Number of maneuvers
TEST_MODE = False				# True for self-test on toy functional, False for real simulation-guided optimization
TEST_OPTIONS = {'xtol':1e-3, 'ftol':1e-4, 'maxiter':3} 	# Optimizer options (test mode)
REAL_OPTIONS = {'xtol':1e-3, 'ftol':1e-3, 'maxiter':100}	# Optimizer options (real mode)
q_min = 0.5					# Minimum allowed value of principal moment of inertia
q_max = 1.5					# Maximum allowed value of principal moment of inertia
I_0 = 10 					# Initial (spherical) moment of inertia

#prog_duration = 20 					# Maneuver duration (absolute time)
#sym_duration = 25 					# Maneuver duration (absolute time)
#base_angvel = 7						# Initial angular velocity of the rotation


prog_duration = 9.0 					# Maneuver duration (absolute time)
sym_duration = 20 					# Maneuver duration (absolute time)
base_angvel = 5						# Initial angular velocity of the rotation


N_timesteps = 1000				# Number of timesteps in <duration>

q_init_values = np.array([1, 1])
q_final_values = np.array([1, 1])


x0 = np.ones(2*N)
#for i in range(2*N):
#	if (i/2==i//2): x0[i] = 0
#	else: x0[i] = np.pi

print(x0)
xc = np.random.rand(3*N) # Known solution (for toy functional only)

delta = np.pi/1000

initial_theta = delta 	  # Initial orientation is along z axis
initial_phi = 0

final_theta = np.pi/2 	  # Final orientation is along x axis
final_phi = 0
"""
final_theta = np.arccos(1/3**0.5) 	  # Final orientation is along y axis
final_phi = np.pi/4
"""

def compute_inertia_profiles(nodal_values, duration, N_t):
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


fun_count = 0
fun_batch = False
# This function computes the value of the functional being optimized based on the simulation run in MercuryDPM.
def real_fun(x):
	import os
	import subprocess
	global source_dir, build_dir, N
	global prog_duration, sym_duration, N_timesteps, fun_count, fun_batch
	global q_max, q_min, I_0

	ret = 0

	
	q_ext_1 = np.ones(N+2)
	q_ext_2 = np.ones(N+2)
	
	q_ext_1[1:N+1] = 4./3.
	q_ext_2[1:N+1] = 3./4.
	

	I_1 = q_ext_2 * I_0
	I_2 = q_ext_1 * I_0
	I_3 = q_ext_1 * q_ext_2 * I_0

	nodal_values = np.hstack((I_1, I_2, I_3))

	print("Nodal values of principal moments of inertia:", nodal_values)
	xx, n_I1, n_I2, n_I3, xc, I1, dI1, I2, dI2, I3, dI3 = compute_inertia_profiles(nodal_values, prog_duration, N_timesteps)

	arr = np.transpose(np.vstack((xc, I1, I2, I3, dI1, dI2, dI3)))

	I1_ext = np.hstack((I1, I1[len(I1)-1] * np.ones((int) (N_timesteps * (sym_duration - prog_duration)/sym_duration))))
	I2_ext = np.hstack((I2, I2[len(I2)-1] * np.ones((int) (N_timesteps * (sym_duration - prog_duration)/sym_duration))))
	I3_ext = np.hstack((I3, I3[len(I3)-1] * np.ones((int) (N_timesteps * (sym_duration - prog_duration)/sym_duration))))
	T_I = np.arange(0, sym_duration, sym_duration / len(I1_ext))
	np.savetxt(build_dir+'/Drivers/Clump/ChangingTOI/opt/inertia_profiles.txt', arr)
	np.savetxt(build_dir+'/Drivers/Clump/ChangingTOI/opt/I_1.txt', I1_ext)
	np.savetxt(build_dir+'/Drivers/Clump/ChangingTOI/opt/I_2.txt', I2_ext)
	np.savetxt(build_dir+'/Drivers/Clump/ChangingTOI/opt/I_3.txt', I3_ext)
	np.savetxt(build_dir+'/Drivers/Clump/ChangingTOI/opt/T_I.txt', T_I)
	np.savetxt(build_dir+'/Drivers/Clump/ChangingTOI/opt/nodal_values.txt', nodal_values)






	# Run Mercury Driver file
	os.chdir(build_dir+'/Drivers/Clump/ChangingTOI')
	command  = [build_dir+'/Drivers/Clump/ChangingTOI/ChangingTOI']
	command.append('-p1')
	command.append(str(prog_duration))
	command.append('-p2')
	command.append(str(sym_duration))
	command.append('-p3')
	command.append(str(base_angvel))

	#for i in range(len(x)):
	#	command.append('-p'+str(i))
	#	command.append(str(x[i]))
	print(command)

	subprocess.run(command)

	f = open(build_dir+'/Drivers/Clump/ChangingTOI/opt/functional.txt', "r")
	ret = float(f.read())

	return ret


def toy_fun(x):
	global xc
	ret = 0
	for m in range(len(x)):
		ret = ret + (x[m] - xc[m])**2
	return ret


def main():
	global source_dir, build_dir, initial_phi, initial_theta, final_phi, final_theta

	# Set up terminal output
	clr = colorClass()
	clr.set_use_colors(True)
	clr.def_colors()

	# Hello tag
	print(clr.END + clr.BOLD + clr.VIOLET + "-------------------------------------------")
	print(clr.END + clr.BOLD + clr.VIOLET + "     MercuryDPM MultiOpt experiment        ")
	print(clr.END + clr.BOLD + clr.VIOLET + "-------------------------------------------")
	
	if (len(sys.argv) == 3):  # Full format: "run source_dir build_dir"
		source_dir = sys.argv[1]
		build_dir = sys.argv[2]

	source_dir = "/home/iostanin/Desktop/MercuryGit/MercurySource"				# Mercury source dir
	build_dir = "/home/iostanin/Desktop/MercuryGit/MercuryBuild"

	print("Source dir: ", source_dir)
	print("Build dir: ", build_dir)
	
	opt_path = build_dir + "/Drivers/Clump/ChangingTOI/opt"
	if ( not os.path.exists( opt_path )):
		os.mkdir( opt_path )
	
	# Let Mercury code know what is the goal orientation
	np.savetxt(build_dir+'/Drivers/Clump/ChangingTOI/opt/init_orientation.txt', (initial_theta, initial_phi), delimiter = ' ')
	np.savetxt(build_dir+'/Drivers/Clump/ChangingTOI/opt/final_orientation.txt', (final_theta, final_phi), delimiter = ' ')

	print(clr.END)

	real_fun(x0)
	
	if (TEST_MODE): print(clr.END + "Single run ended: ")


if __name__ == "__main__":
    sys.exit(main())
