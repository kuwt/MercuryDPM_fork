#Copyright (c) 2013-2023, The MercuryDPM Developers Team. All rights reserved.
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
import numpy as np
from src.color_lib import colorClass
from scipy.optimize import minimize


# Global parameters
source_dir = ""				# Mercury source dir
build_dir = ""				# Mercury Build dir
N = 2					# Number of parameters our functional depends upon
TEST_MODE = True			# True for self-test on toy functional, False for real simulation-guided optimization


x0 = np.ones([N]) # Initial guess for the solution
xc = np.random.rand(N) # Known solution (for toy functional only)
xc = np.array([ 1.57055293e+00, -2.51228303e-08]) # Known solution for the SimpleOpt example


def real_fun(x):
	import subprocess
	global source_dir, build_dir
	ret = 0
	
	command  = [build_dir+'/Drivers/SimpleOpt/Opt']
	
	for i in range(len(x)):
		command.append('-p'+str(i))
		command.append(str(x[i]))
		
	
	#print(command)
	
	subprocess.run(command)	   
	f = open(build_dir+'/Drivers/SimpleOpt/functional.txt', "r")
	ret = float(f.read())
	return ret


def toy_fun(x):   
	global xc
	ret = 0
	for m in range(len(x)):
		ret = ret + (x[m] - xc[m])**2
	return ret


def main():
	global source_dir, build_dir

	# Set up terminal output
	clr = colorClass()
	clr.set_use_colors(True)
	clr.def_colors()

	# Hello tag
	print(clr.END + clr.BOLD + clr.VIOLET + "-------------------------------------------")
	print(clr.END + clr.BOLD + clr.VIOLET + "     MercuryDPM SimpleOpt tool      ")
	print(clr.END + clr.BOLD + clr.VIOLET + "-------------------------------------------")
	
	if (len(sys.argv) == 3):  # Full format: "run source_dir build_dir"
		source_dir = sys.argv[1]
		build_dir = sys.argv[2]

	print("Source dir: ", source_dir)
	print("Build dir: ", build_dir)
	
	res = minimize(real_fun, x0, method='Powell', options = {'xtol':1e-3, 'ftol':1e-4})

	print(res)

	print(clr.END + "Converged to the known solution: ", np.allclose(res["x"], xc))


if __name__ == "__main__":
    sys.exit(main())
