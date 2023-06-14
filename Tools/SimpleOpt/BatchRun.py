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
import numpy as np
import subprocess
from src.color_lib import colorClass


def main():
	
	N = 100 # Number of points in a batch 
	global source_dir, build_dir

	# Set up terminal output
	clr = colorClass()
	clr.set_use_colors(True)
	clr.def_colors()

	# Hello tag
	print(clr.END + clr.BOLD + clr.VIOLET + "-------------------------------------------")
	print(clr.END + clr.BOLD + clr.VIOLET + "     MercuryDPM BatchRun tool      ")
	print(clr.END + clr.BOLD + clr.VIOLET + "-------------------------------------------")
	
	if (len(sys.argv) == 3):  # Full format: "run source_dir build_dir"
		source_dir = sys.argv[1]
		build_dir = sys.argv[2]

	print("Source dir: ", source_dir)
	print("Build dir: ", build_dir)
	
	
	for i in range(100):
		command  = [build_dir + '/Drivers/Clump/Domino']
		command.append('-p0')
		command.append(str(i))
		print(command)
		subprocess.run(command)	   
	
if __name__ == "__main__":
    sys.exit(main())
