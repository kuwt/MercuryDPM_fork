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


# This tool parses the *.ene file and generates paraview-friendly plot data

# usage: python PlotEnergies.py <path_to_the_driver_dir> <name_of_the_driver_file>

import sys
import numpy as np


def main():
	path = sys.argv[1]
	name = sys.argv[2]
	print(path)
	print(name)

	textfile = open(path + name + ".ene", "r")
	content_list = textfile.readlines()


	time = np.zeros(len(content_list))
	gra_ene = np.zeros(len(content_list))

	tra_kin_ene = np.zeros(len(content_list))
	rot_kin_ene = np.zeros(len(content_list))
	kin_ene = np.zeros(len(content_list))

	ela_ene = np.zeros(len(content_list))
	tot_ene = np.zeros(len(content_list))

	ratio = np.zeros(len(content_list))

	for i in range(1,len(content_list)):
		A = content_list[i].split(" ")
		B = []
		for j in range(len(A)):
			if A[j] != "": B.append(A[j])
		time[i] = float(B[0])
		gra_ene[i] = float(B[1])
		tra_kin_ene[i] = float(B[2])
		rot_kin_ene[i] = float(B[3])
		kin_ene[i] = float(B[2]) + float(B[3])
		ela_ene[i] = float(B[4])
		tot_ene[i] = float(B[1]) + float(B[2]) + float(B[3])  + float(B[4])


	tra_int = 0
	rot_int = 0
	for i in range(len(time)):
		tra_int += tra_kin_ene[i]
		rot_int += rot_kin_ene[i]
		if not tra_int == 0: ratio[i] = rot_int / tra_int
	
		 
	np.savetxt(path + "paraview_" + name + "/time.txt", time, delimiter=",")
	np.savetxt(path + "paraview_" + name + "/gra_ene.txt", gra_ene, delimiter=",")
	np.savetxt(path + "paraview_" + name + "/tra_kin_ene.txt", tra_kin_ene, delimiter=",")
	np.savetxt(path + "paraview_" + name + "/rot_kin_ene.txt", rot_kin_ene, delimiter=",")
	np.savetxt(path + "paraview_" + name + "/kin_ene.txt", kin_ene, delimiter=",")
	np.savetxt(path + "paraview_" + name + "/ela_ene.txt", ela_ene, delimiter=",")
	np.savetxt(path + "paraview_" + name + "/tot_ene.txt", tot_ene, delimiter=",")
	np.savetxt(path + "paraview_" + name + "/ratio.txt", ratio, delimiter=",")
	print("Successfully written energy files to paraview dir:", path + "paraview_" + name)


if __name__ == "__main__":
    sys.exit(main())
