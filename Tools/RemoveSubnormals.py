#!/usr/bin/env python3
#Copyright (c) 2013-2020, The MercuryDPM Developers Team. All rights reserved.
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

import sys
import os

'''
Removes subnormal numbers (|x|<1.2e-38) from output files
https://en.cppreference.com/w/cpp/language/types
'''


def main():
    if len(sys.argv) < 2:
        raise Exception("You need to provide a filename as command line argument")
    # get name
    file_name = sys.argv[1]
    print("Input file: %s" % file_name)
    # open
    input = open(file_name, 'r')
    lines = input.readlines()
    input.close()
    words = [line.split() for line in lines]
    for i, wordList in enumerate(words):
        for j, word in enumerate(wordList):
            try:
                if 0 < abs(float(word)) < 1.2e-38:
                    print("subnormal %r %s" % (float(word), word))
                    words[i][j] = "0"
            except:
                pass

    # open
    print("Input file: %s.out" % file_name)
    output = open(file_name+'.out', 'w')
    lines = [" ".join(word) for word in words]
    text = "\n".join(lines)
    output.write(text)

if __name__ == '__main__':
    main()
