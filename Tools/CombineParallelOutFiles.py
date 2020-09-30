#!/usr/bin/env python3
import sys
import os
import operator

'''
Short script to combine out files that are generated in the parallel code, to a single out file with all
particles sorted by time (column one).

Usage: To combine file.out0, file.out1, etc into file.out:
    ./CombineParallelOutFiles.py file.out
'''


def main():
    # checks if the combined outfile already exists
    out_file_name = sys.argv[1]

    #if os.path.isfile(out_file_name):
    #    raise Exception("file " + out_file_name + " already exists; remove first")

    # find the out files of all the cores
    files = sorted([file for file in os.listdir(".") if file.startswith(out_file_name) and not file==out_file_name])

    print("Combining %d out files: " % len(files) + ", ".join(files))

    # For all files
    particles=[]
    header = ""
    for file in files:
        # Open the file and mold it into a list of lists of all the 'words' in the file
        file0 = open(file)
        file0 = file0.readlines()
        # ignore empty files
        if not file0:
            continue
        if not header:
            header=file0[0]
        file0.pop(0)
        file0.pop()
        particles.extend([line.split() for line in file0])
    print("found %d particles (including duplicates)" % len(particles))

    # sort by id
    particles = sorted(particles, key=operator.itemgetter(2))
    # remove duplicates (same id)
    id = -1
    for i in range(len(particles) - 2, -1, -1):
        # check if id and time stamp agrees
        if particles[i][2]==particles[i+1][2]:
            # remove it
            particles.pop(i+1)
    print("found %d particles (after removing duplicates)" % len(particles))

    # sort by time stamp (convert time to float, otherwise it sorts alphabetically)
    for p in particles:
        p[0] = float(p[0])
    particles2 = sorted(particles, key=operator.itemgetter(0))
    for p in particles:
        p[0] = str(p[0])

    # Finally, write everything to the combined out file
    out_file = open(out_file_name, 'w')
    out_file.write("".join(header))
    out_file.write("\n".join([" ".join(p) for p in particles2]))
    print("written %s" % out_file_name)


if __name__ == '__main__':
    main()
