import queue
import sys, os, glob
import subprocess
import shutil, time, datetime
from numpy.ma import ceil

# runs all simulations needed for a given iteration in serial
# - smcTable a parameter file produced by smc.run, usually named smcTable.txt
#   simulations will be run with all parameters specified in this file
# - simDir the name of the directory where the simulation data will be stored
#   most importantly, the data*.txt files should reside in that directory
# - buildDir the directory where the MercuryDPM executables will be built
def runSimulations(smcTable, simDir, buildDir, paramStringTemplate, exeNames, material, samples, onNode=False, nodes=[],
                   cores=0, infrastructure=[], verbose=False):
    # # run make to produce the executables
    # print("Running make in the build directory %s" % buildDir)
    # try:
    #     output = subprocess.check_output('make -j 8 -C ' + buildDir, shell=True)
    # except:
    #     print(output)
    #     raise RuntimeError('Calling make in the buildDirectory failed')

    # make simulation directory absolute
    simDir = os.getcwd() + "/" + simDir
    # build directory is already absolute

    # make the simulation directory if necessary
    if not os.path.exists(simDir):
        if verbose:
            print('Creating directory for simulation output: ' + simDir)
        os.mkdir(simDir)

    # check if params file exists
    if not os.path.exists(smcTable):
        raise RuntimeError('smcTable does not exist')


    # if cores>0 a script is produced, splitting the code into #cores executables
    if cores:
        # todo failure if there are not enough commands per cores
        numCmdTotal = numberOfCommands(smcTable, exeNames)
        numCmd = int(ceil(numCmdTotal/cores))
        cores = int(ceil(numCmdTotal/numCmd))
        print("Number of commands %d, per script %d, number of scripts: %d" % (numCmdTotal, numCmd, cores))
        #override buildDir (only necessary on mercuryCloud, since we cannot run the python script there)
        onMercuryCloud = False
        if onMercuryCloud:
            print("Adjusting script for mercuryCloud")
            buildDir = open('build').read().split()[0]
            print("Copy files to cluster (rsync -avz %s $cloud:) and execute ./run.sh in the sim directory before continuing" % material)
        elif infrastructure == "AWS":
            print("Adjusting script for AWS")
            buildDir = "../../."
        else:
            print("Adjusting script for msm3")
            buildDir = open('build').read().split()[0]
            #buildDir = "~/Code/Lab/build/Drivers/USER/MercuryLab/Buttercup/Calibration"
            print("1) cd %r\n2) Adjust node numbers run.sh.\n3) Execute 'source ./run.sh'" % simDir)
        # set counter into which file the next command is written
        coresCounter = -1
        numCmdCounter = 0
        # write run.sh and create an empty file $coresCounter.sh
        open(simDir + "/run.sh", "w").write("#!/bin/bash\nchmod +x *.sh\n")
        for i in range(cores):
            if onMercuryCloud:
                open(simDir+"/run.sh", "a").write("nohup ./"+str(i)+".sh &\n")
            elif infrastructure == "AWS":
                os.remove(simDir+"/run.sh")
                break
                # open(simDir+"/run.sh", "a").write("nohup ./"+str(i)+".sh 2>&1 &\n")
            else:
                n = nodes[i%len(nodes)]
                if n.isdigit():
                    n = "node%s" % n
                open(simDir + "/run.sh", "a").write("sleep 3\n ssh %s \"cd \"`pwd`\"; nohup ./%d.sh > %s.out\"&\n" % (n,i,i))
        # fix for windows systems
        # cmd = "cd " + simDir + " && sed -i -e 's/\r$//' run.sh"
        # subprocess.Popen(cmd.split(), shell=True)

    # open smcTable
    #print("Opening smcTable: %s" % smcTable)
    file = open(smcTable)
    # ignore header line
    params = file.readline()

    # reading smcTable line by line
    for line in file.readlines():
        # extracting parameters
        params = line.split()[2:]
        # convert list to tuple
        params = tuple(i for i in params)
        # set simulation parameters
        paramString = paramStringTemplate % params
        paramString += "-param _"+material+"_"+"_".join(params)
        #print("Running simulations for: %s" % paramString)

        if cores:
            for executable in exeNames:
                if numCmdCounter == 0:
                    coresCounter += 1
                    open(simDir + str(coresCounter) + ".sh", "w").write("")
                cmd = '%s/%s %s\n' % (buildDir, executable, paramString)
                #print("Writing to %r, %r" % (simDir+str(coresCounter)+".sh", coresCounter))
                open(simDir+str(coresCounter)+".sh", "a").write("#!/bin/sh\n" + cmd)
                numCmdCounter = (numCmdCounter+1)%numCmd
        elif onNode:
            # initialise the freeCores list
            freeCores = getFreeCores(nodes)
            for executable in exeNames:
                cmd = 'cd %s; %s/%s %s' % (simDir, buildDir, executable, paramString)
                # look for new free cores when needed
                while len(freeCores) == 0:
                    print("No free cores, waiting for nodes to free up")
                    # time to wait in seconds between checking for free nodes
                    time.sleep(300)
                    freeCores = getFreeCores(nodes)
                node = freeCores.pop()
                outname = paramString.replace(" ", "")
                sshcmd = "sleep 1\n ssh %s \"%s\"" % (node, cmd + ' &> ' + executable + '_' + outname + '.out &')
                print(sshcmd)
                subprocess.check_output(sshcmd, shell=True)
        else:
            for executable in exeNames:
                cmd = 'cd %s && %s/%s %s' % (simDir, buildDir, executable, paramString)
                if verbose:
                    print("Running in serial: ./%s %s" % (executable, paramString))
                subprocess.check_output(cmd, shell=True)
    if infrastructure == "AWS":
        cwd = os.getcwd()
        startAndWaitForSimulationsToFinish(simDir, samples)
        # Change working directory to simDir to save output files in the right place
        os.chdir(cwd)
    print("Runs finished.")

def runProcess(simDir, id):
    # Output and error logs go to id.out and id.err.
    # When a process is killed/terminated, t.e output file will exist, but is empty.
    cmd = f"{simDir}{id}.sh"
    return subprocess.Popen(cmd)

def startAndWaitForSimulationsToFinish(simDir, samples):

    # Start all processes and add to array.
    processes = []
    cmd = "cd %s && chmod +x *.sh" % simDir
    subprocess.check_output(cmd, shell=True)
    # Change working directory to simDir to save output files in the right place
    os.chdir(simDir)
    for i in range(samples):
        processes.append(runProcess(simDir, i))
        print("Started process with PID %s" % processes[i])
    # Finished processes are removed from the array.
    # Keep looping when there are more than 20% in the array.
    # ratio of simulations to kill when the number of (1-killRatio) * samples has finished.
    killRatio = 0.2
    while len(processes) > killRatio * samples:
        for p in processes:
            # poll() checks if the process is terminated and returns None when this is NOT the case.
            # Remove process from array when finished.
            if p.poll() is not None:
                processes.remove(p)
        # Sleep for a while to only check at certain intervals.
        time.sleep(0.5)

    print(f"At least {(1-killRatio) * 100}% of processes finished.")
    print(f"Killing remaining {len(processes)} processes.")

    # 80% of processes are finished, now kill the remaining processes.
    for p in processes:
        #p.terminate()
        # Since we definitely don't care anymore about the remaining processes, I would use kill.
        p.kill()

    # simulationOutput = []
    # while True:
    #     simulationOutput.clear()
    #     for file in os.listdir(simDir):
    #         if file.endswith(".txt"):
    #             simulationOutput.append(os.path.join(simDir, file))
    #             if len(simulationOutput) >= continueThreshold * samples:
    #                 print("80% of simulations finished, quitting simulations")
    #                 # Kill remaining PIDS HERE to be most efficient
    #                 return

# computes number of commands to be executed
def numberOfCommands(smcTable, exeNames):
    # count lines, ignore header line
    fileCount = open(smcTable)
    numObs = -1
    for line in fileCount.readlines(): numObs += 1
    fileCount.close()
    return len(exeNames) * numObs

def mergeOutputFiles(smcTable, simDir, exeNames, material, verbose=False):
    # open smcTable
    file = open(smcTable)
    # ignore header line
    params = file.readline()
    if verbose:
        print("Writing combined data files in folder %s" % simDir)
    # reading smcTable line by line
    for line in file.readlines():
        # extracting parameters
        params = line.split()[2:]
        # convert list to tuple
        params = tuple(i for i in params)
        # read in output files
        paramStringOut = "_"+material+"_"+"_".join(params)+".txt"
        out=""
        # read files
        for executable in exeNames:
            # skip executables that are labeled -nodata
            if "-nodata" in executable:
                continue
            file = simDir + executable.split()[0] + paramStringOut
            #print("Reading in file: %s" % file)
            # merge the output files into a single file
            try:
                content = open(file).readline()
            except:
                open(file, 'w').write("0 0 0")
                content = "0 0 0"
                # raise Exception(
                # "Directory %s doesn't contain the data files *%s.\nRemove the directory %s and run the simulations again" % (
                # simDir, paramStringOut, simDir))
                # sys.exit(-1)
            out = out + content + " "
        # file that will be written
        outFile = simDir + "data" + paramStringOut
        if verbose:
            print("  data%s: %s" % (paramStringOut, out))
        open(outFile, "w").write(out)

# get number of free cores
def getFreeCores(nodes):
    cwd = os.getcwd()
    freeCores = []
    for i in nodes:
        node = 'node%r' % i
        num = subprocess.check_output('ssh -Y %s python %s/numOfCores.py exit' % (node,cwd), shell=True)
        # use local node
        #num = subprocess.check_output('python %s/numOfCores.py' % (cwd), shell=True)
        print ('%i free nodes in %s' % (int(num), node))
        for i in range(int(num)):
            freeCores.append(node)
    return freeCores
