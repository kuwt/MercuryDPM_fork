#!/usr/bin/env python3
""" Calibrates a MercuryDPM particle model from experimental data.

By default, it calibrates three material parameters (no wall parameters)
 - sliding friction mus
 - rolling friction mur
 - cohesion stiffness kc
The following simulations and experimental results are used for calibration:
 - AngleOfReposeBu
 - DrumBu
 - ShearCellBu

Before you run this code, create a subfolder with the material's name (e.g. "Cornstarch"), and put the data
from experiments into the file PH101/Exp/data.txt. The data file should contain:
 - incipient shear stresses measured at various normal stress in a shear cell
 - angle of repose measured from a heap test
 - dynamic angle of repose measured in a rotating drum
For example, if the shear test was done at four different normal pressures, the file could contain:
    "269 423 609 759 29.7 40.8".

Files:
 - Test/Exp/data.txt contains the experimentally measured data;
   this file has to be in the same format as the simulation data
 - simulateTestCase.py contains the code that produces the simulation data:
   it calls the executables for all the micro-parameter values in a given smcTable.
 - CalibrationTestCase.py this is the main script; execute it to run the calibration procedure
 - Test/Sim* these folders will contain the simulation data; remove old folders before running a calibration
 - smcTable*.tex these files contain all the micro parameters that should be tested in a given iteration
"""

import argparse
from plotResults import *
from simulate import *
import simulateAWS as AWS
from smc import *

# read options
parser = argparse.ArgumentParser()
parser.add_argument("--material", help="Name of the directory with the experimental data", default="", type=str)
parser.add_argument("--iter", help="How many iterations to run", default=2, type=int)
parser.add_argument("--test", help="Runs a simple test instead of a DEM simulation", action='store_true')
parser.add_argument("--analysis", help="Whether to show plots", action='store_true')
parser.add_argument("--overwrite", help="Whether code should overwrite all results of previous run",
                    action='store_true')
parser.add_argument("--verbose", help="Whether smc should be run in verbose mode", action='store_true')
parser.add_argument("--samples", help="Number of samples", default=0, type=int)
parser.add_argument("--cores", help="How many cores we can use", default=0, type=int)
parser.add_argument("--node", help="On which node(s) to run (0=local)", nargs="+", default=[])
parser.add_argument("--infrastructure", help="define the Infrastructure on which you are running",
                    default=[])
args = parser.parse_args()
args.onNode = len(args.node) > 0

# output whether --test has been used
if args.test:
    print("Running a simple test instead of a DEM simulation")

# output whether --node has been specified
if args.onNode:
    print("Running on %r core(s) on node(s) %s" % (args.cores, ", ".join(args.node)))
elif args.cores:
    print("Running on %r cores, splitting simulations among all cores" % args.cores)
else:
    print("Running in serial")

# output whether --material has been specified correctly, i.e. whether calibrateMaterial.txt exists
if not os.path.isfile("calibrate%s.txt" % args.material):
    # find out for which materials a calibrate file exists
    dataFiles = glob.glob('calibrate*.txt')
    materials = [dataFile.replace("calibrate", "").replace(".txt", "") for dataFile in dataFiles]
    # tell user what material values are possible
    sys.exit("Error in interpreting '--material %s':\n File calibrate%s.txt not found.\n Possible arguments: %s."
             % (args.material, args.material, " or ".join(materials)))
print("Material: %s" % args.material)

# read parameters from calibrateMaterial.txt
with open("calibrate%s.txt" % args.material) as paramFile:
    # calibrateMaterial.txt contains sections, each starting with one header line:
    # section 1: fixed material parameters,
    # section 2: variable parameters,
    # section 3: experimental data
    section = 0
    paramNames = []
    paramRanges = {}
    exeNames = []
    obsNames = []  # names of observables
    data = []
    weights = []
    material = args.material.split('_')[0]
    paramString = '-speciesType ' + material + ' '
    # read line by line
    lines = paramFile.readlines()
    # line counter
    c = 0
    while c < len(lines):
        # current line
        line = lines[c]
        c += 1
        if line.startswith('#'):
            section += 1
            continue
        if section == 1:
            # psd's have a special format
            if line.startswith('psd'):
                paramString += '-' + line[:-1] + ' '
                while not lines[c][0].isalpha():
                    line = lines[c]
                    c += 1
                    paramString += line[:-1] + ' '
            else:
                paramString += '-' + line[:-1] + ' '
        elif section == 2:
            val = line.split()
            assert len(val) > 2, "parameter settings should contain name and range: %s" % line
            paramNames.append(val[0])
            paramRanges[val[0]] = [float(val[1]), float(val[2])]
            paramString += '-' + val[0] + ' %s '
        elif section == 3:
            val = line.split()
            assert len(val) > 2, "experimental data should contain exeName, value and weight: %s" % line
            exeNames.append(val[0])
            for i in range(2, len(val), 2):
                data.append(float(val[i - 1]))
                weights.append(float(val[i]))
                obsNames.append(val[0] + (str(int(i / 2)) if (len(val) > 3) else ""))
        else:
            paramString += line + ' '
    # make obs names readable
    par = paramString.split()
    if "-normalStress" in paramString:
        nextCompressiveStress = par.index("-normalStress") + 2
    for i, ob in enumerate(obsNames):
        if ob.startswith("ShearCell"):
            if "-ffc" in paramString:
                obsNames[i] = "FFc"
            else:
                obsNames[i] = "shear stress at %s Pa [Pa]" % par[nextCompressiveStress]
                nextCompressiveStress += 1
        elif ob.startswith("AngleOfRepose"):
            obsNames[i] = "AoR [deg]"
        elif val[0].startswith("Drum"):
            obsNames[i] = "dyn. AoR [deg]"

print("Parameters to be identified: " + ", ".join(paramNames))
print("Parameter ranges: " + ", ".join(str(paramRanges[i]) for i in paramRanges))
print("Executables: " + ", ".join(exeNames))
print("Data values: %s" % data)
print("Weights: %s" % weights)
print("Command line parameters: %s" % paramString)

# check existence of build directory
if args.infrastructure == "AWS":
    buildDir = "../../."
else:
    buildDir = open("build", "r").read().split()[0]
if not os.path.exists(os.path.expanduser(buildDir)):
    raise ValueError("\n Build directory %s does not exist" % os.path.expanduser(buildDir))
print("Build directory for executables: %s" % buildDir)

# deleting content of old output files
if os.path.exists("%sResults.txt" % args.material):
    os.remove("%sResults.txt" % args.material)

# writing data.txt file
if os.path.exists("%s" % args.material):
    shutil.rmtree("%s" % args.material)
if not os.path.exists("%s/Exp" % args.material):
    os.makedirs("%s/Exp" % args.material)
dataFile = "%s/Exp/data.txt" % args.material
open(dataFile, 'w').write(" ".join(str(d) for d in data))

# Covariance 0<sigma<1. Increase if ess is too small
sigma = 1
print("Covariance: %d" % sigma)

# Effective sample size; how much space is covered (should be at least 0.3)
ess = 0.5
print("Goal for effective sample size: %f" % ess)

# Number of Gaussian components
numGaussians = 2
print("Number of Gaussian components: %d" % numGaussians)

# set number of samples
if args.samples == 0:
    args.samples = 10 * len(paramNames)
print("Number of samples: %d" % args.samples)

# set the prior weight (1/numGaussians))
priorWeight = 1.0 / numGaussians
print("Prior weight: %f" % priorWeight)

# dict for analysis data
analysis = []

# iteration of the smc algorithm
for iteration in range(args.iter + 1):

    # define name of output directory, smc table
    outputDir = "%s/Sim%d/" % (args.material, iteration)
    smcTable = 'Tables/smcTable_%s_%i.txt' % (args.material, iteration)

    # delete output directory and its contents
    if args.overwrite and os.path.isdir(outputDir):
        shutil.rmtree(outputDir)

    # if output directory exists, but has no data files, tell the user to run simulations and exit
    # if output directory exists, but has data files, jump to next step
    # else, prepare and do iteration
    if os.path.isdir(outputDir):
        # finish if current iteration directory exists, but has not been executed yet
        dataFiles = glob.glob("%s/*.txt" % outputDir)
        if dataFiles:
            print('\nIteration %d has already been run, continuing with next step' % iteration)
        else:
            print("Iteration %d has been setup, but not yet run. "
                  "Call 'cd %s && source run.sh' to continue the calibration" % (iteration, outputDir))
            raise SystemExit
    else:
        # create smcTable
        print("\nCreating %s. \nOutput will be written to %s." % (smcTable, outputDir))
        if iteration == 0:
            # On the first iteration, create an initial smcTable that does not depend on previous data
            smcTest = smc(sigma, ess, weights, DPMDataDir=outputDir, obsDataFile=dataFile, loadSamples=False,
                          standAlone=True, verbose=args.verbose)
            smcTest.initialize(paramNames, paramRanges, args.samples, numGaussians, priorWeight)
            if not os.path.exists('Tables'):
                os.mkdir('Tables')
            os.system('mv smcTable.txt %s' % smcTable)
        else:
            # On all other iterations, create an smcTable that does depend on the previous data
            # noinspection PyUnboundLocalVariable
            gmm, maxNumComponents = smcTest.resampleParams(caliStep=-1)
            os.system('mv smcTableNew.txt %s' % smcTable)
            # write out gmm properties
            open(smcTable + ".gmm", "w").write("Means %r\n Covariances\n%r" % (gmm.means_, gmm.covariances_))
            np.save(smcTable + ".means", gmm.means_)
            np.save(smcTable + ".covs", gmm.covariances_)
            np.save(smcTable + ".weights", gmm.weights_)

        # run simulations
        print("Preparing simulations. After simulations have finished, rerun this script to continue.")
        if args.infrastructure == "AWS":
            AWS.runSimulations(smcTable, outputDir, buildDir, paramString, exeNames, args.material, args.samples,
                               args.onNode, args.node, args.cores, args.infrastructure, args.verbose)
        else:
            runSimulations(smcTable, outputDir, buildDir, paramString, exeNames, args.material, args.onNode,
                           args.node, args.cores, args.verbose)

        # if run on the cluster, stop after each iteration (because simulations need to finish)
        # if run on AWS we open a pipe (Popen) which waits for the simulations to finish. So the program itself waits
        # for all simulations to finish.
        if args.cores > 0 and args.infrastructure != "AWS":
            sys.exit(0)

    # creating data_*.txt files in outputDir from DEM output files
    # for exe in exeNames:
    #     if not os.path.isfile(exe):
    #         exeNames.remove(exe)
    # if on AWS merge unfinished simulation outputs with "0 0 0"
    if args.infrastructure == "AWS":
        AWS.mergeOutputFiles(smcTable, outputDir, exeNames, material, args.verbose)
    else:
        mergeOutputFiles(smcTable, outputDir, exeNames, material, args.verbose)

    # load parameter table \todo can this be moved to ln 175?
    smcTest = smc(sigma, ess, weights, DPMDataDir=outputDir, obsDataFile=dataFile, loadSamples=True, standAlone=True,
                  simName="data", verbose=args.verbose)
    smcTest.initialize(paramNames, paramRanges, args.samples, numGaussians, priorWeight)

    # run sequential Monte Carlo; returns posterior mean and coefficient of variance
    print("Starting analysis of iteration %r" % iteration)
    smcTest.run(iterNO=iteration)

    # print output data
    print("The optimal parameters for %s in iteration %d are:" % (args.material, iteration))
    posterior = smcTest.getPosterior()
    smcSamples = smcTest.getSmcSamples()[0]
    dpmData = smcTest.DPMData[-1]
    obsData = smcTest.obsData[-1]
    ess = smcTest.getEffectiveSampleSize()
    # mean and averages
    avgValues = smcTest.ips.T[0]
    covariances = smcTest.covs.T[0]
    # best
    bestId = np.argmax(posterior)
    bestValues = smcSamples[bestId]
    bestPosterior = posterior[bestId]
    bestData = dpmData[bestId]
    # compute mean error of best data
    error = 0
    for i in range(len(bestData)):
        error += (bestData[i] / obsData[i] - 1) ** 2
    error = np.sqrt(error / len(bestData))
    bestValuesStr = ', '.join(str(paramNames[i]) + ": " + str(bestValues[i]) for i in range(len(bestValues)))
    bestDataStr = ', '.join("%s: %g (exp: %g)" % (obsNames[i], bestData[i], obsData[i]) for i in range(len(bestData)))
    print("Best params for %s in iteration %d (error %f):\n %s\n data: %s"
          % (args.material, iteration, error, bestValuesStr, bestDataStr))

    print("Means and Covariances:")
    for i in range(len(paramNames)):
        print("%30s = mean %f with covariance %f" % (paramNames[i], avgValues[i], covariances[i]))
    analysis.append(Analysis(smcTest, smcSamples, obsNames, paramNames, posterior, bestId))

    # write output data file
    if args.infrastructure == "AWS":
        with open("%sResults.txt" % args.material, 'a+') as outputFile:

            outputFile.write("Best params for %s in iteration %d (error %f):\n  %s\n  data: %s\n"
                             % (args.material, iteration, error, bestValuesStr, bestDataStr))
            outputFile.write("Means and Covariances:\n")
            for i in range(len(paramNames)):
                outputFile.write("%30s = mean %f with covariance %f\n\n" % (paramNames[i], avgValues[i], covariances[i]))

# plot output data
if args.analysis:
    plotParametersAndObservables(analysis, args.material)
    # plotObservables(smcTest, obsNames)
    # plotPosteriorWeights(paramNames, posterior, smcSamples[0],avgValues)
    # plotCorrelations(smcTest, paramNames, obsNames)

print("Finished Calibration")
