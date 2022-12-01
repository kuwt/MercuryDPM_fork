""" Sequential Monte Carlo: the core of the calibration toolbox
    (https://en.wikipedia.org/wiki/Particle_filter)
"""

import sys, os, glob
import numpy as np
import pickle
from tools import *
#from sklearn import mixture
#from collision import createScene, runCollision, addSimData
from itertools import repeat


class smc:
    """
    SMC base class: sequential Monte Carlo (SMC) filtering
    """

    def __init__(self, sigma, ess, obsWeights, DPMFile='', DPMDataDir='', obsDataFile='', obsCtrl='', simDataKeys='',
                 simName='sim', DPMVersion='yade-batch', scaleWithMax=True, loadSamples=True, skipDEM=True,
                 standAlone=True, verbose=False):
        if verbose:
            print("Called smc.init(sigma=%s, ess=%s, obsWeights=%s, DPMFile=%s, DPMDataDir=%s, obsDataFile=%s, "
                  "obsCtrl=%s, simDataKeys=%s, simName=%s, DPMVersion=%s, scaleWithMax=%s, loadSamples=%s, "
                  "skipDEM=%s, standAlone=%s)"
                  % (sigma, ess, obsWeights, DPMFile, DPMDataDir, obsDataFile, obsCtrl, simDataKeys, simName,
                     DPMVersion, scaleWithMax, loadSamples, skipDEM, standAlone))
        # whether additional debug information is printed
        self.verbose = verbose
        # ?
        self.DPMVersion = DPMVersion
        # ?
        self.DPMFile = DPMFile

        # File in which experimental data is stored (\todo should not be class variable)
        self.obsDataFile = obsDataFile
        # used in getObsDataFromFile to determine how to interpret the data (we always use the default, '') (\todo should not be class variable)
        self.obsCtrl = obsCtrl

        # Directory in which simulation data is stored
        self.DPMDataDir = DPMDataDir
        # ?
        self.numParams = 0
        # ?
        self.numSamples = 0
        # ?
        self.numSteps = 0
        # normalized variance
        self.__sigma = sigma
        # ?
        self.__ess = ess
        # ?
        self.obsWeights = obsWeights
        # ?
        self.simDataKeys = simDataKeys
        # ?
        self.simName = simName

        # obsData: nparray of size [numObs, numSteps], containing the experimental data
        # obsCtrlData: ? (only needed if obsCtrl != '')
        # numObs: Number of observed variables in data file
        # numSteps: Number of steps per variables in data file
        self.obsData, self.obsCtrlData, self.numObs, self.numSteps = self.getObsDataFromFile(obsDataFile, obsCtrl)
        if self.verbose:
            print("Read in obsData=%s, obsCtrlData=%s, numObs=%s, numSteps=%s" % (self.obsData, self.obsCtrlData, self.numObs, self.numSteps))

        # assume all measurements are independent
        self.__obsMatrix = np.identity(self.numObs)
        if self.verbose:
            print("Set obsMatrix=%s"  % self.__obsMatrix)

        # File names of the smcTables \todo (why multiple?)
        self.paramsFiles = []
        # ??
        self.paramNames = []
        # ??
        self.paramRanges = {}
        # ??
        self.smcSamples = []
        # a flag that defines whether Yade is called within python
        self.standAlone = standAlone
        # if run Bayesian filtering on the fly
        if not self.standAlone:
            self.__pool = None
            self.__scenes = None
        # hyper-parameters of Bayesian Gaussian mixture
        self.__maxNumComponents = 0
        self.__priorWeight = 0

        # simulation data: nparray of size [numSteps, numSamples, numObs]
        self.DPMData = None
        # ensemble averages : nparray of size [numParams, numSteps]
        self.ips = None
        # coefficients of variance : nparray of size [numParams, numSteps]
        self.covs = None
        # posterior : nparray of size [numSamples, numSteps]
        self.posterior = None
        # likelihood : nparray of size [numSamples, numSteps]
        self.likelihood = None
        # proposal: nparray of size [numSamples]
        self.proposal = None

        # whether change covariance matrix in time or set it constant (proportional to the maximum observation data)
        self.scaleCovWithMax = scaleWithMax
        # whether load existing parameter samples or create new ones
        self.loadSamples = loadSamples
        # whether run simulations DEM within GrainLearning
        self.skipDEM = skipDEM

    def initialize(self, paramNames, paramRanges, numSamples, maxNumComponents, priorWeight, paramsFile='',
                   proposalFile='', threads=4):
        if self.verbose:
            print("Called smc.initialize(paramNames=%s, paramRanges=%s, numSamples=%s, maxNumComponents=%s, priorWeight=%s, paramsFile=%s, proposalFile=%s, threads=%s)"
                  % (paramNames, paramRanges, numSamples, maxNumComponents, priorWeight, paramsFile, proposalFile, threads))
        # assign parameter names
        self.paramNames = paramNames
        # initialize parameters for Gaussian mixture model
        self.__maxNumComponents = maxNumComponents
        self.__priorWeight = priorWeight
        # initialize sample uniformly if no parameter table available
        if len(self.smcSamples) == 0:
            if self.loadSamples:
                self.numSamples, _ = self.getParamsFromTable(paramsFile, paramNames, paramRanges)
            else:
                self.numSamples = self.getInitParams(paramRanges, numSamples, threads)
        if self.verbose:
            print("Set numParams=%s, paramsFiles=%s, numSamples=%s, smcSamples=%s"
                  % (self.numParams, self.paramsFiles, self.numSamples, self.smcSamples))
            #/todo where is the smcTable.txt file written?
        # initialise simulation data
        self.DPMData = np.zeros([self.numSteps, self.numSamples, self.numObs])
        # identified optimal parameters
        self.ips = np.zeros([self.numParams, self.numSteps])
        self.covs = np.zeros([self.numParams, self.numSteps])
        self.posterior = np.zeros([self.numSamples, self.numSteps])
        self.likelihood = np.zeros([self.numSamples, self.numSteps])
        self.proposal = np.ones([self.numSamples]) / self.numSamples
        if proposalFile != '':
            # load proposal density from file
            self.proposal = self.loadProposalFromFile(proposalFile, -1)
        if self.verbose:
            print("Set proposal=%s" % self.proposal)
        # if run Bayesian filtering on the fly
        if not self.standAlone:
            self.__pool = get_pool(mpi=False, threads=self.numSamples)
            self.__scenes = self.__pool.map(createScene, range(self.numSamples))

    # \todo this function seems never to be used?
    # def getProposalFromSamples(self, iterNO):
    #     if len(self.getSmcSamples()) == 0:
    #         RuntimeError("SMC samples not yet loaded...")
    #     else:
    #         gmm = mixture.BayesianGaussianMixture(n_components=self.__maxNumComponents,
    #                                               weight_concentration_prior=self.__priorWeight, covariance_type='full',
    #                                               tol=1e-5, max_iter=int(1e5), n_init=100)
    #         gmm.fit(self.getSmcSamples()[iterNO])
    #         proposal = np.exp(gmm.score_samples(self.getSmcSamples()[iterNO]))
    #     return proposal / sum(proposal)

    def loadProposalFromFile(self, proposalFile, iterNO):
        if len(self.getSmcSamples()) == 0:
            RuntimeError("SMC samples not yet loaded...")
        else:
            proposalModel = pickle.load(open(proposalFile, 'rb'))
            proposal = np.exp(proposalModel.score_samples(self.getSmcSamples()[iterNO]))
        return proposal / sum(proposal)

    def run(self, iterNO=-1, reverse=False, threads=1):
        if self.verbose:
            print("Called run(iterNO=%s, reverse=%s, threads=%s)" % (iterNO, reverse, threads))
        # if iterating, reload observation data
        if iterNO > 0:
            self.obsData, self.obsCtrlData, self.numObs, self.numSteps = \
                self.getObsDataFromFile(self.obsDataFile, self.obsCtrl)
            if self.verbose:
                print("Read in obsData=%s, obsCtrlData=%s, numObs=%s, numSteps=%s" % (
                self.obsData, self.obsCtrlData, self.numObs, self.numSteps))

        # if use Bayesian filtering as a stand alone tool (data already exist before hand)
        if self.standAlone:
            # if run DEM simulations now, with the new parameter table
            if not self.skipDEM and not self.loadSamples:
                # run DEM simulations in batch.
                raw_input("*** Press confirm if the DPm file name is correct... ***\n" + self.DPMFile
                          + "\nAbout to run DPM in batch mode with " +
                          ' '.join([self.DPMVersion, self.paramsFiles[iterNO], self.DPMFile]))
                os.system(' '.join([self.DPMVersion, self.paramsFiles[iterNO], self.DPMFile]))
                print('All simulations finished')
            # if run DEM simulations outside, with the new parameter table
            elif self.skipDEM and not self.loadSamples:
                print('Leaving GrainLearning... only the parameter table is generated')
                sys.exit()
            # if process the simulation data obtained with the existing parameter table
            else:
                print('Skipping DEM simulations from GrainLearning side, reading in existing simulation data now')
            # get simulation data from DPMDataDir
            DPMDataFiles = glob.glob(os.getcwd() + '/' + self.DPMDataDir + '/' + self.simName + '*')
            DPMDataFiles.sort()
            # if glob.glob(self.DPMDataDir + '/*_*txt*') does not return the list of DPM data file
            while len(DPMDataFiles) == 0:
                self.simName = raw_input("No DEM filename can be found, tell me the simulation name...\n ")
                DPMDataFiles = glob.glob(os.getcwd() + '/' + self.DPMDataDir + '/' + self.simName + '*')
                DPMDataFiles.sort()
            # check if parameter combinations match with the simulation filename.
            for i, f in enumerate(DPMDataFiles):
                # get the file name fore the suffix
                f = f.split('.' + f.split('.')[-1])[0]
                # get parameters from the remaining string
                paramsString = f.split('_')[-self.numParams:]
                # element wise comparison of the parameter vector
                if not (np.equal(np.float64(paramsString), self.getSmcSamples()[-1][i]).all()):
                    raise RuntimeError(
                        "Parameters " + ", ".join(
                            ["%s" % v for v in self.getSmcSamples()[-1][i]]) + " do not match with data file " + f)
            # read simulation data into DPMData and drop keys in obsData
            self.getDPMData(DPMDataFiles)
            # if iteration number is an odd number, reverse the data sequences to ensure continuity
            if reverse:
                self.obsCtrlData = self.obsCtrlData[::-1]
                self.obsData = self.obsData[::-1, :]
                self.DPMData = self.DPMData[::-1, :, :]
                if self.verbose:
                    print("Reversing data sequence")
            # loop over data assimilation steps
            for i in range(self.numSteps):
                self.likelihood[:, i], self.posterior[:, i], \
                self.ips[:, i], self.covs[:, i] = self.recursiveBayesian(i)
            # iterate if effective sample size is too big
            while self.getEffectiveSampleSize()[-1] > self.__ess:
                self.__sigma *= 0.9
                for i in range(self.numSteps):
                    self.likelihood[:, i], self.posterior[:, i], \
                    self.ips[:, i], self.covs[:, i] = self.recursiveBayesian(i)
        # if perform Bayesian filtering while DEM simulations are running
        else:
            exit()
            # parameter list
            paramsList = []
            for i in range(self.numSamples):
                paramsForEach = {}
                for j, name in enumerate(self.paramNames):
                    paramsForEach[name] = self.smcSamples[iterNO][i][j]
                paramsList.append(paramsForEach)
            # pass parameter list to simulation instances FIXME: runCollision is the user-defined DPM script
            simData = self.__pool.map(runCollision, zip(self.__scenes, paramsList, repeat(self.obsCtrlData)))
            self.__pool.close()
            # ~ data = runCollision([self.smc__scenes,paramsList[0]])
            # get observation and simulation data ready for Bayesian filtering
            self.obsData = np.array([self.obsData[name] for name in self.simDataKeys]).transpose()
            for i, data in enumerate(simData):
                for j, name in enumerate(self.simDataKeys):
                    # ~ print len(data[name]),i
                    self.DPMData[:, i, j] = data[name]
            # ~ print np.linalg.norm(data[self.obsCtrl]-self.obsCtrlData)
            # loop over data assimilation steps
            if reverse:
                self.obsCtrlData = self.obsCtrlData[::-1]
                self.obsData = self.obsData[::-1, :]
                self.DPMData = self.DPMData[::-1, :, :]
            for i in range(self.numSteps):
                self.likelihood[:, i], self.posterior[:, i], \
                self.ips[:, i], self.covs[:, i] = self.recursiveBayesian(i)
            # iterate if effective sample size is too big
            while self.getEffectiveSampleSize()[-1] > self.__ess:
                self.__sigma *= 0.9
                for i in range(self.numSteps):
                    self.likelihood[:, i], self.posterior[:, i], \
                    self.ips[:, i], self.covs[:, i] = self.recursiveBayesian(i)
        if self.verbose:
            print("Called recursiveBayesian to set "
                  "effectiveSampleSize=%s, sigma=%s, likelihood=%s, posterior=%s, ips=%s, covs=%s"
                  % (self.getEffectiveSampleSize(), self.__sigma, self.likelihood, self.posterior, self.ips,
                     self.covs))
        return self.ips, self.covs

    def getDPMData(self, DPMDataFiles):
        if self.verbose:
            print("Called getDPMData(DPMDataFiles=%s)" % DPMDataFiles)
        if 0 in self.DPMData.shape:
            raise RuntimeError ("number of Observations, samples or steps undefined!")
        for i, f in enumerate(DPMDataFiles):
            # if do not know the data to control simulation
            if self.obsCtrl == '':
                DPMData = np.genfromtxt(f)
                if len(DPMData.shape) == 0:
                    DPMData = DPMData.reshape([1, 1])
                elif len(DPMData.shape) == 1:
                    DPMData = DPMData.reshape([DPMData.shape[0], 1]).T
                for j in range(self.numObs):
                    self.DPMData[:, i, j] = DPMData[:,j]
            else:
                DPMData = getKeysAndData(f)
                for j, key in enumerate(self.obsData.keys()):
                    self.DPMData[:, i, j] = DPMData[key]
        if self.verbose:
            print("Set DPMData=%s" % self.DPMData)
        # if need to remove the data that controls the simulations
        if self.obsCtrl != '':
            obsData = np.zeros([self.numSteps, self.numObs])
            for j, key in enumerate(self.obsData.keys()):
                obsData[:, j] = self.obsData[key]
            self.obsData = obsData
            if self.verbose:
                print("Set obsData=%s" % self.obsData)

    def recursiveBayesian(self, caliStep, iterNO=-1):
        likelihood = self.getLikelihood(caliStep)
        posterior = self.update(caliStep, likelihood)
        # get ensemble averages and coefficients of variance
        ips = np.zeros(self.numParams)
        covs = np.zeros(self.numParams)
        for i in range(self.numParams):
            # ensemble average
            ips[i] = self.smcSamples[iterNO][:, i].dot(posterior)
            # diagonal variance
            covs[i] = ((self.smcSamples[iterNO][:, i] - ips[i]) ** 2).dot(posterior)
            # get coefficient of variance cov
            covs[i] = np.sqrt(covs[i]) / ips[i]
        return likelihood, posterior, ips, covs

    #calistep: index of output parameter
    def getLikelihood(self, caliStep):
        # state vector y_t = H(x_t)+Sigma_t: the values outputted by DEM
        stateVec = self.DPMData[caliStep, :, :].dot(self.__obsMatrix)
        # obsVec: the values we aim for
        obsVec = self.obsData[caliStep, :]
        # row-wise substraction obsVec[numObs]-stateVec[numSamples,numObs], the error
        vecDiff = obsVec - stateVec
        Sigma = self.getCovMatrix(caliStep, self.obsWeights)
        invSigma = np.linalg.inv(Sigma)
        likelihood = np.zeros(self.numSamples)
        # compute likelihood = exp(-0.5*(y_t-H(x_t))*Sigma_t^{-1}*(y_t-H(x_t)))
        for i in range(self.numSamples):
            power = (vecDiff[i, :]).dot(invSigma.dot(vecDiff[i, :].T))
            likelihood[i] = np.exp(-0.5 * power)
        # regularize likelihood
        likelihood /= np.sum(likelihood)
        return likelihood

    def update(self, caliStep, likelihood):
        # update posterior probability according to Bayes' rule
        posterior = np.zeros(self.numSamples)
        if caliStep == 0:
            posterior = likelihood / self.proposal
        else:
            posterior = self.posterior[:, caliStep - 1] * likelihood
        # regularize likelihood
        posterior /= np.sum(posterior)
        return posterior

    def getCovMatrix(self, caliStep, weights):
        Sigma = np.zeros([self.numObs, self.numObs])
        # scale observation data with normalized variance parameter to get covariance matrix
        for i in range(self.numObs):
            # give smaller weights for better agreement
            if self.scaleCovWithMax:
                Sigma[i, i] = self.__sigma * weights[i] * max(self.obsData[:, i]) ** 2
            else:
                Sigma[i, i] = self.__sigma * weights[i] * self.obsData[caliStep, i] ** 2
        return Sigma

    def getInitParams(self, paramRanges, numSamples, threads):
        if self.verbose:
            print("Called smc.getInitParams(paramRanges=%s, numSamples=%s, threads=%s)" % (paramRanges, numSamples, threads))
        if len(paramRanges.values()) == 0:
            raise RuntimeError(
                "Parameter range is not given. Define the dictionary-type variable paramRanges or loadSamples True")
        self.paramRanges = paramRanges
        self.numParams = len(self.getNames())
        names = self.getNames()
        minsAndMaxs = np.array([paramRanges[key] for key in self.getNames()])
        mins = minsAndMaxs[:, 0]
        maxs = minsAndMaxs[:, 1]
        #print ('Parameters to be identified:', ", ".join(names), '\nMins:', mins, '\nMaxs:', maxs)
        initSmcSamples, initparamsFile = initParamsTable(keys=names, maxs=maxs, mins=mins, num=numSamples,
                                                         threads=threads)
        #print ('Written parameter table name:',initparamsFile )
        self.smcSamples.append(np.array(initSmcSamples))
        self.paramsFiles.append(initparamsFile)
        return numSamples

    def getParamsFromTable(self, paramsFile, names, paramRanges, iterNO=-1):
        self.paramRanges = paramRanges
        self.numParams = len(self.getNames())
        if self.verbose:
            print("Called smc.getParamsFromTable(paramsFile=%s, names=%s, paramRanges=%s, iterNO=%s)"
                % (paramsFile, names, paramRanges, iterNO))
        # if paramsFile exist, read in parameters and append data to self.smcSamples
        if os.path.exists(paramsFile):
            self.paramsFiles.append(paramsFile)
            if self.verbose:
                print("Set paramsFiles=%s" % self.paramsFiles)
            smcSamples = np.genfromtxt(paramsFile, comments='!')[:, -self.numParams:]
            self.smcSamples.append(smcSamples)
            return self.smcSamples[iterNO].shape
        # if cannot find paramsFile, get parameter values from the names of simulation files
        else:
            if self.verbose:
                print ('Cannot find paramsFile... Now get it from the names of simulation files')
            DPMDataFiles = glob.glob(os.getcwd() + '/' + self.DPMDataDir + '/' + self.simName + '*')
            DPMDataFiles.sort()
            while len(DPMDataFiles) == 0:
                self.simName = raw_input("No DEM filename can be found, tell me the simulation name...\n ")
                DPMDataFiles = glob.glob(os.getcwd() + '/' + self.DPMDataDir + '/' + self.simName + '*')
                DPMDataFiles.sort()
            smcSamples = np.zeros([len(DPMDataFiles), len(self.getNames())])
            for i, f in enumerate(DPMDataFiles):
                # get the suffix of the data file name
                suffix = '.' + f.split('.')[-1]
                # get the file name fore the suffix
                f = f.split(suffix)[0]
                # get parameters from the remaining string
                paramsString = f.split('_')[-self.numParams:]
                # if the number of parameters in the string is different from self.numParams, throw a warning
                if len(f.split('_')[2:]) != self.numParams:
                    RuntimeError(
                        "Number of parameters found from the simulation fileName is different from self.numParams\n\
                         Make sure your simulation file is named as 'simName_key_param0_param1_..paramN")
                smcSamples[i, :] = np.float64(paramsString)
            self.smcSamples.append(smcSamples)
            return self.smcSamples[iterNO].shape

    def getObsDataFromFile(self, obsDataFile, obsCtrl):
        if self.verbose:
            print("Called getObsDataFromFile(obsDataFile=%s, obsCtrl=%s)" % (obsDataFile, obsCtrl))
        # if do not know the data to control simulation
        if self.obsCtrl == '':
            keysAndData = np.genfromtxt(obsDataFile)
            # if only one observation data vector exist, reshape it with [numSteps,1]
            if len(keysAndData.shape) == 0:
                keysAndData = keysAndData.reshape(1, 1).T
            elif len(keysAndData.shape) == 1:
                keysAndData = keysAndData.reshape([keysAndData.shape[0], 1]).T
            return keysAndData, None, keysAndData.shape[1], keysAndData.shape[0]
        else:
            keysAndData = getKeysAndData(obsDataFile)
            # separate obsCtrl for controlling simulations from obsData
            obsCtrlData = keysAndData.pop(obsCtrl)
            numSteps = len(obsCtrlData)
            numObs = len(keysAndData.keys())
            return keysAndData, obsCtrlData, numObs, numSteps

    def resampleParams(self, caliStep, thread=4, iterNO=-1):
        if self.verbose:
            print("Called resampleParams(caliStep=%s, thread=%s, iterNO=%s)" % (caliStep, thread, iterNO))
        names = self.getNames()
        smcSamples = self.smcSamples[iterNO]
        numSamples = self.numSamples
        # posterior at caliStep is used as the proposal distribution
        proposal = self.posterior[:, caliStep]
        newSmcSamples, newparamsFile, gmm, maxNumComponents = \
            resampledParamsTable(keys=names, smcSamples=smcSamples, proposal=proposal, num=numSamples, thread=thread,
                                 maxNumComponents=self.__maxNumComponents, priorWeight=self.__priorWeight,verbose=self.verbose)
        #print('Parameter table name',newparamsFile)
        self.smcSamples.append(newSmcSamples)
        self.paramsFiles.append(newparamsFile)
        if self.verbose:
            print("Set smcSamples=%s, self.paramsFiles=%s, gmm=%s, maxNumComponents=%s" % (self.smcSamples, self.paramsFiles, gmm, maxNumComponents))
        return gmm, maxNumComponents

    def getPosterior(self):
        return self.posterior

    def getSmcSamples(self):
        return self.smcSamples

    def getNumSteps(self):
        return self.numSteps

    def getEffectiveSampleSize(self):
        nEff = 1. / sum(self.getPosterior() ** 2)
        return nEff / self.numSamples

    def getNames(self):
        return self.paramNames

    def getObsData(self):
        return np.hstack((self.obsCtrlData.reshape(self.numSteps, 1), self.obsData))

    def trainGMMinTime(self, maxNumComponents, iterNO=-1):
        gmmList = []
        smcSamples = self.smcSamples[iterNO]
        for i in range(self.numSteps):
            print ('Train DP mixture at time %i...' % i)
            posterior = self.posterior[:, i]
            gmmList.append(getGMMFromPosterior(smcSamples, posterior, maxNumComponents))
        return gmmList

    def removeDegeneracy(self, caliStep=-1):
        effIDs = np.where(self.posterior[:, caliStep] < 0.5)[0]
        self.proposal = self.proposal[effIDs, :]
        self.likelihood = self.likelihood[effIDs, :]
        self.posterior = self.posterior[effIDs, :]
        self.smcSamples[0] = self.smcSamples[0][effIDs, :]
        self.DPMData = self.DPMData[:, effIDs, :]
        self.numSamples = len(effIDs)
        for i in range(self.numSteps):
            self.likelihood[:, i], self.posterior[:, i], \
            self.ips[:, i], self.covs[:, i] = self.recursiveBayesian(i, self.proposal[:, i])

    def writeBayeStatsToFile(self, reverse):
        np.savetxt(self.DPMDataDir + '/particle.txt', self.getSmcSamples()[0])
        np.savetxt(self.DPMDataDir + '/IP.txt', self.ips[:, ::(-1) ** reverse].T)
        np.savetxt(self.DPMDataDir + '/weight.txt', self.getPosterior()[:, ::(-1) ** reverse])
