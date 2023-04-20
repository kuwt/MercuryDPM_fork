import glob
import os

import pandas as pd
import numpy as np
from DotDict import DotDict

# loadstatistics (version 1.0) by Thomas Weinhart and conversion to Python by Timo Plath
#
# loads a MercuryDPM stat (or MDCLR data) file into a matlab struct
#
# it accepts a single filename, a cell of filenames, or a string that will
# be inserted into ls
#
# Usages:
# data=loadstatistics('chuteDemo.stat');
# data=loadstatistics('*.stat');
# data=loadstatistics({'chuteDemo.stat','hopperDemo.stat'})

# BEWARE This python code currently works for .stat files and .ene files.
# reading in .data or .restart files is not possible at the moment. If anyone knows a solution he is welcome to
# implement it.s

def load_statistics(filenames,opt=None):
    if opt is None:
        opt = {}
    if isinstance(filenames,list):
        data = [load_file(filename,opt) for filename in filenames]
    else:
        if '?' not in filenames and '*' not in filenames:
            data = load_file(filenames,opt)
        else:
            data = load_statistics(glob.glob(filenames),opt)
    # Transform to a DotDict dictionary to be able to easily access the dictionary by dot.notation
    data = DotDict(data)
    return data

def load_file(filename,opt):
    print('loading data from ' + filename)
    if filename[-5:] == '.data':
        data = load_data_file(filename)
    else:
        data = load_stat_file(filename)
    if not isinstance(data,list):
        data = make_readable(data)
        if 'basic' not in opt:
            data = get_standard_variables(data,opt)
    else:
        for i in range(len(data)):
            data[i] = make_readable(data[i])
            if 'basic' not in opt:
                data[i] = get_standard_variables(data[i],opt)
    # # Transpose arrays in the data dictionary; when converting this file to Python
    # # it swapped indices compared to Matlab
    # for key, value in data.items():
    #     if isinstance(value,np.ndarray):
    #         if len(value.shape) > 1:
    #             if value.shape[1] == value.shape[0]:
    #                 data[key] = value.transpose()
    return data


def load_stat_file(filename):
    data = {}
    data['name'] = filename[:filename.find('.stat')]

    rawdata = np.genfromtxt(filename, skip_header=3)

    if rawdata.ndim == 0:
        return data

    if rawdata.shape[1] == 26:
        # for old .stat files
        print('WARNING: outdated stat file')
        rawdata = np.hstack(
            (rawdata[:, :14], np.zeros((rawdata.shape[0], 3)), rawdata[:, 15:], np.zeros((rawdata.shape[0], 10))))
        data['variable_names'] = ['VolumeFraction', 'Density',
            'MomentumX', 'MomentumY', 'MomentumZ',
            'MomentumFluxXX', 'MomentumFluxXY', 'MomentumFluxXZ', 'MomentumFluxYY', 'MomentumFluxYZ', 'MomentumFluxZZ',
            'EnergyFluxX', 'EnergyFluxY', 'EnergyFluxZ',
            'NormalStressXX', 'NormalStressXY', 'NormalStressXZ', 'NormalStressYY', 'NormalStressYZ', 'NormalStressZZ',
            'TangentialStressXX', 'TangentialStressXY', 'TangentialStressXZ', 'TangentialStressYY',
            'TangentialStressYZ', 'TangentialStressZZ',
            'FabricXX', 'FabricXY', 'FabricXZ', 'FabricYY', 'FabricYZ', 'FabricZZ',
            'CollisionalHeatFluxX', 'CollisionalHeatFluxY', 'CollisionalHeatFluxZ',
            'Dissipation',]
    else:
        variable_names = np.genfromtxt(filename, skip_header=0, max_rows=1, dtype=str)
        data['variable_names'] = variable_names

    # also allows input from the restart and ene files if they exist
    data = read_restart(filename, data)
    data = read_ene(filename, data)

    text = np.genfromtxt(filename, skip_header=1, max_rows=1, dtype=str)
    data['text'] = text

    index_time = np.where(np.isnan(rawdata[:, 2]))[0]
    index_time = np.append(index_time, rawdata.shape[0])
    data['time'] = np.genfromtxt(filename, skip_header=2, max_rows=1)
    data['w'] = np.genfromtxt(filename, skip_header=1, max_rows=1, usecols=1)

    data['coordinates'] = rawdata[:index_time[0], :3]
    doGradient = any(['doGradient' in t for t in text])
    doVariance = any(['doVariance' in t for t in text])

    if len(index_time) == 1:
        data['variables'] = rawdata[:index_time[0], 3:]
    elif doVariance and len(index_time) == 2:
        data['variables'] = rawdata[:index_time[0], 3:]
        data['variance'] = rawdata[index_time[0] + 1:index_time[1], 3:]
    elif doGradient and len(index_time) == 4:
        data['variables'] = rawdata[:index_time[0], 3:]
        data['gradx'] = rawdata[index_time[0] + 1:index_time[1], 3:]
        data['grady'] = rawdata[index_time[1] + 1:index_time[2], 3:]
        data['gradz'] = rawdata[index_time[2] + 1:index_time[3], 3:]
    elif not doGradient:
        dataTemplate = data
        data = [dataTemplate]
        data[0]['variables'] = rawdata[:index_time[0], 3:]
        for i in range(1, len(index_time) - 1):
            data.append(dataTemplate)
            data[i]['time'] = rawdata[index_time[i - 1], :2]
            data[i]['variables'] = rawdata[index_time[i - 1] + 1:index_time[i], 3:]
        print('multiple time steps ({}); creating cell output'.format(len(index_time) - 1))
    else:
        dataTemplate = data
        data = [dataTemplate]
        data[0]['variables'] = rawdata[:index_time[0], 3:]
        data[0]['gradx'] = rawdata[index_time[0] + 1:index_time[1], 3:]
        data[0]['grady'] = rawdata[index_time[1] + 1:index_time[2], 3:]
        data[0]['gradz'] = rawdata[index_time[2] + 1:index_time[3], 3:]
        for i in range(8, len(index_time), 4):
            data.append(dataTemplate)
            data[i / 4]['time'] = rawdata[index_time[i - 4], :2]
            data[i / 4]['variables'] = rawdata[index_time[i - 4] + 1:index_time[i - 3], 3:]
            data[i / 4]['gradx'] = rawdata[index_time[i - 3] + 1:index_time[i - 2], 3:]
            data[i / 4]['grady'] = rawdata[index_time[i - 2] + 1:index_time[i - 1], 3:]
            data[i / 4]['gradz'] = rawdata[index_time[i - 1] + 1:index_time[i], 3:]
        print('multiple time steps ({}); creating cell output'.format(len(index_time) / 4))

    return data

def load_data_file(filename):
    data = {}
    data['name'] = filename[:-5]
    rawdata = pd.read_csv(filename, delim_whitespace=True, skiprows=17)
    data['variable_names'] = ['VolumeFraction', 'Density', 'MomentumX', 'MomentumY', 'MomentumZ', 'MomentumFluxXX', 'MomentumFluxXY', 'MomentumFluxXZ', 'MomentumFluxYY', 'MomentumFluxYZ', 'MomentumFluxZZ', 'EnergyFluxX', 'EnergyFluxY', 'EnergyFluxZ', 'NormalStressXX', 'NormalStressXY', 'NormalStressXZ', 'NormalStressYY', 'NormalStressYZ', 'NormalStressZZ', 'TangentialStressXX', 'TangentialStressXY', 'TangentialStressXZ', 'TangentialStressYY', 'TangentialStressYZ', 'TangentialStressZZ', 'FabricXX', 'FabricXY', 'FabricXZ', 'FabricYY', 'FabricYZ', 'FabricZZ', 'CollisionalHeatFluxX', 'CollisionalHeatFluxY', 'CollisionalHeatFluxZ', 'Dissipation']
    # data['text'] = rawdata.iloc[0,4] + str(rawdata.iloc[0,8]*100)
    rho = rawdata.iloc[1,4]
    data['time'] = 0
    data['coordinates'] = rawdata.iloc[:,3:6]
    VolumeFraction = np.maximum(rawdata.iloc[:,10], 1e-200*np.ones(rawdata.shape[0]))
    data['variables'] = pd.DataFrame(np.concatenate((VolumeFraction,
                                                     VolumeFraction*rho,
                                                     rawdata.iloc[:,51:53].multiply(VolumeFraction.values[:,np.newaxis]),
                                                     rawdata.iloc[:,[41,42,43,45,46,49]],
                                                     np.nan*np.ones((rawdata.shape[0],3)),
                                                     rawdata.iloc[:,[21,22,23,25,26,29]],
                                                     rawdata.iloc[:,[31,32,33,35,36,39]],
                                                     rawdata.iloc[:,[71,72,73,75,76,79]],
                                                     np.nan*np.ones((rawdata.shape[0],4))), axis=1))
    return data

# def load_data_file(filename):
#     data = {}
#     data['name'] = filename[:-5]
#
#     # load raw data from file, with the first three lines as header files (such
#     # that matlab recognizes the array size)
#     rawdata = np.loadtxt(filename,delimiter=' ', skiprows=17)
#
#     # reads variable names from line 1 and text output
#     # from line 2
#     data['variable_names'] = ['VolumeFraction', 'Density', 'MomentumX', 'MomentumY', 'MomentumZ', 'MomentumFluxXX',
#                               'MomentumFluxXY', 'MomentumFluxXZ', 'MomentumFluxYY', 'MomentumFluxYZ', 'MomentumFluxZZ',
#                               'EnergyFluxX', 'EnergyFluxY', 'EnergyFluxZ', 'NormalStressXX', 'NormalStressXY',
#                               'NormalStressXZ', 'NormalStressYY', 'NormalStressYZ', 'NormalStressZZ',
#                               'TangentialStressXX', 'TangentialStressXY', 'TangentialStressXZ', 'TangentialStressYY',
#                               'TangentialStressYZ', 'TangentialStressZZ', 'FabricXX', 'FabricXY', 'FabricXZ',
#                               'FabricYY', 'FabricYZ', 'FabricZZ', 'CollisionalHeatFluxX', 'CollisionalHeatFluxY',
#                               'CollisionalHeatFluxZ', 'Dissipation']
#     with open(filename, 'r') as f:
#         text = f.readlines()
#     data['text'] = text[7][4:7] + str(float(text[7][8:12]) * 100)
#     rho_text = text[10].split()
#     rho = float(rho_text[-1])
#
#     # writes the raw data into bits of data for each timestep
#     # index_time = [find(isnan(rawdata.data(:,2))); size(rawdata.data,1)+1];
#     data['time'] = 0
#     data['coordinates'] = rawdata[:, 3:6]
#     VolumeFraction = np.maximum(rawdata[:, 9], 1e-200 * np.ones((rawdata.shape[0], 1)))
#     data['variables'] = [VolumeFraction, VolumeFraction * rho, rawdata.data[:, 51:53] * VolumeFraction[:, 1, 1,
#                                                                                          1] * rho,
#                       rawdata.data[:, 41:43, 45, 46, 49], np.nan(np.size(rawdata.data, 1), 3),
#                       rawdata.data[:, 21:23, 25, 26, 29], rawdata.data[:, 31:33, 35, 36, 39],
#                       rawdata.data[:, 71:73, 75, 76, 79], np.nan(np.size(rawdata.data, 1), 4)]
#     return data

def make_readable(data):
    data['nx'] = int(len(data['coordinates'][:, 0]) / int(sum(data['coordinates'][:, 0] == data['coordinates'][0, 0])))
    data['ny'] = int(len(data['coordinates'][:, 1]) / int(sum(data['coordinates'][:, 1] == data['coordinates'][0, 1])))
    data['nz'] = int(len(data['coordinates'][:, 2]) / int(sum(data['coordinates'][:, 2] == data['coordinates'][0, 2])))
    shape = np.shape(np.squeeze(np.zeros((data['nz'], data['ny'], data['nx']))))
    if (np.prod(shape) != np.shape(data['coordinates'])[0]):
        print('Warning: cannot put xyz values on mesh')
        shape = (np.shape(data['coordinates'])[0], 1)
    data['x'] = np.reshape(data['coordinates'][:, 0], shape)
    data['y'] = np.reshape(data['coordinates'][:, 1], shape)
    data['z'] = np.reshape(data['coordinates'][:, 2], shape)
    del data['coordinates']

    # the array 'variable' is expanded into smaller arrays, each containing one
    # variable only (f.e. Density), with the field names defined in
    # variable_names.
    # If available, variance, gradx, grady, gradz are also expanded the same way

    if 'variance' in data:
        for i in range(len(data.variable_names)):
            data[data['variable_names'][i].split()[0] + '_var'] = np.reshape(data['variance'][:, i], shape)
        del data['variance']

    if 'gradx' in data:
        for i in range(len(data['variable_names'])):
            data[data['variable_names'][i].split()[0] + '_dx'] = np.reshape(data['gradx'][:, i], shape)
            data[data['variable_names'][i].split()[0] + '_dy'] = np.reshape(data['grady'][:, i], shape)
            data[data['variable_names'][i].split()[0] + '_dz'] = np.reshape(data['gradz'][:, i], shape)
        del data['gradx']
        del data['grady']
        del data['gradz']

    for i in range(len(data['variable_names'])):
        data[data['variable_names'][i].split()[0]] = np.reshape(data['variables'][:, i], shape)
    del data['variables']
    del data['variable_names']

    if 'variance' in data:
        data['VelocityX_var'] = (data['MomentumX'] / data['Density)']) ** 2 * (
                    data['MomentumX_var'] / data['MomentumX'] ** 2 + data['Density_var'] / data['Density'] ** 2)
        data['VelocityY_var'] = (data['MomentumY'] / data['Density']) ** 2 * (
                    data['MomentumY_var'] / data['MomentumY'] ** 2 + data['Density_var'] / data['Density'] ** 2)
        data['VelocityZ_var'] = (data['MomentumZ'] / data['Density']) ** 2 * (
                    data['MomentumZ_var'] / data['MomentumZ'] ** 2 + data['Density_var'] / data['Density'] ** 2)

    if 'Nu' in data:
        data['VolumeFraction'] = data['Nu']
        del data['Nu']
    return data


def get_standard_variables(data, opt):
    data['MomentumX'][np.isnan(data['MomentumX'])] = 0
    data['MomentumY'][np.isnan(data['MomentumY'])] = 0
    data['MomentumZ'][np.isnan(data['MomentumZ'])] = 0

    data['VelocityX'] = data['MomentumX'] / data['Density']
    data['VelocityY'] = data['MomentumY'] / data['Density']
    data['VelocityZ'] = data['MomentumZ'] / data['Density']

    data['VelocityX'][np.isnan(data['VelocityX'])] = 0
    data['VelocityY'][np.isnan(data['VelocityY'])] = 0
    data['VelocityZ'][np.isnan(data['VelocityZ'])] = 0

    data['Temperature'] = (data['MomentumFluxXX'] + data['MomentumFluxYY'] + data['MomentumFluxZZ']) / data[
        'Density'] - (data['VelocityX'] ** 2 + data['VelocityY'] ** 2 + data['VelocityZ'] ** 2) / 3
    data['Temperature'][np.isnan(data['Temperature'])] = 0

    data['Coordinationnumber'] = (data['FabricXX'] + data['FabricYY'] + data['FabricZZ']) / data['VolumeFraction']
    data['Coordinationnumber'][np.isnan(data['Coordinationnumber'])] = 0

    if 'basic' in opt:
        return data

    data['TractionX'] = data['NormalTractionX'] + data['TangentialTractionX']
    data['TractionY'] = data['NormalTractionY'] + data['TangentialTractionY']
    data['TractionZ'] = data['NormalTractionZ'] + data['TangentialTractionZ']
    data['KineticStressXX'] = data['MomentumFluxXX'] - data['Density'] * data['VelocityX'] * data['VelocityX']
    data['KineticStressXY'] = data['MomentumFluxXY'] - data['Density'] * data['VelocityX'] * data['VelocityY']
    data['KineticStressXZ'] = data['MomentumFluxXZ'] - data['Density'] * data['VelocityX'] * data['VelocityZ']
    data['KineticStressYY'] = data['MomentumFluxYY'] - data['Density'] * data['VelocityY'] * data['VelocityY']
    data['KineticStressYZ'] = data['MomentumFluxYZ'] - data['Density'] * data['VelocityY'] * data['VelocityZ']
    data['KineticStressZZ'] = data['MomentumFluxZZ'] - data['Density'] * data['VelocityZ'] * data['VelocityZ']

    data['ContactStressXX'] = data['NormalStressXX'] + data['TangentialStressXX']
    data['ContactStressXY'] = data['NormalStressXY'] + data['TangentialStressXY']
    data['ContactStressXZ'] = data['NormalStressXZ'] + data['TangentialStressXZ']
    data['ContactStressYX'] = data['NormalStressYX'] + data['TangentialStressYX']
    data['ContactStressYY'] = data['NormalStressYY'] + data['TangentialStressYY']
    data['ContactStressYZ'] = data['NormalStressYZ'] + data['TangentialStressYZ']
    data['ContactStressZX'] = data['NormalStressZX'] + data['TangentialStressZX']
    data['ContactStressZY'] = data['NormalStressZY'] + data['TangentialStressZY']
    data['ContactStressZZ'] = data['NormalStressZZ'] + data['TangentialStressZZ']

    data['StressXX'] = data['NormalStressXX'] + data['TangentialStressXX'] + data['KineticStressXX']
    data['StressXY'] = data['NormalStressXY'] + data['TangentialStressXY'] + data['KineticStressXY']
    data['StressXZ'] = data['NormalStressXZ'] + data['TangentialStressXZ'] + data['KineticStressXZ']
    data['StressYX'] = data['NormalStressYX'] + data['TangentialStressYX'] + data['KineticStressXY']
    data['StressYY'] = data['NormalStressYY'] + data['TangentialStressYY'] + data['KineticStressYY']
    data['StressYZ'] = data['NormalStressYZ'] + data['TangentialStressYZ'] + data['KineticStressYZ']
    data['StressZX'] = data['NormalStressZX'] + data['TangentialStressZX'] + data['KineticStressXZ']
    data['StressZY'] = data['NormalStressZY'] + data['TangentialStressZY'] + data['KineticStressYZ']
    data['StressZZ'] = data['NormalStressZZ'] + data['TangentialStressZZ'] + data['KineticStressZZ']
    data['Pressure'] = (data['StressXX'] + data['StressYY'] + data['StressZZ']) / 3

    data['KineticPressure'] = (data['KineticStressXX'] + data['StressYY'] + data['StressZZ']) / 3

    # macroscopic friction coefficient, or deviator-stress ratio sD1:=(S1-S3)/(2*p)  or sD2 (to be discussed)
    data['MacroFrictionCoefficient'] = np.zeros(data['x'].shape)
    data['Sd'] = np.zeros(data['x'].shape)
    data['SdStar'] = np.zeros(data['x'].shape)
    for i in range(len(data['x'])):
        for j in range(len(data['x'])):
            Stress = np.array([[data['StressXX'][j,i], data['StressXY'][j,i], data['StressXZ'][j,i]],
                               [data['StressXY'][j,i], data['StressYY'][j,i], data['StressYZ'][j,i]],
                               [data['StressXZ'][j,i], data['StressYZ'][j,i], data['StressZZ'][j,i]]])
            if np.array_equal(Stress, np.zeros((3, 3))) or np.sum(np.isnan(Stress)) > 0:
                data['MacroFrictionCoefficient'][j] = np.nan
                data['Sd'][j] = np.nan
                data['SdStar'][j] = np.nan
            else:
                d1, v1 = np.linalg.eig(Stress)
                v, d = dsort_ev(v1, np.diag(d1))
                d = np.diag(d)
                p = np.sum(d) / 3
                data['MacroFrictionCoefficient'][j] = -data['StressXZ'][j] / data['StressZZ'][j]
                data['Sd'][j] = (d[0] - d[2]) / (2 * p)
                data['SdStar'][j] = np.sqrt(((d[0] - d[2]) ** 2 + (d[1] - d[2]) ** 2 + (d[0] - d[1]) ** 2) / (6 * p ** 2))

    if 'VolumeFraction_dz' in data:
        for i in 'XYZ':
            for k in 'xyz':
                vector = i + '_d' + k
                data['Velocity' + vector] = (data['Momentum' + vector]*data['Density']-data['Momentum' + i]*data[
                    'Density_d' + k])/data['Density']**2
                data['Velocity' + vector][np.isnan(data['Velocity' + vector])] = 0
        for i in 'XYZ':
            for j in 'XYZ':
                for k in 'xyz':
                    tensor = i + j + '_d' + k
                    if i<j:
                        symtensor = tensor
                    else:
                        symtensor = tensor[2] + tensor[1] + tensor[3:]
                    data['KineticStress' + tensor] = data['MomentumFlux' + symtensor] - data['Momentum' + i + '_d' + k]*data['Velocity' + j] - data['Momentum' + j + '_d' + k]*data['Velocity' + i]
                    data['Stress' + tensor] = data['NormalStress' + tensor] + data['TangentialStress' + tensor] + data['KineticStress' + tensor]
    if data['nz']>1 and data['nx']==1 and data['ny']==1:
        data = get_depth_variables(data)
        if 'ParticleDensity' in data:
            data['MomentumEquationsRemainder'], data['MomentumEquationsMaximum'] = get_momentum_equation(data)
    return data

def get_depth_variables(data):
    dz = np.diff(data['z'][0:2, 0])
    nz = data['z'].shape[0]
    data['ShearRate'] = np.diff(data['VelocityX']) / dz
    data['ShearRate'] = np.append(data['ShearRate'], data['ShearRate'][-1])
    if 'd' in data:
        data['InertialParamter'] = data['ShearRate'] * data['d'] / np.sqrt(data['StressZZ'] / data['ParticleDensity'])
    data['FlowVolumeFraction'] = np.mean(data['VolumeFraction'])
    data['FlowDensity'] = np.mean(data['Density'], 1)
    data['FlowMomentumX'] = np.mean(data['MomentumX'], 1)
    data['FlowMomentumY'] = np.mean(data['MomentumY'], 1)
    data['FlowMomentumZ'] = np.mean(data['MomentumZ'], 1)
    data['FlowVelocityX'] = data['FlowMomentumX'] / data['FlowDensity']
    data['FlowVelocityY'] = data['FlowMomentumY'] / data['FlowDensity']
    data['FlowVelocityZ'] = data['FlowMomentumZ'] / data['FlowDensity']
    data['BottomFrictionCoefficient'] = -data['StressXZ'][0, :] / data['StressZZ'][0, :]
    if not (data['nx'] == 1 and data['ny'] == 1):
        ContactStressZZ = data['NormalStressZZ'] + data['TangentialStressZZ']
        data['FlowHeight'] = np.max(data['z'] * (ContactStressZZ != 0), 1)
    else:
        kappa = 0.02
        IntDensity = np.cumsum(data['Density'])
        Base = np.min(data['z'][IntDensity >= kappa * np.max(IntDensity)])
        Surface = np.max(data['z'][IntDensity <= (1 - kappa) * np.max(IntDensity)])
        data['FlowHeight'] = (Surface - Base) / (1 - 2 * kappa)
        data['Base'] = Base - kappa * data['FlowHeight']
        data['Surface'] = Surface + kappa * data['FlowHeight']
        FlowHeightIndex = data['z'] > data['Base'] and data['z'] < data['Base'] + data['FlowHeight']
        if not FlowHeightIndex.any() and 'Gravity' in data:
            data['Froude'] = np.linalg.norm([data['FlowVelocityX'], data['FlowVelocityY'], data['FlowVelocityZ']]) / np.sqrt(
                data['FlowHeight'] * (-data['Gravity'][2]))
    return data

def dsort_ev(V, D):
    Ds=D
    Vs=V
    tmpvec = np.empty(3)
    if (Ds[0,0])<(Ds[1,1]):
        tmp=Ds[0,0]
        Ds[0,0]=Ds[1,1]
        Ds[1,1]=tmp
        tmpvec[:]=Vs[:,0]
        Vs[:,0]=Vs[:,1]
        Vs[:,1]=tmpvec[:]
    if (Ds[1,1])<(Ds[2,2]):
        tmp=Ds[1,1]
        Ds[1,1]=Ds[2,2]
        Ds[2,2]=tmp
        tmpvec[:] = Vs[:, 1]
        Vs[:, 1] = Vs[:, 2]
        Vs[:, 2] = tmpvec[:]
    if (Ds[0,0])<(Ds[1,1]):
        tmp=Ds[0,0]
        Ds[0,0]=Ds[1,1]
        Ds[1,1]=tmp
        tmpvec[:]=Vs[:,0]
        Vs[:,0]=Vs[:,1]
        Vs[:,1]=tmpvec[:]
    if (Ds[1,1])<(Ds[2,2]):
        tmp=Ds[1,1]
        Ds[1,1]=Ds[2,2]
        Ds[2,2]=tmp
        tmpvec[:]=Vs[:,1]
        Vs[:,1]=Vs[:,2]
        Vs[:,2]=tmpvec[:]
    return Vs, Ds


def get_momentum_equation(data):
    if 'VolumeFraction_dx' in data:
        NablaMomentumX = data['ParticleDensity'][0] * np.array(
            [data['MomentumX_dx'], data['MomentumX_dy'], data['MomentumX_dz']])
        NablaMomentumY = data['ParticleDensity'][0] * np.array(
            [data['MomentumY_dx'], data['MomentumY_dy'], data['MomentumY_dz']])
        NablaMomentumZ = data['ParticleDensity'][0] * np.array(
            [data['MomentumZ_dx'], data['MomentumZ_dy'], data['MomentumZ_dz']])
        NablaDotStress = np.array([data['StressXX_dx'] + data['StressXY_dy'] + data['StressXZ_dz'],
                                   data['StressYX_dx'] + data['StressYY_dy'] + data['StressYZ_dz'],
                                   data['StressZX_dx'] + data['StressZY_dy'] + data['StressZZ_dz']])
    else:
        print('estimating gradient!')
        NablaMomentumX = nabla(data['ParticleDensity'][0] * data['VolumeFraction'] * data['VelocityX'], data['x'],
                               data['y'], data['z'])
        NablaMomentumY = nabla(data['ParticleDensity'][0] * data['VolumeFraction'] * data['VelocityY'], data['x'],
                               data['y'], data['z'])
        NablaMomentumZ = nabla(data['ParticleDensity'][0] * data['VolumeFraction'] * data['VelocityZ'], data['x'],
                               data['y'], data['z'])

        NablaDotStress = nabla(np.array(
            [data['StressXX'], data['StressXY'], data['StressXZ'], data['StressYY'], data['StressYZ'],
             data['StressZZ']]), data['x'], data['y'], data['z'])

    VelocityDotNablaMomentum = np.array([np.sum(np.array([data['VelocityX'], data['VelocityX'],
                                                          data['VelocityZ']]) * NablaMomentumX, axis=0),
                                         np.sum(np.array([data['VelocityX'], data['VelocityX'],
                                                          data['VelocityZ']]) * NablaMomentumY, axis=0),
                                         np.sum(np.array([data['VelocityX'], data['VelocityX'],
                                                          data['VelocityZ']]) * NablaMomentumZ, axis=0)])

    DensityGravity = data['ParticleDensity'][0] * data['VolumeFraction'] * data['Gravity']

    Traction = np.array([data['TractionX'], data['TractionY'], data['TractionZ']])

    remainder = VelocityDotNablaMomentum - DensityGravity + NablaDotStress + Traction
    maximum = np.max(np.array([np.max(VelocityDotNablaMomentum), np.max(DensityGravity), np.max(NablaDotStress)]))
    return remainder, maximum

def read_ene(statname,data):
    #load ene data
    filename = statname[:-5] + '.ene'
    if not os.path.isfile(filename):
        #if unsuccessful, load ene data with last part of filename removed
        dots = statname.count('.')
        if dots > 1:
            splitfilename = os.path.splitext(filename)
            for dot in range(1,dots):
                splitfilename = os.path.splitext(splitfilename[0])
            filename = splitfilename[0] + '.ene'
        if not os.path.isfile(filename):
            print(filename + ' not found')
        else:
            #if unsuccessful, give up
            print(statname[:-5] + '.ene not found, using ' + filename + ' instead')
    if os.path.isfile(filename):
        # load raw data from file, with the first three lines as header files
        # (such that matlab recognizes the array size)
        rawdata = np.genfromtxt(filename,skip_header=1)
        # if tabs are used as delimiters
        if rawdata.shape[1] != 8:
            rawdata = np.genfromtxt(filename,delimiter='\t',skip_header=1)
        if rawdata.shape[1] != 8:
            print(filename + ' not readable')
            return
        data['Ene'] = dict()
        data['Ene']['Time'] = rawdata[:,0]
        data['Ene']['Gra'] = rawdata[:,1]
        data['Ene']['Kin'] = rawdata[:,2]
        data['Ene']['Rot'] = rawdata[:,3]
        data['Ene']['Ela'] = rawdata[:,4]
        data['Ene']['ComX'] = rawdata[:,5]
        data['Ene']['ComY'] = rawdata[:,6]
        data['Ene']['ComZ'] = rawdata[:,7]
    return data

def read_restart(statname,data):
    #load restart data
    filename = statname.replace('.stat','.restart')
    if not os.path.isfile(filename):
    # if name.restart could not be opened, it maybe because name has an appendix
    # (e.g. problem.X.stat tries to load problem.X.restart).
    # Thus, try loading restart data with last part of name removed
    # (e.g. load problem.restart)
    # if unsuccessful, load ene data with last part of filename removed
        dots = statname.count('.')
        if dots > 1:
            splitfilename = os.path.splitext(filename)
            for dot in range(1, dots):
                splitfilename = os.path.splitext(splitfilename[0])
            filename = splitfilename[0] + '.restart'
            if not os.path.isfile(filename):
                print(filename + ' not found')
                return data
    if os.path.isfile(filename):
        # read the file into a numpy array
        with open(filename, 'r') as f:
            lines = f.readlines()
        tdata = []
        for line in lines:
            splitLines = line.split()
            tdata.extend(splitLines)
        # parse the file into a dictionary
        #this assumes monodispersed particles
        if (len(tdata[0])==1):
            #old restart files
            data.ParticleDensity = float(tdata[28])
            data.d = 2*float(tdata[27])
            data.Gravity = np.array([float(tdata[1]),float(tdata[2]),float(tdata[3])])
            data.N = float(tdata[-4])
            data.Domain=np.array([float(tdata[4]),float(tdata[5]),float(tdata[6]),float(tdata[7]),float(tdata[8]),float(tdata[9])])
        else:
            #split tdata into multiple strings
            #new restart version
            i = [j for j, x in enumerate(tdata) if x == "density"]
            if i==[]: i = [j for j, x in enumerate(tdata) if x == "rho"]
            if len(i)>1: print('multiple species detected')
            if i!=[]: data["ParticleDensity"] = np.array([tdata[j+1] for j in i]).astype(float)
            i = [j for j, x in enumerate(tdata) if x == 'slidingFrictionCoefficient']
            if i==[]: i = [j for j, x in enumerate(tdata) if x == 'mu']
            if i!=[]: data["Mu"] = np.array([tdata[j+1] for j in i]).astype(float)
            i = [j for j, x in enumerate(tdata) if x == 'gravity']
            if i!=[]: data["Gravity"] = np.array([[tdata[j+1],tdata[j+2],tdata[j+3]] for j in i]).astype(float)
            i = [j for j, x in enumerate(tdata) if x == 'Particles']
            if i!=[]: data["N"] = float(np.array([tdata[j+1] for j in i]))
            i = [j for j, x in enumerate(tdata) if x == 'xMin']
            if i==[]: i = [j for j, x in enumerate(tdata) if x == 'xmin']
            if i!=[]: data["Domain"] = np.array([[tdata[j+1],tdata[j+3],tdata[j+5],tdata[j+7],tdata[j+9],tdata[j+11]]
                                                 for j in i]).astype(float)
            i = [j for j, x in enumerate(tdata) if x == 'FixedParticleRadius']
            if i!=[]: data["FixedParticleRadius"] = np.array([tdata[j+1] for j in i]).astype(float)
            i = [j for j, x in enumerate(tdata) if x == 'MinInflowParticleRadius']
            if i!=[]:
                data["InflowParticleRadius"] = np.array([[tdata[j+1],tdata[j+3]] for j in i]).astype(float)
                data["d"] = sum(data["InflowParticleRadius"])
            else:
                i = [j for j, x in enumerate(tdata) if x == 'Particles']
                if i!=[]:
                    if ([tdata[j+1] for j in i] == 0):
                        data["d"] = np.nan
                    else:
                        k = [min([j for j, x in enumerate([tdata[v+2:] for v in i][0]) if x == 'radius'])]
                        data["d"] = 2 * float(np.array([tdata[j+v+3] for j,v in zip(i,k)]))
        if 'd' in data: data["ParticleVolume"] = np.pi/6*data["d"]**3
        if 'Domain' in data: data["DomainVolume"] = np.prod(data["Domain"][0,1::2].astype(float)-data["Domain"][0,::2].astype(float))
        if 'Gravity' in data: data["ChuteAngle"] = round(np.arctan(-data["Gravity"][0,0].astype(float)/data[
            "Gravity"][0,2].astype(float))*400)/400
        return data

def nabla(variable, x, y, z):
    # stores the grid dimensions
    n = [len(x) / sum(x[0] == x), len(y) / sum(y[0] == y), len(z) / sum(z[0] == z)]

    # stores coordinates in grid
    X = np.zeros(n[::-1])
    X[:] = x
    Y = np.zeros(n[::-1])
    Y[:] = y
    Z = np.zeros(n[::-1])
    Z[:] = z
    V = np.zeros(n[::-1])

    # stores differentials (small number if X/Y/Z is constant)
    dX = np.maximum(1e-60, (X[:, :, 2:] - X[:, :, :-1]))
    dY = np.maximum(1e-60, (Y[:, 2:, :] - Y[:, :-1, :]))
    dZ = np.maximum(1e-60, (Z[2:, :, :] - Z[:-1, :, :]))

    if variable.shape[1] == 1:  # scalar variable
        V[:] = variable
        dXV = (V[:, :, 2:] - V[:, :, :-1]) / dX
        dYV = (V[:, 2:, :] - V[:, :-1, :]) / dY
        dZV = (V[2:, :, :] - V[:-1, :, :]) / dZ

        NablaVariable = np.vstack((dXV.flatten(), dYV.flatten(), dZV.flatten())).T
    elif variable.shape[1] == 3:  # vector
        V[:] = variable[:, 0]
        dXV = (V[:, :, 2:] - V[:, :, :-1]) / dX
        V[:] = variable[:, 1]
        dYV = (V[:, 2:, :] - V[:, :-1, :]) / dY
        V[:] = variable[:, 2]
        dZV = (V[2:, :, :] - V[:-1, :, :]) / dZ

        NablaVariable = dXV.flatten() + dYV.flatten() + dZV.flatten()
    elif variable.shape[1] == 6:  # symmetric matrix
        V[:] = variable[:, 0]
        dXVXX = (V[:, :, 2:] - V[:, :, :-1]) / dX
        V[:] = variable[:, 1]
        dXVYX = (V[:, :, 2:] - V[:, :, :-1]) / dX
        V[:] = variable[:, 2]
        dXVZX = (V[:, :, 2:] - V[:, :, :-1]) / dX

        V[:] = variable[:, 1]
        dYVXY = (V[:, 2:, :] - V[:, :-1, :]) / dY
        V[:] = variable[:, 3]
        dYVYY = (V[:, 2:, :] - V[:, :-1, :]) / dY
        V[:] = variable[:, 4]
        dYVZY = (V[:, 2:, :] - V[:, :-1, :]) / dY

        V[:] = variable[:, 2]
        dZVXZ = (V[2:, :, :] - V[:-1, :, :]) / dZ
        V[:] = variable[:, 4]
        dZVYZ = (V[2:, :, :] - V[:-1, :, :]) / dZ
        V[:] = variable[:, 5]
        dZVZZ = (V[2:, :, :] - V[:-1, :, :]) / dZ

        NablaVariable = np.vstack((dXVXX.flatten() + dYVXY.flatten() + dZVXZ.flatten(),
                                   dXVYX.flatten() + dYVYY.flatten() + dZVYZ.flatten(),
                                   dXVZX.flatten() + dYVZY.flatten() + dZVZZ.flatten())).T
    else:
        print('Error')

    return NablaVariable


def nablaCenter(variable, x, y, z):
    # stores the grid dimensions
    n = [len(x), len(y), len(z)]

    # stores coordinates in grid
    X = np.zeros(n[::-1])
    X.flat = x
    Y = np.zeros(n[::-1])
    Y.flat = y
    Z = np.zeros(n[::-1])
    Z.flat = z
    V = np.zeros(n[::-1])

    # stores differentials (small number if X/Y/Z is constant)
    dX = np.maximum(1e-60, (X[:, :, 1:] - X[:, :, :-1]))
    dY = np.maximum(1e-60, (Y[:, 1:] - Y[:, :-1]))
    dZ = np.maximum(1e-60, (Z[1:] - Z[:-1]))

    if variable.shape[1] == 1:  # scalar variable
        V.flat = variable
        dXV = (V[:, :, 1:] - V[:, :, :-1]) / dX
        dYV = (V[:, 1:] - V[:, :-1]) / dY
        dZV = (V[1:] - V[:-1]) / dZ

        NablaVariable = np.vstack((dXV.flat, dYV.flat, dZV.flat)).T
    elif variable.shape[1] == 3:  # vector
        V.flat = variable[:, 0]
        dXV = (V[:, :, 1:] - V[:, :, :-1]) / dX
        V.flat = variable[:, 1]
        dYV = (V[:, 1:] - V[:, :-1]) / dY
        V.flat = variable[:, 2]
        dZV = (V[1:] - V[:-1]) / dZ

        NablaVariable = dXV.flat + dYV.flat + dZV.flat
    elif variable.shape[1] == 6:  # symmetric matrix
        V.flat = variable[:, 0]
        dXVXX = (V[:, :, 1:] - V[:, :, :-1]) / dX
        V.flat = variable[:, 1]
        dXVYX = (V[:, :, 1:] - V[:, :, :-1]) / dX
        V.flat = variable[:, 2]
        dXVZX = (V[:, :, 1:] - V[:, :, :-1]) / dX

        V.flat = variable[:, 1]
        dYVXY = (V[:, 1:] - V[:, :-1]) / dY
        V.flat = variable[:, 3]
        dYVYY = (V[:, 1:] - V[:, :-1]) / dY
        V.flat = variable[:, 4]
        dYVZY = (V[:, 1:] - V[:, :-1]) / dY

        V.flat = variable[:, 2]
        dZVXZ = (V[1:] - V[:-1]) / dZ
        V.flat = variable[:, 4]
        dZVYZ = (V[1:] - V[:-1]) / dZ
        V.flat = variable[:, 5]
        dZVZZ = (V[1:] - V[:-1]) / dZ

        NablaVariable = np.vstack((dXVXX.flat + dYVXY.flat + dZVXZ.flat, dXVYX.flat + dYVYY.flat + dZVYZ.flat,
                                   dXVZX.flat + dYVZY.flat + dZVZZ.flat)).T
    else:
        print('Error')
    return NablaVariable