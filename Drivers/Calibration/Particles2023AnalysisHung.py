"""!
@brief This is the calibration script from the BSc of Q.H. Nguyen https://purl.utwente.nl/essays/91991
This script contains functions for data analysis, model training, and evaluation.
It includes functions for plotting data, training neural network models, random forest models,
and generating correlation matrices.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from keras import layers, Sequential
from sklearn.model_selection import train_test_split
from tensorflow import keras
import os
import tensorflow_addons as tfa
import pandas as pd
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.kernel_ridge import KernelRidge
import joblib
from sklearn.metrics import mean_absolute_error
from keras_visualizer import visualizer
from sklearn.model_selection import cross_val_score
from sklearn import metrics
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms
import mpl_scatter_density
from matplotlib.colors import LinearSegmentedColormap

params = {'lines.linewidth': 1, 'axes.labelsize': 12, 'font.size': 12, 'legend.fontsize': 9,
          'xtick.labelsize': 9, 'ytick.labelsize': 9, 'text.usetex': True, 'font.family': 'serif',
          'legend.edgecolor': 'k', 'legend.fancybox': False}
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
matplotlib.rcParams.update(params)

def plot_result_GL(gl_file):
    """!
    @brief Plot the GL data.

    This function reads the data from a given file and plots the data.

    @param[in] gl_file Path to the GL data file.
    """
    data = np.loadtxt(gl_file, delimiter=";", encoding='utf-8-sig')
    angle = data[:, 4]
    res_coef = data[:, 0]
    sliding = data[:, 1]
    rolling = data[:, 2]
    bond = data[:, 3]

    fig, axes = plt.subplots(1, 4, dpi=100, figsize=(7, 5), sharey=True, sharex=False)
    axes[0].scatter(res_coef[0], angle[0], c='b')
    color = ['r', 'g', 'y', 'k']
    for i in range(1, 4, 1):
        for j in range(4):
            axes[i].scatter(data[:, i][j], angle[j], c=color[i])
            axes[i].scatter(data[:, i][j], angle[j], c=color[i])
            axes[i].scatter(data[:, i][j], angle[j], c=color[i])
            axes[i].scatter(data[:, i][j], angle[j], c=color[i])

    plt.show()
    # axes[1].scatter(sliding[1], angle[1], c='r')
    # axes[1].scatter(sliding[2], angle[2], c='r')
    # axes[1].scatter(sliding[3], angle[3], c='r')
    # axes[1].scatter(sliding[4], angle[4], c='r')


# def plot_data_distribution(file_name, save):
#     data = np.loadtxt(file_name, skiprows=1, delimiter=",")
#     data_fig_name = file_name.split(".")[0] + "_dataDistribution"
#     angle_fig_name = file_name.split(".")[0] + "_angleDistribtion"
#     restitution_coefficient = data[:, 0]
#     sliding_friction = data[:, 1]
#     rolling_friction = data[:, 2]
#     bond_number = data[:, 3]
#     angle = data[:, 4]
#
#     fig, axes = plt.subplots(1, 5, dpi=100, figsize=(7, 5), sharey=True, sharex=False)
#     axes[0].hist(restitution_coefficient, bins=10, color="red", alpha=0.5, edgecolor='black')
#     axes[1].hist(sliding_friction, bins=10, color="green", alpha=0.5, edgecolor='black')
#     axes[2].hist(rolling_friction, bins=10, color="yellow", alpha=0.5, edgecolor='black')
#     axes[3].hist(bond_number, bins=10, color="blue", alpha=0.5, edgecolor='black')
#     axes[4].hist(angle, bins=10, color="pink", alpha=0.5, edgecolor='black')
#
#     axes[0].set(xlabel="Restitution Coefficient")
#     axes[1].set(xlabel="Sliding Friction")
#     axes[2].set(xlabel="Rolling Friction")
#     axes[3].set(xlabel="Bond Number")
#     axes[4].set(xlabel="Static AoR")
#     # plt.savefig(data_fig_name)

def predict_model(file_name, num_nodes, num_layer):
    """!
    @brief Train and evaluate a neural network model for angle of repose prediction.

    This function trains a neural network model using the given data and evaluates its performance.

    @param[in] file_name Path to the input data file.
    @param[in] num_nodes Number of nodes in the hidden layers of the neural network.
    @param[in] num_layer Number of hidden layers in the neural network.

    @return The evaluation result of the test data set.
    """
    data = np.loadtxt(file_name, skiprows=1, delimiter=",")
    angle = data[:, 4] / 90
    micro_param = data[:, 0:4]
    x_train, x_test, y_train, y_test = train_test_split(micro_param, angle, test_size=0.15)

    model_relu = Sequential() # create a sequential model.
    model_relu.add(layers.Dense(num_nodes, input_dim=4, activation='relu')) # input layer with 4 microparameters.
    for _ in range(num_layer - 2): # add hidden layers with hidden nodes
        model_relu.add(layers.Dense(num_nodes, input_dim=num_nodes, activation='relu'))
    model_relu.add(layers.Dense(1, input_dim=num_nodes, activation='linear'))

    model_relu.compile(optimizer='Adam', loss="mean_absolute_error", metrics=tfa.metrics.RSquare()) # compile model/
    fit_result = model_relu.fit(x_train, y_train, epochs=50, batch_size=32, verbose=0, validation_split=0.05)
    val_result = model_relu.evaluate(x_test, y_test)
    print(fit_result.history['val_loss'][-1], fit_result.history['val_r_square'][-1])
    if val_result[0] < 0.02 and fit_result.history['val_loss'][-1] < 0.02:
        model_name = file_name + "_" + str(num_layer) + "-" + str(num_nodes) + "-" + str(val_result[0]) + \
                     "-" + str(val_result[1])
        if "Sand_bew" in file_name:
            model_relu.save("model/Sand/" + model_name)
        elif "Eskal" in file_name:
            model_relu.save("model/Eskal/" + model_name)
    return val_result #



def random_forest(file_name):
    """!
    @brief Train and evaluate a random forest model for angle of repose prediction.

    This function trains a random forest model using the given data and evaluates its performance.

    @param[in] file_name Path to the input data file.

    @return The mean absolute error of the model.
    """
    data = np.loadtxt(file_name, skiprows=1, delimiter=",")
    angle = data[:, 4] / 90
    micro_param = data[:, 0:4]
    while True:
        x_train, x_test, y_train, y_test = train_test_split(micro_param, angle, test_size=0.2)
        regressor = RandomForestRegressor(random_state=0, criterion="absolute_error") # random forest in scikit learn.
        regressor.fit(x_train, y_train)

        y_pred = regressor.predict(x_test)
        error = mean_absolute_error(y_test, y_pred)
        cross_val = abs(np.mean(cross_val_score(regressor, x_train, y_train, cv=5, scoring='neg_mean_absolute_error')))
        print(cross_val)
        if error < 0.012:
            print("dumping")
            joblib.dump(regressor, "model/random_forest_{}.joblib".format(file_name))
            break

    return error



def plot_4d(file_name):
    """!
    @brief Plot a 3D scatter plot of the training data.

    This function creates a 3D scatter plot for the first three dimensions of the data,
    using angle of repose for color coding.

    @param[in] file_name Path to the input data file.
    """
    data = np.loadtxt(file_name, skiprows=1, delimiter=",")
    # np.meshgrid
    X, Y, Z = (data[:, 0], data[:, 1], data[:, 2])
    T = data[:, 4]
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    # Creating plot
    ax.scatter3D(X, Y, Z, c=T, alpha=0.7, marker='.')
    plt.show()


def plot_relationship(file_name, save):
    """!
    @brief Plot the relationship between input parameters and the angle of repose.

    This function creates scatter plots to visualize the relationship between each input parameter
    and the angle of repose.

    @param[in] file_name Path to the input data file.
    @param[in] save Path to save the generated plot.
    """
    data = np.loadtxt(file_name, skiprows=1, delimiter=",")
    # data_fig_name = file_name.split(".")[0] + "_dataDistribution"
    # angle_fig_name = file_name.split(".")[0] + "_angleDistribtion"
    fig, axes = plt.subplots(1, 4, dpi=100, figsize=(12, 3))
    name = ['Restitution coefficient', 'Sliding friction', 'Rolling friction', 'Bond number', 'Angle of repose']
    for i in range(4):
        if i == 0:
            axes[i].set_ylabel(name[4])
        axes[i].scatter(data[:, i], data[:, 4])
        axes[i].axhline(y=33, xmin=0, xmax=1, c='red' , )
                        # label='Experimental value = $33^{\circ}$')
        axes[i].set_xlabel(name[i])
        # axes[i].legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
        #                mode="expand", borderaxespad=0, ncol=3)
        axes[i].set_ylim((20, 35))

    plt.tight_layout()
    plt.savefig(save, bbox_inches="tight")


def get_correlation(file_name):
    """
    @brief Calculate the correlation between input parameters and the angle of repose.

    This function calculates the correlation matrix between input parameters and the angle of repose.

    @param[in] file_name Path to the input data file.

    @return The correlation matrix.
    """
    pd_data = pd.read_csv(file_name)
    pd_data.iloc[:, 4] = pd_data.iloc[:, 4] / 90
    return pd_data.corr().iloc[[4], 0:4]


def write_calibration_script(var: str, material: str):
    """!
    @brief Generate a calibration script for simulation.

    This function generates a calibration script for simulation based on the provided parameters.

    @param[in] var String containing the parameters for the simulation.
    @param[in] material The material type for the simulation.

    @return The generated calibration script.
    """
    param_list = var.split(" ")
    param_str = ""
    for param in param_list:
        param_str += param + "_"

    if material == "Sand":
        density = 2653
        psd_str = "cumulative volume diameter 3e-4 0.0621 4.25e-4 0.24 5e-4 0.55 6e-4 1"
    elif material == "Eskal":
        density = 2737
        psd_str = "cumulative volume diameter 9.7e-5 0.1 1.38e-4 0.5 1.94e-4 0.9"
    else:
        return

    return "srun --ntasks=1 " \
           "/home/s2096307/Calibration/cmake-build-release-remote-host/Drivers/Calibration/CalibrationHeap " \
           "-speciesType {} -species LinearViscoelasticFrictionReversibleAdhesiveSpecies -density {} " \
           "-psd {} -collisionTime 0.000067997 -torsionFriction 0 " \
           "-normalStress 1000 800 600 400 200 -restitutionCoefficient {} -slidingFriction {} " \
           "-rollingFriction {} -bondNumber {} " \
           "-param {} &".format(material, density, psd_str,
                                param_list[0], param_list[1], param_list[2], param_list[3], param_str)


# def remake_file(file_name):
#     fin = open(file_name, "rt")
#     # read file contents to string
#     data = fin.read()
#     # replace all occurrences of the required string
#     data = data.replace('[', '')
#     data = data.replace(']', '')
#     # close the input file
#     fin.close()
#     # open the input file in write mode
#     fin = open(file_name, "wt")
#     # overrite the input file with the resulting data
#     fin.write(data)
#     # close the file
#     fin.close()


# def plot_density(file_name, save):
#     data = np.loadtxt(file_name, skiprows=1, delimiter=",")
#     fig, axes = plt.subplots(1, 4, dpi=100, figsize=(12, 3))
#     for i in range(4):
#         h = axes[i].hist2d(data[:, i], data[:, 4])
#         axes[i].axhline(y=33, xmin=0, xmax=1, c='red')
#         # , label='Experimental value = $33^{\circ}$')
#
#
#     fig.colorbar(h[3])
#     plt.tight_layout()



file_name_eskal = "Eskal150NN.csv"
file_name_sand = "SandNN.csv"

plot_relationship("valid_combinations/nn_sand.csv", 'figures/nn_sand.png')
# remake_file('rf_sand_predict_old.csv')
# plot_density("rf_eskal_predict.csv", 'rfsand_predict.png')

plt.show()

# fig, axes = plt.subplots(4, 1, dpi=100, figsize=(7, 5), sharey=False, sharex=True)
# for i in range(4):
#     file = "Eskal150GL" + str(i) + ".csv"
#     corr = get_correlation(file)
#     sns.heatmap(corr,  annot=True, ax=axes[i],cbar=False)
#
#     # sns.heatmap(corr_data_sand,  annot=True, ax=axes[1], cbar=False)
#     # sns.heatmap(corr_data_eskal_GL,  annot=True, ax=axes[2], cbar=False)
#
#
# # plot_data_distribution(file_name)
# # plot_relationship(file_name)
# plt.show()


file_name = "SandNN.csv"

# JOBLIB IS USED TO SAVE THE RF MODEL, WHILE NN MODELS ARE SAVED BY THE KERAS LIBRARY.

# rfm = joblib.load("model/random_forest_eskal.joblib")

# print(rfm.predict([[0.9105937, 0.0691358, 0.9256, 0.202332]]) * 90)

# if "Eskal" in file_name:
#     model_data = open("EskalNNerror.csv", 'w+')
#     model_list = os.listdir('model/Eskal')
# elif "Sand" in file_name:
#     model_data = open("SandNNerror.csv", 'w+')
#     model_list = os.listdir('model/Sand')

# visualizer(model_relu, format='png', view=True)


for i in range(2, 15):
    for j in range(5, 15):
        res = predict_model(file_name, j, i)
        out = str(i) + ',' + str(j) + ',' + str(res) + "\n"

for _ in range(20):
    res = predict_model(file_name, num_nodes=13, num_layer=14)

