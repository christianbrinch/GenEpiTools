# -*- coding: utf-8 -*-
""" Tools for doing machine learning

"""

__author__ = "Christian Brinch"
__copyright__ = "Copyright 2019"
__credits__ = ["Christian Brinch"]
__license__ = "AFL 3.0"
__version__ = "0.1"
__maintainer__ = "Christian Brinch"
__email__ = "cbri@food.dtu.dk"

import somoclu
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial as ss
from genepitools.tf_som import SelfOrganizingMap


def get_tfsom(data, mapsize, epochs=10):
    # if type(data) == 'pandas.core.frame.DataFrame':
    data = np.array([list(data[i]) for i in data.columns]).astype(np.float32)

    graph = tf.Graph()
    with graph.as_default():
        session = tf.Session(config=tf.ConfigProto(
                             allow_soft_placement=True,
                             log_device_placement=False))
        dataset = tf.data.Dataset.from_tensor_slices(data.astype(np.float32))
        dataset = dataset.repeat()
        dataset = dataset.batch(128)
        iterator = dataset.make_one_shot_iterator()
        next_element = iterator.get_next()

        som = SelfOrganizingMap(m=mapsize[1],
                                n=mapsize[0],
                                dim=np.shape(data)[1],
                                max_epochs=epochs,
                                gpus=1,
                                session=session,
                                graph=graph,
                                input_tensor=next_element,
                                batch_size=128,
                                initial_learning_rate=0.1)

        session.run([tf.global_variables_initializer()])

        som.train(num_inputs=np.shape(data)[0])
    return np.reshape(som.output_weights, (mapsize[1], mapsize[0], np.shape(data)[1]))


def get_somoclu(data, mapsize, epochs=10):
    som = somoclu.Somoclu(*mapsize)
    # if type(data) == 'pandas.core.frame.DataFrame':
    data = np.array([list(data[i]) for i in data.columns]).astype(np.float32)
    som.train(data=data, epochs=epochs)

    '''
    codebook = np.reshape(som.codebook, (mapsize[0]*mapsize[1], np.shape(data)[0],))
    QE = []
    TE = 0
    dist = []
    for sample in data.columns:
        for idx, w in enumerate(codebook):
            dist.append(np.linalg.norm(w-data[sample]))

    min_idx = dist.index(min(dist))
    QE.append(dist.pop(min_idx))
    if np.abs(min_idx-dist.index(min(dist))) > 1:
        TE += 1

    print "QE = ", 1./(np.shape(data)[1]) * sum(QE)
    print "TE = ", 1./(np.shape(data)[1]) * TE
    '''
    return som.codebook


def plot_som(som, data, labels=None):
    data = np.array([list(data[i]) for i in data.columns]).astype(np.float32)
    umatrix = _get_umatrix(som)
    fig = plt.figure()
    plt.imshow(umatrix, cmap='Spectral_r', origin='lower')

    weights = som.reshape(-1, som.shape[-1])

    for i, sample in enumerate(data):
        d = 1e6
        cell = 0
        for idx, w in enumerate(weights):
            dist = np.linalg.norm(w-sample)
            if dist < d:
                d = dist
                cell = idx

        if labels:
            plt.scatter([cell % np.shape(som)[1]], [cell // np.shape(som)[1]],
                        color=labels[i])
        else:
            plt.scatter([cell % np.shape(som)[1]], [cell // np.shape(som)[1]],
                        color='white')  # , color=colors[i])

    '''
    points = som.bmus
    cenx, ceny = mapsize[1]/2., mapsize[0]/2.
    rad = 1.1 * np.sqrt(pow(mapsize[0], 2)+pow(mapsize[1], 2))/2.
    nump = 100
    for i in range(nump):
        angle = i * 2*np.pi/nump
        points = np.concatenate(
            (points, [[rad*np.cos(angle)+cenx, rad*np.sin(angle)+ceny]]), axis=0)

    vor = ss.Voronoi(points)
    for r, p in enumerate(vor.point_region):
        region = vor.regions[p]
        if not -1 in region and r < len(som.bmus):
            polygon = [vor.vertices[i] for i in region]
            plt.fill(*zip(*polygon), alpha=0.4, color=labels[r])

    plt.xlim(0, mapsize[1]-1)
    plt.ylim(0, mapsize[0]-1)
    '''


def _get_umatrix(weights):
    m, n, _ = np.shape(weights)
    weights = np.reshape(weights, (m*n, -1))

    umatrix = np.zeros((m * n, 1))
    neuron_locs = list()
    for i in range(m):
        for j in range(n):
            neuron_locs.append(np.array([i, j]))
    neuron_distmat = ss.distance_matrix(neuron_locs, neuron_locs)

    for i in range(m * n):
        # Get the indices of the units which neighbor i
        # Change this to `< 2` if you want to include diagonal neighbors
        neighbor_idxs = neuron_distmat[i] <= 1

        neighbor_weights = weights[neighbor_idxs]
        # Get the average distance between unit i and all of its neighbors
        # Expand dims to broadcast to each of the neighbors
        umatrix[i] = ss.distance_matrix(np.expand_dims(weights[i], 0), neighbor_weights).mean()

    return umatrix.reshape((m, n))
