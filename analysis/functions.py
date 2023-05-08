# Defining Analysis Functions

import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Dihedral
from MDAnalysis.analysis.dihedrals import Ramachandran
from MDAnalysis.core.topologyattrs import Atomnames
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from shapely.geometry import Polygon, Point
import sys
import csv
import os
import datetime


def smooth_diff(RG, kernel_size):

    data = np.absolute(np.gradient(RG, 1))
    kernel = np.ones(kernel_size) / kernel_size
    smooth_diff = np.convolve(data, kernel, mode="same")

    return smooth_diff

def smooth_data(data, kernel_size):

    kernel = np.ones(kernel_size) / kernel_size
    smooth_data = np.convolve(data, kernel, mode="same")

    return smooth_data

