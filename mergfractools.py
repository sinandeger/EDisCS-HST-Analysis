'''
Created on Feb 23, 2016

@author: sinandeger
'''

import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd


def masscomp_clust(BV_zclust, LVbest_zclust):
    h = 0.7

    lgmlv = -0.734 + (1.404 * (BV_zclust + 0.084))

    Mstar_zclust = np.power(10, lgmlv) * LVbest_zclust * 1e10 / np.power(h, 2)

    K = np.log10(Mstar_zclust)

    return K


"""masscomp_clust finds the stellar masses for cluster members"""


def masscomp_field(BV_ZS, LVBEST_ZS, BV_ZP, LVBEST_ZP):
    h = 0.7

    if BV_ZS != -999:

        lgmlv = -0.734 + (1.404 * (BV_ZS + 0.084))

        Mstar_zclust = np.power(10, lgmlv) * LVBEST_ZS * 1e10 / np.power(h, 2)

        K = np.log10(Mstar_zclust)

    else:

        lgmlv = -0.734 + (1.404 * (BV_ZP + 0.084))

        Mstar_zclust = np.power(10, lgmlv) * LVBEST_ZP * 1e10 / np.power(h, 2)

        K = np.log10(Mstar_zclust)

    return K


"""masscomp_field finds the stellar masses for cluster members"""


def pickstar(STARFLAG, FLUX_RADIUS, MAG_AUTO):
    sizecut = 0.7 + 3 * (0.04)

    if MAG_AUTO < 22:

        if STARFLAG == 1 or FLUX_RADIUS <= sizecut:
            return 1

    if MAG_AUTO > 22:

        if FLUX_RADIUS < sizecut:
            return 1


"""pickstar will choose stars accordingly and return 1. If you want to exclude stars, have pickstar != 1"""


def PhotMember(WMIN, SPECMEMB, MEMBFLAG):
    if WMIN > 0.3:

        if SPECMEMB == 1 | (SPECMEMB != 0 & MEMBFLAG == 1):
            return 1


"""PhotMember will pick photometric members, require PhotMember = 1 to pick"""

"""Use PhotMember_pandas with pandas calls"""


def PhotMember_pandas(SPECMEMB, MEMBFLAG):
    return np.where((SPECMEMB == 1) | ((SPECMEMB != 0) & (MEMBFLAG == 1)))


def eliminate_star_pandas(STARFLAG, FLUX_RADIUS, MAG_AUTO):
    sizecut = 0.7 + 3 * 0.04

    return np.where(
        ((MAG_AUTO < 22) & ((STARFLAG == 0) | (FLUX_RADIUS >= sizecut))) | ((MAG_AUTO >= 22) & (FLUX_RADIUS > sizecut)))


def PhotField(WMIN, SPECMEMB, MEMBFLAG):
    if WMIN > 0.3:

        if SPECMEMB == 0 | (SPECMEMB != 1 & MEMBFLAG == 0):
            return 1


"""PhotField will pick field objects, require PhotField == 1 to pick"""


def R200_indegrees(R200_inMPC, scale_at_z_in_kpc_arcsec):
    R200_inarcsec = (1000.0 * R200_inMPC) / (scale_at_z_in_kpc_arcsec)
    R200_in_degrees = R200_inarcsec / 3600.0

    return R200_in_degrees


"""Careful with the units for the function R200_indegrees. Needs R200 in MPC, scale at kpc/" to work correctly"""


def bootstrap_2_arrays(A1, A2, A3, iter_num):
    array_ind = len(A1)

    B1 = []
    B2 = []
    B3 = []

    for j in range(iter_num):
        Rand_int = random.randint(0, array_ind - 1)

        B1.append(A1[Rand_int])
        B2.append(A2[Rand_int])
        B3.append(A3[Rand_int])

    return B1, B2, B3


"""bootstrap_GM20 will pick 2 arrays, A1 & A2, and a number of iteration. It will return 2 arrays, with length"""
"""equal to iter_num. Their elements are random picks of the elements of A1 and A2"""
"""Update: A3 will pick index numbers and the function will also return an array of index numbers corresponding"""
"""to the G-M20 values"""


def bootstrap_1_array(A, array_len):

    #array_ind = len(A)

    B = np.zeros(array_len, dtype=int)

    for j in range(array_len):
        Rand_int = random.randint(0, array_len - 1)

        B[j] = A[Rand_int]

    return B

#
# def above_line(X_coord, Y_coord, Y_intercept, line_slope):
#     point_val = (np.abs(Y_coord) - Y_intercept) / np.abs(X_coord)
#
#     if point_val > line_slope:
#
#         return 1
#
#     else:
#
#         return 0
#
# """ If above_line = 1, point is above the line"""
#
# def below_line(X_coord, Y_coord, Y_intercept, line_slope):
#     point_val = (np.abs(Y_coord) - Y_intercept) / np.abs(X_coord)
#
#     if point_val > line_slope:
#
#         return 0
#
#     else:
#
#         return 1

""" If below_line == 1, point is below the line """

"""Vectorized implementation of above line"""

def pandas_above_line(X_coord, Y_coord, Y_intercept, line_slope):

    y_one = Y_intercept
    y_two = 1.5

    x_one = (y_one-Y_intercept)/line_slope
    x_two = (y_two-Y_intercept)/line_slope

    point_val = (X_coord - x_one)*(y_two - y_one) - (Y_coord - y_one)*(x_two - x_one)

    return np.where(point_val > 0.0)

""" If pandas_above_line = 1, point is above the line"""

"""Use the below function to check if a single entry is above the line"""

def pandas_above_line_postage(X_coord, Y_coord, Y_intercept, line_slope):

    y_one = Y_intercept
    y_two = 1.5

    x_one = (y_one-Y_intercept)/line_slope
    x_two = (y_two-Y_intercept)/line_slope

    point_val = (X_coord - x_one)*(y_two - y_one) - (Y_coord - y_one)*(x_two - x_one)

    if point_val > 0.0:

        return 1

    else:

        return  0

"""##########################"""


def pandas_below_line(X_coord, Y_coord, Y_intercept, line_slope):

    y_one = Y_intercept
    y_two = 1.5

    x_one = (y_one-Y_intercept)/line_slope
    x_two = (y_two-Y_intercept)/line_slope

    point_val = (X_coord - x_one)*(y_two - y_one) - (Y_coord - y_one)*(x_two - x_one)

    return np.where(point_val < 0.0)

""" If pandas_below_line = 1, point is above the line"""

""" If below_line == 1, point is below the line """

#
# def pandas_above_line(X_coord, Y_coord, Y_intercept, line_slope):
#     point_val = (np.abs(Y_coord) - Y_intercept) / np.abs(X_coord)
#
#     return np.where((point_val > line_slope))
#
# """ If pandas_above_line = 1, point is above the line"""
#
#
# def pandas_below_line(X_coord, Y_coord, Y_intercept, line_slope):
#     point_val = (np.abs(Y_coord) - Y_intercept) / np.abs(X_coord)
#
#     return np.where((point_val < line_slope))
#
# """ If pandas_below_line = 1, point is above the line"""


def conf_interval(A1, conf_limit):
    A2 = sorted(A1)
    A3 = []

    length = len(A2)

    for j in range(int(((100.0 - conf_limit) / 200.0) * length),
                   int(length - (((100.0 - conf_limit) / 200.0) * length))):
        A3.append(A2[j])

    return A3, min(A3), max(A3)


""" Below are the functions specific for use with pandas. They are written to work with vectors instead of array elements """


def pandas_pick_spec_field(membflag, cluster_z, galaxy_z, qual):

    return np.where(
        (membflag == 0) & (0.4 < galaxy_z) & (galaxy_z < 0.8) & (qual == -999) & (galaxy_z < cluster_z + 0.2) &
        (galaxy_z > cluster_z - 0.2))


def PhotField_pandas(SPECMEMB, MEMBFLAG, cluster_z, galaxy_spec_z, galaxy_phot_z, qual):

    return np.where(((SPECMEMB == 0) | ((SPECMEMB != 1) & (MEMBFLAG == 0)))
                    & (((galaxy_spec_z != -999) & (galaxy_spec_z < cluster_z + 0.2) & (galaxy_spec_z > cluster_z - 0.2) & (0.4 < galaxy_spec_z) & (galaxy_spec_z < 0.8) & (qual == -999))
                    | ((galaxy_spec_z == -999) & (galaxy_phot_z < 0.8) & (galaxy_phot_z > 0.5))))

# def PhotField_pandas(SPECMEMB, MEMBFLAG, cluster_z, galaxy_spec_z, galaxy_phot_z):
#
#     return np.where(((SPECMEMB == 0) | ((SPECMEMB != 1) & (MEMBFLAG == 0)))
#                     & (((galaxy_spec_z != -999) & (galaxy_spec_z < cluster_z + 0.1) & (galaxy_spec_z > cluster_z - 0.1))
#                     | ((galaxy_spec_z == -999) & (galaxy_phot_z < cluster_z + 0.1) & (galaxy_phot_z > cluster_z - 0.1))))


#def pandas_phot_field_stellarmass(BV_ZS, LVBEST_ZS, BV_ZP, LVBEST_ZP):

def pandas_phot_field_stellarmass(BV_ZP, LVBEST_ZP):

    h = 0.7

    lgmlv = -0.734 + (1.404 * (BV_ZP + 0.084))

    Mstar_zclust = np.power(10, lgmlv) * LVBEST_ZP * 1e10 / np.power(h, 2)

    K = np.log10(Mstar_zclust)

    return K


def pandas_spec_field_stellarmass(BV_ZS, LVBEST_ZS):

    h = 0.7

    lgmlv = -0.734 + (1.404 * (BV_ZS + 0.084))

    Mstar_zclust = np.power(10, lgmlv) * LVBEST_ZS * 1e10 / np.power(h, 2)

    K = np.log10(Mstar_zclust)

    return K

# def above_line_new_method(X_coord, Y_coord, Y_intercept, line_slope):
#     point_on_line_val = line_slope*X_coord + Y_intercept
#
#     if point_on_line_val < Y_coord:
#
#         return 1
#
#     else:
#
#         return 0



"""Create region files from pandas dataframes. Use "text={...}" to add text to selected objects"""

def pandas_region_file(RA, DEC, arcsec, name):

    reg_file = open(name + '.reg', 'w')

    print >> reg_file, '# Region file format: DS9 version 4.1\n',\
        'global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n',\
        'fk5\n'

    for k in range(len(RA)):

        print >>reg_file, 'circle(', RA.iloc[k], ',', DEC.iloc[k], ',', str(arcsec) + '"', ')'

    reg_file.close()
