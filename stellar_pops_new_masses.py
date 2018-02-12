import sys
import os

sys.path.insert(0, '/home/sinandeger/PycharmProjects/SinanDeger/Projects/SinanLibrary/')

from mergfractools import *
from structure_select import *

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, NullFormatter
import matplotlib.patches as patches

from sklearn import svm
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

import xgboost as xgb
from xgboost import XGBClassifier, plot_tree

import time
import aplpy

import scipy
from scipy import stats

import scipy.interpolate

"""########################"""

specz_y_intercept = -0.97
specz_slope = -0.87

""""#######################"""

PSF_FWHM = 0.18

"""Below we select the spectroscopic galaxy clusters and the field"""

structures = spec_sample_select()

cl1040 = structures['cl1040']
cl1054_1146 = structures['cl1054_1146']
cl1054_1245 = structures['cl1054_1245']
cl1138_cluster = structures['cl1138_cluster']
cl1138a_cluster = structures['cl1138a_cluster']
cl1227_cluster = structures['cl1227_cluster']
cl1227a_cluster = structures['cl1227a_cluster']
cl1216 = structures['cl1216']
cl1232 = structures['cl1232']
cl1354_cluster = structures['cl1354_cluster']
cl1354a_cluster = structures['cl1354a_cluster']

"""Below we select the galaxy groups"""

cl1037_group = structures['cl1037_group']
cl1040a_group = structures['cl1040a_group']
cl1040b_group = structures['cl1040b_group']
cl1054_1146a_group = structures['cl1054_1146a_group']
cl1054_1245a_group = structures['cl1054_1245a_group']
cl1103a_group = structures['cl1103a_group']
cl1103b_group = structures['cl1103b_group']

"Below we select the field"

spec_field = structures['Field']

spec_cls = [cl1040, cl1054_1146, cl1054_1245, cl1138_cluster, cl1138a_cluster, cl1227_cluster, cl1216, cl1232, cl1354_cluster, cl1354a_cluster]

grp = [cl1037_group, cl1040a_group, cl1040b_group, cl1054_1146a_group, cl1054_1245a_group, cl1103a_group, cl1103b_group]

groups = pd.concat(grp)

spec_clusters = pd.concat(spec_cls)

spec_str = [groups, spec_clusters, spec_field]
spec_entire = pd.concat(spec_str)

"""Spec structures split into visual TIM and undisturbed"""

spec_clusters_tim = spec_clusters.where(spec_clusters['class_D'] != '0').dropna()
spec_clusters_undisturbed = spec_clusters.where(spec_clusters['class_D'] == '0').dropna()

groups_tim = groups.where(groups['class_D'] != '0').dropna()
groups_undisturbed = groups.where(groups['class_D'] == '0').dropna()

spec_field_tim = spec_field.where(spec_field['class_D'] != '0').dropna()
spec_field_undisturbed = spec_field.where(spec_field['class_D'] == '0').dropna()

"""U-V vs V-J plots"""

ax = plt.gca()

xminorLocator = plt.MultipleLocator(0.1)
xmajorLocator = plt.MultipleLocator(0.5)
yminorLocator = plt.MultipleLocator(0.1)
ymajorLocator = plt.MultipleLocator(0.5)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.yaxis.set_major_locator(ymajorLocator)

ax.xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%0.1f'))

ax.spines['bottom'].set_linewidth(2.0)
ax.spines['left'].set_linewidth(2.0)
ax.spines['top'].set_linewidth(2.0)
ax.spines['right'].set_linewidth(2.0)

plt.tick_params(which='both', width=1.4, labelsize=14)
plt.tick_params(which='major', length=8)
plt.tick_params(which='minor', length=5)

csfont = {'fontname': 'serif'}
hfont = {'fontname': 'monospace'}

ax.text(0.52, 2.35, '0.4 < z < 0.8', fontsize=16, **csfont)
ax.text(0.6, 2.0, 'Quiescent', fontsize=16, **csfont)
ax.text(1.4, 0.9, 'Star Forming', fontsize=16, **csfont)

plt.xlabel('V-J', labelpad=-3, fontsize=16, **csfont)
plt.ylabel('U-V', labelpad=3, fontsize=16, **csfont)

plt.ylim([0.5, 2.5])
plt.xlim([0.5, 2.0])

"""Clusters+Groups"""

# plt.scatter(spec_clusters_undisturbed['VJ_ZCLUST'],
#             spec_clusters_undisturbed['UV_ZCLUST'], color='darkslateblue', marker='*', label='Undisturbed')
#
# plt.scatter(groups_undisturbed['VJ_ZCLUST'],
#             groups_undisturbed['UV_ZCLUST'], color='darkslateblue', marker='*', label='_nolabel_')
#
# plt.scatter(spec_clusters['VJ_ZCLUST'].iloc[pandas_above_line(spec_clusters_tim['M20'],
#                                                               spec_clusters_tim['G'],
#                                                               specz_y_intercept, specz_slope)],
#             spec_clusters['UV_ZCLUST'].iloc[pandas_above_line(spec_clusters_tim['M20'],
#                                                               spec_clusters_tim['G'],
#                                                               specz_y_intercept, specz_slope)],
#             color='crimson', marker='o', label='$G-M_{20}$ TIM')
#
#
# plt.scatter(groups['VJ_ZCLUST'].iloc[pandas_above_line(groups_tim['M20'],
#                                                               groups_tim['G'],
#                                                               specz_y_intercept, specz_slope)],
#             groups['UV_ZCLUST'].iloc[pandas_above_line(groups_tim['M20'],
#                                                               groups_tim['G'],
#                                                               specz_y_intercept, specz_slope)],
#             color='crimson', marker='o', label='_nolabel_')

"""Field"""

plt.scatter(spec_field_undisturbed['VJ_ZS'],
            spec_field_undisturbed['UV_ZS'], color='darkslateblue', marker='*', label='Undisturbed')

plt.scatter(spec_field['VJ_ZS'].iloc[pandas_above_line(spec_field_tim['M20'],
                                                       spec_field_tim['G'],
                                                              specz_y_intercept, specz_slope)],
            spec_field['UV_ZS'].iloc[pandas_above_line(spec_field_tim['M20'],
                                                       spec_field_tim['G'],
                                                       specz_y_intercept, specz_slope)],
            color='crimson', marker='o', label='$G-M_{20}$ TIM')

plt.legend(loc='lower right')

plt.savefig('UVJ-Field-only.pdf', format='pdf')
plt.show()

