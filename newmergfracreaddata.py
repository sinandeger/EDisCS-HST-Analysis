'''
Created on Feb 21, 2016

@author: sinandeger
'''

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

import seaborn as sns
from matplotlib import rc

sns.set_color_codes(palette='deep')

#sns.palplot(sns.color_palette("Paired", 10))
sns.set(font_scale=1.5)
sns.set_style("ticks")

rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

from mergfractools import * 
from plot_it import *

data = ascii.read('/home/sinandeger/Research/MasterCatalog/index+MasterCatalog-GroupCat-match-FillNa')
hst_data = ascii.read('/home/sinandeger/Research/MasterCatalog/HST-ACS-clusters-only-Spectral')

f1 = open('StellarMass-z.txt','w')
f2 = open('merger_picks.txt','w')

length = len(data)
print length

length_hst = len(hst_data)
print length_hst

G_merger = []
M20_merger = []

G_non = []
M20_non = []

G_str = []
M20_str = []
Index_no_str = []

PSF_FWHM = 0.18

"""Slope and y-intercept values of our line"""

# line_slope = 0.28
# Y_intercept = 0.03

line_slope = 0.32
Y_intercept = 0.0

my_dpi = 96 * 1.1

count = 0
merger_above = 0
ntot_above = 0

"""Below configuration selects group members"""

#(data['membflag_1'][i] == '1a' or
# data['membflag_1'][i] == '1b' or
# data['membflag_1'][i] == '1A' or
# data['membflag_1'][i] == '1B' or
# data['qual'][i] == 1):

################

"""Below configuration selects field"""

#data['membflag_1'][i] == '0' and data['qual'][i] != 1 and data['qual'][i] != 0
#data['ZSPEC'][i] - 0.2 <= data['ZCLUST'][i] <= data['ZSPEC'][i] + 0.2 and
#data['ZSPEC'][i] < 0.9 and

#data['membflag_1'][i] == '0' and
#data['qual'][i] == -999 and
#data['ZSPEC'][i] - 0.2 <= data['ZCLUST'][i] <= data['ZSPEC'][i] + 0.2 and
#data['ZSPEC'][i] < 0.65 and
#data['ZSPEC'][i] > 0.4 and


################

#     if (data['membflag_1'][i] == '1a' or
#         data['membflag_1'][i] == '1b' or
#         data['membflag_1'][i] == '1A' or
#         data['membflag_1'][i] == '1B' or
#         data['qual'][i] == 1):
#         
#         print data['membflag_1'][i]
#         
#         print >>f1, masscomp_clust(data['BV_ZCLUST'][i],data['LVBEST_ZCLUST'][i])

        #data['qual'][i] != 1 and

"""Below config picks all cluster members"""

# if (data['membflag_1'][i] != '0' and
#             data['qual'][i] == -999 and
#             data['WMIN'][i] > 0.3 and
#             data['FLAG'][i] == 0 and
#             data['<S/N>'][i] > 2.5 and
#             data['Rpet_ell'][i] > 2 * PSF_FWHM and
#             data['ITOT'][i] < 24.5 and
#             data['CLUSTNAME_1'][i] != 'cl1037-1243' and
#             data['CLUSTNAME_1'][i] != 'cl1103-1245' and
#             data['ZSPEC'][i] > 0.67 and
#             10.4 < masscomp_clust(data['BV_ZCLUST'][i], data['LVBEST_ZCLUST'][i]) < 12):

#data['CLUSTNAME_1'][i] == 'cl1354-1230' and
#data['membflag_1'][i] == '1a' and
#data['qual'][i] == 1 and
#data['ZSPEC'][i] <= 0.67 and
#data['CLUSTNAME_1'][i] != 'cl1037-1243' and
        

#data['CLUSTNAME_1'][i] != 'cl1232-1250' and
#data['CLUSTNAME_1'][i] != 'cl1227-1138' and
#        data['qual'][i] == -999 and

clust_index = 0
itera = 0

for i in range(length):

    if ((data['membflag_1'][i] != '0' or data['qual'][i] == 1) and
        data['WMIN'][i] > 0.3 and
        data['FLAG'][i] == 0 and
        data['<S/N>'][i] > 2.5 and
        data['Rpet_ell'][i] > 2*PSF_FWHM and
        data['ITOT'][i] < 24.5 and
        10.4 < masscomp_clust(data['BV_ZCLUST'][i],data['LVBEST_ZCLUST'][i]) < 12):
          
        print >>f1, data['CLUSTNAME_1'][i], masscomp_clust(data['BV_ZCLUST'][i], data['LVBEST_ZCLUST'][i]), data['ZCLUST'][i], data['ZSPEC'][i], data['membflag_1'][i]
          
        G_str.append(data['G'][i])
        M20_str.append(data['M20'][i])
        Index_no_str.append(data['index_no'][i])

        count += 1
        
        # if data['CLUSTNAME_1'][i] == hst_data['CLUSTNAME'][clust_index]:
        #
        #     Del_RA = data['RA_2'][i] - hst_data['BCG_RA'][clust_index]
        #     Del_Dec = data['DEC_2'][i] - hst_data['BCG_Dec'][clust_index]
        #
        #     Coord = np.sqrt(np.power(Del_RA,2)+np.power(Del_Dec,2))
        #
        #     R200 = R200_indegrees(hst_data['R200'][clust_index], hst_data['Scale_at_z'][clust_index])
        #
        #     #print R200, Coord, clust_index, i, data['class_D'][i], data['CLUSTNAME_1'][i], hst_data['CLUSTNAME'][clust_index]
        #
        #     if Coord <= 0.5*R200:
        #
        #         count += 1
        #
        #         G_str.append(data['G'][i])
        #         M20_str.append(data['M20'][i])
        #         Index_no_str.append(data['index_no'][i])
                
        if data['class_D'][i] == '0':

            G_non.append(float(data['G'][i]))
            M20_non.append(float(data['M20'][i]))

            if above_line(data['M20'][i], data['G'][i], Y_intercept, line_slope) == 1:

                ntot_above += 1

        if data['class_D'][i] != '0':

            G_merger.append(float(data['G'][i]))
            M20_merger.append(float(data['M20'][i]))

            if above_line(data['M20'][i], data['G'][i], Y_intercept, line_slope) == 1:

                merger_above += 1
                ntot_above += 1

        # else:
        #
        #     clust_index += 1
        #     print clust_index
        #
        #     if clust_index == length_hst:
        #
        #         clust_index = 0
        #
        #         itera += 1
        #         print 'iteration', itera
        #
        #         if itera == 3:
        #
        #             continue
        #             print 'Code break'
        #
        #             itera = 0
        
        
f1.close()
print len(G_str)

len_boot = len(G_str)
 
iter_num = 10000
 
Merg_fracs = []
 
for z in range(iter_num):
 
    Index_no_boot = np.zeros(len_boot)
         
    n_merg_above = 0.0
     
    Index_no_boot = bootstrap_1_array(Index_no_str,len_boot)
     
    for w in range(len_boot):
 
            index_temp = Index_no_boot[w]-1
            #print index_temp
                  
            if (above_line(data['M20'][index_temp], data['G'][index_temp], Y_intercept, line_slope) == 1 and 
                data['class_D'][index_temp] != '0'):
                      
                n_merg_above += 1.0
 
    print >>f2, n_merg_above
     
    MF = float(n_merg_above/len_boot)
     
    Merg_fracs.append(MF)
 
 
# print Merg_fracs
print len(Merg_fracs)
 
plt.hist(Merg_fracs)
 
Conf_lim_merg_frac, lower_limit, upper_limit = conf_interval(Merg_fracs, 68)
 
#print Conf_lim_merg_frac
print len(Conf_lim_merg_frac)
print 'median:',np.median(Conf_lim_merg_frac),'mean:',np.mean(Conf_lim_merg_frac),'st dev:',np.std(Conf_lim_merg_frac)
print "lower limit:", lower_limit, "upper limit:", upper_limit
 
plt.hist(Conf_lim_merg_frac)
plt.xlabel('Merger Fraction')
plt.ylabel('N')
 
plt.show()

G2_merger = np.zeros(len(G_merger),dtype = float)
M20_2_merger = np.zeros(len(M20_merger),dtype = float)

G2_merger = G_merger
M20_2_merger = M20_merger

G2_non = np.zeros(len(G_non),dtype = float)
M20_2_non = np.zeros(len(M20_non),dtype = float)

G2_non = G_non
M20_2_non = M20_non

print 'Total # of objects in plot:', count
print 'Total # of mergers above our line', merger_above
print 'Total # of objects above our line', ntot_above

#print G2_merger
#print M20_2_merger

ax = plt.gca()

xminorLocator = plt.MultipleLocator(0.2)
xmajorLocator = plt.MultipleLocator(1.0)
yminorLocator = plt.MultipleLocator(0.05)
ymajorLocator = plt.MultipleLocator(0.1)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)
ax.yaxis.set_major_locator(ymajorLocator)

ax.xaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%0.2f'))

ax.spines['bottom'].set_linewidth(2.0)
ax.spines['left'].set_linewidth(2.0)
ax.spines['top'].set_linewidth(2.0)
ax.spines['right'].set_linewidth(2.0)

plt.tick_params(which='both', width=2.0, labelsize=16)
plt.tick_params(which='major', length=8)
plt.tick_params(which='minor', length=5)

csfont = {'fontname':'Serif'}
hfont = {'fontname':'Serif'}

#plt.title('$Field$',fontsize=20,**hfont)
#plt.text('Field', loc='upper left', frameon=True)

plt.xlabel('$M_{20}$',labelpad=5,fontsize=20, **csfont)
plt.ylabel('G',labelpad=20, fontsize=23, rotation=0, **csfont)

# plot_it(M20_2_merger,G2_merger,'$M_{20}$','$G$','AAS cl1040 - Mass Complete','rs')
# plot_it(M20_2_non,G2_non,'$M_{20}$','$G$','AAS cl1040 - Mass Complete','bo')

p3, = plt.plot([0,-3],[0.33,0.75],'g', lw=1.5)
p4, = plt.plot([0,-3],[0.0,0.96],'r', lw=1.5)

p2, = plt.plot(M20_2_non,G2_non,'bo')
p1, = plt.plot(M20_2_merger,G2_merger,'rs')

plt.legend((p1,p2,p3,p4), ("Merger", "Undisturbed", "Lotz Line", "Our Line"), loc='lower right', frameon=True)


plt.xlim([-0.5,-3])
plt.ylim([0.3,0.7])

plt.savefig('Clusters-Field.pdf', format='pdf', dpi=my_dpi)

plt.show()
