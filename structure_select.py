import numpy as np
import random
import pandas as pd


def PhotMember_pandas(SPECMEMB, MEMBFLAG):
    return np.where((SPECMEMB == 1) | ((SPECMEMB != 0) & (MEMBFLAG == 1)))

def pandas_pick_spec_field(membflag, cluster_z, galaxy_z, qual):

    return np.where(
        (membflag == 0) & (0.4 < galaxy_z) & (galaxy_z < 0.8) & (qual == -999) & (galaxy_z < cluster_z + 0.2) &
        (galaxy_z > cluster_z - 0.2))


"""Here we define structures we use for TIM analysis of EDisCS sample"""

# = pd.read_csv('/home/sinandeger/Research/MasterCatalog/Spectroscopic-Catalog(Sinan-Edit).csv')

spec_sample_df = pd.read_csv('/home/sinandeger/Research/Moustakas-EDisCS-Catalog/At_spec_z/Mous_at_spec_z-Sinan_cat_match.csv')

PSF_FWHM = 0.18

def spec_sample_select():

    cl1037_group = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                        (spec_sample_df['FLAG'] == 0) &
                                        (spec_sample_df['<S/N>'] > 2.5) &
                                        (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                        (spec_sample_df['ITOT'] < 24.5) &
                                        (spec_sample_df['clustname_2'] == 'cl1037.9-1243') &
                                        (spec_sample_df['qual'] == 1) &
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    cl1040 = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                  (spec_sample_df['FLAG'] == 0) &
                                  (spec_sample_df['<S/N>'] > 2.5) &
                                  (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                  (spec_sample_df['ITOT'] < 24.5) &
                                  (spec_sample_df['CLUSTFULLNAME'] == 'cl1040.7-1155') &
                                  (spec_sample_df['SPECMEMB'] == 1) &
                                  (spec_sample_df['qual'] == -999)&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    cl1040a_group = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                         (spec_sample_df['FLAG'] == 0) &
                                         (spec_sample_df['<S/N>'] > 2.5) &
                                         (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                         (spec_sample_df['ITOT'] < 24.5) &
                                         (spec_sample_df['clustname_2'] == 'cl1040.7-1155a') &
                                         (spec_sample_df['ZSPEC'] < 0.7) &
                                         (spec_sample_df['qual'] == 1)&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    cl1040b_group = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                         (spec_sample_df['FLAG'] == 0) &
                                         (spec_sample_df['<S/N>'] > 2.5) &
                                         (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                         (spec_sample_df['ITOT'] < 24.5) &
                                         (spec_sample_df['clustname_2'] == 'cl1040.7-1155b') &
                                         (spec_sample_df['ZSPEC'] > 0.71) &
                                         (spec_sample_df['qual'] == 1)&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    cl1054_1146 = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                       (spec_sample_df['FLAG'] == 0) &
                                       (spec_sample_df['<S/N>'] > 2.5) &
                                       (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                       (spec_sample_df['ITOT'] < 24.5) &
                                       (spec_sample_df['CLUSTFULLNAME'] == 'cl1054.4-1146') &
                                       (spec_sample_df['SPECMEMB'] == 1) &
                                       (spec_sample_df['qual'] == -999)&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    cl1054_1146a_group = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                              (spec_sample_df['FLAG'] == 0) &
                                              (spec_sample_df['<S/N>'] > 2.5) &
                                              (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                              (spec_sample_df['ITOT'] < 24.5) &
                                              (spec_sample_df['clustname_2'] == 'cl1054.4-1146a') &
                                              (spec_sample_df['qual'] == 1)&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    cl1054_1245 = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                       (spec_sample_df['FLAG'] == 0) &
                                       (spec_sample_df['<S/N>'] > 2.5) &
                                       (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                       (spec_sample_df['ITOT'] < 24.5) &
                                       (spec_sample_df['CLUSTFULLNAME'] == 'cl1054.7-1245') &
                                       (spec_sample_df['SPECMEMB'] == 1) &
                                       (spec_sample_df['qual'] == -999)&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    cl1054_1245a_group = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                              (spec_sample_df['FLAG'] == 0) &
                                              (spec_sample_df['<S/N>'] > 2.5) &
                                              (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                              (spec_sample_df['ITOT'] < 24.5) &
                                              (spec_sample_df['clustname_2'] == 'cl1054.7-1245a') &
                                              (spec_sample_df['qual'] == 1)&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    cl1103a_group = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                         (spec_sample_df['FLAG'] == 0) &
                                         (spec_sample_df['<S/N>'] > 2.5) &
                                         (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                         (spec_sample_df['ITOT'] < 24.5) &
                                         (spec_sample_df['clustname_2'] == 'cl1103.7-1245a') &
                                         (spec_sample_df['qual'] == 1)&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    cl1103b_group = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                         (spec_sample_df['FLAG'] == 0) &
                                         (spec_sample_df['<S/N>'] > 2.5) &
                                         (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                         (spec_sample_df['ITOT'] < 24.5) &
                                         (spec_sample_df['clustname_2'] == 'cl1103.7-1245b') &
                                         (spec_sample_df['qual'] == 1)&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    cl1138_cluster = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                          (spec_sample_df['FLAG'] == 0) &
                                          (spec_sample_df['<S/N>'] > 2.5) &
                                          (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                          (spec_sample_df['ITOT'] < 24.5) &
                                          (spec_sample_df['CLUSTFULLNAME'] == 'cl1138.2-1133') &
                                          (spec_sample_df['membflag_1'] == '1')&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    cl1138a_cluster = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                           (spec_sample_df['FLAG'] == 0) &
                                           (spec_sample_df['<S/N>'] > 2.5) &
                                           (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                           (spec_sample_df['ITOT'] < 24.5) &
                                           (spec_sample_df['CLUSTFULLNAME'] == 'cl1138.2-1133') &
                                           (spec_sample_df['membflag_1'] == '1a')&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    cl1227_cluster = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                          (spec_sample_df['FLAG'] == 0) &
                                          (spec_sample_df['<S/N>'] > 2.5) &
                                          (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                          (spec_sample_df['ITOT'] < 24.5) &
                                          (spec_sample_df['CLUSTFULLNAME'] == 'cl1227.9-1138') &
                                          (spec_sample_df['membflag_1'] == '1')&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    cl1227a_cluster = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                           (spec_sample_df['FLAG'] == 0) &
                                           (spec_sample_df['<S/N>'] > 2.5) &
                                           (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                           (spec_sample_df['ITOT'] < 24.5) &
                                           (spec_sample_df['CLUSTFULLNAME'] == 'cl1227.9-1138') &
                                           (spec_sample_df['membflag_1'] == '1a')&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    cl1216 = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                  (spec_sample_df['FLAG'] == 0) &
                                  (spec_sample_df['<S/N>'] > 2.5) &
                                  (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                  (spec_sample_df['ITOT'] < 24.5) &
                                  (spec_sample_df['CLUSTFULLNAME'] == 'cl1216.8-1201') &
                                  (spec_sample_df['SPECMEMB'] == 1)&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    cl1232 = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                  (spec_sample_df['FLAG'] == 0) &
                                  (spec_sample_df['<S/N>'] > 2.5) &
                                  (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                  (spec_sample_df['ITOT'] < 24.5) &
                                  (spec_sample_df['CLUSTFULLNAME'] == 'cl1232.5-1250') &
                                  (spec_sample_df['SPECMEMB'] == 1)&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    cl1354_cluster = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                          (spec_sample_df['FLAG'] == 0) &
                                          (spec_sample_df['<S/N>'] > 2.5) &
                                          (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                          (spec_sample_df['ITOT'] < 24.5) &
                                          (spec_sample_df['CLUSTFULLNAME'] == 'cl1354.2-1230') &
                                          (spec_sample_df['membflag_1'] == '1')&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    cl1354a_cluster = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                           (spec_sample_df['FLAG'] == 0) &
                                           (spec_sample_df['<S/N>'] > 2.5) &
                                           (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                           (spec_sample_df['ITOT'] < 24.5) &
                                           (spec_sample_df['CLUSTFULLNAME'] == 'cl1354.2-1230') &
                                           (spec_sample_df['membflag_1'] == '1a')&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()

    """Field"""

    spec_sample_field_base = spec_sample_df.where((spec_sample_df['WMIN'] > 0.3) &
                                                  (spec_sample_df['FLAG'] == 0) &
                                                  (spec_sample_df['<S/N>'] > 2.5) &
                                                  (spec_sample_df['Rpet_ell'] > 2 * PSF_FWHM) &
                                                  (spec_sample_df['ITOT'] < 24.5)&
                                        (spec_sample_df['MSTAR'] >= 10.4)).dropna()
    
    spec_sample_field_data = spec_sample_field_base.iloc[pandas_pick_spec_field(spec_sample_field_base['SPECMEMB'],
                                                                                spec_sample_field_base['ZCLUST'],
                                                                                spec_sample_field_base['ZSPEC'],
                                                                                spec_sample_field_base['qual'])]


    spec_structures = {"cl1037_group": cl1037_group,
                       "cl1040": cl1040,
                       "cl1040a_group": cl1040a_group,
                       "cl1040b_group": cl1040b_group,
                       "cl1054_1146": cl1054_1146,
                       "cl1054_1146a_group": cl1054_1146a_group,
                       "cl1054_1245": cl1054_1245,
                       "cl1054_1245a_group": cl1054_1245a_group,
                       "cl1103a_group": cl1103a_group,
                       "cl1103b_group": cl1103b_group,
                       "cl1138_cluster": cl1138_cluster,
                       "cl1138a_cluster": cl1138a_cluster,
                       "cl1216": cl1216,
                       "cl1227_cluster": cl1227_cluster,
                       "cl1227a_cluster": cl1227a_cluster,
                       "cl1232": cl1232,
                       "cl1354_cluster": cl1354_cluster,
                       "cl1354a_cluster": cl1354a_cluster,
                       "Field": spec_sample_field_data}

    return spec_structures


"""Pick phot+spec sample cluster members"""

phot_sample = pd.read_csv('/home/sinandeger/Research/MasterCatalog/AfterEdit-Photometric-Catalog(Sinan-Edit)+qual-column-FillNA.csv')

mous_phot_at_clustz = pd.read_csv('/home/sinandeger/Research/Moustakas-EDisCS-Catalog/Phot-at-Cluster-z/allclust-moustakas-phot-at-clustz-match.csv')

mous_phot_for_field = pd.read_csv('/home/sinandeger/Research/Moustakas-EDisCS-Catalog/Phot_at_Phot_z/Mous-at-Phot-z-Match-EDisCS.csv')

def phot_sample_select():

    phot_sample_members = phot_sample.iloc[PhotMember_pandas(np.array(phot_sample['SPECMEMB']), np.array(phot_sample['MEMBFLAG']))]

    phot_sample_members_data = phot_sample_members.where((phot_sample_members['WMIN'] > 0.3) &
                                                 (phot_sample_members['FLAG'] == 0) &
                                                 (phot_sample_members['<S/N>'] > 2.5) &
                                                 (phot_sample_members['Rpet_ell'] > 2 * PSF_FWHM) &
                                                 (phot_sample_members['ITOT'] < 24.5) &
                                                 (phot_sample_members['CLUSTNAME'] != 'cl1037-1243') &
                                                 (phot_sample_members['CLUSTNAME'] != 'cl1103-1245')).dropna()

    len(phot_sample_members_data)

    cl1040_phot = phot_sample_members_data.where(phot_sample_members_data['CLUSTNAME'] == 'cl1040-1155').dropna()
    cl1054_1146_phot = phot_sample_members_data.where(phot_sample_members_data['CLUSTNAME'] == 'cl1054-1146').dropna()
    cl1054_1245_phot = phot_sample_members_data.where(phot_sample_members_data['CLUSTNAME'] == 'cl1054-1245').dropna()
    cl1216_phot = phot_sample_members_data.where(phot_sample_members_data['CLUSTNAME'] == 'cl1216-1201').dropna()
    cl1227_phot = phot_sample_members_data.where(phot_sample_members_data['CLUSTNAME'] == 'cl1227-1138').dropna()
    cl1232_phot = phot_sample_members_data.where(phot_sample_members_data['CLUSTNAME'] == 'cl1232-1250').dropna()
    cl1354_phot = phot_sample_members_data.where(phot_sample_members_data['CLUSTNAME'] == 'cl1354-1230').dropna()

    cl1354a_spec_subtract_df = cl1354_phot.where((cl1354_phot['ZSPEC'] < 0.6) &
                                                 (cl1354_phot['ZSPEC'] > 0.5)).dropna()

    cl1354_phot = cl1354_phot[~cl1354_phot.isin(cl1354a_spec_subtract_df)].dropna()

    phot_structures = {"cl1040_phot": cl1040_phot,
                       "cl1054_1146_phot": cl1054_1146_phot,
                       "cl1054_1245_phot": cl1054_1245_phot,
                       "cl1216_phot": cl1216_phot,
                       "cl1227_phot": cl1227_phot,
                       "cl1232_phot": cl1232_phot,
                       "cl1354_phot": cl1354_phot}

    return phot_structures