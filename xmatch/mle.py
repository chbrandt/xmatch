# -*- coding:utf-8 -*-

'''
Maximum Likelihood Estimator

The MLE method, learnt from [LaMassa et al. 2016], takes into account
the distance between the primary source and the ancillary objects,
their astrometic errors and magnitude distribution. Ancillary objects
will be the sources found within the search radius in the second catalog.

In MLE it is defined the Reliability value, `R`, which is a way to
distinguish between true counterparts and chance detections. The greater
the `R`, greater are the chances to be a true counterpart. `R` is computed
as follows:
$$
    R = \frac{LR}{\Sum_i LR_i + (1 - Q)}
$$

`Q` is the ratio of primary sources that present ancillary objects within
the search radius [DOUBT-1] divided by the total number of (primary) sources.

`LR` is the likelihood ratio, it is assigned to each ancillary object
within the search radius and it is the probability that one of those
objects is the true counterpart divided by the probability that a background
source is there by chance.
$$
    LR = \frac{q(m) f(r)}{n(m)}
$$

$q(m)$ is the expected normalized magnitude distribution of counterparts
within the radius search. $q(m)$ is estimated by subtracting the histogram
of the sources inside the the radius from the histogram of background
sources; the histograms are, each, normalized by the respective search
areas. Background objects are taken from an annulus around the primary
source, where the internal radius should be larger than the search radius.

$f(r)$ is the probability distribution of the astrometric errors; it is
modeled as a two-dimensional Gaussian distribution where the primary
ancillary positional errors are added in quadrature.

$n(m)$ is the normalized magnitude distribution of background sources.


LaMassa et al. 2016: http://doi.org/10.3847/0004-637X/817/2/172

DOUBT-1: are we talking about the SAME searched radius --i.e, intersecting
regions-- or in general, throughout the catalog/sky?!

'''
import logging
import astropy
import pandas
import sklearn

import warnings
warnings.filterwarnings(action='ignore')

import booq
import math

def nearest_neighbors(cat,N=2):
    '''
    Compute the distances from each entry to its 'N' nearest neighbours
    '''
    from sklearn.neighbors import NearestNeighbors

    X = cat[['RA','DEC']].values
    nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(X)
    distances, indices = nbrs.kneighbors(X)
    distances = distances.transpose()
    return distances,indices

def mle(catalog_A, catalog_B, columns_A=None, columns_B=None,
        feature_column=None,
        ancillary_radius=None,
        background_radii=None,
        separation_unit='arcsec',
        inner_join=False):
    '''

    Input:
     - catalog_[AB]     : ~pandas.DataFrame
                          DFs providing at least columns 'RA','DEC','ID'
                          "ra" and "dec" are assumed to be in degres if not Astropy' Skycoords
                          (see columns_[AB] for mapping column names)
     - columns_[AB]     : dict mapping 'ra','dec','id' and 'pos_err' columns
                          In case catalog(s) have different column names for 'RA','DEC','ID'
                          e.g, {'ra':'RA', 'dec':'Dec', 'id':'ObjID'}
     - feature_column   : string
                          Column name from 'catalog_B' to be used as MLE estimator
     - ancillary_radius : float or Astropy' Angle
                          Radius to search for ancillary object in "B",
                          if a float is given, 'degree' is the default unit
     - background_radii : float or [float,float] or Astropy' Angle(s)
                          If only a value is given it is assumed to be the outer radius
                          and the inner radius would that of 'ancillary_radius'
     - separation_unit  : string ['arcsec','arcmin','degree']
                          unit to express the distance between matches
     - inner_join       : boolean
                          True for output table with only matching entries
                          False will output the equivalente of a 'left' join

    Output:
     - matched_catalog : ~pandas.DataFrame
                         Column 'R', 'LR', ...
    '''
    from astropy.coordinates import Angle
    try:
        rs = ancillary_radius.to(separation_unit)
    except:
        rs = Angle(ancillary_radius,separation_unit)

    if not background_radii:
        ri,ro = define_background_annulus(search_radius)
    else:
        try:
            ri,ro = background_radii
            if ri>ro:
                _r = ro
                ro = ri
                ri = _r
        except:
            ro = background_radii
            ri = rs

    if not all([isinstance(ri,Angle),isinstance(ro,Angle)]):
        ro = Angle(ro,'degree')
        ri = Angle(ri,'degree')

    # Alias
    feature = feature_column

    from booq.pipelines.xmatch import xmatch
    # Notice that the following 'ancillary' and 'background'
    # matched catalogs have the same index.
    # 'xmatch' returns the index of 'catalog_A'.

    # Ancillary match
    match_ancillary = xmatch(catalog_A, catalog_B,
                             columns_A=columns_A,
                             columns_B=columns_B,
                             method='gc',radius=rs)
    # remove non-matching entries
    non_matches = match_ancillary.loc[:,('AB','separation')].isnull()
    match_ancillary = match_ancillary.loc[~non_matches]

    ancillary_sample = flatten_matches(match_ancillary,
                                        ('A',columns_A['id']),
                                        ('B',columns_B['id']))

    ancillary_sample.index.names = [columns_A['id'],columns_B['id']]
    ancillary_sample.reset_index(inplace=True)

    N_anc = len(ancillary_sample[columns_A['id']].unique())
    Area_radius = math.pi * rs**2
    Area_radius = Area_radius.to('arcmin2')
    Area_total_ancillary = Area_radius * N_anc
    del N_anc

    ancillary_sample = ancillary_sample.merge(catalog_A[[columns_A['id'],
                                                         columns_A['pos_err']]],
                                              on=columns_A['id'])

    # Background match
    match_background = xmatch(catalog_A, catalog_B,
                              columns_A=columns_A,
                              columns_B=columns_B,
                              method='gc',radius=ro)
    # remove non-matching entries
    # This sample is composed by all matched objects distant more then 'ri'
    match_background = match_background.loc[match_ancillary.index]

    background_sample = flatten_matches(match_background,
                                        ('A',columns_A['id']),
                                        ('B',columns_B['id']))

    background_sample.index.names = [columns_A['id'],columns_B['id']]
    background_sample.reset_index(inplace=True)

    # Keep only objects beyond "inner-radius"
    background_sample = background_sample.loc[background_sample['distances'] > ri]

    N_bkg = len(background_sample[columns_A['id']].unique())
    Area_annulus = math.pi * (ro**2 - ri**2)
    Area_annulus = Area_annulus.to('arcmin2')
    Area_total_background = Area_annulus * N_bkg
    del N_bkg

    background_sample = background_sample.merge(catalog_A[[columns_A['id'],columns_A['pos_err']]],on=columns_A['id'])

    e_pos = ancillary_sample['e_Pos'].mean()
    func_fr = define_fr(e_pos,e_pos)


    ancillary_sample = ancillary_sample.merge(catalog_B[[columns_B['id'],feature]],on=columns_B['id'])
    background_sample = background_sample.merge(catalog_B[[columns_B['id'],feature]],on=columns_B['id'])

    func_nm = surface_density(background_sample[feature],
                              Area_total_background, False)

    func_qm = define_qm(ancillary_sample[feature],
                        background_sample[feature],
                        Area_total_ancillary,
                        Area_total_background)

    func_LR = define_LR( func_qm=func_qm.evaluate,
                    func_fr=func_fr,
                    func_nm=func_nm.evaluate)


    # Alias
    df = ancillary_sample

    lr_col = 'LR_{}'.format(feature)
    df[lr_col] = df.apply(lambda x:func_LR(x['distances'],x[feature]), axis=1)

    Q = 1
    def reliability(lr,Q=1):
        return lr/(sum(lr)+(1-Q))

    r_col = 'R_{}'.format(feature)
    df[r_col] = df.groupby(columns_A['id'])[lr_col].apply(reliability)

    from pandas import DataFrame
    def collapse_table(group):
        out = DataFrame(columns=['LR','duplicates','LRs'])
        out.index.name = columns_B['id']
        if len(group)==1:
            out['Reliability'] = group[r_col]
            out['LR'] = group[lr_col]
            out['duplicates'] = None
            out['duplicates_LR'] = None
            out['duplicates_R'] = None
            out.index = group[columns_B['id']]
            return out
        ind_max = group[r_col] == group[r_col].max()
        out['Reliability'] = group.loc[ind_max,r_col]
        out['LR'] = group.loc[ind_max,lr_col]
        out['duplicates'] = ';'.join([ str(r) for r in group.loc[~ind_max,columns_B['id']] ])
        out['duplicates_LR'] = ';'.join([ str(r) for r in group.loc[~ind_max,lr_col] ])
        out['duplicates_R'] = ';'.join([ str(r) for r in group.loc[~ind_max,r_col] ])
        out.index = group.loc[ind_max,columns_B['id']]
        return out
    df_r = df.groupby(columns_A['id']).apply(collapse_table)

    matching_columns = ['Reliability','LR','duplicates','duplicates_LR','duplicates_R']

    # output_minimal=True
    if inner_join is True:
        A_matched = catalog_A.set_index(columns_A['id']).loc[df_r.index.droplevel(1)].reset_index()
        # print(len(A_matched))
        # A_matched.head()
        B_matched = catalog_B.set_index(columns_B['id']).loc[df_r.index.droplevel(0)].reset_index()
        # print(len(B_matched))
        # B_matched.head()
        AB_matched = df_r[matching_columns].reset_index(drop=True)
        # print(len(AB_matched))
        # AB_matched.head()

    else:
        # A_matched = catalog_A[['Seq','RAdeg','DEdeg','e_Pos']].reset_index(drop=True)
        A_matched = catalog_A[list(columns_A.values())].reset_index(drop=True)
        A_matched.index = A_matched[columns_A['id']]

        B_matched = catalog_B.set_index(columns_B['id']).loc[df_r.index.droplevel(0)].reset_index()
        B_matched.index = catalog_A.set_index(columns_A['id']).loc[df_r.index.droplevel(1)].reset_index()[columns_A['id']]

        AB_matched = df_r[matching_columns].reset_index(level=1,drop=True)

    from pandas import concat
    df = concat([ A_matched, B_matched, AB_matched ], axis=1,
                keys=['A','B','AB_'+feature])

    for col in df.columns:
        ind = df[col].isnull()
        df.loc[ind,col] = None
    return df



# Surface brightness density distribution estimation
#
# Remember, the prior is based on the surface brightness
# distribution, assuming they are different between
# ancillary and background samples.
#

from .fit import surface_density,compile_interpolation_function

def flatten_matches(match_table,column_id_A,column_id_B):
    from pandas import DataFrame
    # from numpy import nan
    def expand_duplicates(group):
        '''
        Expand duplicates/distances ;-separated strings in multiple entries
        '''
        assert len(group)==1
        df_row = group.iloc[0]

        assert ('AB','duplicates') in df_row.index
        assert ('AB','distances') in df_row.index

        dups = [ str(df_row[column_id_B]) ]
        seps = [ float(df_row[('AB','separation')]) ]
        if df_row[('AB','duplicates')] is not None:
            dups.extend( str(df_row[('AB','duplicates')]).split(';') )
            _dists_str = str(df_row[('AB','distances')]).split(';')
            seps.extend([float(sep) for sep in _dists_str])

        df = DataFrame({'distances':seps}, index=dups)
        # Maybe change the index name?
        df.index.name = 'duplicates'

        # Notice that the nearest objects from column ('B','"id"') is not included; conciously
        return df
    # Grouping by the Object-ID --which are unique-- will produce length-1 groups;
    # The goal is to the Object-ID (A/Seq) as the returned DataFrame's index.
    return match_table.groupby(by=[column_id_A]).apply(expand_duplicates)

def build_features_catalog(sample_table,catalog,column_id):
    '''
    Return inner join from sample_table' index and catalog' column_id

    'sample_table' has a Pandas' MultiIndex, level 1 is used
    to match 'catalog' primary-key.
    '''
    # from pandas import merge
    # remove eventual duplicated entries (objects metched by multiple targets)
    bgs = sample_table.copy()
    bgs.index = sample_table.index.droplevel(0).astype(catalog[column_id].dtype)
    bgs = bgs.loc[~bgs.index.duplicated()]

    features_table = bgs.merge(catalog,
                                left_index=True,
                                right_on=column_id,
                                sort=False,
                                how='inner')

    features_table.set_index(column_id,inplace=True)
    return features_table

def define_qm(feature_data_ancillary, feature_data_background,
              surface_area_ancillary, surface_area_background):
    assert feature_data_ancillary.ndim == 1
    assert feature_data_background.ndim == 1

    from numpy import histogram,diff

    veca = feature_data_ancillary
    Area_total_ancillary = surface_area_ancillary
    vecb = feature_data_ancillary
    Area_total_background = surface_area_background

    _min = min(veca.min(),vecb.min())
    _max = min(veca.max(),vecb.max())

    ha,ba = histogram(veca, bins=50, range=(_min,_max))
    ha = ha/ Area_total_ancillary

    hb,bb = histogram(vecb, bins=50, range=(_min,_max))
    hb = hb/ Area_total_background

    x = diff(ba)/2 + ba[:-1]

    h_diff = ha - hb
    # hist_diff[col] = h_diff
    # x_diff[col] = x
    return compile_interpolation_function(x, h_diff/h_diff.sum(),normal=True)

def define_qms(acs_magnitudes,bkg_magnitudes,columns,
                Area_total_ancillary,Area_total_background):
    '''
    Return compiled versions of q(m) for each of 'columns'
    '''
    func_qm = {}
    for col in columns:#['YAPERMAG3','J_1APERMAG3','HAPERMAG3','KAPERMAG3']:
        func_qm[col] = define_qm(acs_magnitudes[col], bkg_magnitudes[col],
                         Area_total_ancillary, Area_total_background)
    return func_qm

# Spatial distribution
#
def define_fr(sigma_x,sigma_y):
    '''
    Return a compiled gaussian f(r) from 'sigma_[xy]'
    '''
    from cmath import sqrt,pi,exp
    sigma = sqrt(sigma_x**2 + sigma_y**2)
    f_2_pi_sigma = 1/(2 * pi * sigma)
    f_2_sigma_2 = -1/(2 * sigma**2)
    def g(r):
        gaussian = f_2_pi_sigma * exp( r**2 * f_2_sigma_2 )
        return gaussian.real
    return g

# Likelihood Reliability
#
def define_LR(func_qm, func_fr, func_nm):
    '''
    Return compiled Likelihood Reliability f(r,m)
    '''
    def LR(r,m, qm=func_qm, fr=func_fr, nm=func_nm):
        return qm(m)*fr(r)/nm(m)
    return LR
