from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd
import pyvo
import matplotlib.pyplot as plt
from astropy.io import fits

pd.set_option('display.max_columns', None)
service = pyvo.dal.TAPService("https://almascience.eso.org/tap")


def query_coordinate_list(service, ra_list, dec_list, radius):
    """Runs a single query for an entire list of coordinates. Query only 100 or less coordinates in one go.

       service        pyvo TAPService instance
       ra_list        List of Right Ascension values in decimal degrees
       dec            List of Declination values in decimal degrees
       radius         Search radius in decimal degrees

       returns        pandas table
    """

    conditionlist = [f"INTERSECTS(CIRCLE('ICRS',{c[0]},{c[1]},{radius}),s_region)=1" for c in zip(ra_list, dec_list)]

    query = f"""
            SELECT *
            FROM ivoa.obscore
            WHERE {" OR ".join(conditionlist)}"""

    return service.search(query).to_table().to_pandas()


def query_alma_output(ra_list, dec_list, radius):
    result = pd.concat([query_coordinate_list(service, ra_list[i:i + 100], dec_list[i:i + 100], radius) for i in
                        range(0, len(ra_list), 100)], ignore_index=True)
    print(f"Number of ALMA rows returned for the {len(ra_list)} catalog entries: {len(result)}")
    final_result = result.head(len(result))
    return final_result


def match_to_catalog(ra_list, dec_list, ra_cat_list, dec_cat_list, redshift, identifier):
    """

    :param ra_list: List of Right Ascension values in decimal degrees (Input)
    :param dec_list:  List of Declination values in decimal degrees (Input)
    :param ra_cat_list: List of Right Ascension values in decimal degrees from the ALMA catalog
    :param dec_cat_list:  List of Declination values in decimal degrees from the ALMA catalog
    :param redshift: Redshift value from the input catalog
    :return: cross match for nthclosest neighbour
    """
    sdss_cat = SkyCoord(ra=ra_list * u.degree, dec=dec_list * u.degree)  # Input catalog in astropy units.
    alma_cat = SkyCoord(ra=ra_cat_list * u.degree, dec=dec_cat_list * u.degree)
    # The base catalog in which to search for matches.
    RA_S = []
    DEC_S = []
    REDSHIFT = []
    IDENTIFIER = []
    D2D = []
    CLOSEST = []
    RA_A = []
    DEC_A = []
    ACCESS_URL = []
    PROPOSAL_ID = []
    TARGET_NAME = []
    FOV = []
    RESOLUTION = []
    T_EXPTIME = []
    BAND_LIST = []
    SENSITIVITY_10KMS = []
    CONT_SENSITIVITY_BANDWIDTH = []
    BANDWIDTH = []
    SPATIAL_RESOLUTION = []
    FREQUENCY_SUPPORT = []
    FREQUENCY = []
    VELOCITY_RESOLUTION = []
    SCIENCE_KEYWORD = []

    for closest in range(1, 180):
        idx, d2d, d3d = sdss_cat.match_to_catalog_sky(alma_cat, nthneighbor=closest)
        # matches of this coordinate in a set of input catalog coordinates.
        max_sep = 10.0 * u.arcsec  # create a mask
        sep_constraint = d2d < max_sep
        sdss_matches = sdss_cat[sep_constraint]
        alma_matches = alma_cat[idx[sep_constraint]]
        D2D = D2D + d2d.arcsec[sep_constraint].tolist()
        # angular distance
        RA_S = RA_S + sdss_matches.ra.deg.tolist()
        DEC_S = DEC_S + sdss_matches.dec.deg.tolist()
        REDSHIFT = REDSHIFT + redshift[sep_constraint].tolist()
        IDENTIFIER = IDENTIFIER + identifier[sep_constraint].tolist()
        # SDSS Input
        CLOSEST = CLOSEST + [closest for kk in range(len(sdss_matches))]
        # nth closest match
        RA_A = RA_A + alma_matches.ra.deg.tolist()
        DEC_A = DEC_A + alma_matches.dec.deg.tolist()
        ACCESS_URL = ACCESS_URL + catalog_query_output['access_url'][idx[sep_constraint]].values.tolist()
        PROPOSAL_ID = PROPOSAL_ID + catalog_query_output['proposal_id'][idx[sep_constraint]].values.tolist()
        TARGET_NAME = TARGET_NAME + catalog_query_output['target_name'][idx[sep_constraint]].values.tolist()
        FOV = FOV + catalog_query_output['s_fov'][idx[sep_constraint]].values.tolist()
        RESOLUTION = RESOLUTION + catalog_query_output['s_resolution'][idx[sep_constraint]].values.tolist()
        T_EXPTIME = T_EXPTIME + catalog_query_output['t_exptime'][idx[sep_constraint]].values.tolist()
        BAND_LIST = BAND_LIST + catalog_query_output['band_list'][idx[sep_constraint]].values.tolist()
        SENSITIVITY_10KMS = SENSITIVITY_10KMS + catalog_query_output['sensitivity_10kms'][
            idx[sep_constraint]].values.tolist()
        CONT_SENSITIVITY_BANDWIDTH = CONT_SENSITIVITY_BANDWIDTH + catalog_query_output['cont_sensitivity_bandwidth'][
            idx[sep_constraint]].values.tolist()
        BANDWIDTH = BANDWIDTH + catalog_query_output['bandwidth'][idx[sep_constraint]].values.tolist()
        SPATIAL_RESOLUTION = SPATIAL_RESOLUTION + catalog_query_output['spatial_resolution'][
            idx[sep_constraint]].values.tolist()
        FREQUENCY_SUPPORT = FREQUENCY_SUPPORT + catalog_query_output['frequency_support'][
            idx[sep_constraint]].values.tolist()
        FREQUENCY = FREQUENCY + catalog_query_output['frequency'][idx[sep_constraint]].values.tolist()
        VELOCITY_RESOLUTION = VELOCITY_RESOLUTION + catalog_query_output['velocity_resolution'][
            idx[sep_constraint]].values.tolist()
        SCIENCE_KEYWORD = SCIENCE_KEYWORD + catalog_query_output['science_keyword'][idx[sep_constraint]].values.tolist()
        # ALMA Input

    catalog_alma_match_catalog = pd.DataFrame(
            {'Closest': CLOSEST, 'ra_sdss': RA_S, 'dec_sdss': DEC_S, 'ra_alma': RA_A, 'dec_alma': DEC_A,
             'angular_distance': D2D,
             'redshift': REDSHIFT, 'SDSS_identifier': IDENTIFIER, 'access_url': ACCESS_URL, 'proposalid': PROPOSAL_ID, 'target_name': TARGET_NAME,
             's_fov': FOV, 's_resolution': RESOLUTION, 't_exptime': T_EXPTIME,
             'band_list': BAND_LIST, 'sensitivity_10kms': SENSITIVITY_10KMS,
             'cont_sensitivity_bandwidth': CONT_SENSITIVITY_BANDWIDTH, 'bandwidth': BANDWIDTH,
             'spatial_resolution': SPATIAL_RESOLUTION, 'frequency': FREQUENCY, 'velocity_resolution': VELOCITY_RESOLUTION,
	     'science_key': SCIENCE_KEYWORD, 'frequency_support': FREQUENCY_SUPPORT})

    return catalog_alma_match_catalog


def search_around_coordinate(ra_list, dec_list, ra_cat_list, dec_cat_list, redshift, identifier):
    """

    :param ra_list:  List of Right Ascension values in decimal degrees  (Input)
    :param dec_list:  List of Declination values in decimal degrees (Input)
    :param ra_cat_list: List of Right Ascension values in decimal degrees from the ALMA catalog
    :param dec_cat_list: List of Declination values in decimal degrees from the ALMA catalog
    :param redshift: Redshift value from the input catalog
    :return: cross match within 10arcsec radius for all objects
    """
    sdss_cat = SkyCoord(ra=ra_list * u.degree, dec=dec_list * u.degree)  # Input catalog in astropy units.
    alma_cat = SkyCoord(ra=ra_cat_list * u.degree, dec=dec_cat_list * u.degree)
    #  The base catalog in which to  search for matches.
    idxcatalog, idxc, d2d, d3d = sdss_cat.search_around_sky(alma_cat, 10 * u.arcsec)  # Searches for all coordinates
    # in this object around a supplied set of points within 10arcsec on-sky separation.

    ra_sdss = sdss_cat[idxc].ra.deg
    dec_sdss = sdss_cat[idxc].dec.deg
    redshift_sdss = redshift[idxc]
    sdss_identifier = identifier[idxc]
    # SDSS Input
    ra_alma = alma_cat[idxcatalog].ra.deg
    dec_alma = alma_cat[idxcatalog].dec.deg
    angluar_distance = d2d.arcsec
    access_url = catalog_query_output['access_url'][idxcatalog].values
    proposal_id = catalog_query_output['proposal_id'][idxcatalog].values
    target_name = catalog_query_output['target_name'][idxcatalog].values
    s_fov = catalog_query_output['s_fov'][idxcatalog].values
    s_resolution = catalog_query_output['s_resolution'][idxcatalog].values
    t_exptime = catalog_query_output['t_exptime'][idxcatalog].values
    band_list = catalog_query_output['band_list'][idxcatalog].values
    sensitivity_10kms = catalog_query_output['sensitivity_10kms'][idxcatalog].values
    cont_sensitivity_bandwidth = catalog_query_output['cont_sensitivity_bandwidth'][idxcatalog].values
    bandwidth = catalog_query_output['bandwidth'][idxcatalog].values
    spatial_resolution = catalog_query_output['spatial_resolution'][idxcatalog].values
    frequency_support = catalog_query_output['frequency_support'][idxcatalog].values
    frequency = catalog_query_output['frequency'][idxcatalog].values
    velocity_resolution = catalog_query_output['velocity_resolution'][idxcatalog].values
    science_keyword = catalog_query_output['science_keyword'][idxcatalog].values

    Skycoord_search_around = pd.DataFrame(
        {'Closest': 0, 'ra_sdss': ra_sdss, 'dec_sdss': dec_sdss, 'ra_alma': ra_alma, 'dec_alma': dec_alma,
         'angular_distance': angluar_distance, 'redshift': redshift_sdss, 'SDSS_identifier': sdss_identifier, 'access_url': access_url,
         'proposalid': proposal_id, 'target_name': target_name,
         's_fov': s_fov, 's_resolution': s_resolution, 't_exptime': t_exptime,
         'band_list': band_list, 'sensitivity_10kms': sensitivity_10kms,
         'cont_sensitivity_bandwidth': cont_sensitivity_bandwidth, 'bandwidth': bandwidth,
         'spatial_resolution': spatial_resolution, 'frequency': frequency, 'velocity_resolution': velocity_resolution, 
	 'science_key': science_keyword, 'frequency_support': frequency_support})

    return Skycoord_search_around


def final_match(match_catalog, search_around_catalog):
    """

    :param match_catalog:  nearest on-sky matches of the coordinate in a set of catalog coordinates.
    :param search_around_catalog: all coordinates cross matched  around a supplied set of points within a 10arcsec on-sky separation.
    :return: merged catalog and duplicate removed catalogs
    """
    match_catalog = match_catalog.append(search_around_catalog, ignore_index=True)

    catalog_duplicate_removed_keeping_frequency = match_catalog.drop_duplicates(
        subset=['ra_alma', 'dec_alma', 'frequency'], ignore_index=True)
    catalog_duplicate_removed = match_catalog.drop_duplicates(subset=['ra_alma', 'dec_alma'],
                                                                   ignore_index=True)

    return match_catalog, catalog_duplicate_removed_keeping_frequency, catalog_duplicate_removed


def alma_plots(ra_list, dec_list, ra_cat_list, dec_cat_list):
    """

    :param ra_list:  List of Right Ascension values in decimal degrees  (Input)
    :param dec_list:  List of Declination values in decimal degrees (Input)
    :param ra_cat_list: List of Right Ascension values in decimal degrees from the ALMA catalog
    :param dec_cat_list: List of Declination values in decimal degrees from the ALMA catalog
    :return:
    """

    plt.rcParams["figure.figsize"] = (15, 10)

    plt.plot(ra_list, dec_list, ls='', marker='.', ms=1, fillstyle='none', label='all')
    plt.plot(ra_cat_list, dec_cat_list, ls='', marker='o', ms=12, fillstyle='none',
             label='observed by ALMA')
    plt.legend(loc='best')

    return


if __name__ == '__main__':
    DR16_file = fits.open('/data/beegfs/astro-storage/groups/walter/hegde/alma_query/DR16Q_v4.fits')
    ra_list = DR16_file[1].data['RA']
    dec_list = DR16_file[1].data['DEC']
    SDSS_identifier = DR16_file[1].data['THING_ID']
    redshift = DR16_file[1].data['Z']

    catalog_query_output = query_alma_output(ra_list, dec_list, 10 / 3600)

    ra_cat = catalog_query_output['s_ra'].values
    dec_cat = catalog_query_output['s_dec'].values

    catalog_match = match_to_catalog(ra_list, dec_list, ra_cat, dec_cat, redshift, SDSS_identifier)
    coordinate_search = search_around_coordinate(ra_list, dec_list, ra_cat, dec_cat, redshift, SDSS_identifier)
    merged_catalog, catalog_with_frequency, catalog_without_frequency = final_match(coordinate_search, catalog_match)
   
    catalog_query_output.to_csv('ALMA_query_catalog.csv', sep='\t')
    catalog_match.to_csv('match_to_catalog_nthclosest_neighbor.csv', sep='\t')
    coordinate_search.to_csv("search_around_catalog_all.csv", sep='\t')
    merged_catalog.to_csv("merged_catalog_for_all_coordinates.csv", sep='\t')
    catalog_with_frequency.to_csv("duplicate_frequency_ra_dec_removed.csv", sep='\t')
    catalog_without_frequency.to_csv('duplicate_ra_dec_removed.csv', sep='\t')

    print(f" Query in ALMA for input catalog : {len(catalog_query_output)}")
    print(f" Cross matched for n closest neighbor: {len(catalog_match)}")
    print(
        f" Cross matched  around a supplied set of points within a given on-sky separation  : {len(coordinate_search)}")
    print(f" Cross matched catalogs after merging : {len(merged_catalog)}")
    print(f" Duplicate ra, dec and frequency removed :  {len(catalog_with_frequency)}")
    print(f" Duplicate ra and dec  removed :  {len(catalog_without_frequency)}")
  
    plots = alma_plots(ra_list, dec_list, ra_cat, dec_cat)
    plt.savefig("SDSS_ALMA_ra_dec.pdf")
