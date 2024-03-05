import geopandas as gpd
gpd.options.use_pygeos = False
import pandas as pd
import numpy as np
import maup
import argparse
import os
from gerrychain import Graph
from gerrychain.tree import recursive_tree_part
import sklearn.cluster as skc
import libpysal as lp
from tqdm import tqdm

import shapely
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)

import warnings; warnings.filterwarnings('ignore', 'GeoSeries.isna', UserWarning)

def make_assignment(args, shp, num_dists):
    if args.districting_method == 'agglom':
        rook_graph = lp.weights.Rook.from_dataframe(shp)
        clusters = skc.AgglomerativeClustering(n_clusters=num_dists, connectivity=rook_graph.sparse, linkage='ward').fit(np.array(shp.TOTPOPdp).reshape(-1, 1))
        return shp.assign(DISTRICT=clusters.labels_ + 1)
    elif args.districting_method == 'gerrychain':
        graph = Graph.from_geodataframe(shp)
        totpop = shp['TOTPOPdp'].sum()
        print(f'Creating {num_dists} districts')
        assignments = recursive_tree_part(graph, range(1, num_dists + 1), totpop/num_dists, 'TOTPOPdp', args.pop_dev)
        x = list(assignments.keys())
        x.sort()
        shp['DISTRICT'] = [assignments[i] for i in x]
        return shp
    elif args.districting_method == 'districts':
        if args.level == 'block':
            districts = pd.read_csv(args.bef, dtype={'BLOCKID': object})
            districts = districts.rename(columns={'BLOCKID': 'GEOID10'})

            sorted_dict = {x:(ind+1) for ind,x in enumerate(sorted(districts.DISTRICT.unique()))}
            districts = districts[['GEOID10', 'DISTRICT']]
            districts.DISTRICT = districts.DISTRICT.map(sorted_dict)

            shp = shp.merge(districts, on='GEOID10')

            return shp
        else:
            districts = gpd.read_file('zip://' + args.dist_path)
            districts.to_crs(shp.crs, inplace=True)
            assignment = maup.assign(shp, districts)

            assignment = assignment.rank(method='dense').astype(int)
            return shp.assign(DISTRICT=assignment)
    else:
        print(f'Districting method {args.districting_method} not implemented')
        exit(1)
    return None

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('shp_path', type=str, help='Path to TIGER shapefile')
    parser.add_argument('dist_path', type=str, help='Path to district shapefile')
    parser.add_argument('dp_path', type=str, help='Path to DP csv')
    parser.add_argument('save_path', type=str, help='Path to save output (.shp)')
    parser.add_argument('level', type=str, choices=['block','bg','tract'], help='Geographic level')
    parser.add_argument('--boundary_name', type=str, default=None, help='Name of boundary of place to clip data to (must also provide path to places')
    parser.add_argument('--place', type=str, default=None, help='Path to TIGER shapefile of Census-designated places')
    parser.add_argument('--force', action='store_true', default=False, help='Run without checking if save_path exists')
    parser.add_argument('--lang', choices=['julia', 'python'], default='julia', help='Language used for ensemble generation')
    parser.add_argument('--detailed', action='store_true', default=False, help='Include larger set of racial groups and housing data')
    parser.add_argument('--districting_method', choices=['gerrychain', 'agglom', 'districts'], default='gerrychain', help='Which method to use for districting')
    parser.add_argument('--pop_dev', default=0.1, type=float, help='Maxs deviation from the mean district population')
    parser.add_argument('--num_districts', type=int, default=None, help='Number of districts to create')
    parser.add_argument('--county', type=str, default=None, help='Only export data for this county')
    parser.add_argument('--all_counties', action='store_true', default=False, help='Create separate shapefiles for each county')
    parser.add_argument('--pop_target', type=int, default=None, help='Target population for each district')
    parser.add_argument('--no_districts', action='store_true', default=False, help='Do not add district column to output')
    parser.add_argument('--bef', type=str, default=None, help='Path to block equivalency file')
    parser.add_argument('--save_graph', action='store_true', default=False, help='Save graph to file')
    # parser.add_argument('--hh_index', action='store_true', default=False, help='Include extra demographic information for computing hh index')
    args = parser.parse_args()

    if os.path.exists(args.save_path) and not args.force:
        print(f'{args.save_path} already exists, continuing...')
        exit()

    # Load TIGER shapefile, add GISJOIN attribute
    print('Loading TIGER data')
    shp = gpd.read_file(args.shp_path)
    
    if args.county:
        shp = shp[shp['COUNTYFP10']==args.county]

    if args.level == 'block':
        shp['gisjoin'] = 'G' + shp['STATEFP10'] + '0' + shp['COUNTYFP10'] + '0' + shp['TRACTCE10'] + shp['BLOCKCE10']
    if args.level == 'bg':
        shp['gisjoin'] = 'G' + shp['STATEFP10'] + '0' + shp['COUNTYFP10'] + '0' + shp['TRACTCE10'] + shp['BLKGRPCE10']
    elif args.level == 'tract':
        shp['gisjoin'] = 'G' + shp['STATEFP10'] + '0' + shp['COUNTYFP10'] + '0' + shp['TRACTCE10']
    shp = shp.to_crs('epsg:26918')
    # Filter out districts that are all water
    # shp = shp[shp['ALAND10'] > 0]

    # If --place, load places shapefile and assign places to TIGER
    if args.place:
        print('Loading Census-designated places')
        place = gpd.read_file(args.place)
        place = place.to_crs('epsg:26918')
        # TODO: Temporary, only chooses large cities and towns
        if args.boundary_name:
            print(f'Clipping data to {args.boundary_name}')
            boundary = place[place['NAME10'] == args.boundary_name]
            shp = shp.clip(boundary)
        else:
            place = place[place['ALAND10'] > 20000000]
            place_assignment = maup.assign(shp, place)
            place_assignment[place_assignment.isna()] = max(place_assignment) + 1
            shp['PLACE'] = place_assignment

    # Load DP csv file
    print('Loading DP data')
    dp = pd.read_csv(args.dp_path)

    dp_col_names = {'H73001_dp':'TOTPOPdp', 'H73005_dp':'NH_WHITEdp', 'H73006_dp':'NH_BLACKdp', 'H73007_dp':'NH_AMINdp', 'H73008_dp':'NH_ASIANdp', 'H73009_dp':'NH_NHPIdp', 'H73010_dp':'NH_OTHERdp', 'H73011_dp':'NH_2MOREdp', 'H73002_dp':'HISPdp', 'H75001_dp':'VAPdp', 'H75002_dp':'HVAPdp', 'H75005_dp':'WVAPdp', 'H75006_dp':'BVAPdp', 'H75007_dp':'AMINVAPdp', 'H75008_dp':'ASIANVAPdp', 'H75009_dp':'NHPIVAPdp', 'H75010_dp':'OTHERVAPdp', 'H75011_dp':'2MOREVAPdp'}#, 'IFE001_dp':'UNITSdp', 'IFE002_dp':'OCCdp', 'IFE003_dp':'VACANTdp'}
    sf_col_names = {'H73001_sf':'TOTPOPsf', 'H73005_sf':'NH_WHITEsf', 'H73006_sf':'NH_BLACKsf', 'H73007_sf':'NH_AMINsf', 'H73008_sf':'NH_ASIANsf', 'H73009_sf':'NH_NHPIsf', 'H73010_sf':'NH_OTHERsf', 'H73011_sf':'NH_2MOREsf', 'H73002_sf':'HISPsf', 'H75001_sf':'VAPsf', 'H75002_sf':'HVAPsf', 'H75005_sf':'WVAPsf', 'H75006_sf':'BVAPsf', 'H75007_sf':'AMINVAPsf', 'H75008_sf':'ASIANVAPsf', 'H75009_sf':'NHPIVAPsf', 'H75010_sf':'OTHERVAPsf', 'H75011_sf':'2MOREVAPsf'}#, 'IFE001_sf':'UNITSsf', 'IFE002_sf':'OCCsf', 'IFE003_sf':'VACANTsf'}

    dp = dp[['gisjoin'] + list(dp_col_names.keys()) + list(sf_col_names.keys())]
    dp = dp.rename(columns=dp_col_names).rename(columns=sf_col_names)

    # Merge shapefile with dp data
    print('Merging tables')
    shp = shp.merge(dp, on='gisjoin')

    if args.place:
        indices = ['COUNTYFP10', 'PLACE', 'geometry', 'ALAND10', 'GEOID10']
    else:
        indices = ['COUNTYFP10', 'geometry', 'ALAND10', 'GEOID10']
    shp2 = shp.loc[:,indices]

    # Add complement values
    print('Calculating complement values')
    dts = ['sf', 'dp']
    for dt in dts:
        pop = 'TOTPOP' + dt
        demo_pop = ['NH_WHITE', 'NH_BLACK', 'HISP']
        if args.detailed:
            demo_pop += ['NH_AMIN', 'NH_ASIAN', 'NH_NHPI']
        demo_pop = [d + dt for d in demo_pop]
        vap = 'VAP' + dt
        demo_vap = ['WVAP', 'BVAP', 'HVAP']
        if args.detailed:
            demo_vap += ['AMINVAP', 'ASIANVAP', 'NHPIVAP']
        demo_vap = [d + dt for d in demo_vap]

        shp2[[pop, vap] + demo_pop + demo_vap] = shp[[pop, vap] + demo_pop + demo_vap]

        new_pop = ['NWHITE', 'NBLACK', 'NHISP']
        if args.detailed:
            new_pop += ['NAMIN', 'NASIAN', 'NNHPI']
        new_pop = [d + dt for d in new_pop]
        new_vap = ['NWVAP', 'NBVAP', 'NHVAP']
        if args.detailed:
            new_vap += ['NAMINVAP', 'NASIAVAP', 'NNHPIVAP']
        new_vap = [d + dt for d in new_vap]

        for i in range(len(demo_pop)):
            shp2.loc[:,new_pop[i]] = shp2[pop] - shp2[demo_pop[i]]
        for i in range(len(demo_vap)):
            shp2.loc[:,new_vap[i]] = shp2[vap] - shp2[demo_vap[i]]

    # Write final shapefile
    if args.all_counties:
        print('Writing shapefiles')
        counties = shp2['COUNTYFP10'].unique()
        for county in tqdm(counties):
            county_fname = args.save_path[:-4] + '_' + str(county) + '.shp'
            county_shp = shp2[shp2['COUNTYFP10'] == county]

            if args.num_districts:
                num_dists = args.num_districts
            elif args.pop_target:
                tot_pop = county_shp['TOTPOPdp'].sum()
                num_dists = int(tot_pop / args.pop_target)
            else:
                districts = gpd.read_file('zip://' + args.dist_path)
                num_dists = len(districts)

            county_shp = make_assignment(args, county_shp, num_dists)

            county_shp.to_file(county_fname)
            if args.save_graph:
                county_graph = Graph.from_geodataframe(county_shp)
                county_graph.to_json(county_fname[:-4] + '.json')
    else:
        print('Writing shapefile')
        print(shp2.columns)

        if args.num_districts:
            num_dists = args.num_districts
        elif args.pop_target:
            tot_pop = shp2['TOTPOPdp'].sum()
            num_dists = int(tot_pop / args.pop_target)
        else:
            districts = gpd.read_file('zip://' + args.dist_path)
            num_dists = len(districts)

        shp2 = make_assignment(args, shp2, num_dists)

        shp2.to_file(args.save_path)
        if args.save_graph:
            print('Writing graph')
            graph = Graph.from_geodataframe(shp2)
            graph.to_json(args.save_path[:-4] + '_graph.json')
