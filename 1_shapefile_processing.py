#probablistic WWR
import shapefile
import os
import pandas as pd
from shapely.geometry import shape, Polygon, LineString, MultiLineString, LinearRing
from shapely.wkt import dumps, loads
from shapely.ops import unary_union
import random as rd
from ast import literal_eval

import geopandas as gpd


ROOT_DIR = os.path.abspath(os.path.dirname(__file__))
SHAPEFILES_DIR = os.path.join(ROOT_DIR, 'shapefiles')
file_name = 'amd_res' #shapefile name
constructions_file = 'constructions_amd.csv'

tolerance = 0.1

def preprocessing():
    # Convert shapefile to the Pandas DataFrame
    df = shp_to_df()
    df = resolve_intersecting_polygons(df)
    df = remove_duplicated_coordinates(df)
    df = polygon_topology_check(
        df, 'AMD_initial_touching', 'AMD_initial_intersect',
        'AMD_initial_perimeter', 'AMD_initial_exposed', 'AMD_initial_partition')
    df = polygon_tolerance(df)
    simplify_polygon_no = 0
    for row in df.itertuples():
        if row.AMD_polygon_simplify:
            simplify_polygon_no += 1
    print('Number of polygons to be simplified: %s' % simplify_polygon_no)
    if simplify_polygon_no > 0:
        df_removed = pd.DataFrame()
        df, _ = polygon_simplification(df, df_removed, simplify_polygon_no)
    df = polygon_topology_check(
        df, 'AMD_simplified_touching', 'AMD_simplified_intersect',
        'AMD_simplified_perimeter', 'AMD_simplified_exposed',
        'AMD_simplified_partition')

    df = collinear_exterior(df)

    df = polygon_topology_check(
        df, 'AMD_collinear_touching', 'AMD_collinear_intersect',
        'AMD_collinear_perimeter', 'AMD_collinear_exposed',
        'AMD_collinear_partition')
    #df = catchment_area(df)
    df = osgb_3D_out(df)
    df = expand_floors(df)
    #df = built_island(df)
    return df


def shp_to_df():

    def polygon_testing(name, p):
        # Internal function which test the polygon
        if not p.is_valid:
            print('%s polygon is not valid' % name)
        if p.exterior.is_ccw:
            print('%s polygon outer ring is not clock-wise' % name)
        if p.interiors:
            for inner_ring in p.interiors:
                if not inner_ring.is_ccw:
                    print('%s polygon inner ring is not counter clock-wise'
                          % name)
        return

    def bi_adj(df):
        gdf = gpd.GeoDataFrame(df, geometry='polygon_original')
        polygon_union = gdf.polygon_original.unary_union

        lst = list()
        for i, item in enumerate(polygon_union):
            lst.append([i+1, item])
        new_df = pd.DataFrame(lst, columns =['id', 'geom'])
        new_gdf = gpd.GeoDataFrame(new_df, geometry='geom')

        for index, row in gdf.iterrows():
            bi = new_gdf[new_gdf.geom.intersects(row['polygon_original'])].id.tolist()
            gdf.at[index, 'bi'] = 'bi_{}'.format(str(bi[0]))
            adj = gdf[gdf.polygon_original.touches(row['polygon_original'])].osgb.tolist()
            gdf.at[index, "adjacent"] = str(adj)
        return gdf

    # Load the shapefile from files provided by Sophie. Get the filename from
    # the configuration file
    sf = shapefile.Reader(os.path.join(SHAPEFILES_DIR, file_name))
    # Load shapes and records stored in the shapefile
    shapes = sf.shapes()
    records = sf.records()
    # Empty list to store shape items which will be converted to the Pandas df
    shapes_to_pd = list()
    # Loop through shapes:
    for i, shp in enumerate(shapes):
        # Check if shape record exists
        if shp.points:
            # Extract from records: building name, height, minimum/maximum Z,
            # perimeter, area
            #bld_name = records[i][1]
            #height = round(float(records[i][6]), 1)
            #z_min = round(float(records[i][7]), 1)
            #z_max = round(float(records[i][8]), 1)
            #perimeter = round(float(records[i][9]), 3)
            #area = round(float(records[i][32]), 3)
            #print(area)
            new_bld_id = records[i][7]
            age = records[i][4]
            cond = records[i][29]
            #print(new_bld_id)
            floors = round(float(records[i][13]),1)
            #print(floors)
            use = records[i][10]
            area = records[i][0]
            #print(use)
            #print (age)
            #print(cond)
            # Create a polygon
            polygon = shape(shp.__geo_interface__)
            #print(polygon)
            # Test the polygon for validity and coordinates direction
            polygon_testing(new_bld_id, polygon)
            # Calculate the area of a polygon
            area_calculated = round(polygon.area, 3)
            #print(area,area_calculated)# row in the df with: building name, polygon, area,
            # unique z elements, and maximum z
            #shapes_to_pd.append([bld_name, height, z_min, z_max, perimeter,
                                 #area, dumps(polygon, rounding_precision=2),
                                 #area_calculated])
            # shapes_to_pd.append([new_bld_id, use, dumps(polygon,
            #                      rounding_precision=2),age, cond, floors, area])
            shapes_to_pd.append([new_bld_id, use, polygon, age, cond, floors, area])
    # DataFrame labels
    labels = ['osgb','use', 'polygon_original', 'age','condition', 'floors','area_from_file']
    # Create a DataFrame and save it to the csv file
    df = pd.DataFrame.from_records(shapes_to_pd, columns=labels)

    df = bi_adj(df)

    df['polygon_original'] = df['polygon_original'].apply(lambda x: dumps(x, rounding_precision=2))
    #df.to_csv(os.path.join(ROOT_DIR, 'amd.csv'), index=False)
    #df.to_csv(os.path.join(ROOT_DIR, 'amd_IK.csv'), index=False)
    #print(df.head())
    # Return the DataFrame
    return df

def resolve_intersecting_polygons(df):
    # Extract building_name-polygon pairs from the DataFrame
    bld_polygon_pair = variable_polygon_pairs(df, 'osgb',
                                                  'polygon_original')
    # List for intersecting polygons
    intersecting_polygons_pairs = list()
    for i, (bld, polygon) in enumerate(bld_polygon_pair[:-1]):
        for adj_bld, adj_polygon in bld_polygon_pair[i+1:]:
            if polygon.intersects(adj_polygon):
                if not polygon.touches(adj_polygon):
                    intersecting_polygons_pairs.append([bld, adj_bld])
                    print('%s polygon intersects %s polygon' % (bld, adj_bld))
    if intersecting_polygons_pairs:
        for pair in intersecting_polygons_pairs:
            polygon = loads(df.loc[df['osgb'] == pair[0], 'polygon_original'].values[0])
            intersect_polygon = loads(df.loc[df['osgb'] == pair[1], 'polygon_original'].values[0])
            new_polygon = polygon.difference(intersect_polygon)
            new_intersect_polygon = intersect_polygon.difference(new_polygon)
            df.loc[df['osgb'] == pair[0], 'polygon_original'] = dumps(
                new_polygon, rounding_precision=2)
            df.loc[df['osgb'] == pair[1], 'polygon_original'] = dumps(
                new_intersect_polygon, rounding_precision=2)
    return df

def expand_floors(df):
    #create a separate entry for each floor, assign mixed use to commercial for ground floor
    #and residential above
    df= df.loc[df.index.repeat(df.floors)]
    df['floor_no'] = df.groupby(['osgb']).cumcount()+1
    mask1 = (df['use'] == 'Mixed Use_res') & (df['floor_no']==1)
    df['use'] = df['use'].mask(mask1,'Ret_Sales')
    mask2 = (df['use'] == 'Mixed Use_res') & (df['floor_no']>1)
    df['use'] = df['use'].mask(mask2,'Residential')
    mask3 = (df['use'] == 'Mixed Use_comm') & (df['floor_no']==1)
    df['use'] = df['use'].mask(mask3,'Ret_Sales')
    mask4 = (df['use'] == 'Mixed Use_comm') & (df['floor_no']>1)
    df['use'] = df['use'].mask(mask4,'Commercial')
    mask5 = (df['use'] == 'Mixed Use_ed') & (df['floor_no']==1)
    df['use'] = df['use'].mask(mask5,'Ret_Sales')
    mask6 = (df['use'] == 'Mixed Use_ed') & (df['floor_no']>1)
    df['use'] = df['use'].mask(mask6,'Education')
    mask7 = (df['use'] == 'Mixed Use_res') & (df['floor_no']==1)
    df['use'] = df['use'].mask(mask7,'Ret_Sales')
    mask8 = (df['use'] == 'Mixed Use_hotel') & (df['floor_no']>1)
    df['use'] = df['use'].mask(mask8,'Hotel')
    return(df)

def catchment_area(df):
    # Load the catchment area radius from the configuration file
    radius = config.getfloat('shading', 'catchment_radius')
    osgb_polygon_pairs = variable_polygon_pairs(df, 'osgb', 'AMD_polygon')
    for osgb, polygon in osgb_polygon_pairs:
        polygon_convex = polygon.convex_hull
        ca = polygon_convex.buffer(radius)

        osgbs_in_ca = list()
        for sur_osgb, sur_polygon in osgb_polygon_pairs:
            if sur_osgb != osgb:
                if sur_polygon.within(ca) or (sur_polygon.intersects(ca) and (
                        sur_polygon.intersection(ca).geom_type not in ['Point'])):
                    osgbs_in_ca.append(sur_osgb)
        df.loc[df['osgb'] == osgb, 'osgbs_in_ca'] = str(osgbs_in_ca)
    return df

def osgb_3D_out(df):
    df['height']= df.floors*3
    df.to_csv(os.path.join(ROOT_DIR, '1_preprocessed_osgb_3D_blocks.csv'), index=False)
    return df




def preprocessing_constructions():

    def constructions():

        # List of available construction sets
        construction_sets = ['improved', 'good', 'basic']
        # List of pre-defined construction elements
        construction_elements = [
            'AMD_wall', 'AMD_roof', 'AMD_ground_floor',
            'AMD_ceiling', 'AMD_floor', 'AMD_partition',
            'AMD_glazing']
        # Dictionary of construction element names for each construction set
        improved = {construction_elements[0]: 'improved_wall',
                            construction_elements[1]: 'improved_flat_roof',
                            construction_elements[2]: 'improved_solid_ground_floor',
                            construction_elements[3]: 'ceiling_improved',
                            construction_elements[4]: 'ceiling_improved_inverse',
                            construction_elements[5]: 'partition',
                            construction_elements[6]: 'improved_glazing'}
        good = {construction_elements[0]: 'good_wall',
                        construction_elements[1]: 'good_flat_roof',
                        construction_elements[2]: 'good_solid_ground_floor',
                        construction_elements[3]: 'ceiling_good',
                        construction_elements[4]: 'ceiling_good_inverse',
                        construction_elements[5]: 'partition',
                        construction_elements[6]: 'good_glazing'}
        basic = {construction_elements[0]: 'basic_wall',
                        construction_elements[1]: 'basic_flat_roof',
                        construction_elements[2]: 'basic_solid_ground_floor',
                        construction_elements[3]: 'ceiling_basic',
                        construction_elements[4]: 'ceiling_basic_inverse',
                        construction_elements[5]: 'partition',
                        construction_elements[6]: 'basic_glazing'}
        construction_dict = dict()
        for item in construction_sets:
            construction_dict[item] = eval(item)
        return construction_sets, construction_elements, construction_dict


    def age_related_constructions(row, construction_sets,
                             construction_elements, construction_dict, ):
        if (row['Use'])=='Residential':
            row['AMD_glazing_ratio'] = int(10 * rd.randint(1,3))
        elif (row['Use'])=='Commercial':
            row['AMD_glazing_ratio'] = int(10 * rd.randint(3,6))
        elif (row['Use'])=='Health':
            row['AMD_glazing_ratio'] = int(10 * rd.randint(3,6))
        elif (row['Use'])=='Hotel':
            row['AMD_glazing_ratio'] = int(10 * rd.randint(4,7))		
        elif (row['Use'])=='Education':
            row['AMD_glazing_ratio'] = int(10 * rd.randint(2,4))	
        else:
            row['AMD_glazing_ratio'] = int(10 * rd.randint(2, 3))
        #print(int(row['age']))
        if int(row['age']) <= 20:
            for item in construction_elements:
                row[item] = construction_dict['improved'][item]
        elif 20 < int(row['age']) <= 40:
            for item in construction_elements:
                row[item] = construction_dict['good'][item]		
        else:
            for item in construction_elements:
                row[item] = construction_dict['basic'][item]
        return row

    df = pd.read_csv(os.path.join(ROOT_DIR, constructions_file))
    const_sets, const_elements, const_dict = constructions()
    # Populate construction elements with random constructions unless data
    # provided externally
    df = df.apply(age_related_constructions, args=(
        const_sets, const_elements, const_dict,), axis=1)
    df.to_csv(os.path.join(
        ROOT_DIR, constructions_file[:-4] + '_populated.csv'), index=False)
    return df

def main():
    shapes_df = preprocessing()
    const_df = preprocessing_constructions()
    print(shapes_df.head)
    df = pd.merge(shapes_df,
     const_df, how='left',
                  on=['osgb'])
    df.to_csv(os.path.join(ROOT_DIR, 'amd_all.csv'), index=False)

def remove_duplicated_coordinates(df):
    '''
    Function which removes duplicated coordinates from Polygons within Pandas
    DataFrame
    '''
    # Loop through the DataFrame row-by-row
    for row in df.itertuples():
        # Load the original polygon as a Shapely object
        polygon = loads(row.polygon_original)
        # Extract exterior coordinates and remove duplicates
        ext_ring_coords = list(polygon.exterior.coords)
        ext_ring_no_dup = remove_duplicated_coords_from_list(ext_ring_coords)
        # Empty list for interior ring coordinates
        int_ring_no_dup_list = []
        # If there are holes in Polygon, loop through the list of holes
        if polygon.interiors:
            for item in polygon.interiors:
                # Extract inner coordinates and remove duplicates
                item_coords = list(item.coords)
                int_ring_no_dup = remove_duplicated_coords_from_list(
                    item_coords)
                # Append a list of holes
                int_ring_no_dup_list.append(int_ring_no_dup)
        # Create a new Polygon from coordinates without duplicates
        polygon_no_dup = Polygon(ext_ring_no_dup, int_ring_no_dup_list)
        # Save new Polygon with a two decimal places
        df.loc[row.Index, 'AMD_polygon_no_dup'] = dumps(polygon_no_dup,
                                                       rounding_precision=2)
    # Create a new column with a copy of Polygons without duplicated coords
    df['AMD_polygon'] = df['AMD_polygon_no_dup']
    # Return the DataFrame
    return df

def polygon_topology_check(df, touching, intersect, perimeter,
                           exposed, partition):
    '''
    Function which checks the relationship between polygons (separated,
    touching, intersected). Intersected polygons are not permitted. It also
    calculates the basic polygon parameters such as total perimeter length,
    exposed perimeter length and partition perimeter length.
    '''
    # Create five columns either empty or with empty lists in them
    df[touching] = '[]'
    df[intersect] = '[]'
    df[perimeter] = ''
    df[exposed] = ''
    df[partition] = ''

    # Create the list of osgb-polygon pairs
    osgb_polygon_pairs = variable_polygon_pairs(df, 'osgb', 'AMD_polygon')
    # Loop through osgb_polygon_pairs list
    for osgb, osgb_polygon in osgb_polygon_pairs:
        # Empty lists for touching and intersecting polygons
        osgb_touching = []
        osgb_intersect = []
        # Loop through osgb_polygon_pairs list to identify adjacent osgb's
        for osgb_adj, adj_polygon in osgb_polygon_pairs:
            # Do comparison only for unique osgb's
            if osgb_adj != osgb:
                # Determine if two polygons are touching (shapely) in more
                # than a single point
                if osgb_polygon.touches(adj_polygon) and (
                    osgb_polygon.intersection(
                        adj_polygon).geom_type not in ['Point']):
                    # Add touching polygon osgb to list
                    osgb_touching.append(osgb_adj)
                # Determine if two polygons are intersecting (shapely) in
                # more than a single point
                if osgb_polygon.intersects(adj_polygon) and (
                    osgb_polygon.intersection(
                        adj_polygon).geom_type not in ['Point']):
                    # Add intersecting polygon osgb to list
                    osgb_intersect.append(osgb_adj)

        # Save lists of touching and intersecting polygons to DataFrame
        df.loc[df['osgb'] == osgb, touching] = str(osgb_touching)
        df.loc[df['osgb'] == osgb, intersect] = str(osgb_intersect)
        # Check for intersection polygons and issue warning (all should
        # touch only)
        if len(osgb_touching) < len(osgb_intersect):
            difference = list(set(osgb_intersect) - set(osgb_touching))
            print('***WARNING: OSGB %s intersects following polygon(s): %s'
                  % (osgb, difference))

        # Convert outer/inner rings to LineString/MultiLineString
        outer_ring = LineString(osgb_polygon.exterior)
        inner_ring = MultiLineString(osgb_polygon.interiors)
        # Create the union of outer ring with inner ring
        total_ring = unary_union((outer_ring, inner_ring))
        # Calculate and save the polygon perimeter length
        polygon_perimeter = round(osgb_polygon.length, 3)
        df.loc[df['osgb'] == osgb, perimeter] = round(polygon_perimeter, 3)

        # If there are touching polygons
        if osgb_touching:
            # Loop through the list of touching polygons
            for item in osgb_touching:
                # Load the polygon as shapely object
                item_polygon = loads(
                    df.loc[df['osgb'] == item, 'AMD_polygon'].values[0])
                # Subtract the touching elements from the total ring
                total_ring -= osgb_polygon.intersection(item_polygon)
            # Exposed length is the remaining while partition length is
            # difference between perimeter length  and exposed length
            exposed_length = round(total_ring.length, 3)
            partition_length = polygon_perimeter - exposed_length
            # Save exposed length and partition length
            df.loc[df['osgb'] == osgb, exposed] = round(exposed_length, 3)
            df.loc[df['osgb'] == osgb, partition] = round(partition_length, 3)
        else:
            # If no touching elements, then exposed length is equal to
            # perimeter length while partition length is 0
            df.loc[df['osgb'] == osgb, exposed] = round(polygon_perimeter, 3)
            df.loc[df['osgb'] == osgb, partition] = 0

    # Delete the intersect column
    df = df.drop([intersect], axis=1)
    # Return the DataFrame
    return df

def polygon_tolerance(df):
    '''
    Function which checks whether polygon simplification is required or not.
    Simplification is required if the difference between two consecutive
    Polygon coordinates is less than a tolerance previously defined (default
    0.1m)
    '''

    def distance_within_tolerance(coords_list, tolerance):
        '''
        Internal function which loops through the list of coordinates and check
        the distance between two consecutive coordinates. If the distance is
        less than a tolerance the looping stops and True is returned. Otherwise
        the False is returned.
        '''
        for i, coord in enumerate(coords_list[:-1]):
            first = coords_list[i]
            second = coords_list[i + 1]
            distance = LineString([first, second]).length
            if distance < tolerance:
                return True
        return False
    # Loop through the DataFrame row-by-row
    for row in df.itertuples():
        # Load the polygon as shapely object
        polygon = loads(row.AMD_polygon)
        # Extract exterior ring coordinates and check the tolerance
        ext_ring_coords = list(polygon.exterior.coords)
        simplify_required = distance_within_tolerance(ext_ring_coords,
                                                      tolerance)
        # Create a new column with a polygon simplification boolean
        df.loc[row.Index, 'AMD_polygon_simplify'] = simplify_required
        # If simplification is required skip to the next row, otherwise loop
        # through interior rings if they exist
        if simplify_required:
            continue
        if polygon.interiors:
            # For hole within polygon check the tolerance within coordinates
            for item in polygon.interiors:
                item_coords = list(item.coords)
                simplify_required = distance_within_tolerance(item_coords,
                                                              tolerance)
                # If simplification is required, update the df column
                if simplify_required:
                    df.loc[row.Index,
                           'AMD_polygon_simplify'] = simplify_required
                    break
    # Return the DataFrame
    return df

def polygon_simplification(df, df_removed, simplify_polygon_no):
    '''
    Function which simplifies polygons. During simplification procedure some
    polygons initial excluded from simplification might be affected and get to
    the shape to require simplification. That's the reason of checking
    tolerance after simplification procedure and repeats simplification if
    there is a need for.
    '''

    def polygon_buffer(df):
        '''
        Internal function which searches for invalid polygons; buffer them; and
        than check for removed coordinates after buffering in order to update
        affected adjacent polygons (if any)
        '''

        def remove_buffered_coordinates(osgb, coords, new, removed):
            '''
            Internal function which removes coordinates based on the
            items in removed and new lists
            '''
            # Loop through the list of removed coordinates
            for r_c in removed:
                # For each removed coordinate check whether coordinate exists
                # in the list of coordinates. If it exists, replace it
                # with the closest coordinate from the list of new coordinates
                for i, coord in enumerate(coords):
                    if coord == r_c:
                        minimum_dist = LineString([r_c, new[0]]).length
                        replacement_coord = new[0]
                        for n_c in new:
                            dist = LineString([r_c, n_c]).length
                            if dist < minimum_dist:
                                minimum_dist = dist
                                replacement_coord = n_c
                        coords[i] = replacement_coord
                        print('Affected adjacent %s polygon adjusted' % osgb)
            # Remove duplicated coordinates
            coords = remove_duplicated_coords_from_list(coords)
            # Return coordinates
            return coords

        # Loop through the DataFrame and check each polygon for validity
        for row in df.itertuples():
            polygon = loads(row.AMD_polygon)
            # Buffer invalid polygons and check adjacent polygons for update
            if not polygon.is_valid:
                # Buffer a polygon and report the osgb of buffered polygon
                new_polygon = polygon.buffer(0)
                print('Self-intersection of %s polygon buffered' % row.osgb)
                # Save buffered polygon back to the DataFrame
                df.loc[row.Index, 'AMD_polygon'] = dumps(new_polygon,
                                                        rounding_precision=2)
                # Buffered polygon with two decimal places
                new_polygon = loads(dumps(new_polygon, rounding_precision=2))
                # Check if polygon has adjacent polygons
                osgb_touching = literal_eval(row.AMD_initial_touching)
                # If there are adjacent polygons, check if they need update
                if osgb_touching:
                    # List of new coordinates as a result of buffering
                    new_coords = list(set(list(new_polygon.exterior.coords)) - set(list(polygon.exterior.coords)))
                    # List of dropped coordinates as a result of buffering
                    removed_coords = list(set(list(polygon.exterior.coords)) - set(list(new_polygon.exterior.coords)))
                    if new_coords:
                        # If there are new coordinates loop through the list of
                        # adjacent polygons; load a polygon as Shapely object;
                        # extract outer ring; remove/replace buffered
                        # coordinates from outer ring (if any); save updated
                        # polygon back to the DataFrame
                        for t in osgb_touching:
                            t_polygon = loads(
                                df.loc[df['osgb'] == t,
                                       'AMD_polygon'].values[0])
                            t_polygon_coords = list(t_polygon.exterior.coords)
                            t_polygon_coords = remove_buffered_coordinates(t,
                                t_polygon_coords, new_coords, removed_coords)
                            new_t_polygon = Polygon(t_polygon_coords,
                                                    t_polygon.interiors)
                            df.loc[df['osgb'] == t, 'AMD_polygon'] = dumps(
                                new_t_polygon, rounding_precision=2)
        # Return the DataFrame
        return df

    while simplify_polygon_no > 0:
        # Determine polygons within holes
        df = polygon_within_hole(df)
        # Simplify polygons
        df = polygon_simplify(df)
        # Remove records which Polygons are invalid (less than 3 coordinates)
        dropped = df.loc[df['AMD_polygon'].isin([False])].reset_index(drop=True)
        df_removed = df_removed.append(dropped)
        df = df.loc[~df['AMD_polygon'].isin([False])].reset_index(drop=True)
        # In extremely rare occasions the simplification results in invalid
        # polygons due to self-intersection of outer ring. This can be sorted
        # in most cases with buffering polygon while preserving surrounding
        # topology
        df = polygon_buffer(df)
        # Check whether further simplification has to be conducted
        df = polygon_tolerance(df)
        simplify_polygon_no = 0
        for row in df.itertuples():
            if row.AMD_polygon_simplify:
                simplify_polygon_no += 1
    # Remove two columns from the DataFrame
    df = df.drop(['AMD_polygon_simplify', 'AMD_polygon_within_hole'], axis=1)
    # Return DataFrames
    return df, df_removed

def collinear_exterior(df):
    '''
    Function which searches for collinear points and creates the polygon which
    is used for horizontal surfaces (roof/ceiling/floor). Also checks for
    non-convex horizontal surfaces
    '''

    def collinear_points_list(objects_list):
        '''
        Function which creates a list of collinear points
        '''
        coll_list = []
        # For MultiLineStrins get collinear points for each LineSting
        if objects_list.geom_type in ['MultiLineString', 'GeometryCollection']:
            for item in objects_list:
                coll_points = coollinear_points(list(item.coords))
                if coll_points:
                    coll_list.append(coll_points)
        # Obtain collinear point for LineSting
        elif objects_list.geom_type == 'LineString':
            coll_points = coollinear_points(list(objects_list.coords))
            if coll_points:
                coll_list.append(coll_points)
        # Generate the list of collinear points
        collinear_points_list = []
        for item in coll_list:
            for i in item:
                collinear_points_list.append(i)
        # Return the list of collinear points
        return collinear_points_list
    def coollinear_points(coord_list):
        '''
        Internal function which searches for collinear points and return a list
        of collinear points
        '''

        def check_collinearity(f, m, l):
            '''
            Internal sub function which checks collinearity of three
            consecutive coordinate pairs by checking the area of polygon formed
            # of these three consecutive coordinate pairs
            '''
            if Polygon([f, m, l]).area <= 1e-9:
                return True
            else:
                return False

        # Empty list for removed coordinates
        removed_coll = []
        if len(coord_list) >= 3:
            # Take three consecutive coordinate pairs and check collinearity
            for i in range(len(coord_list) - 2):
                first = coord_list[i]
                middle = coord_list[i + 1]
                last = coord_list[i + 2]
                if check_collinearity(first, middle, last):
                    removed_coll.append(coord_list[i + 1])
        # Return the list of collinear points
        return removed_coll

    def update_polygon(polygon, points_to_remove):
        '''
        Function which updates a polygon by removing collinear points (provided
        in the list "points_to_remove") from outer and inner rings
        '''
        outer_ring = list(polygon.exterior.coords)
        new_outer_ring = remove_items_from_list(outer_ring, points_to_remove)
        new_inner_rings = []
        if polygon.interiors:
            for item in polygon.interiors:
                inner_ring = list(item.coords)
                new_inner_ring = remove_items_from_list(inner_ring,
                                                        points_to_remove)
                new_inner_rings.append(new_inner_ring)
        new_polygon = Polygon(new_outer_ring, new_inner_rings)
        # Return the updated polygon
        return new_polygon

    def remove_items_from_list(coords, items):
        '''
        Function which removes items from the coordinates list. In case the
        first and last coordinates are they are removed, append the list of
        coordinates with the new first coordinate
        '''
        if coords[0] == coords[-1]:
            for i in items:
                coords = [x for x in coords if x != i]
            if coords[0] != coords[-1]:
                coords.append(coords[0])
        else:
            for i in items:
                coords = [x for x in coords if x != i]
        # Return the coordinates list
        return coords

    def update_exposed(exposed_ring, points_to_remove):
        '''
        Function which removes collinear points from the exposed parts of a
        polygon
        '''
        # In the case of MultiLineSting check each LineString in order to
        # remove collinear points (stored in the points_to_remove list)
        if exposed_ring.geom_type == 'MultiLineString':
            new_ms = []
            for item in exposed_ring:
                new_item = remove_items_from_list(list(item.coords),
                                                  points_to_remove)
                # If the new LineString has more than one coordinate than
                # append the list of new MultiLineStirngs
                if len(new_item) > 1:
                    new_ms.append(new_item)
            # If MultiLineString list has more than one item create a
            # MultiLineStrng Shapely object. If there is only one item create
            # a LineString Shapely object. If there are no items than create
            # geometryCollection Empty Shapely object.
            if len(new_ms) > 1:
                new_exposed_ring = MultiLineString(new_ms)
            elif len(new_ms) == 1:
                new_exposed_ring = LineString(new_ms[0])
            else:
                new_exposed_ring = loads('GEOMETRYCOLLECTION EMPTY')
        elif exposed_ring.geom_type == 'LineString':
            new_item = remove_items_from_list(list(exposed_ring.coords),
                                              points_to_remove)
            if len(new_item) > 1:
                new_exposed_ring = LineString(new_item)
            else:
                new_exposed_ring = loads('GEOMETRYCOLLECTION EMPTY')
        else:
            new_exposed_ring = loads('GEOMETRYCOLLECTION EMPTY')
        # Return the new exposed ring
        return new_exposed_ring

    def remove_coollinear_points_horizontal(polygon):
        '''
        Function which removes collinear points from polygon to be used in
        creation of horizontal surfaces. This can be different polygon from
        initial polygon in cases where two collinear points are placed on
        partition wall and exposed wall.
        '''
        coll_list = []
        # Convert outer/inner rings to LineString/MultiLineString
        o_r = LineString(polygon.exterior)
        i_r = MultiLineString(polygon.interiors)
        # Create the union of outer ring with inner ring
        t_t = unary_union((o_r, i_r))
        # For MultiLineStrins get collinear points for each LineSting
        if t_t.geom_type == 'MultiLineString':
            for item in t_t:
                coords = list(item.coords)
                coords.append(coords[1])
                coll_points = coollinear_points(coords)
                if coll_points:
                    coll_list.append(coll_points)
        # Obtain collinear points for LineString
        elif t_t.geom_type == 'LineString':
            coords = list(t_t.coords)
            coords.append(coords[1])
            coll_points = coollinear_points(coords)
            if coll_points:
                coll_list.append(coll_points)
        # Returned the cleaned coordinates list
        collinear_points_list = []
        for item in coll_list:
            for i in item:
                collinear_points_list.append(i)
        # Create a new polygon by removing collinear points
        new_polygon = update_polygon(polygon, collinear_points_list)
        # Return the new polygon
        return new_polygon

    def triangulate(polygon):
        '''
        Internal function which checks for non-convex polygons
        '''
        # If there is a hole within polygon it is automatically non-convex
        if polygon.interiors:
            return True
        else:
            # Extract the exterior ring and append it with the second point
            outer = polygon.exterior
            coords = outer.coords
            coords_ext = list(coords)
            coords_ext.append(coords[1])
            # Loop through three consecutive coordinate pairs and check ccw
            for i in range(len(coords_ext) - 2):
                first = coords_ext[i]
                middle = coords_ext[i + 1]
                last = coords_ext[i + 2]
                if LinearRing([first, middle, last]).is_ccw:
                    return True
            return False

    # ************** Main body of the function ****************
    # Loop through the DataFrame row-by-row to remove collinear points from
    # partitions
    for row in df.itertuples():
        # Get the touching polygons list
        if 'AMD_simplified_touching' in df.columns:
            osgb_touching = literal_eval(row.AMD_simplified_touching)
        else:
            osgb_touching = literal_eval(row.AMD_initial_touching)
        # Load the polygon as Shapely object and get the osgb number
        polygon = loads(row.AMD_polygon)
        osgb = row.osgb
        # If there are touching polygons than check for collinear points on
        # partition walls only
        if osgb_touching:
            # Loop through the list of touching polygons and load touching
            # polygon as Shapely object
            for t in osgb_touching:
                t_polygon = loads(
                    df.loc[df['osgb'] == t, 'AMD_polygon'].values[0])
                # Partition wall (intersecting polygon with touching polygon)
                partition = polygon.intersection(t_polygon)
                # If partition wall is MultiLineString try to merge LineStrings
                if partition.geom_type == 'MultiLineString':
                    partition = linemerge(partition)
                # Check for collinear points
                partition_collinear_points = collinear_points_list(partition)
                # if there are collinear points update both polygons
                if partition_collinear_points:
                    polygon = update_polygon(
                        polygon, partition_collinear_points)
                    t_polygon = update_polygon(
                        t_polygon, partition_collinear_points)
                    # Save updated polygons to the DataFrame
                    df.loc[df['osgb'] == t, 'AMD_polygon'] = dumps(
                        t_polygon, rounding_precision=2)
                    df.loc[df['osgb'] == osgb, 'AMD_polygon'] = dumps(
                        polygon, rounding_precision=2)

    # Loop through the DataFrame row-by-row to remove collinear points from
    # exposed walls
    for row in df.itertuples():
        # Get the touching polygons list
        if 'AMD_simplified_touching' in df.columns:
            osgb_touching = literal_eval(row.AMD_simplified_touching)
        else:
            osgb_touching = literal_eval(row.AMD_initial_touching)
        # Load the polygon as Shapely object and get the osgb number
        polygon = loads(row.AMD_polygon)
        osgb = row.osgb
        # If there are touching polygons than for the initial exposed surfaces
        # from outer and inner coordinates
        if osgb_touching:
            # Convert outer/inner rings to LineString/MultiLineString
            outer_ring = LineString(polygon.exterior)
            inner_ring = MultiLineString(polygon.interiors)
            # Create the union of outer ring with inner ring
            exposed = unary_union((outer_ring, inner_ring))
            # Loop through the list of touching polygons and a subtract
            # partition walls from initial exposed surfaces
            for t in osgb_touching:
                t_polygon = loads(
                    df.loc[df['osgb'] == t, 'AMD_polygon'].values[0])
                exposed -= polygon.intersection(t_polygon)
            # Clean a total ring from collinear points
            exposed_collinear_points = collinear_points_list(exposed)
            # if there are collinear points update both exposed surfaces and
            # polygon
            if exposed_collinear_points:
                exposed = update_exposed(exposed, exposed_collinear_points)
                polygon = update_polygon(polygon, exposed_collinear_points)
            # Remove collinear points from polygon used for horizontal surfaces
            horizontal = remove_coollinear_points_horizontal(polygon)
        else:
            # If there are no touching polygons all surfaces are exposed
            polygon = remove_coollinear_points_horizontal(polygon)
            horizontal = polygon
            # Convert outer/inner rings to LineString/MultiLineString
            outer_ring = LineString(polygon.exterior)
            inner_ring = MultiLineString(polygon.interiors)
            # Create the union of outer ring with inner ring
            exposed = unary_union((outer_ring, inner_ring))

        # Create a new column for info on non-convex polygons
        polygon_triangulation = triangulate(polygon)
        # Store data on cleaned polygons, exposed walls, horizontal surfaces
        # and triangulation in the DataFrame
        df.loc[df['osgb'] == osgb, 'AMD_polygon_exposed_wall'] = dumps(
            exposed, rounding_precision=2)
        df.loc[df['osgb'] == osgb, 'AMD_polygon'] = dumps(
            polygon, rounding_precision=2)
        df.loc[df['osgb'] == osgb, 'AMD_polygon_horizontal'] = dumps(
            horizontal, rounding_precision=2)
        df.loc[df['osgb'] == osgb, 'AMD_triangulation'] = polygon_triangulation
    # Return the DataFrame
    return df

def variable_polygon_pairs(df, variable_value, polygon_value):
    '''
    Function which converts two DataFrame columns (second is polygon) into list
    of variable-polygon pairs (polygon is loaded with Shapely)
    '''
    variable_polygon_values = df[[variable_value, polygon_value]].values
    variable_polygon = []
    for variable, polygon in variable_polygon_values:
        polygon = loads(polygon)
        variable_polygon.append([variable, polygon])
    return variable_polygon

def remove_duplicated_coords_from_list(coords_list):
    '''
    Function which removes duplicated coordinates from the list of coordinates
    '''
    # Empty list for cleaned list of coordinates
    coords_list_no_dup = []
    if coords_list[0] == coords_list[-1]:
        # If first and last coordinates are the same (Polygon/LinearRing),
        # remove the last coordinate.
        del coords_list[-1]
        # Compare lengths of the list of coordinates with the list of unique
        # elements. If lengths differ than loop through the list of coordinates
        # and append coords_list_no_dup with unique elements. Otherwise, there
        # are no duplicates
        if len(coords_list) != len(set(coords_list)):
            [coords_list_no_dup.append(x)
                for x in coords_list if x not in coords_list_no_dup]
        else:
            coords_list_no_dup = coords_list
        # Append the list without duplicates with the first coordinate
        coords_list_no_dup.append(coords_list_no_dup[0])
    else:
        # Compare lengths of the list of coordinates with the list of unique
        # elements. If lengths differ than loop through the list of coordinates
        # and append coords_list_no_dup with unique elements. Otherwise, there
        # are no duplicates
        if len(coords_list) != len(set(coords_list)):
            [coords_list_no_dup.append(x)
                for x in coords_list if x not in coords_list_no_dup]
        else:
            coords_list_no_dup = coords_list
    # Return the list of coordinates without duplicates
    return coords_list_no_dup

def polygon_within_hole(df):
    '''
    Function searches for "polygon" within other "polygon with hole".
    Used in polygon simplification
    '''
    # Loop through the DataFrame row-by-row
    for row in df.itertuples():
        # Load the Polygon as Shapely object
        osgb_polygon = loads(row.AMD_polygon)
        # Empty list to store polygons within other polygons
        osgb_within_hole = []
        # Check if Polygon has the hole
        if osgb_polygon.interiors:
            # Extract polygons osgb
            osgb = row.osgb
            # Create the list of osgb-polygon pairs
            osgb_polygon_pairs = variable_polygon_pairs(df, 'osgb', 'AMD_polygon')
            # Loop through the holes
            for item in osgb_polygon.interiors:
                # Convert hole to the Polygon by reversing coordinates
                item_polygon = Polygon(item.coords[::-1])
                # Loop through the list of osgbs
                for osgb_adj, adj_polygon in osgb_polygon_pairs:
                    # For other polygons load the Polygon as Shapely object
                    if osgb_adj != osgb:
                        # If other Polygon is within the hole and touches hole
                        if item_polygon.contains(adj_polygon) and osgb_polygon.touches(adj_polygon):
                            # Append the list of polygons within other polygons
                            # with hole and issue warning if polygon has hole
                            osgb_within_hole.append(osgb_adj)
                            if adj_polygon.interiors:
                                print('***WARNING: OSGB with the hole (%s) WITHIN other OSGB with the hole (%s)' % (osgb_adj, osgb))
        # Create a new column in the DataFrame for a list of polygons
        df.loc[row.Index, 'AMD_polygon_within_hole'] = str(osgb_within_hole)
    # Return DataFrame
    return df

def polygon_simplify(df):
    '''
    Internal function which simplifies polygons while preserving topology
    '''
    def touching_poly(df, osgb, polygon, osgb_list, osgb_touching):
        '''
        Internal function to determine touching polygons including polygons
        which touch in the single point
        '''
        # Loop through the list of polygons
        for t in osgb_list:
            if t != osgb:
                # if potential touching polygon exists then check if it touches
                if df.loc[df['osgb'] == t, 'AMD_polygon'].values[0]:
                    t_polygon = loads(df.loc[df['osgb'] == t,
                                             'AMD_polygon'].values[0])
                    if polygon.touches(t_polygon):
                        # Append the list of touching polygons
                        osgb_touching.append(t)
        # Return the list of touching polygons
        return osgb_touching

    def polygon_simplifying(polygon, tolerance, df, osgb, osgb_touching):
        '''
        Internal function which does the main simplification work
        '''

        def coords_cleaning(coords, remove_leave_pairs, tolerance):
            '''
            Internal function which applies radial distance method to remove
            coordinates which are within tolerance distance
            '''
            def radial_dist_simplify(coords, tolerance):
                '''
                Internal function which checks two consecutive coordinates for
                tolerance distance. Breaks as soon as it finds the first one.
                '''

                def remove_item_from_list(coords, item):
                    '''
                    Internal function which remove duplicated coordinates
                    '''
                    if coords[0] == coords[-1]:
                        coords = [x for x in coords if x != item]
                        if coords[0] != coords[-1]:
                            coords.append(coords[0])
                    else:
                        coords = [x for x in coords if x != item]
                    # Return coordinates
                    return coords

                # Empty list for remove-leave-pair
                remove_leave_pair = []
                # Loop through the coordinates list
                for i, coord in enumerate(coords[:-1]):
                    # First and second coordinate and their distance
                    first = coords[i]
                    second = coords[i + 1]
                    distance = LineString([first, second]).length
                    # If  distance is smaller than tolerance
                    if distance <= tolerance:
                        # For all coordinates except the last one remove second
                        # coordinate while keeping the first one
                        if i < (len(coords) - 2):
                            coord_remove = second
                            coord_leave = first
                        # In the case of the last coordinate remove the first
                        # one and keep the second
                        else:
                            coord_remove = first
                            coord_leave = second
                        # Make a list of remove-leave coordinate pair and
                        # append list of remove-leave-pairs
                        remove_leave_pair = [coord_remove, coord_leave]
                        coords = remove_item_from_list(coords, coord_remove)
                        # Return cleaned coordinates and remove-leave pairs
                        return coords, remove_leave_pair
                # Return cleaned coordinates and remove-leave pairs
                return coords, remove_leave_pair

            # Initial coordinates length increased by 1 for looping purposes
            coords_lenght = len(coords) + 1
            # Apply radial simplification as long as there are more than 3
            # coordinates and returned list of coordinates is smaller than
            # initial list of coordinates (this means a coordinate is removed
            # due to simplification)
            while (len(coords) < coords_lenght) and (len(coords) > 3):
                coords_lenght = len(coords)
                # Apply simplification and return coordinates list and
                # remove-leave-pair
                coords, r_l_pair = radial_dist_simplify(coords, tolerance)
                # If remove-leave-pair exists append remove-leave-pair list
                if r_l_pair:
                    remove_leave_pairs.append(r_l_pair)
            # Return cleaned coordinates and remove-leave-pair list
            return coords, remove_leave_pairs

        def simplified_coords(polygon, ring_position,
                              remove_leave_pairs, tolerance):
            '''
            Internal function which checks the distance between consecutive
            coordinates and if distance is less than tolerance then merge these
            two consecutive points into one
            '''
            # Apply coords cleaning function to the outer ring coordinates list
            if ring_position == 'outer':
                coords = list(polygon.exterior.coords)
                coords, remove_leave_pairs = coords_cleaning(
                    coords, remove_leave_pairs, tolerance)
                # Return cleaned coordinates and remove-leave-pairs list
                return coords, remove_leave_pairs
            # Apply coords cleaning function to each of inner rings
            elif ring_position == 'inner':
                inner_coords_list = []
                for inner in polygon.interiors:
                    coords = list(inner.coords)
                    coords, remove_leave_pairs = coords_cleaning(
                        coords, remove_leave_pairs, tolerance)
                    # If list of cleaned coordinated has more than 3 items
                    # append inner_coords list
                    if len(coords) > 3:
                        inner_coords_list.append(coords)
                # Return inner coords list and remove-leave-pairs list
                return inner_coords_list, remove_leave_pairs

        def remove_cleaned_coordinates(coords, remove_leave_pairs):
            '''
            Internal function which removes coordinates based on the
            remove-leave-pairs
            '''
            # Loop through the list of remove-leave-pairs
            for pair in remove_leave_pairs:
                # For each remove-leave pair check whether remove coordinate
                # exists in the list of coordinates. If it exists, replace it
                # with leave coordinate
                for i, coord in enumerate(coords):
                    if coord == pair[0]:
                        coords[i] = pair[1]
            # Remove duplicated coordinates
            coords = remove_duplicated_coords_from_list(coords)
            # Return coordinates
            return coords

        def simplification_affects_inner_ring(adjacent_polygon,
                                              polygon, remove_leave_pairs):
            '''
            Internal function which adjust coordinates of the adjacent polygon
            hole in cases where simplified polygon is within the the adjacent
            polygon hole
            '''
            # Check only inner coordinates
            adjacent_inner_coords_list = []
            # Loop through the list of inner coordinates convert holes to
            # Polygons
            for inner in adjacent_polygon.interiors:
                inner_coords = list(inner.coords)
                inner_polygon = Polygon(inner_coords)
                # If hole polygon contains simplified polygon
                if inner_polygon.contains(polygon):
                    # apply cleaning coordinates function
                    inner_coords = remove_cleaned_coordinates(
                        inner_coords, remove_leave_pairs)
                adjacent_inner_coords_list.append(inner_coords)
            # Form the adjacent polygon from exterior ring and updated inner
            # rings
            adjacent_polygon = Polygon(adjacent_polygon.exterior,
                                       adjacent_inner_coords_list)
            # Return the adjacent polygon
            return adjacent_polygon

        def ccw_iterior_ring(coords):
            '''
            Function which checks the direction of coordinates (different for
            holes and outer ring). Return boolean
            '''
            if not LinearRing(coords).is_ccw:
                coords = coords[::-1]
            return coords

        def not_valid_polygons(osgb, polygon, df, outer_coords,
                               polygon_within_hole, remove_leave_pairs,
                               tolerance):
            '''
            Function which does polygon correction if polygon becomes irregular
            after simplification. Irregular polygons have outer ring
            intersecting the inner ring
            '''

            def remove_hole_if_inner_is_removed(df, inner_polygon,
                                                polygon_within_hole):
                '''
                Internal function which remove the polygons from hole if hole
                is removed as a result of simplification
                '''
                # Loop through the list of polygons within hole
                for p in polygon_within_hole:
                    # Check whether the polygon string exists
                    p_polygon = df.loc[df['osgb'] == p, 'AMD_polygon'].values[0]
                    if p_polygon:
                        # Load polygon as Shapely object
                        p_polygon = loads(p_polygon)
                        # Check if hole contains the polygon
                        if inner_polygon.contains(p_polygon):
                            # If yes, remove the polygon (set to False)
                            df.loc[df['osgb'] == p, 'AMD_polygon'] = False
                            # Check if removed polygon has its own polygons
                            # within holes. If yes, remove them too.
                            p_polygon_within_hole = literal_eval(
                                df.loc[df['osgb'] == p,
                                       'AMD_polygon_within_hole'].values[0])
                            if p_polygon_within_hole:
                                df = remove_holes(df, p_polygon_within_hole)
                # Return the DataFrame
                return df

            def p_is_not_valid(p):
                '''
                Function which examines the polygon validity. If inner rings
                touches outer ring, polygon is valid; if they intersect,
                polygon is not valid
                '''
                ex = LinearRing(p.exterior)
                for i, inner in enumerate(p.interiors):
                    in_i = Polygon(inner.coords)
                    if not ex.touches(in_i):
                        if ex.intersects(in_i):
                            return True
                return False

            # ************ Main body of this internal function ****************
            # Check if polygon is valid
            if p_is_not_valid(polygon):
                # If not create an eroded polygon from outer ring
                eroded_outer = Polygon(outer_coords).buffer(-tolerance)
                # List to store updated inner rings to form valid polygon
                eroded_inner_list = []
                # Loop through the list of inner rings
                for inner in polygon.interiors:
                    # Convert hole to the polygon
                    inner_polygon = Polygon(inner.coords)
                    # If hole can fit to the eroded polygon than append the
                    # list of inner holes
                    if inner_polygon.within(eroded_outer):
                        eroded_inner_list.append(inner)
                    else:
                        # Check if eroded polygon intersect the hole. If not
                        # that means hole is outside the eroded polygon and it
                        # has to be removed altogether with eventual polygons
                        # within the hole
                        if eroded_outer.intersects(inner_polygon):
                            # Hole new polygon is an intersection of the eroded
                            # polygon and the hole
                            new_inner_polygon = eroded_outer.intersection(
                                inner_polygon)
                            new_inner_coords = list(
                                new_inner_polygon.exterior.coords)
                            # Check if the hole new polygon has more than three
                            # coordinates. If not than remove it altogether
                            # with eventual polygons within the hole
                            if len(new_inner_coords) > 3:
                                # Check if there are polygons within hole
                                if polygon_within_hole:
                                    # Loop through polygons within hole list
                                    for p in polygon_within_hole:
                                        # Extract polygon string (if any)
                                        p_polygon = df.loc[
                                            df['osgb'] == p,
                                            'AMD_polygon'].values[0]
                                        # If there is polygon string
                                        if p_polygon:
                                            # Load the polygon and check if the
                                            # original hole contains the loaded
                                            # polygon
                                            p_polygon = loads(p_polygon)
                                            if inner_polygon.contains(p_polygon):
                                                # Extract the list of outer coordinates and clean them (remove close points)
                                                p_outer_coords = list(p_polygon.exterior.coords)
                                                p_outer_coords = remove_cleaned_coordinates(p_outer_coords, remove_leave_pairs)
                                                # Create the new polygon within hole by intersecting the hole new polygon
                                                # with cleaned polygon within hole
                                                new_p_polygon = Polygon(new_inner_coords).intersection(Polygon(p_outer_coords))
                                                # Create the difference between the hole new polygon and old polygon within hole outer ring
                                                new_inner_diff = Polygon(new_inner_coords).difference(Polygon(p_outer_coords))
                                                # updated hole new polygon is made as union of the inner difference
                                                # and the polygon within hole outer ring
                                                united_inner_polygon = unary_union((new_inner_diff, new_p_polygon))
                                                # Check for further simplification of the hole new polygon
                                                remove_leave_pairs_inner = []
                                                new_inner_coords, remove_leave_pairs_inner = simplified_coords(united_inner_polygon, 'outer', remove_leave_pairs_inner, tolerance)
                                                # If simplification appeared than clean coordinates of the new polygon within the hole new polygon
                                                p_outer_coords = list(new_p_polygon.exterior.coords)
                                                p_outer_coords = remove_cleaned_coordinates(p_outer_coords, remove_leave_pairs_inner)
                                                # Check if new polygon within hole has more than three coordinates
                                                if len(p_outer_coords) > 3:
                                                    # If yes than check if the new polygon within the hole has hole
                                                    if p_polygon.interiors:
                                                        # If yes create the new polygon within hole from new outer ring and holes
                                                        print('***WARNING - Simplification: OSGB with the hole (%s) WITHIN other OSGB with the hole (%s)' % (p, osgb))
                                                        new_p_polygon = Polygon(p_outer_coords, p_polygon.interiors)
                                                        # Check for the validity of this new polygon within hole
                                                        p_polygon_within_hole = literal_eval(df.loc[df['osgb'] == p, 'AMD_polygon_within_hole'].values[0])
                                                        mock_list = []
                                                        df, new_p_polygon = not_valid_polygons(p, new_p_polygon, df, p_outer_coords, p_polygon_within_hole, mock_list, tolerance)
                                                    else:
                                                        # If no holes than new polygon within hole is formed form outer ring only
                                                        new_p_polygon = Polygon(p_outer_coords)
                                                    # Save it to the DataFrame
                                                    df.loc[df['osgb'] == p, 'AMD_polygon'] = dumps(
                                                        new_p_polygon, rounding_precision=2)
                                                # Remove polygon within the hole
                                                else:
                                                    df.loc[df['osgb'] == p, 'AMD_polygon'] = False
                                # If there are no polygons within hole than
                                # check if the hole new polygon need to be
                                # further simplified
                                else:
                                    mock_list = []
                                    new_inner_coords, _ = simplified_coords(
                                        new_inner_polygon, 'outer',
                                        mock_list, tolerance)
                                # Check if the hole new polygon has more than
                                # three coordinates after simplification
                                if len(new_inner_coords) > 3:
                                    # If yes, than convert the coordinates
                                    # direction and append the list of holes
                                    new_inner_coords = ccw_iterior_ring(
                                        new_inner_coords)
                                    eroded_inner_list.append(new_inner_coords)
                            # Remove polygons within hole (if any)
                            elif polygon_within_hole:
                                df = remove_hole_if_inner_is_removed(
                                    df, inner_polygon, polygon_within_hole)
                        # Remove polygons within hole (if any)
                        elif polygon_within_hole:
                            df = remove_hole_if_inner_is_removed(
                                df, inner_polygon, polygon_within_hole)
                # Create a valid polygon from outer coordinates and a list of
                # updated holes
                polygon = Polygon(outer_coords, eroded_inner_list)
            # Return the DataFrame and updated polygon (if update occurred)
            return df, polygon

        def remove_holes(df, polygon_within_hole):
            '''
            Internal function which removes polygon within holes if the hole
            does not exist anymore after simplification
            '''
            # If there are polygons within holes than loop through the list
            if polygon_within_hole:
                for p in polygon_within_hole:
                    # Load the polygon within hole string and check if that
                    # polygon has hole and another polygon within that hole
                    p_polygon = df.loc[df['osgb'] == p, 'AMD_polygon'].values[0]
                    p_polygon_within_hole = literal_eval(
                        df.loc[df['osgb'] == p,
                               'AMD_polygon_within_hole'].values[0])
                    # In case of having polygon within polygon with hole than
                    # apply this function again to remove that polygon (False)
                    if p_polygon and p_polygon_within_hole:
                        df = remove_holes(df, p_polygon_within_hole)
                    # Remove polygon (False)
                    df.loc[df['osgb'] == p, 'AMD_polygon'] = False
            # Return the DataFrame
            return df

        # *********** Body of the polygon simplifying internal function ******
        # Load the list of polygons within holes (if any)
        polygon_within_hole = literal_eval(
            df.loc[df['osgb'] == osgb, 'AMD_polygon_within_hole'].values[0])
        rlp = []  # empty list to store remove-leave coordinate pairs
        # List of simplified outer coordinates and list of remove-leave
        # coordinate pairs
        outer_coords, rlp = simplified_coords(polygon, 'outer', rlp, tolerance)
        # Check if simplified coordinates list has more than three coordinates.
        # 4 coordinates form a triangle. If not, collapse polygon (False) and
        # (if required) remove polygons within holes
        if len(outer_coords) > 3:
            # Conduct simplification on holes if they exist
            if polygon.interiors:
                # Simplified hole coordinates list and appended list of
                # remove-leave coordinate pairs
                inner_coords_list, rlp = simplified_coords(polygon, 'inner',
                                                           rlp, tolerance)
                # Form the new polygon from simplified lists of outer
                # coordinates and inner coordinates
                new_polygon = Polygon(outer_coords, inner_coords_list)
                # In case of completely despairing inner coordinates delete
                # polygons within  holes (if any)
                if not inner_coords_list:
                    df = remove_holes(df, polygon_within_hole)
                # Check the polygon validity. In some cases polygon becomes
                # invalid (hole intersects outer ring) due to simplification
                df, new_polygon = not_valid_polygons(osgb, new_polygon, df,
                                                     outer_coords,
                                                     polygon_within_hole,
                                                     rlp, tolerance)
            # If there are no inner holes, form the new polygon from only
            # simplified list of outer coordinates and save it back to the
            # DataFrame with rounding precision of 2 decimal places
            else:
                new_polygon = Polygon(outer_coords)
            df.loc[df['osgb'] == osgb,
                   'AMD_polygon'] = dumps(new_polygon, rounding_precision=2)
        # Save 'False' new polygon and remove polygons from holes (if any)
        else:
            df.loc[df['osgb'] == osgb, 'AMD_polygon'] = False
            df = remove_holes(df, polygon_within_hole)

        # If there are removed points and touching polygons than adjust
        # coordinates in touching polygons
        if rlp and osgb_touching:
            # Loop through the list of touching polygons
            for t in osgb_touching:
                # Check if touching polygon string exists (it can be False too)
                t_polygon = df.loc[df['osgb'] == t, 'AMD_polygon'].values[0]
                if t_polygon:
                    # Load polygons within holes (if any)
                    t_polygon_within_hole = literal_eval(
                        df.loc[df['osgb'] == t,
                               'AMD_polygon_within_hole'].values[0])
                    # Load touching polygon as Shapely object
                    t_polygon = loads(t_polygon)
                    # Load osgb_polygon string (it can be False too)
                    osgb_polygon = df.loc[df['osgb'] == osgb,
                                          'AMD_polygon'].values[0]
                    # Check if osbg is a polygon within the hole of a touching
                    # polygon and osgb polygon still exists after polygon
                    # simplification than update the hole of a touching polygon
                    if osgb_polygon and t_polygon_within_hole and (osgb in t_polygon_within_hole):
                        t_polygon = simplification_affects_inner_ring(
                            t_polygon, polygon, rlp)
                    # List the outer ring coordinates of a touching polygon
                    t_outer_coords = list(t_polygon.exterior.coords)
                    # Adjust coordinates by checking the remove-leave pairs
                    t_outer_coords = remove_cleaned_coordinates(
                        t_outer_coords, rlp)
                    # Check if simplified outer coordinates list has more than
                    # three coordinates. If not, collapse polygon (False) and
                    # (if required) remove polygons within holes
                    if len(t_outer_coords) > 3:
                        if t_polygon.interiors:
                            # Form the new touching polygon from lists of outer
                            # coordinates and inner coordinates
                            t_polygon = Polygon(t_outer_coords,
                                                t_polygon.interiors)
                            # Check the polygon validity.
                            df, t_polygon = not_valid_polygons(
                                t, t_polygon, df, t_outer_coords,
                                t_polygon_within_hole, rlp, tolerance)
                        # If there are no inner holes, form the new polygon
                        # from only simplified list of outer coordinates and
                        # save it back to the DataFrame
                        else:
                            t_polygon = Polygon(t_outer_coords)
                        df.loc[df['osgb'] == t, 'AMD_polygon'] = dumps(
                            t_polygon, rounding_precision=2)
                    # Save 'False' new touching polygon and remove polygons
                    # from holes (if any)
                    else:
                        df.loc[df['osgb'] == t, 'AMD_polygon'] = False
                        df = remove_holes(df, t_polygon_within_hole)
        # Return the DataFrame
        return df

    # ***************** Body of the function *********************
    # Create the list of unique polygons (osgb)
    osgb_list = df['osgb'].unique().tolist()
    # Loop through the list of polygons and identify ones which needs
    # simplification
    for osgb in osgb_list:
        if df.loc[df['osgb'] == osgb, 'AMD_polygon_simplify'].values[0]:
            # Empty list to hold touching polygons (Note: Including
            # polygons which touches in the single point. That is the
            # reason for not taking polygons from the already existing list
            # of touching polygons)
            osgb_touching = []
            # Load the polygon sting (It can be False if polygon collapsed
            # due to previous simplification)
            polygon = df.loc[df['osgb'] == osgb, 'AMD_polygon'].values[0]
            if polygon:
                # Load the polygon as Shapely object
                osgb_polygon = loads(polygon)
                # Get the touching polygons
                osgb_touching = touching_poly(df, osgb, osgb_polygon,
                                              osgb_list, osgb_touching)
                # Conduct polygon simplification
                df = polygon_simplifying(osgb_polygon, tolerance, df,
                                         osgb, osgb_touching)
    # Return DataFrame
    return df

if __name__ == '__main__':
    main()
