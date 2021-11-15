# Based on IK's 2_idf_geometry.py
# Takes basic settings file as an argument and produces IDF files, saved in directory for
# parallel processing on legion
# Does not include random assignment of ages and activities per IK script but uses data from
# preprocessing output files
import os
from time import strftime, localtime, time
import sys
import platform
from eppy.modeleditor import IDF
import pandas as pd
from shapely.wkt import loads
from shapely.ops import unary_union
from ast import literal_eval
from shapely.geometry import LineString, MultiLineString
import numpy as np
import math
from random import randrange, uniform


ROOT_DIR = os.path.abspath(os.path.dirname(__file__))
EP_DIR = os.path.join(ROOT_DIR, 'EnergyPlus')
IDF_DIR = os.path.join(ROOT_DIR, 'idf_files')
ep_basic_settings = os.path.join(ROOT_DIR, 'amd_basic_settings.idf')
input_data = os.path.join(ROOT_DIR, 'amd_all.CSV')
osgb3D = os.path.join(ROOT_DIR, '1_preprocessed_osgb_3D_blocks.csv')

# Do not place window if the wall width is less than this number
min_avail_width_for_window = 1
# Do not place window if partially exposed external wall is less than this number % of zone height
min_avail_height = 80
f2f = 3
# Include surrounding polygons within this radius [m]
radius = 30


def main():
    start = time()
   # print(strftime('%H:%M:%S', localtime()), '- idf development start time')
    # sys.stdout.flush()
    # Find the computer's operating system and set path to E+ idd file
    system = platform.system().lower()
    if system in ['windows', 'linux', 'darwin']:
        iddfile = os.path.join(EP_DIR, 'ep8.9_{}/Energy+.idd'.format(system))
    IDF.setiddname(iddfile)

    # Load input data (preprocessing outputs) and find unique built islands
    osgb3D_df = pd.read_csv(osgb3D)
    osgb3D_list = variable_polygon_pairs(osgb3D_df, 'osgb', 'AMD_polygon')
    df = pd.read_csv(input_data)
    # bi_list = df['osgb'].unique().tolist()
    # print(bi_list[:1])
    bi_list = df['bi'].unique().tolist()

    # Create EnergyPlus idf file for each built islands
    for bi in bi_list:
        pst = time()
        zone_use_dict = dict()
        # bi_df = df[df['osgb'] == bi]
        bi_df = df[df['bi'] == bi]
        idf = IDF(ep_basic_settings)

        # Change the name filed of the building object
        building_object = idf.idfobjects['BUILDING'][0]
        building_object.Name = bi

        # Move all objects towards origins
        origin = loads(bi_df['AMD_polygon'].iloc[0])
        origin = list(origin.exterior.coords[0])
        origin.append(0)

        # Built island polygon - convex and buffer
        osgb_in_bi = bi_df['osgb'].unique().tolist()
        polygons_in_bi = bi_df['AMD_polygon'].unique().tolist()
        polygons_in_bi = [loads(p) for p in polygons_in_bi]
        bi_polygon = unary_union(polygons_in_bi).convex_hull
        buffered_zone = bi_polygon.buffer(radius)
        # Surrounding context converted to shading elements
        osgbs_in_buffered_zone = list()
        for osgb, polygon in osgb3D_list:
            if osgb not in osgb_in_bi:
                if polygon.intersects(buffered_zone):
                    osgbs_in_buffered_zone.append(osgb)
        surrounding_df = osgb3D_df.loc[
            osgb3D_df['osgb'].isin(osgbs_in_buffered_zone)]
        surrounding_df.apply(shading_volumes,
                             args=(osgb3D_df, idf, origin,), axis=1)
        # Adiabatic volumes converted to shading elements
        """ shading_osgb = bi_df.loc[bi_df['AMD_adiabatic_volume_only'],
                                 'osgb'].unique().tolist()
        shading_df = osgb3D_df.loc[osgb3D_df['osgb'].isin(shading_osgb)]
        shading_df.apply(shading_volumes,
                         args=(osgb3D_df, idf, origin,), axis=1) """
        osgb_list = bi_df['osgb'].unique().tolist()

        for osgb in osgb_list:
            floor_use = bi_df.loc[bi_df['osgb'] == osgb][[
                'floor_no', 'use']].values.tolist()
            # print(non_dom_floors)
            #floor_use = [item for item in non_dom_floors if not np.isnan(item[0])]
            row = bi_df[bi_df['osgb'] == osgb].iloc[0]
            polygon = loads(row.AMD_polygon)
            hor_polygon = row.AMD_polygon_horizontal
            hor_poly_coord_dict = polygon_coordinates_dictionary(hor_polygon)
            horiz_surf_coord = horizontal_surface_coordinates(
                hor_poly_coord_dict, origin)
            ext_surf_polygon = loads(row.AMD_polygon_exposed_wall)
            ext_surf_coord = surface_coordinates(ext_surf_polygon, origin)
            adj_osgb_list = literal_eval(row.AMD_collinear_touching)
            # print(adj_osgb_list)
            glazing_ratio = row.AMD_glazing_ratio



            for floor_no, use in floor_use:  # non_dom_floors:
                if floor_no == 1:
                    if row.floors == 1:
                        zone_name = '{}_floor_{}'.format(row.osgb, floor_no)
                        zone_use_dict[zone_name] = use
                        zone_floor_h = (int(floor_no)-1) * f2f
                        space_below_floor = 'Ground'
                        zone_ceiling_h = int(floor_no) * f2f
                        space_above_floor = 'Outdoors'

                        idf.newidfobject('ZONE', Name=zone_name)

                        floor_const = row.AMD_ground_floor
                        floor(idf, zone_name, space_below_floor,
                              horiz_surf_coord, zone_floor_h, floor_const)
                        roof_const = row.AMD_roof
                        roof_ceiling(idf, zone_name, space_above_floor,
                                     horiz_surf_coord, zone_ceiling_h, roof_const)
                        glazing_const = row.AMD_glazing
                        zone_height = zone_ceiling_h - zone_floor_h
                        wall_const = row.AMD_wall
                        external_walls(
                            idf, zone_name, floor_no, ext_surf_coord,
                            zone_ceiling_h, zone_floor_h, zone_height,
                            min_avail_height, min_avail_width_for_window,
                            wall_const, glazing_const, glazing_ratio)
                    else:
                        zone_name = '{}_floor_{}'.format(row.osgb, floor_no)
                        zone_use_dict[zone_name] = use
                        zone_floor_h = (int(floor_no)-1) * f2f
                        space_below_floor = 'Ground'
                        zone_ceiling_h = int(floor_no) * f2f
                        space_above_floor = '{}_floor_{}'.format(
                            row.osgb, (floor_no + 1))

                        idf.newidfobject('ZONE', Name=zone_name)

                        floor_const = row.AMD_ground_floor
                        floor(idf, zone_name, space_below_floor,
                              horiz_surf_coord, zone_floor_h, floor_const)
                        roof_const = row.AMD_ceiling
                        roof_ceiling(idf, zone_name, space_above_floor,
                                     horiz_surf_coord, zone_ceiling_h,
                                     roof_const)
                        glazing_const = row.AMD_glazing
                        zone_height = zone_ceiling_h - zone_floor_h
                        wall_const = row.AMD_wall
                        external_walls(
                            idf, zone_name, floor_no, ext_surf_coord,
                            zone_ceiling_h, zone_floor_h, zone_height,
                            min_avail_height, min_avail_width_for_window,
                            wall_const, glazing_const, glazing_ratio)

                elif floor_no == row.floors:
                    zone_name = '{}_floor_{}'.format(row.osgb, floor_no)
                    zone_use_dict[zone_name] = use
                    zone_floor_h = (int(floor_no)-1) * f2f
                    space_below_floor = '{}_floor_{}'.format(
                        row.osgb, (floor_no - 1))
                    zone_ceiling_h = int(floor_no) * f2f
                    space_above_floor = 'Outdoors'

                    idf.newidfobject('ZONE', Name=zone_name)

                    floor_const = row.AMD_floor
                    floor(idf, zone_name, space_below_floor,
                          horiz_surf_coord, zone_floor_h, floor_const)
                    roof_const = row.AMD_roof
                    roof_ceiling(idf, zone_name, space_above_floor,
                                 horiz_surf_coord, zone_ceiling_h, roof_const)
                    glazing_const = row.AMD_glazing
                    zone_height = zone_ceiling_h - zone_floor_h
                    wall_const = row.AMD_wall
                    external_walls(
                        idf, zone_name, floor_no, ext_surf_coord,
                        zone_ceiling_h, zone_floor_h, zone_height,
                        min_avail_height, min_avail_width_for_window,
                        wall_const, glazing_const, glazing_ratio)

                else:
                    zone_name = '{}_floor_{}'.format(row.osgb, floor_no)
                    zone_use_dict[zone_name] = use
                    zone_floor_h = (int(floor_no)-1) * f2f
                    space_below_floor = '{}_floor_{}'.format(
                        row.osgb, (floor_no - 1))
                    zone_ceiling_h = floor_no * f2f
                    space_above_floor = '{}_floor_{}'.format(
                        row.osgb, (floor_no + 1))
                    idf.newidfobject('ZONE', Name=zone_name)

                    floor_const = row.AMD_floor
                    floor(idf, zone_name, space_below_floor,
                          horiz_surf_coord, zone_floor_h, floor_const)
                    roof_const = row.AMD_ceiling
                    roof_ceiling(idf, zone_name, space_above_floor,
                                 horiz_surf_coord, zone_ceiling_h,
                                 roof_const)
                    glazing_const = row.AMD_glazing
                    zone_height = zone_ceiling_h - zone_floor_h
                    wall_const = row.AMD_wall
                    external_walls(
                        idf, zone_name, floor_no, ext_surf_coord,
                        zone_ceiling_h, zone_floor_h, zone_height,
                        min_avail_height, min_avail_width_for_window,
                        wall_const, glazing_const, glazing_ratio)

        use_list = list(set(list(zone_use_dict.values())))
        for use in use_list:
            zone_list = list()
            for key, value in zone_use_dict.items():
                if value == use:
                    zone_list.append(key)
            idf.newidfobject('ZONELIST', Name=use)
            objects = idf.idfobjects['ZONELIST'][-1]
            for i, zone in enumerate(zone_list):
                exec('objects.Zone_%s_Name = zone' % (i + 1))  

        objects_to_delete = list()
        for obj in ['PEOPLE', 'LIGHTS', 'ELECTRICEQUIPMENT',
                    'ZONEINFILTRATION:DESIGNFLOWRATE',
                    'ZONEVENTILATION:DESIGNFLOWRATE',
                    'ZONECONTROL:THERMOSTAT']:
            objects = idf.idfobjects[obj]
            for item in objects:
                if item.Zone_or_ZoneList_Name not in use_list:
                    objects_to_delete.append(item)

        for item in objects_to_delete:
            idf.removeidfobject(item)

        all_zones = list(zone_use_dict.keys())
        for zone in all_zones:
            system_name = '{}_HVAC'.format(zone)
            eq_name = '{}_Eq'.format(zone)
            supp_air_node = '{}_supply'.format(zone)
            air_node = '{}_air_node'.format(zone)
            ret_air_node = '{}_return'.format(zone)

            idf.newidfobject('ZONEHVAC:IDEALLOADSAIRSYSTEM',
                             Name=system_name,
                             Zone_Supply_Air_Node_Name=supp_air_node,
                             Dehumidification_Control_Type='None')

            idf.newidfobject('ZONEHVAC:EQUIPMENTLIST',
                             Name=eq_name,
                             Zone_Equipment_1_Object_Type='ZONEHVAC:IDEALLOADSAIRSYSTEM',
                             Zone_Equipment_1_Name=system_name,
                             Zone_Equipment_1_Cooling_Sequence=1,
                             Zone_Equipment_1_Heating_or_NoLoad_Sequence=1)

            idf.newidfobject('ZONEHVAC:EQUIPMENTCONNECTIONS',
                             Zone_Name=zone,
                             Zone_Conditioning_Equipment_List_Name=eq_name,
                             Zone_Air_Inlet_Node_or_NodeList_Name=supp_air_node,
                             Zone_Air_Node_Name=air_node,
                             Zone_Return_Air_Node_or_NodeList_Name=ret_air_node)

        ep_bs = str(ep_basic_settings).strip("\basic_settings\.idf")
        idf.saveas(os.path.join(IDF_DIR, '{}.idf'.format(bi)))
        pt('{}.idf file created in:'.format(bi), pst)

    print(strftime('%H:%M:%S', localtime()), '- idf development finish time')
    sys.stdout.flush()
    pt('##### idf development completed in:'.format(bi), start)

# END OF MAIN  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def pt(printout, pst):
    pft = time()
    process_time = pft - pst
    if process_time <= 60:
        unit = 'sec'
    elif process_time <= 3600:
        process_time = process_time / 60
        unit = 'min'
    else:
        process_time = process_time / 3600
        unit = 'hr'
    loctime = strftime('%H:%M:%S', localtime())
    print('{0} - {1} {2:.2f} {3}'.format(loctime, printout, process_time, unit))
    sys.stdout.flush()
    return

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def variable_polygon_pairs(df, variable_value, polygon_value):
    variable_polygon_values = df[[variable_value, polygon_value]].values
    variable_polygon = list()
    for variable, polygon in variable_polygon_values:
        polygon = loads(polygon)
        variable_polygon.append([variable, polygon])
    return variable_polygon

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def shading_volumes(row, df, idf, origin):
    osgb, polygon = row.osgb, loads(row.AMD_polygon)
    hor_polygon = row.AMD_polygon_horizontal
    hor_poly_coord_dict = polygon_coordinates_dictionary(hor_polygon)
    adj_osgb_list = literal_eval(str(row.AMD_collinear_touching))
    ext_surf_polygon = loads(row.AMD_polygon_exposed_wall)
    ext_surf_coord = surface_coordinates(ext_surf_polygon, origin)
    horiz_surf_coord = horizontal_surface_coordinates(
        hor_poly_coord_dict, origin)
    zone_floor_h = 0
    zone_ceiling_h = row.height
    adiabatic_roof(idf, osgb, horiz_surf_coord, zone_ceiling_h)
    adiabatic_wall_name = 'AdiabaticWall'
    adiabatic_external_walls(idf, osgb, ext_surf_coord, zone_ceiling_h,
                             zone_floor_h, adiabatic_wall_name,
                             adj_osgb_list, df, polygon, origin)
    return

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def surface_coordinates(polygon, origin):
    coordinates_list = list()
    if polygon.geom_type in ['MultiLineString', 'GeometryCollection']:
        for item in polygon:
            if not item.geom_type == 'Point':
                coordinates_list.append(item.coords)
    elif polygon.geom_type == 'LineString':
        coordinates_list.append(polygon.coords)
    coordinates = coordinates_move_origin(coordinates_list, origin)
    return coordinates

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def coordinates_move_origin(coordinates_list, origin):
    coordinates_list_moved_origin = list()
    for coordinates in coordinates_list:
        coordinates_moved_origin = list()
        for ordered_pair in coordinates:
            ordered_pair_moved_origin = [i - j for i, j in zip(ordered_pair,
                                                               origin)]
            ordered_pair_moved_origin = [
                round(coord, 2) for coord in ordered_pair_moved_origin]
            coordinates_moved_origin.append(ordered_pair_moved_origin)
        coordinates_list_moved_origin.append(coordinates_moved_origin)
    return coordinates_list_moved_origin

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def adiabatic_external_walls(idf, polygon_name, perimeter_surface_coordinates,
                             ceiling_height, floor_height, wall_name,
                             adjacent_polygons_list, df,
                             polygon_shapely, origin):

    def adiabatic_walls(idf, polygon_name, perimeter_surface_coordinates,
                        ceiling_height, floor_height, wall_name):
        ceiling_coordinates = coordinates_add_height(
            ceiling_height, perimeter_surface_coordinates)
        floor_coordinates = coordinates_add_height(
            floor_height, perimeter_surface_coordinates)
        for n, item in enumerate(ceiling_coordinates):
            ceil_coord = ceiling_coordinates[n]
            floor_coord = floor_coordinates[n]
            for i, p in enumerate(ceil_coord[:-1]):
                wcc = wall_centre_coordinate(
                    ceil_coord[i + 1], ceil_coord[i], floor_coord[i])
                surface_name = polygon_name + '_' + wall_name + '_' + wcc
                coordinates = idf_wall_coordinates(i, ceil_coord, floor_coord)
                shading_building_detailed(idf, surface_name, coordinates)
        return

    adiabatic_walls(idf, polygon_name, perimeter_surface_coordinates,
                    ceiling_height, floor_height, wall_name)
    if adjacent_polygons_list:
        for polygon in adjacent_polygons_list:
            adjacent_polygon_df = df.loc[df['osgb'] == polygon]
            adjacent_polygon = adjacent_polygon_df['AMD_polygon'].iloc[0]
            adjacent_polygon = loads(adjacent_polygon)
            part_wall_polygon = polygon_shapely.intersection(adjacent_polygon)
            ajd_wall_parti_surf_coord = surface_coordinates(part_wall_polygon,
                                                            origin)
            adjacent_height = adjacent_polygon_df['height'].iloc[0]
            if ceiling_height > adjacent_height:
                ceil_h = ceiling_height
                if floor_height > adjacent_height:
                    floor_h = floor_height
                else:
                    floor_h = adjacent_height
                adiabatic_walls(idf, polygon_name, ajd_wall_parti_surf_coord,
                                ceil_h, floor_h, wall_name)
    return

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def adiabatic_roof(idf, polygon_name, horizontal_surface_coordinates,
                   ceiling_height):
    ceiling_coordinates_list = coordinates_add_height(
        ceiling_height, horizontal_surface_coordinates)
    coordinates_list = idf_ceiling_coordinates_list(ceiling_coordinates_list)
    surface_name = polygon_name + '_AdiabaticRoof'
    for coordinates in coordinates_list:
        shading_building_detailed(idf, surface_name, coordinates)
    return

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def coordinates_add_height(height, coordinates_list):
    coordinates_with_height = []
    for coordinates in coordinates_list:
        ordered_pair_with_height = []
        for op in coordinates:
            op_with_height = op + [round(height, 2)]
            ordered_pair_with_height.append(op_with_height)
        coordinates_with_height.append(ordered_pair_with_height)
    return coordinates_with_height

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def wall_centre_coordinate(ceil_1, ceil_0, floor_0):
    hor = [ceil_1, ceil_0]
    ver = [ceil_0, floor_0]
    hor = [item[:-1] for item in hor]
    ver = [item[1:] for item in ver]
    hor = LineString(hor).centroid
    hor = hor.coords[0]
    ver = LineString(ver).centroid
    ver = ver.coords[0]
    wcc = hor + (ver[-1],)
    wcc = ['%.2f' % item for item in wcc]
    wcc = '(' + '_'.join(wcc) + ')'
    return wcc

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def idf_wall_coordinates(i, ceiling_coordinates, floor_coordinates):
    idf_wall_coordinates = [ceiling_coordinates[i + 1],
                            floor_coordinates[i + 1],
                            floor_coordinates[i],
                            ceiling_coordinates[i]]
    return idf_wall_coordinates

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def shading_building_detailed(idf, surface_name, coordinates):
    idf.newidfobject(
        'Shading:Building:Detailed'.upper(),
        Name=surface_name)
    objects = idf.idfobjects['Shading:Building:Detailed'.upper()][-1]
    for i, ordered_pair in enumerate(coordinates):
        exec('objects.Vertex_{}_Xcoordinate = ordered_pair[0]'.format(i + 1))
        exec('objects.Vertex_{}_Ycoordinate = ordered_pair[1]'.format(i + 1))
        exec('objects.Vertex_{}_Zcoordinate = ordered_pair[2]'.format(i + 1))
    return

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def polygon_coordinates_dictionary(polygon):
    polygon = loads(polygon)
    polygon_coordinates_dict = dict()
    polygon_coordinates_dict['outer_ring'] = [polygon.exterior.coords]
    if polygon.interiors:
        polygon_coordinates_dict['inner_rings'] = list()
        for item in polygon.interiors:
            polygon_coordinates_dict['inner_rings'].append(item.coords)
    return polygon_coordinates_dict

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def horizontal_surface_coordinates(coordinates_dictionary, origin):
    if len(coordinates_dictionary) == 1:
        coordinates_list = coordinates_dictionary['outer_ring']
    else:
        coordinates_list = polygon_with_holes(coordinates_dictionary)
    coordinates = coordinates_move_origin(coordinates_list, origin)
    return coordinates

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def polygon_with_holes(coordinates_dictionary):

    def polygon_with_holes_coordinates(coordinates, outer_ring_linestring,
                                       interceptors_op_dict):
        first_op = outer_ring_linestring.coords[0]
        if first_op in interceptors_op_dict:
            coordinates.append(first_op)
            holes = interceptors_op_dict[first_op]
            for hole in holes:
                if hole.geom_type == 'LineString':
                    for coord in hole.coords:
                        coordinates.append(coord)
                if hole.geom_type == 'MultiLineString':
                    for coord in hole[1].coords:
                        coordinates.append(coord)
                    for coord in hole[0].coords[1:]:
                        coordinates.append(coord)
                coordinates.append(first_op)
            for coord in outer_ring_linestring.coords[1:-1]:
                coordinates.append(coord)
        else:
            for coord in outer_ring_linestring.coords[:-1]:
                coordinates.append(coord)
        return coordinates

    def dist_two_points(p1, p2):
        distance = math.sqrt(math.pow(
            (p1[0] - p2[0]), 2) + math.pow((p1[1] - p2[1]), 2))
        return distance

    def inner_string(inner_ring_coordinates, irop):
        inner_coordinates = list(inner_ring_coordinates)

        for i, coord in enumerate(inner_coordinates):
            if coord == irop:
                split_position = i
                break
        if split_position == 0:
            inner_linestring = LineString(inner_coordinates)
        else:
            first = inner_coordinates[:(split_position + 1)]
            last = inner_coordinates[split_position:]
            inner_linestring = MultiLineString([first, last])
        return inner_linestring

    outer_ring_coordinates = coordinates_dictionary['outer_ring']
    interceptors_op_dict = dict()
    oi_min_linestring_list = list()
    for inner_ring_coordinates in coordinates_dictionary['inner_rings']:
        for i, orop in enumerate(outer_ring_coordinates[0][:-1]):
            for j, irop in enumerate(inner_ring_coordinates[:-1]):
                if i == 0 and j == 0:
                    oi_min_distance = dist_two_points(orop, irop)
                    if oi_min_distance > 0.015:
                        oi_min_linestring = LineString([orop, irop])
                        inner_linestring = inner_string(inner_ring_coordinates,
                                                        irop)
                    else:
                        oi_min_distance = 1e9
                else:
                    distance = dist_two_points(orop, irop)
                    if (distance < oi_min_distance) and (distance > 0.015):
                        oi_min_distance = distance
                        oi_min_linestring = LineString([orop, irop])
                        inner_linestring = inner_string(inner_ring_coordinates,
                                                        irop)

        if oi_min_linestring.coords[0] in interceptors_op_dict:
            interceptors_op_dict[oi_min_linestring.coords[0]].append(
                inner_linestring)
        else:
            interceptors_op_dict[
                oi_min_linestring.coords[0]] = [inner_linestring]
        oi_min_linestring_list.append(oi_min_linestring)
    outer_ring = MultiLineString(coordinates_dictionary['outer_ring'])
    for oi_min_linestring in oi_min_linestring_list:
        outer_ring = outer_ring.difference(oi_min_linestring)
    coordinates = list()
    if outer_ring.geom_type == 'LineString':
        start_end_op = outer_ring.coords[0]
        coordinates = polygon_with_holes_coordinates(coordinates, outer_ring,
                                                     interceptors_op_dict)
    elif outer_ring.geom_type == 'MultiLineString':
        for i, linestring in enumerate(outer_ring):
            if i == 0:
                start_end_op = linestring.coords[0]
            coordinates = polygon_with_holes_coordinates(coordinates,
                                                         linestring,
                                                         interceptors_op_dict)
    coordinates.append(start_end_op)
    return [coordinates]

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def floor(idf, zone_name, space_below_floor, horizontal_surface_coordinates,
          floor_height, ground_floor_const):

    def idf_floor_coordinates_list(floor_coordinates_list):
        idf_floor_coordinates_list = list()
        for floor_coordinates in floor_coordinates_list:
            idf_floor_coordinates = floor_coordinates[:-1]
            idf_floor_coordinates_list.append(idf_floor_coordinates)
        return idf_floor_coordinates_list

    surface_name = zone_name + '_Floor'
    surface_type = 'Floor'
    sun_exposure = 'NoSun'
    wind_exposure = 'NoWind'
    if space_below_floor == 'Ground':
        outside_boundary_condition = space_below_floor
        outside_boundary_condition_object = ''
    else:
        outside_boundary_condition = 'Surface'
        outside_boundary_condition_object = space_below_floor + '_Ceiling'
    floor_coordinates_list = coordinates_add_height(
        floor_height, horizontal_surface_coordinates)
    coordinates_list = idf_floor_coordinates_list(floor_coordinates_list)
    for coordinates in coordinates_list:
        building_surface_detailed(idf, surface_name, surface_type,
                                  ground_floor_const, zone_name,
                                  outside_boundary_condition,
                                  outside_boundary_condition_object,
                                  sun_exposure, wind_exposure, coordinates)
    return

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def building_surface_detailed(idf, surface_name, surface_type,
                              construction_name, zone_name,
                              outside_boundary_condition,
                              outside_boundary_condition_object,
                              sun_exposure, wind_exposure, coordinates):
    idf.newidfobject(
        'BuildingSurface:Detailed'.upper(),
        Name=surface_name,
        Surface_Type=surface_type,
        Construction_Name=construction_name,
        Zone_Name=zone_name,
        Outside_Boundary_Condition=outside_boundary_condition,
        Outside_Boundary_Condition_Object=outside_boundary_condition_object,
        Sun_Exposure=sun_exposure,
        Wind_Exposure=wind_exposure)
    objects = idf.idfobjects['BuildingSurface:Detailed'.upper()][-1]
    for i, ordered_pair in enumerate(coordinates):
        exec('objects.Vertex_%s_Xcoordinate = ordered_pair[0]' % (i + 1))
        exec('objects.Vertex_%s_Ycoordinate = ordered_pair[1]' % (i + 1))
        exec('objects.Vertex_%s_Zcoordinate = ordered_pair[2]' % (i + 1))
    return

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def roof_ceiling(idf, zone_name, space_above_floor,
                 horizontal_surface_coordinates, ceiling_height, roof_const):
    if space_above_floor == 'Outdoors':
        surface_name = zone_name + '_Roof'
        surface_type = 'Roof'
        sun_exposure = 'SunExposed'
        wind_exposure = 'WindExposed'
        outside_boundary_condition = space_above_floor
        outside_boundary_condition_object = ''
    else:
        surface_name = zone_name + '_Ceiling'
        surface_type = 'Ceiling'
        sun_exposure = 'NoSun'
        wind_exposure = 'NoWind'
        outside_boundary_condition = 'Surface'
        outside_boundary_condition_object = space_above_floor + '_Floor'
    ceililng_coordinates_list = coordinates_add_height(
        ceiling_height, horizontal_surface_coordinates)
    coordinates_list = idf_ceiling_coordinates_list(ceililng_coordinates_list)
    for coordinates in coordinates_list:
        building_surface_detailed(idf, surface_name, surface_type,
                                  roof_const, zone_name,
                                  outside_boundary_condition,
                                  outside_boundary_condition_object,
                                  sun_exposure, wind_exposure, coordinates)
    return

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def idf_ceiling_coordinates_list(ceiling_coordinates_list):
    idf_ceiling_coordinates_list = []
    for ceiling_coordinates in ceiling_coordinates_list:
        ceiling_coordinates = ceiling_coordinates[:-1]
        ceiling_coordinates = list(reversed(ceiling_coordinates))
        idf_ceiling_coordinates_list.append(ceiling_coordinates)
    return idf_ceiling_coordinates_list

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def external_walls(idf, zone_name, floor_number,
                   vertical_surface_coordinates, ceiling_height,
                   floor_height, zone_height, min_avail_height,
                   min_window_width, wall_const, glazing_const, glazing_ratio):

    def wall_width_height(i, ceil_coord, floor_coord):
        ulc = ceil_coord[i + 1]
        blc = floor_coord[i + 1]
        brc = floor_coord[i]
        w = math.sqrt(math.pow(
            (brc[0] - blc[0]), 2) + math.pow((brc[1] - blc[1]), 2) +
            math.pow((brc[2] - blc[2]), 2))
        h = math.sqrt(math.pow(
            (ulc[0] - blc[0]), 2) + math.pow((ulc[1] - blc[1]), 2) +
            math.pow((ulc[2] - blc[2]), 2))
        return w, h

    def window(idf, surface_name, construction_name, building_surface_name,
               starting_x_coordinate, starting_z_coordinate, length, height):
        idf.newidfobject(
            'Window'.upper(),
            Name=surface_name,
            Construction_Name=construction_name,
            Building_Surface_Name=building_surface_name,
            Starting_X_Coordinate=starting_x_coordinate,
            Starting_Z_Coordinate=starting_z_coordinate,
            Length=length,
            Height=height)
        return

    surface_type = 'Wall'
    outside_boundary_condition_object = ''
    sun_exposure = 'SunExposed'
    wind_exposure = 'WindExposed'
    outside_boundary_condition = 'Outdoors'
    ceiling_coordinates = coordinates_add_height(ceiling_height,
                                                 vertical_surface_coordinates)
    floor_coordinates = coordinates_add_height(floor_height,
                                               vertical_surface_coordinates)
    for n, _ in enumerate(ceiling_coordinates):
        ceil_coord = ceiling_coordinates[n]
        floor_coord = floor_coordinates[n]
        for i, _ in enumerate(ceil_coord[:-1]):
            wcc = wall_centre_coordinate(
                ceil_coord[i + 1], ceil_coord[i], floor_coord[i])
            surface_name = zone_name + '_Wall_' + wcc
            coordinates = idf_wall_coordinates(i, ceil_coord, floor_coord)
            building_surface_detailed(idf, surface_name, surface_type,
                                      wall_const, zone_name,
                                      outside_boundary_condition,
                                      outside_boundary_condition_object,
                                      sun_exposure, wind_exposure, coordinates)
            w, h = wall_width_height(i, ceil_coord, floor_coord)
            if (w >= min_window_width) and (
                    h >= (min_avail_height * zone_height / 100) - 0.001):
                wl = w * math.sqrt(glazing_ratio / 100)
                wh = h * math.sqrt(glazing_ratio / 100)
                x = (w - wl) / 2
                z = (h - wh) / 2
                win_surface_name = surface_name + '_Window'
                building_surface_name = surface_name
                starting_x_coordinate = '%.2f' % x
                starting_z_coordinate = '%.2f' % z
                win_length = '%.2f' % wl
                win_height = '%.2f' % wh
                window(idf, win_surface_name, glazing_const,
                       building_surface_name, starting_x_coordinate,
                       starting_z_coordinate, win_length, win_height)
    return

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def partition_walls(idf, zone_name, adj_osgb, vertical_surface_coordinates,
                    ceiling_height, floor_height, partition_const):
    surface_type = 'Wall'
    sun_exposure = 'NoSun'
    wind_exposure = 'NoWind'
    opposite_zone = adj_osgb
    outside_boundary_condition = 'Adiabatic'
    obco = ''
    ceiling_coordinates = coordinates_add_height(ceiling_height,
                                                 vertical_surface_coordinates)
    floor_coordinates = coordinates_add_height(floor_height,
                                               vertical_surface_coordinates)
    for n, item in enumerate(ceiling_coordinates):
        ceil_coord = ceiling_coordinates[n]
        floor_coord = floor_coordinates[n]
        for i, p in enumerate(ceil_coord[:-1]):
            wcc = wall_centre_coordinate(
                ceil_coord[i + 1], ceil_coord[i], floor_coord[i])
            surface_name = zone_name + '_Part_' + opposite_zone + '_' + wcc
            coordinates = idf_wall_coordinates(
                i, ceil_coord, floor_coord)
            building_surface_detailed(idf, surface_name, surface_type,
                                      partition_const, zone_name,
                                      outside_boundary_condition, obco,
                                      sun_exposure, wind_exposure, coordinates)
    return

# END OF FUNCTION  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


if __name__ == '__main__':
    main()
