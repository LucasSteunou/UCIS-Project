import h3
import math as m
import heapq
import matplotlib.pyplot as plt
import shapely
from matplotlib import patches
import geopandas as gpd
from shapely.geometry import Polygon, Point, LineString
import numpy as np
import networkx as nx
import pickle
import matplotlib.patheffects as path_effects
from path_error import PathError, WeightValueError
import time as t


def distance_h3(id_h3, end_h3):
    r = 6378137  # Rayon de la Terre en mètres
    dep = h3.h3_to_geo(id_h3)
    end = h3.h3_to_geo(end_h3)
    lat1, lon1 = m.radians(dep[0]), m.radians(dep[1])
    lat2, lon2 = m.radians(end[0]), m.radians(end[1])
    return r * m.acos(m.sin(lat1) * m.sin(lat2) + m.cos(lat1) * m.cos(lat2) * m.cos(lon1 - lon2))


def reconstruct_path(came_from, current_node):
    path = []
    longueur = 0
    while current_node:
        path.append(current_node)
        if came_from.get(current_node):
            longueur += distance_h3(current_node, came_from.get(current_node))
        current_node = came_from.get(current_node)
    return [path[::-1], longueur]


def neighbors(node):
    v = h3.k_ring(node)
    v.remove(node)
    return list(v)


def is_in_buffed_uas(neighbor, poly, res, margin):
    return is_in_polygon(h3.h3_to_geo(neighbor), poly) or distance_point_to_polygon(h3.h3_to_geo(neighbor),
                                                                                    poly) < margin + cote_hexa(res)


def cote_hexa(res):
    return m.sqrt(h3.hex_area(res, unit='m^2') * 2 / (3 * m.sqrt(3)))


def distance_point_to_polygon(point, polygon):
    point_geom = Point(point)
    polygon_geom = Polygon(polygon)
    return 6378137 * m.sin(2 * m.pi / 360 * point_geom.distance(polygon_geom))


def is_in_polygon(point_geo, polygon):
    polygon = Polygon(polygon)
    point_geo = Point(point_geo)
    return polygon.contains(point_geo)


def is_in_same_direction(dep, end, way, alpha):  # alpha = angle maxi pour être considéré comme dans le bon sens
    traj_array = np.array(h3.h3_to_geo(end)) - np.array(h3.h3_to_geo(dep))
    direction_array = np.array(way[1]) - np.array(way[0])
    return m.acos(np.dot(traj_array, direction_array) / (
            np.linalg.norm(traj_array) * np.linalg.norm(direction_array))) <= m.radians(alpha)


def is_in_corridor(hexa, corridors):
    is_in = False
    for corridor in corridors:
        if is_in_polygon(h3.h3_to_geo(hexa), corridor):
            is_in = True
            break
    return is_in


def segment_in_buffed_uas(dep_h3, end_h3, uas, margin):
    dep, end = h3.h3_to_geo(dep_h3), h3.h3_to_geo(end_h3)
    segment = LineString([dep, end])
    is_in = False
    for poly in uas:
        polygon = Polygon(poly)
        if polygon.crosses(segment) or 6378137 * m.sin(2 * m.pi / 360 * segment.distance(polygon)) < margin:
            is_in = True
            break
    return is_in


def trace_line(P1, P2):  # Plot uniquement
    n = 1000
    p1, p2 = h3.h3_to_geo(P1), h3.h3_to_geo(P2)
    if p1[1] == p2[1]:
        return [[p1[1] for _ in range(n + 1)], [p1[0] + (p2[0] - p1[0]) * i / n for i in range(n + 1)]]
    else:
        if p1[1] < p2[1]:
            p1, p2 = p2, p1
        a = (p1[0] - p2[0]) / (p1[1] - p2[1])
        X = [p1[1] + (p2[1] - p1[1]) * i / n for i in range(n + 1)]
        Y = [p2[0] - a * p2[1] + a * i for i in X]
        return [Y, X]


def invert_coordinates(UAS):
    return [[[coordinate[1], coordinate[0]] for coordinate in coordinates] for coordinates in UAS]


def invert_single_uas(uas):
    return [[coordinate[1], coordinate[0]] for coordinate in uas]


class Astar:
    def __init__(self, dep_geo, end_geo, res=9):
        Graphe = open('Visualisation/Graphe', 'rb')
        Posi = open('Visualisation/Posi', 'rb')
        self.graphe_toulouse = pickle.load(Graphe)
        self.position_toulouse = pickle.load(Posi)
        self.dep_geo = dep_geo
        self.end_geo = end_geo
        self.dep_h3 = h3.geo_to_h3(dep_geo[0], dep_geo[1], res)
        self.end_h3 = h3.geo_to_h3(end_geo[0], end_geo[1], res)
        self.res = res
        self.restricted_uas_margin = 30
        self.corridors_uas = []
        self.corridors_ways = []
        self.restricted_uas = []
        self.restricted_uas_names = []
        self.smoothed_path_length = None
        self.path_length = None
        self.length_in_corridor = None
        self.ideal_length = distance_h3(self.dep_h3, self.end_h3)
        self.corridor_weight = 0.8  # Poids attribué au corridor dans l'heuristique, entre 0 et 1, si >1 on tourne à
        # l'infini dans le corridor (plus opti = 0.95?)
        self.min_angle_in_corridor = 90
        self.smoothed_path = None
        self.astar_path = None
        self.max_number_of_studied_paths = 50  # Entier >1, Plus grand -> longueur chemin plus courte. Plus petit -> temps
        # de calcul réduit
        fig, ax = plt.subplots()
        self.ax = ax

    def set_margin(self, margin):
        self.restricted_uas_margin = margin

    def set_dep_end(self, dep_geo, end_geo):
        self.dep_h3 = h3.geo_to_h3(dep_geo)
        self.end_h3 = h3.geo_to_h3(end_geo)

    def set_corridor_weight(self, corridor_weight):
        if corridor_weight > 1 or corridor_weight < 0:
            raise WeightValueError
        else:
            self.corridor_weight = corridor_weight

    def set_min_angle_corridor(self, min_angle_corridor):
        self.min_angle_in_corridor = min_angle_corridor

    def set_max_number_of_studied_paths(self, max_number_of_studied_paths):
        self.max_number_of_studied_paths = max_number_of_studied_paths

    def get_smoothed_path(self):
        if not self.smoothed_path:
            self.compute_smoothed_path()
        return self.smoothed_path

    def get_smoothed_path_length(self):
        if not self.smoothed_path_length:
            self.compute_smoothed_path()
        return self.smoothed_path_length

    def get_astar_path(self):
        if not self.astar_path:
            self.compute_astar_path()
        return self.astar_path

    def get_dep_end(self, coord_type):
        if coord_type == 'h3':
            return [self.dep_h3, self.end_h3]
        elif coord_type == 'geo':
            return [self.dep_geo, self.end_geo]
        else:
            raise ValueError("coord_type must be 'h3' or 'geo'")

    def add_corridor_uas(self, poly, way):  # way = [coord_dep, coord_end], vecteur définissant le sens du flux
        self.corridors_uas.append(poly)
        self.corridors_ways.append(way)

    def add_restricted_uas(self, poly, name):
        self.restricted_uas.append(poly)
        self.restricted_uas_names.append(name)

    def compute_astar_path(self):
        open_list = []
        came_from = {}
        g_score = {self.dep_h3: 0}
        f_score = {self.dep_h3: distance_h3(self.dep_h3, self.end_h3)}
        heapq.heappush(open_list, (f_score[self.dep_h3], self.dep_h3))
        while open_list:
            _, current = heapq.heappop(open_list)
            open_list = open_list[:self.max_number_of_studied_paths]
            if current == self.end_h3:
                path = reconstruct_path(came_from, current)
                total_distance = path[1]
                self.path_length = total_distance
                return path[0]

            for neighbor in neighbors(current):
                is_in_uas = False
                for restricted_area in self.restricted_uas:
                    if is_in_buffed_uas(neighbor, restricted_area, self.res, self.restricted_uas_margin):
                        is_in_uas = True
                        break
                if is_in_uas:
                    continue
                wrong_way = False
                for corridor, way in zip(self.corridors_uas, self.corridors_ways):
                    if is_in_polygon(h3.h3_to_geo(neighbor), corridor) and not is_in_same_direction(current, neighbor,
                                                                                                    way,
                                                                                                    self.min_angle_in_corridor):
                        wrong_way = True
                        break
                if wrong_way:
                    continue

                coeff_corridor_current = 0
                if is_in_corridor(neighbor, self.corridors_uas):
                    coeff_corridor_current = self.corridor_weight

                tentative_g_score = g_score[current] + distance_h3(current, neighbor) * (1 - coeff_corridor_current)

                if neighbor not in g_score or tentative_g_score < g_score[neighbor]:
                    came_from[neighbor] = current
                    g_score[neighbor] = tentative_g_score
                    f_score[neighbor] = tentative_g_score + distance_h3(neighbor, self.end_h3)
                    heapq.heappush(open_list, (f_score[neighbor], neighbor))

        raise PathError

    def compute_smoothed_path(self):
        t_init = t.time()
        path = self.compute_astar_path()
        print('temps astar: ', t.time() - t_init)
        t_init = t.time()
        id_current, init = 0, path[0]
        length = 0
        length_in_corridor = 0
        opti_path = [self.dep_geo]
        while id_current != len(path) - 1:
            if is_in_corridor(path[id_current], self.corridors_uas):
                opti_path.append(h3.h3_to_geo(path[id_current]))
                length += distance_h3(init, path[id_current])
                init = path[id_current]
                while is_in_corridor(path[id_current], self.corridors_uas) and id_current < len(path) - 2:
                    id_current += 1
                opti_path.append(h3.h3_to_geo(path[id_current - 1]))
                opti_path.append(h3.h3_to_geo(path[id_current]))
                length_in_corridor += distance_h3(init, path[id_current])
                length += distance_h3(init, path[id_current])
                init = path[id_current]
            else:
                if segment_in_buffed_uas(init, path[id_current], self.restricted_uas, self.restricted_uas_margin):
                    opti_path.append(h3.h3_to_geo(path[id_current]))
                    length += distance_h3(init, path[id_current])
                    init = path[id_current]
                id_current += 1
        opti_path.append(h3.h3_to_geo(path[id_current]))
        length += distance_h3(init, path[id_current])
        self.length_in_corridor = length_in_corridor
        self.smoothed_path_length = length
        self.smoothed_path = opti_path
        print('temps lissage: ', t.time() - t_init)

    def compare_length(self):
        if not self.path_length or not self.smoothed_path_length:
            self.compute_astar_path()
            self.compute_smoothed_path()
        print(
            f'Longueur avec h3 : {self.path_length} mètres, soit {round((self.path_length - self.ideal_length) / self.ideal_length * 100, 2)}% de plus que le vol d\'oiseau')
        print(
            f'Longueur opti : {self.smoothed_path_length} mètres, soit {round((self.smoothed_path_length - self.ideal_length) / self.ideal_length * 100, 2)}% de plus que le vol d\'oiseau')
        print(
            f'Pourcentage du trajet dans un corridor : {round(self.length_in_corridor / self.smoothed_path_length * 100, 2)}%')

    def plot_smoothed_path(self):
        path = self.compute_astar_path()
        id_current, init = 0, path[0]
        while id_current != len(path) - 1:
            if is_in_corridor(path[id_current], self.corridors_uas):
                plt.plot(trace_line(init, path[id_current])[1], trace_line(init, path[id_current])[0], color='yellow',
                         linewidth=3)
                init = path[id_current]
                while is_in_corridor(path[id_current], self.corridors_uas) and id_current < len(path) - 2:
                    id_current += 1
                plt.plot(trace_line(init, path[id_current - 1])[1], trace_line(init, path[id_current - 1])[0],
                         color='yellow', linewidth=3)
                plt.plot(trace_line(path[id_current - 1], path[id_current])[1],
                         trace_line(path[id_current - 1], path[id_current])[0], color='yellow', linewidth=3)
                init = path[id_current]
            else:
                if segment_in_buffed_uas(init, path[id_current + 1], self.restricted_uas, self.restricted_uas_margin):
                    plt.plot(trace_line(init, path[id_current])[1], trace_line(init, path[id_current])[0],
                             color='yellow', linewidth=3)
                    init = path[id_current]
                id_current += 1
        plt.plot(trace_line(init, path[id_current])[1], trace_line(init, path[id_current])[0], color='yellow',
                 linewidth=3)

    def plot_toulouse(self):
        nx.draw_networkx_edges(self.graphe_toulouse, self.position_toulouse)

    def plot_uas(self):
        uas = invert_coordinates(self.restricted_uas)
        for polygon, name in zip(uas, self.restricted_uas_names):
            polygon_shape = patches.Polygon(polygon, closed=True, fill=True, edgecolor='red', facecolor='red',
                                            alpha=0.6)
            self.ax.add_patch(polygon_shape)
            centroid = Polygon(polygon).centroid
            text = self.ax.text(centroid.x, centroid.y, name, horizontalalignment='center', verticalalignment='center',
                                fontsize=10, color='white', fontweight='bold')
            text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

    def plot_corridors(self):
        corridors = invert_coordinates(self.corridors_uas)
        for corridor, way in zip(corridors, self.corridors_ways):
            polygon_shape = patches.Polygon(corridor, closed=True, fill=True, edgecolor='green', facecolor='green',
                                            alpha=0.6)
            self.ax.add_patch(polygon_shape)
            self.ax.quiver(way[0][1], way[0][0], way[1][1] - way[0][1], way[1][0] - way[0][0], color='black',
                           linewidth=4, angles='xy', scale=1, alpha=0.9)

    def plot_hexa_path(self):
        polygons = []
        for h3_index in self.compute_astar_path():
            hex_boundary = h3.h3_to_geo_boundary(h3_index, geo_json=True)
            polygons.append(Polygon(hex_boundary))
        gdf = gpd.GeoDataFrame({'geometry': polygons})
        gdf.plot(ax=self.ax, color='blue', edgecolor='black')
        self.ax.quiver(self.dep_geo[1], self.dep_geo[0], self.end_geo[1] - self.dep_geo[1],
                       self.end_geo[0] - self.dep_geo[0], color='green', linewidth=4, angles='xy', scale_units='xy',
                       scale=1)

    def show_all(self):
        self.plot_uas()
        self.plot_toulouse()
        self.plot_corridors()
        self.plot_hexa_path()
        self.plot_smoothed_path()
        plt.show()


if __name__ == '__main__':
    t0 = t.time()
    # path_finder = Astar([43.545221, 1.46423], [43.6328342937271, 1.3614776566658])
    path_finder = Astar([43.58, 1.5], [43.6328342937271, 1.3614776566658864])
    path_finder.add_restricted_uas([[43.6077, 1.4069], [43.5930, 1.3882], [43.5823, 1.4048], [43.5992, 1.4205]],
                                   'Aéroport')
    path_finder.add_restricted_uas([[43.5794, 1.4558], [43.5705, 1.4330], [43.5865, 1.4334]], 'Héliport')
    path_finder.add_corridor_uas([[43.5874, 1.4169], [43.6117, 1.4351], [43.6078, 1.4442], [43.5817, 1.427]],
                                 [[43.5874, 1.4169], [43.6117, 1.4351]])
    path_finder.add_corridor_uas([[43.6661, 1.3779], [43.6276, 1.4114], [43.6326, 1.4183], [43.6706, 1.3854]],
                                 [[43.6661, 1.3779], [43.6276, 1.4114]])
    path_finder.add_restricted_uas(
        invert_single_uas([[1.4285628927779612, 43.608731901218505], [1.4294792329866937, 43.60109830223689],
                           [1.4368072243693746, 43.59180840900919], [1.4532873331101825, 43.59645130226747],
                           [1.4523800960861308, 43.604412512263224], [1.440477511539541, 43.61204936076652], ]),
        'Hypercentre')
    path_finder.add_restricted_uas(
        invert_single_uas([[1.4028904968766653, 43.615375653647305], [1.419394474708895, 43.61404741799933],
                           [1.4203024649028748, 43.64628833508161]]), 'Parc')
    path_finder.add_restricted_uas(
        invert_single_uas([[1.46243695142482, 43.58982021236878], [1.4514521700985483, 43.579209093732175],
                           [1.4720357090032508, 43.56960063949242], [1.4747645768079565, 43.59512357439195]]),
        'Manifestation')
    path_finder.set_corridor_weight(0.999)
    path_finder.set_max_number_of_studied_paths(2000)
    path_finder.set_min_angle_corridor(55)
    print('init :', t.time() - t0)
    t0 = t.time()
    path_finder.get_smoothed_path()
    print('Temps à opti < 0.03 :', t.time() - t0)
    print("-" * 25 + "\n")
    path_finder.compare_length()




