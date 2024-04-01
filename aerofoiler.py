# import matplotlib.pyplot as plt
#
#
# # class Aerofoil:
# #     def init(self, name, chord, upper, lower):
# #         self.
#
# def bernstein(n, i, t):
#     return [binomial_coefficient(n, i) * (t_i ** i) * ((1 - t_i) ** (n - i)) for t_i in t]
#
#
# def binomial_coefficient(n, k):
#     if k < 0 or k > n:
#         return 0
#     result = 1
#     for i in range(min(k, n - k)):
#         result = result * (n - i) // (i + 1)
#     return result
#
#
# class Camber_Profile:
#     def __init__(self, p):
#         self.p = p
#         self.camberline = self.create_camberpoints()
#
#     def create_camberpoints(self, n=101):
#         t = [i / 100.0 for i in range(n)]
#         n = len(self.p) - 1
#         result = [(0.0, 0.0) for _ in range(len(t))]
#         for i in range(n + 1):
#             b = bernstein(n, i, t)
#             result = [[r[0] + self.p[i][0] * b_j, r[1] + self.p[i][1] * b_j] for r, b_j in
#                       zip(result, b)]
#         return result
#
#     def show_camberline(self):
#         fig, ax = plt.subplots(nrows=1, ncols=1)
#         ax.plot([i[0] for i in self.camberline], [i[1] for i in self.camberline], 'g-', label="camberline")
#         ax.plot([i[0] for i in self.p], [i[1] for i in self.p], 'ro-', markersize=10, markerfacecolor="None",
#                 label="Ctrl Points")
#         ax.grid()
#         ax.axis('equal')
#         plt.show()
#
#
# class Profile():
#     def __init__(self, camberline, ppoints):
#         self.ppoints = points
#         self.camberpoints = camberline
#         self.normals = self.get_normals()
#         self.controlpoints = self.get_controlpoints();
#
#     def get_normals(self):
#         # for i in range(len(self.ppoints)):
#         # get x,y, from s
#
#         # get tangent from s
#
#         # dot product to get unit vector normal
#
#         print("getting normals")
#         return
#
#     def get_controlpoints(self):
#         print("getting points")
#         return
#
#     def show_profile(self):
#         print("show profile")
#         # fig, ax = plt.subplots(nrows=1, ncols=1)
#         # ax.plot([i[0] for i in self.camberline], [i[1] for i in self.camberline], 'g-', label="camberline")
#         # ax.plot([i[0] for i in self.p], [i[1] for i in self.p], 'ro-',markersize=10, markerfacecolor="None", label="Ctrl Points")
#         # ax.grid()
#         # ax.axis('equal')
#         # plt.show()
#
#
# if __name__ == "__main__":
#     print("#### Aerofoil modder ######")
#     points = [[0, 0], [0.25, 0.25], [0.75, 0.1], [1, 0]]
#     a = Camber_Profile(points)
#     # a.show_camberline()
#     ppoints = [[0, 0.1], [0.5, 0.1], [1, 0.1], [45, 0.1]]
#     upper = Profile(a.camberline, ppoints)
#     upper.show_profile()

import matplotlib.pyplot as plt
from matplotlib import cycler
import os
import numpy as np
from scipy import interpolate
import math as m


# Read in aerofoil

# class for aerofoil
# scaled 0-1
# zero incidence
# camber line


class Aerofoil:
    def __init__(self, name, filename, chord=1, te_thickness=0.01, npoints=201):
        self.name = name
        self.chord = 1
        self.npoints = npoints
        self.raw_data = self.read_in_aerofoil(filename)
        self.points = self.raw_data
        self.te_thickness = te_thickness
        self.bspline = self.get_bspline_knots(self.raw_data)
        self.centre_s = self.get_le_s()
        self.s_dist = self.create_sdist()
        self.point = self.bspline_interp(self.bspline, self.s_dist)
        self.chord_knots = self.get_bspline_knots(self.get_chordline())
        self.chord_points = self.bspline_interp(self.chord_knots, np.linspace(0, 1, 100))
        self.set_scale_and_te(chord)

    def get_chord_length(self):
        xmin = min(self.points[:, 0])
        xmax = max(self.points[:, 0])
        return (xmax - xmin)

    def set_scale_and_te(self, chord):
        c_chord = self.get_chord_length()
        self.points = self.scale_aerofoil(chord / c_chord)
        self.bspline = self.get_bspline_knots(self.point)
        self.chord_knots = self.get_bspline_knots(self.get_chordline())
        self.chord_points = self.bspline_interp(self.chord_knots, np.linspace(0, 1, 100))
        nloops = 10
        for i in range(nloops):
            s_te, x_te, te_thick = self.find_te_thickness(self.te_thickness)
            self.bspline = self.get_bspline_knots(self.points)
            c_chord = x_te
            # print("Current chord is ",c_chord)
            # print("TE thickness is ",te_thick)
            # print("scaling factor is ",chord/c_chord)
            self.points = self.scale_aerofoil(chord / c_chord)
            self.bspline = self.get_bspline_knots(self.points)
            self.chord_knots = self.get_bspline_knots(self.get_chordline())
            self.chord_points = self.bspline_interp(self.chord_knots, np.linspace(0, 1, 100))
            if i == nloops - 1:
                s_te, x_te, te_thick = self.find_te_thickness(self.te_thickness, True)
            else:
                s_te, x_te, te_thick = self.find_te_thickness(self.te_thickness)
            # print("Current chord is ", c_chord)
            # print("TE thickness is ", te_thick)
            # print("scaling factor is ", chord / c_chord)
        self.points = self.truncate_points(s_te)
        self.bspline = self.get_bspline_knots(self.point)
        self.centre_s = self.get_le_s()
        self.s_dist = self.create_sdist(self.point.shape[0])
        self.points = self.bspline_interp(self.bspline, self.s_dist)
        self.chord_knots = self.get_bspline_knots(self.get_chordline())
        self.chord_points = self.bspline_interp(self.chord_knots, np.linspace(0, 1, 100))
        return

    def get_te_thickness(self):
        upper_y = self.point[-1, 1]
        lower_y = self.point[0, 1]
        return (upper_y - lower_y)

    def truncate_points(self, s_te):
        upper_x, upper_y = self.get_surface_point(self.bspline, s_te, self.centre_s, 1)
        lower_x, lower_y = self.get_surface_point(self.bspline, s_te, self.centre_s, -1)
        points = self.points
        i_dels = []
        for i in range(len(self.points)):
            if self.points[i, 0] > upper_x:
                i_dels.append(i)
        points = np.delete(self.points, i_dels, 0)
        # add te_points
        xx = np.asarray(lower_x)
        yy = np.asarray(lower_y)
        out = np.transpose((xx, yy))
        points = np.insert(points, 0, out, 0)
        xx = np.asarray(upper_x)
        yy = np.asarray(upper_y)
        out = np.transpose((xx, yy))
        points = np.insert(points, points.shape[0], out, 0)
        return points

    def find_te_thickness(self, thick, add_feature=False):
        s_test = 0.75
        i = 0
        tol = 0.0000001
        sstep = -0.05
        sdir = -1
        error = 10
        itmax = 100
        while abs(error) > tol:
            # s_test=np.linspace(0,1,10)
            upper_x, upper_y = self.get_surface_point(self.bspline, s_test, self.centre_s, 1)
            lower_x, lower_y = self.get_surface_point(self.bspline, s_test, self.centre_s, -1)
            error = abs(upper_y - lower_y) - thick
            # print(i, s_test, error, abs(upper_y-lower_y), thick)
            if error < 0 and abs(error) > tol:
                s_test = s_test + sstep
                if sdir == -1:
                    sstep = sstep / 2
                    sdir = 1
            elif error > 0 and abs(error) > tol:
                s_test = s_test - sstep
                if sdir == 1:
                    sstep = sstep / 2
                    sdir = -1
            if i == itmax:
                error = 0
            i += 1
        # print(s_test,upper_y,lower_y,abs(upper_y-lower_y), thick)
        # if add_feature:
        #     x = [upper_x, lower_x]
        #     y = [upper_y, lower_y]
        #     self.features.append(np.transpose((x, y)))
        return s_test, upper_x, abs(upper_y - lower_y)

    def get_surface_point(self, knots, s, centre_s, side):
        s_target = interpolate.splev(s, self.chord_knots)
        if side == 1:
            max_s = 1
            min_s = centre_s
            sstep = -0.05
            sdir = -1
        elif side == -1:
            max_s = centre_s
            min_s = 0
            sstep = 0.05
            sdir = 1
        else:
            max_s = 1
            min_s = 0
        s_test = (max_s + min_s) / 2
        itmax = 100
        error = 10
        i = 0
        tol = 0.0000001
        while abs(error) > tol:
            # s_test=np.linspace(0,1,10)
            out = interpolate.splev(s_test, self.bspline)
            error = s_target[0] - out[0]
            # print(i, s_test, error, out[0], s_target[0])
            if error < 0:
                s_test = s_test + sstep
                if sdir == -1:
                    sstep = sstep / 2
                    sdir = 1
            elif error > 0:
                s_test = s_test - sstep
                if sdir == 1:
                    sstep = sstep / 2
                    sdir = -1
            if i == itmax:
                error = 0
            i += 1
        # print(s,s_target[0],out[0],out[1])
        return out[0], out[1]

    def get_chordline(self):
        npoints = self.points.shape[0]
        x = []
        y = []
        # order .. LE to TE
        x.append(self.points[int((npoints - 1) / 2), 0])
        y.append(self.points[int((npoints - 1) / 2), 1])
        for i in range(int((npoints - 1) / 2) - 1, -1, -1):
            j = (npoints - 1) - i
            x.append((self.points[i, 0] + self.points[j, 0]) / 2)
            y.append((self.points[i, 1] + self.points[j, 1]) / 2)
        xx = np.asarray(x)
        yy = np.asarray(y)
        out = np.transpose((xx, yy))
        return out

    def scale_aerofoil(self, scale):
        out_points = self.points * scale
        return out_points

    def show_aerofoil(self):
        # plt.style.use('dark_background')
        # colors = cycler('color',
        #                 ['#EE6666', '#3388BB', '#9988DD',
        #                  '#EECC55', '#88BB44', '#FFBBBB'])
        # plt.rc('axes', facecolor='black', edgecolor='none',
        #        axisbelow=True, grid=True, prop_cycle=colors)
        # plt.rc('grid', color='gray', linestyle='solid')
        # plt.rc('xtick', direction='out', color='gray')
        # plt.rc('ytick', direction='out', color='gray')
        # plt.rc('patch', edgecolor='#E6E6E6')
        # plt.rc('lines', linewidth=2)

        fig, ax = plt.subplots(nrows=1, ncols=1)
        ax.plot(self.points[:, 0], self.points[:, 1], 'g-', label=self.name)
        ax.plot(self.chord_points[:, 0], self.chord_points[:, 1], 'r--', label=self.name)
        ax.grid()
        ax.axis('equal')
        title_string = self.name + ', npoints=' + str(self.points.shape[0])
        ax.set_title(title_string)
        plt.show()

    def print_aerofoil(self):
        print(self.name)
        print(self.points.shape[0])
        for i in range(self.points.shape[0]):
            print(f'{self.points[i][0]:.3f},{self.points[i][1]:.3f}')

    def read_in_aerofoil(self, filename):
        if os.path.isfile(filename):
            print("- raw data file found")
        else:
            print("- filename not found")
            return 0
        raw = np.loadtxt(filename, skiprows=1)
        npoints = len(raw)
        print("Reading aerofoil, npoints =", npoints)
        x = []
        y = []
        for i in range(npoints):
            x.append(raw[i][0])
            y.append(raw[i][1])
        raw_data = np.array(list(zip(x, y)))
        nn = int(raw_data.shape[1] / 2)
        ysum1 = 0
        ysum2 = 0
        for i in range(nn):
            ysum1 += raw_data[1][i]
        for i in np.arange(nn, raw_data.shape[1], 1):
            ysum2 += raw_data[1][i]
        if ysum1 > ysum2:
            print("re-ordering")
            raw_data = raw_data[::-1]
        if ysum1 < ysum2:
            print("ordering correct")

        return raw_data

    def get_bspline_knots(self, points):
        tck, u = interpolate.splprep((points[:, 0], points[:, 1]), s=0, k=3)
        return tck

    def bspline_interp(self, knots, u_vec):
        out = interpolate.splev(u_vec, knots)
        x = out[0]
        y = out[1]
        data = np.transpose((x, y))
        return data

    def get_le_s(self, tol=0.000000001):
        itmax = 10
        s_test = 0.5
        error = 10
        sstep = 0.001
        sdir = -1
        i = 0
        while abs(error) > tol:
            out = interpolate.splev(s_test, self.bspline)
            error = out[0] + out[1]
            if error < 0:
                s_test = s_test + sstep
                if sdir == -1:
                    sstep = sstep / 2
                    sdir = 1
            elif error > 0:
                s_test = s_test - sstep
                if sdir == 1:
                    sstep = sstep / 2
                    sdir = -1
            if i == itmax:
                error = 0
            i += 1
        # print(f"LE s located at {s_test:.3f}")
        return s_test

    def create_sdist(self, mids=0.1, tes=0.1, showplot=False):
        npoints = self.npoints
        uu = np.linspace(0, 1, num=npoints)
        plist = [(0, 0), (tes, tes / 2), (0.5 - mids, self.centre_s - 0.05), (0.5, self.centre_s),
                 (0.5 + mids, self.centre_s + 0.05),
                 (1 - tes, 1 - (tes / 2)), (1.0, 1.0)]
        ctr = np.array(plist)
        x = ctr[:, 0]
        y = ctr[:, 1]
        ll = len(x)
        t = np.linspace(0, 1, ll - 2, endpoint=True)
        t = np.append([0, 0, 0], t)
        t = np.append(t, [1, 1, 1])
        tck = [t, [x, y], 3]
        u = interpolate.splev(uu, tck)

        if showplot:
            zy = np.zeros(len(u[0]))
            plt.figure()
            plt.plot(u[0], u[1], 'ro-')
            plt.plot(x, y, 'bo--')
            plt.plot(zy, u[1], 'ko')
            plt.show()
        return u[1]


def calc_curvature(data):
    x = data[0]
    y = data[1]
    npoints = len(x)
    angle = []
    c = []
    s = []
    ss = []
    stotal = 0
    for i in range(npoints - 1):
        dy = y[i + 1] - y[i]
        dx = x[i + 1] - x[i]
        aa = np.degrees(m.atan2(dy, dx))
        if aa > 90:
            if aa < 180:
                aa = aa - 360
        angle.append(aa)
        s.append(m.sqrt(dx ** 2 + dy ** 2))
        stotal += s[i]
    for i in range(npoints - 2):
        c.append((angle[i + 1] - angle[i]))
    for j in range(len(s)):
        ss.append(s[j] / stotal)
    return c, angle, stotal


if __name__ == "__main__":
    a1 = Aerofoil('PW51_1', 'PW51.dat', chord=280, te_thickness=0, npoints=201)
    a1.show_aerofoil()
    a2 = Aerofoil('HT_12_1', 'ht12.dat', chord=280, te_thickness=0, npoints=201)
    a2.show_aerofoil()

# class for a baspline foil
# parameters
# camber line - LE anagle , TE angle, midpoint
# upper profile - using camberline normal, TE angle, TE thicknes
# lower profile - using camberline normal, TE angle, TE thickness

# Plot and interation of foil
#

# Optimiser to match control points of bspline foil to minimise error to aerofoil


# App
# custom tkinter
# Read in multiple aerofoils
# create a bspline foil
# turn bspline foilf in aerofoil
# save out aerofoil
