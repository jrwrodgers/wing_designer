import matplotlib.pyplot as plt
import math
import numpy as np
from Bezier import Bezier
from aerofoil_tools import *

def spine(x1 ,y1 ,x2 ,y2):
    # y=mx+c
    # ToDo consider bspline
    m = (y2 -y1 ) /(x2 -x1)
    c = y1 -( m *x1)
    return m ,c

def rotate_points(points,angle,axis_p):
    ox=axis_p[0]
    oy=axis_p[1]
    points_x=points[:,0]
    points_y=points[:,1]
    tt=math.radians(angle)
    points_xr=[ox+math.cos(tt)*(px-ox)-math.sin(tt)*(py-oy) for px,py in zip(points_x,points_y)]
    points_yr=[ox+math.sin(tt)*(px-ox)+math.cos(tt)*(py-oy) for px,py in zip(points_x,points_y)]
    points_r=[[x,y] for x,y in zip(points_xr,points_yr)]
    return np.array(points_r)


def get_bspline(s,p):
    point = Bezier.Curve([s], np.array(p))
    return point

class Wing():
    def __init__(self,c_chord,t_chord,hspan,t_chord_x,sections,twist,dihedral):
        self.c_chord=c_chord
        self.t_chord=t_chord
        self.hspan=hspan
        self.t_chord_x=t_chord_x
        self.nsecs=1
        self.sections=sections
        self.m_le ,self.c_le = spine(0 ,0  ,self.t_chord_x,self.hspan) #def spine(x1 ,y1 ,x2 ,y2)
        self.m_te, self.c_te = spine(self.c_chord,0 ,(self.t_chord_x +self.t_chord), self.hspan) #def spine(x1 ,y1 ,x2 ,y2)
        self.taper_ratio = self.t_chord / self.c_chord
        self.area = (self.hspan * self.c_chord * (1 + self.taper_ratio)) / 1000000 #mm span to meters area
        self.le_sweep = math.atan(self.t_chord_x / self.hspan)
        self.te_sweep = math.atan(((self.t_chord_x + self.t_chord) - self.c_chord) / self.hspan)
        self.wing_sweep = math.atan((self.t_chord_x + (0.25 * self.t_chord) - (self.c_chord * 0.25)) / self.hspan)
        self.twist=twist
        self.dihedral=dihedral

    def __str__(self):
        print(f"Wing Area = {self.area}")
        print(f"Taper Ratio = {self.taper_ratio}")
        print(f"Wing Sweep = {math.degrees(self.wing_sweep):.2f}")
        return "Done"


    def plot_wing(self):

        for s in range(len(self.sections)):
            print (f"Plotting section {s}")
            nj=20
            y_stations = [self.sections[s].ib_span_fraction*self.hspan + (j / (nj - 1)) * self.sections[s].ob_span_fraction*self.hspan for j in range(nj)]
            for j in range(len(y_stations)):
                c2_frac=j/(nj-1)
                c1_frac=1-(c2_frac)
                b_chord = np.add(self.sections[s].ib_section.points * c1_frac,self.sections[s].ob_section.points * c2_frac )
                #rotate
                angle=get_bspline(c2_frac, self.twist)
                print(angle[:,1])
                br_chord=rotate_points(b_chord,angle[:,1],[0,0])
                #dihedral


                le_offset=((c2_frac*self.hspan -self.c_le )/self.m_le)
                te_offset=((c2_frac*self.hspan -self.c_te )/self.m_te)
                scale=te_offset-le_offset
                br_chord = br_chord * scale
                br_chord[:, 0] = br_chord[:, 0] + le_offset

                #print j,c2_frac,angle)
                #dihedral
                plt.plot(br_chord[:, 0],br_chord[:, 1], '-')
                #plt.plot(br_chord[:, 0], br_chord[:, 1], 'r-')
                #plt.plot(le_offset, y_stations[j], te_offset,y_stations[j], 'g-')

        # nj =30
        # y_stations =[( j /(nj -1) ) *hspan for j in range(nj)]
        # # print(j_stations)
        # x_stations =[]
        #
        # for j in y_stations:
        #     le_x =m_le *j + c_le # or return from a spine function
        #     te_x =m_te *j + c_te # or return from a spine function
        #     chord =le_x -te_x
        #     x_stations.append([le_x -(chord * i /(ni -1)) for i in range(ni)])
        #
        # # print(x_stations)
        #
        # for j in range(nj -1):
        #     for i in range(ni -1):
        #         p1 = [y_stations[j] ,x_stations[j][i]]
        #         p2 = [y_stations[j] ,x_stations[j][ i +1]]
        #         p3 = [y_stations[j + 1], x_stations[j + 1][i + 1]]
        #         p4 = [y_stations[ j +1] ,x_stations[ j +1][i]]
        #         plt.plot([p1[0] ,p2[0]] ,[p1[1] ,p2[1]] ,'b', linestyle="-")
        #         plt.plot([p4[0], p1[0]], [p4[1], p1[1]], 'b', linestyle="-")
        #         # reflect
        #         plt.plot([-p1[0] ,-p2[0]] ,[p1[1] ,p2[1]] ,'b', linestyle="-")
        #         plt.plot([-p4[0], -p1[0]], [p4[1], p1[1]], 'b', linestyle="-")
        #         if i== ni - 2:
        #             plt.plot([p2[0], p3[0]], [p2[1], p3[1]], 'b', linestyle="-")
        #             # reflect
        #             plt.plot([-p2[0], -p3[0]], [p2[1], p3[1]], 'b', linestyle="-")
        #         if j == nj - 2:
        #             plt.plot([p3[0], p4[0]], [p3[1], p4[1]], 'b', linestyle="-")
        #             # reflect
        #             plt.plot([-p3[0], -p4[0]], [p3[1], p4[1]], 'b', linestyle="-")

        #ax = fig.add_subplot(111, projection='3d')

        plt.axis('equal')
        plt.grid()
        plt.show()


class Section():
    def __init__(self,ib_section,ob_section,ib_span_fraction,ob_span_fraction):
        self.ib_section=ib_section
        self.ob_section=ob_section
        self.ib_span_fraction=ib_span_fraction
        self.ob_span_fraction=ob_span_fraction


def main():
    print("Wing Designer")
    a1 = Aerofoil('PW51_1', 'PW51.dat', 1, 0.05, 51)
    a2 = Aerofoil('HT12_1', 'ht12.dat', 1, 0.05, 51)
    sections = []
    sections.append(Section(a1, a1, 0, 0.4))
    sections.append(Section(a1, a2,0.4,0.75))
    sections.append(Section(a2, a2, 0.75, 1.0))
    twist=[[0,0],[0.2,0],[0.8,10],[1,10]]
    dihedral = 0
    wing = Wing(400,164,1000,500,sections,twist,dihedral)
    wing.plot_wing()

if __name__ == "__main__":
    main()