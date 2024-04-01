import matplotlib.pyplot as plt
import math
import numpy as np

def spine(x1 ,y1 ,x2 ,y2):
    # y=mx+c
    # ToDo consider bspline
    m = (y2 -y1 ) /(x2 -x1)
    c = y1 -( m *x1)
    return m ,c


# Central AOA  0
# Tip AOA  -4

# Area 1      m2
# Area 2      m2
# Area 3      m2
#
# Section 1 Span   400
# Section 2 Span   750
# Section 1 AOA    0
# Section 2 AOA    -2.3
# Flap IB  100
# Flap OB  400
# Flap IB Chord    0.8
# Flap OB Chord    0.8
# Elevon IB Chord
# Elevon OB Chord
# Spar %Chord  0.3

print("Wing Designer")
# Input params
# ToDo make this a json read/write file
c_chord = 400
t_chord = 164
hspan = 1000
t_chord_y = 500

#Segments
##IB_Section
#OB_Section
#IB_Span_fraction
#OB_Span_fraction

# Stats
# For linear properties
# ToDo consider integration functions to get Area and MAC
# https://www.sciencedirect.com/topics/engineering/mean-aerodynamic-chord
taper_ratio =t_chord /c_chord
area= (hspan *c_chord *( 1 +taper_ratio) ) /1000000
le_sweep =math.atan(t_chord_y /hspan)
te_sweep =math.atan(((t_chord_y +t_chord ) -c_chord ) /hspan)
wing_sweep =math.atan((t_chord_y +(0.25 *t_chord ) -(c_chord *0.25) ) /hspan)
#


# Output properties
print(f"taper ration = {taper_ratio:.2f}")
print(f"area = {area:.2f} m2")
print(f"le_sweep = {math.degrees(le_sweep):.2f} deg")
print(f"te_sweep = {math.degrees(te_sweep):.2f} deg")
print(f"wing_sweep = {math.degrees(wing_sweep):.2f} deg")


# Create Profiles
# Make LE Profile
m_le ,c_le = spine(0 ,0 ,1000 ,-500)

# Make TE profile
m_te ,c_te = spine(0 ,-400 ,1000 ,(-500 -164))

# Make Spar Profile



# Make Chord Stations
#nj =30
#ni =10
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
# plt.axis('equal')
# plt.grid()
# plt.show()


