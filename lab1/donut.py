import gmsh
import sys

gmsh.initialize()
gmsh.model.add("donut")


lc = 0.04
R = 2.0
r = 1
W = 0.04

pi = 3.14

p_0 = gmsh.model.geo.addPoint(0, 0, 0, lc)
p_1 = gmsh.model.geo.addPoint(R - r, 0, 0, lc)
p_2 = gmsh.model.geo.addPoint(R, 0, r, lc)
p_3 = gmsh.model.geo.addPoint(R, 0, 0, lc)
p_4 = gmsh.model.geo.addPoint(R, 0, -r, lc)
p_5 = gmsh.model.geo.addPoint(R + r, 0, 0, lc)

arc_1 = gmsh.model.geo.addCircleArc(p_1, p_3, p_2)
arc_2 = gmsh.model.geo.addCircleArc(p_2, p_3, p_5)
arc_3 = gmsh.model.geo.addCircleArc(p_5, p_3, p_4)
arc_4 = gmsh.model.geo.addCircleArc(p_4, p_3, p_1)

c_loop_1 = gmsh.model.geo.addCurveLoop([arc_1, arc_2, arc_3, arc_4])
oval = gmsh.model.geo.copy([(1, arc_1), (1, arc_2), (1, arc_3), (1, arc_4)])
gmsh.model.geo.dilate(oval, R, 0, 0, 1 - W/r, 1 - W/r, 1 - W/r)

c_loop_2 = gmsh.model.geo.addCurveLoop([oval[0][1], oval[1][1], oval[2][1], oval[3][1]])
s1 = gmsh.model.geo.addPlaneSurface([c_loop_1, -c_loop_2])

gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(2)


rot1 = gmsh.model.geo.revolve([(2, s1)], 0, 0, 0, 0, 0, 1, pi /2)
rot2 = gmsh.model.geo.revolve([(2, rot1[0][1])], 0, 0, 0, 0, 0, 1, pi / 2)
rot3 = gmsh.model.geo.revolve([(2, rot2[0][1])], 0, 0, 0, 0, 0, 1, pi / 2)
rot4 = gmsh.model.geo.revolve([(2, rot3[0][1])], 0, 0, 0, 0, 0, 1, pi / 2)

gmsh.model.geo.addPhysicalGroup(3, [rot1[1][1], rot2[1][1], rot3[1][1], rot4[1][1]], 1)
gmsh.model.geo.synchronize()
gmsh.model.mesh.generate(3)

gmsh.write("lab1/donut.msh")


if '-nopopup' not in sys.argv:
    gmsh.fltk.run()
    

gmsh.finalize()