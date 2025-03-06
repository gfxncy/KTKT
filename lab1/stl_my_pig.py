import gmsh
import math
import os
import sys

gmsh.initialize()


path = os.path.dirname(os.path.abspath(__file__))
gmsh.merge(os.path.join(path, '16456_CoinBank_Pig_V2.stl'))

gmsh.option.setNumber("Mesh.MeshSizeMax", 0.5)
angle = 40


forceParametrizablePatches = True

includeBoundary = True

curveAngle = 180

gmsh.model.mesh.classifySurfaces(angle * math.pi / 180., includeBoundary,
                                 forceParametrizablePatches,
                                 curveAngle * math.pi / 180.)


gmsh.model.mesh.createGeometry()

s = gmsh.model.getEntities(2)
l = gmsh.model.geo.addSurfaceLoop([s[i][1] for i in range(len(s))])
gmsh.model.geo.addVolume([l])

gmsh.model.geo.synchronize()

f = gmsh.model.mesh.field.add("MathEval")
gmsh.model.mesh.field.setString(f, "F", "4")
gmsh.model.mesh.field.setAsBackgroundMesh(f)

gmsh.model.mesh.generate(4)
gmsh.write('lab1/pig.msh')

# Launch the GUI to see the results:
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()