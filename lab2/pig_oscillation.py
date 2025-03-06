import gmsh
import math
import os
import sys
import numpy as np
import vtk

#ПОЛЯ НА СЕТКЕ

class CalcMesh:

    # Конструктор сетки, полученной из stl-файла
    def __init__(self, nodes_coords, tetrs_points):

        self.nodes = np.array([nodes_coords[0::3],nodes_coords[1::3],nodes_coords[2::3]])
        self.smth = np.power(self.nodes[0, :], 2) + np.power(self.nodes[1, :], 2)
        self.velocity = np.zeros(shape=(3, int(len(nodes_coords) / 3)), dtype=np.double)

        self.tetrs = np.array([tetrs_points[0::4],tetrs_points[1::4],tetrs_points[2::4],tetrs_points[3::4]])
        self.tetrs -= 1

        #Запишем начальное положение
        self.initial_nodes = 1 * self.nodes

    #Метод поступательного движения
    def move_translation(self, dx, dy, dz):
        self.nodes[0] += dx
        self.nodes[1] += dy
        self.nodes[2] += dz

    #Метод вращательного движения вокруг центра масс
    def move_rotation(self, phi_x, phi_y, phi_z):
        x_cm = np.mean(self.nodes[0])
        y_cm = np.mean(self.nodes[1])
        z_cm = np.mean(self.nodes[2])
        self.nodes[0] += phi_y * (self.nodes[2] - z_cm) - phi_z * (self.nodes[1] - y_cm)
        self.nodes[1] += phi_z * (self.nodes[0] - x_cm) - phi_x * (self.nodes[2] - z_cm)
        self.nodes[2] += phi_x * (self.nodes[1] - y_cm) - phi_y * (self.nodes[0] - x_cm)

    #Распространение поперечной волны в виде гауссова пакета вдоль оси x
    def move_gaussian_wave_packet(self, t, v_ph, width):
        f_t = 1 / np.exp((1 / 2 / width ** 2) * np.power(self.initial_nodes[0] - t * v_ph + v_ph * 40, 2))
        self.nodes[1] = self.initial_nodes[1] * (1 + f_t)
        self.nodes[2] = self.initial_nodes[2] * (1 + f_t)
    
    # Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    def move(self, tau, n):
        omega = 2 * np.pi / (40 * 0.01)
        self.velocity[2] += 4 * np.cos(0.05 * np.power(self.nodes[0] - self.nodes[1], 2)) * np.cos(omega * n * tau)
        self.nodes += self.velocity * tau

    # Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    def snapshot(self, snap_number):
        unstructuredGrid = vtk.vtkUnstructuredGrid()
        points = vtk.vtkPoints()

        smth = vtk.vtkDoubleArray()
        smth.SetName("smth")

        vel = vtk.vtkDoubleArray()
        vel.SetNumberOfComponents(3)
        vel.SetName("vel")

        for i in range(0, len(self.nodes[0])):
            points.InsertNextPoint(self.nodes[0,i], self.nodes[1,i], self.nodes[2,i])
            smth.InsertNextValue(self.smth[i])
            vel.InsertNextTuple((self.velocity[0,i], self.velocity[1,i], self.velocity[2,i]))

        unstructuredGrid.SetPoints(points)

        unstructuredGrid.GetPointData().AddArray(smth)
        unstructuredGrid.GetPointData().AddArray(vel)

        for i in range(0, len(self.tetrs[0])):
            tetr = vtk.vtkTetra()
            for j in range(0, 4):
                tetr.GetPointIds().SetId(j, self.tetrs[j,i])
            unstructuredGrid.InsertNextCell(tetr.GetCellType(), tetr.GetPointIds())

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetInputDataObject(unstructuredGrid)
        writer.SetFileName("lab2/pig oscillations/pig_oscillations-" + str(snap_number) + ".vtu")
        writer.Write()


#ПОСТРОЕНИЕ СЕТКИ

gmsh.initialize()

path = os.path.dirname(os.path.abspath(__file__))
gmsh.merge(os.path.join(path, '16456_CoinBank_Pig_V2.stl'))

gmsh.option.setNumber("Mesh.MeshSizeMax", 1)
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

gmsh.model.mesh.generate(3)

# РЕАЛИЗАЦИЯ ДВИЖЕНИЙ
nodeTags, nodesCoord, parametricCoord = gmsh.model.mesh.getNodes()

GMSH_TETR_CODE = 4
tetrsNodesTags = None
elementTypes, elementTags, elementNodeTags = gmsh.model.mesh.getElements()
for i in range(0, len(elementTypes)):
    if elementTypes[i] != GMSH_TETR_CODE:
        continue
    tetrsNodesTags = elementNodeTags[i]

if tetrsNodesTags is None:
    print("Can not find tetra data. Exiting.")
    gmsh.finalize()
    exit(-2)

print("The model has %d nodes and %d tetrs" % (len(nodeTags), len(tetrsNodesTags) / 4))

mesh = CalcMesh(nodesCoord, tetrsNodesTags)

mesh.snapshot(0)
for i in range(1, 120):
    mesh.move_gaussian_wave_packet(i, 0.2, 10)
    mesh.snapshot(i)


gmsh.finalize()