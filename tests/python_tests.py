import vtk
import numpy as np

x_min = 0
x_max = 20
v_max = 20

def test_file_writer_output(path):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(path)
    reader.Update()

    pdata = reader.GetOutput()

    for i in range(pdata.GetNumberOfPoints()):
        v = 0
        if (pdata.GetPointData().GetArray('Boundary').GetValue(i) == 0):
            for j in range(2):
                assert (pdata.GetPoint(i)[j] >= x_min and pdata.GetPoint(i)[j] <= x_max)
                v += pdata.GetPointData().GetArray('Velocity').GetTuple(i)[j] ** 2
            assert np.sqrt(v) < v_max
            assert (pdata.GetPointData().GetArray('Pressure').GetValue(i) < 200000)
            assert (pdata.GetPointData().GetArray('Density').GetValue(i) > 100 and pdata.GetPointData().GetArray('Density').GetValue(i) < 2000)

for i in range(299):
    path = "./step_" + str(i) + ".vtp"
    test_file_writer_output(path)
    print ("step %d passed" % i)

    
