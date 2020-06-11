import vtk
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as anim

# Get info from vtp data structure
def vtp_to_dataframe(filepath):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(filepath)
    reader.Update()

    pdata = reader.GetOutput()

    p_x, p_y, density, velocity, boundary = [], [], [], [], []
    for i in range(pdata.GetNumberOfPoints()):
        p_x.append(pdata.GetPoint(i)[0])    # Extract x locations of points and append
        p_y.append(pdata.GetPoint(i)[1])    # Extract y locations of points and append
        v_x = pdata.GetPointData().GetArray('Velocity').GetTuple(i)[0]  # Extract x velocity locations of points
        v_y = pdata.GetPointData().GetArray('Velocity').GetTuple(i)[1]  # Extract y velocity locations of points
        # Add to list
        velocity.append(np.sqrt(v_x**2 + v_y**2))
        density.append(pdata.GetPointData().GetArray('Density').GetValue(i))
        boundary.append(pdata.GetPointData().GetArray('Boundary').GetValue(i))

    return pd.DataFrame({'p_x':p_x, 'p_y':p_y, 'density':density, 'velocity':velocity, 'boundary':boundary})

X, Y, step = [], [], []

# Get inventory length (number of images to generate)
num_items = len(os.listdir('./tests/'))

print("Generating density animations")
imgs_dens = []
fig = plt.figure()
for i in range(num_items):
    # Get pandas dataframe from vtp data structure
    A = vtp_to_dataframe("./tests/step_" + str(i) + ".vtp")
    # Isolate boundary
    boundary = A[(A.boundary == 1)]
    # Begin plotting data for this info
    plot_title = plt.title("Density animation for %i iterations" % i)
    plot_liquid_dens = plt.scatter(A['p_x'], A['p_y'], c=A['density'], s = 0.05)
    plot_boundaries = plt.scatter(boundary['p_x'], boundary['p_y'], c='red', s = 0.1)
    # Only consider points that aren't the boundary and are above a certain density
    A = A[(A.boundary == 0) & (A.density > 800)]
    # And subsequently find the highest point
    A = A[A['p_y'].isin([A['p_y'].max()])]
    # And append to max-point x and y coordinate list    
    X.append(A['p_x'].max())
    Y.append(A['p_y'].max())
    # Then plot it
    plot_highest = plt.scatter(X[-1], Y[-1], c="green", s=1)
    # Append list of paintings to master artist lists
    imgs_dens.append([plot_liquid_dens, plot_boundaries, plot_highest, plot_title])

print("Rendering mp4")
plt.close(fig)  # prevent final frame plot from showing up inline below
ani_dens = anim.ArtistAnimation(fig, imgs_dens, interval=50, blit=True)
ani_dens.save('wave_density.mp4')

print("Generating velocity animations")
imgs_vel = []
fig = plt.figure()
for i in range(num_items):
    # Get pandas dataframe from vtp data structure
    A = vtp_to_dataframe("./tests/step_" + str(i) + ".vtp")
    # Isolate boundary
    boundary = A[(A.boundary == 1)]
    # Begin plotting data for this info
    plot_title = plt.title("Velocity animation for %i iterations" % i)
    plot_liquid_vel = plt.scatter(A['p_x'], A['p_y'], c=A['velocity'], s = 0.05)
    plot_boundaries = plt.scatter(boundary['p_x'], boundary['p_y'], c='red', s = 0.1)
    # Only consider points that aren't the boundary and are above a certain density
    A = A[(A.boundary == 0) & (A.density > 800)]
    # And subsequently find the highest point
    A = A[A['p_y'].isin([A['p_y'].max()])]
    # And append to max-point x and y coordinate list     X.append(A['p_x'].max())
    Y.append(A['p_y'].max())
    # Then plot it
    plot_highest = plt.scatter(X[-1], Y[-1], c="green", s=1)
    # Append list of paintings to master artist lists
    imgs_vel.append([plot_liquid_vel, plot_boundaries, plot_highest, plot_title])


print("Rendering mp4")
plt.close(fig)  # prevent final frame plot from showing up inline below
ani_vel = anim.ArtistAnimation(fig, imgs_vel, interval=50, blit=True)
ani_vel.save('wave_velocity.mp4')




