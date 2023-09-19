#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import vtk

colors = vtk.vtkNamedColors()

filename = 'MASK (copy).vtp'

    # Read all the data from the file
reader = vtk.vtkXMLPolyDataReader()
reader.SetFileName(filename)
reader.Update()

    # Visualize
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(reader.GetOutputPort())

actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetColor(colors.GetColor3d('Yellow'))

renderer = vtk.vtkRenderer()
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

renderer.AddActor(actor)
renderer.SetBackground(colors.GetColor3d('Gray'))
renderer.GetActiveCamera().Pitch(90)
renderer.GetActiveCamera().SetViewUp(0, 0, 1)
renderer.ResetCamera()

renderWindow.SetSize(600, 600)
renderWindow.Render()
renderWindow.SetWindowName('ReadPolyData')
renderWindowInteractor.Start()


# In[2]:


import sys
import vtk

colors = vtk.vtkNamedColors()

filename = 'FinalGeometry.vtu'

    # Read all the data from the file
reader = vtk.vtkXMLUnstructuredGridReader() #VTU
reader.SetFileName(filename)
reader.Update()
output = reader.GetOutput()

    # Visualize
mapper = vtk.vtkDataSetMapper()
mapper.SetInputData(output)

actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetColor(colors.GetColor3d('Yellow'))

renderer = vtk.vtkRenderer()
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

renderer.AddActor(actor)
renderer.SetBackground(colors.GetColor3d('DarkOliveGreen'))
renderer.GetActiveCamera().Pitch(90)
renderer.GetActiveCamera().SetViewUp(0, 0, 1)
renderer.ResetCamera()

renderWindow.SetSize(600, 600)
renderWindow.Render()
renderWindow.SetWindowName('ReadPolyData')
renderWindowInteractor.Start()


# In[ ]:




