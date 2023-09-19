#!/usr/bin/env python
# coding: utf-8

# In[1]:


import vtk

colors = vtk.vtkNamedColors()

filename1 = 'DISP_SUBTRAIDO.stl'
reader1 = vtk.vtkSTLReader()
reader1.SetFileName(filename1)
reader1.Update()

poly1 = reader1.GetOutputPort()

filtro1 = vtk.vtkGeometryFilter()
filtro1.SetInputConnection(poly1)
filtro1.Update()

polydata1 = filtro1.GetOutput()

Mass1 = vtk.vtkMassProperties()
Mass1.SetInputConnection(filtro1.GetOutputPort())
Mass1.Update()

Vol_1 = Mass1.GetVolume()
Vol_1 = round(Vol_1, 4)

#print('Volume 1 =', Vol_1)

volume = str(Vol_1)

file = open('VOLUME.txt', 'w')
file.write(volume)
file.close()


# In[ ]:




