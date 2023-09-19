#!/usr/bin/env python
# coding: utf-8

# In[6]:


import sys
import vtk

filename = 'DISP_TESTE.vtu'
    
reader = vtk.vtkXMLUnstructuredGridReader() #VTU
reader.SetFileName(filename)

surface_filter = vtk.vtkDataSetSurfaceFilter()
surface_filter.SetInputConnection(reader.GetOutputPort())

triangle_filter = vtk.vtkTriangleFilter()
triangle_filter.SetInputConnection(surface_filter.GetOutputPort())

writer = vtk.vtkSTLWriter()
writer.SetFileName('DISP_TESTE.stl')
writer.SetInputConnection(triangle_filter.GetOutputPort())
writer.Write()     

