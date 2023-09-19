#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
import vtk

filename = 'ARQ_FINAL.vtp'
    
reader = vtk.vtkXMLPolyDataReader() #VTp
reader.SetFileName(filename)

surface_filter = vtk.vtkDataSetSurfaceFilter()
surface_filter.SetInputConnection(reader.GetOutputPort())

triangle_filter = vtk.vtkTriangleFilter()
triangle_filter.SetInputConnection(surface_filter.GetOutputPort())

writer = vtk.vtkSTLWriter()
writer.SetFileName('MASK.stl')
writer.SetInputConnection(triangle_filter.GetOutputPort())
writer.Write()     

