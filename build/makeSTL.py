import vtk

# Carregar o arquivo VTU
reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("DISP_TESTE.vtu")
reader.Update()

# Criar uma caixa delimitadora para a área desejada
bounding_box = vtk.vtkBox()
bounding_box.SetBounds(-30, 30, -15, 15, -1000, 1000)  # Ajuste conforme necessário

# Filtro de recorte para a área delimitada pela caixa
clip_box = vtk.vtkClipDataSet()
clip_box.SetInputData(reader.GetOutput())
clip_box.SetClipFunction(bounding_box)
clip_box.GenerateClippedOutputOn()
clip_box.SetInsideOut(True)
clip_box.Update()

# Converter a malha para polígono
triangulator = vtk.vtkDataSetTriangleFilter()
triangulator.SetInputData(clip_box.GetOutput())
triangulator.Update()

# Converter de UnstructuredGrid para PolyData
geometry_filter = vtk.vtkGeometryFilter()
geometry_filter.SetInputData(triangulator.GetOutput())
geometry_filter.Update()


# Salvar o resultado em um novo arquivo STL
writer = vtk.vtkSTLWriter()
writer.SetFileName("DISP_TESTE.stl")
writer.SetInputData(geometry_filter.GetOutput())
writer.Write()
