using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements
using BIOSSFemDriver


FemMesh = RegularMesh(_2D,Quad4,Difusion,[20, 20], [0.1, 0.1])


### test for fuctions: get_node_indexes and get_node_coordinates

for i in 1: FemMesh.NofNodes
   indexes = get_node_indexes(i, FemMesh)
   value = eval(Meta.parse("FemMesh.Nodes"*"$indexes"))
   if value != i
      println(" Nodo $i incorrecto; index = $indexes, valor = $value")
   end
   coordinates = get_node_coordinates(i, FemMesh)

end
