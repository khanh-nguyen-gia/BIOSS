using Revise

using BIOSSFemCore
using BIOSSFemMaterial

material1 = MultiDifusion(2, [1., 2.], _2D)
display(material1.DMat)

getSpecies(material1)

material2 = MultiDifusion{2}(_2D)
display(material2.DMat)

mat1 = Difusion(_3D)
display(mat1.DMat)

mat3 = Elastic(_3D)
display(mat3.DMat)

mat1 = Elastic3D(_3D)
display(mat1.DMat)
