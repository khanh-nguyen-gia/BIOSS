using Revise

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements

EjemploMaterial = CellVegf(_2D)

display(EjemploMaterial.DMat)

EjemploMaterial = CellVegf(_2D, 5., 6.)

display(EjemploMaterial.DMat)


EjemploMaterial2 = CellVegf2(_2D)

display(EjemploMaterial2.DMat)
display(EjemploMaterial2.VEGF)

EjemploMaterial3 = VEGF(_2D)
display(EjemploMaterial3.DMat)

EjemploMaterial4 = Fibronectin(_2D)
display(EjemploMaterial4.DMat)
