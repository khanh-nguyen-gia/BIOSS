using Revise
#using Plots

using BIOSSFemCore
using BIOSSFemMaterial
using BIOSSFemElements

EjemploEndothelium = Endothelium(_2D)
EjemploTAF = TAF(_2D)
EjemploHypoxCell = HypoxCell(1, [0., 0.], [1, 2], [1, 2, 3, 4])
HypoxCellDeactivation!([1], [EjemploHypoxCell])

EjemploTipCell = TipCell(1,_2D)

println("Tip Cell succesfully created,
        ID= $(EjemploTipCell.ID)
        BirthDate = $(EjemploTipCell.BirthDate)
        Genealogy = $(EjemploTipCell.Genealogy); mother vessel: $(mothervessel(EjemploTipCell))
        coordinates: $(EjemploTipCell.coords)
        Position: Element $(EjemploTipCell.Element), Nodes $(EjemploTipCell.Nodes), distances: $(EjemploTipCell.distances)
        velocity: $(EjemploTipCell.velocity)
        radius: $(EjemploTipCell.radius)
        Density = $(EjemploTipCell.rho)")

EjemploTipCell2 = TipCell(2, 0.4, [0, 1], _2D)
EjemploTipCell3 = TipCell(3, 14.5, [0], _2D)
EjemploTipCell4 = TipCell(4, 1.3, [0, 1, 2], _2D)

mothervessel(EjemploTipCell)
mothervessel(EjemploTipCell2)

List = TipCellList(_2D)
#in case we want to use quad9, quad8, etc elemts: previous function would have to be modified
TC_vector =[EjemploTipCell, EjemploTipCell2, EjemploTipCell3, EjemploTipCell4]

expandTipCellList!(List, TC_vector)

mothervessel(4, List)
TipCellDeactivation!(5, List)
TipCellDeactivation!([5], List)
TipCellDeactivation!([4, 5], List)
TipCellDeactivation!([2, 3], List)
TipCellDeactivation!(1, List)
TipCellDeactivation!([1], List)
TipCellDeactivation!(1, List)
