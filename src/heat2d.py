import numpy as np
import variables as var
import math
import FE_subroutines as FE


for e in range(var.nel):
    FE.assembly(e)

FE.src_flux(var.neq, var.nbe, var.ngp)

d = FE.solvedr(var.neq, var.nd)
    
print(d)
