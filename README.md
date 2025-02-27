# 2D-2-order-corner-state
Calculate 2D second-order corner state in WannierTools.

cp the ek_ribbon.f90 and ham_twoborder.f90 in the wanniertools-2.7.1/src folder and recompile.

BBH-example can be tested and you should open the WireBand_calc = T in the wt.in.

This is more useful for systems with chiral symmetry when calculating specific materials, because their corner states are pinned at zero energy. Otherwise you need to modify line 103 in ek_ribbon.
