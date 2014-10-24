#!/bin/csh

bin/Choose_any_df pot/PJM11_convenient.Tpot df/iso_GaiaSchoolMexico.df test.J 100000
bin/Create_DF_tori pot/PJM11_convenient.Tpot test.J test.list 1000 > & ! test.J_out
bin/Sample_list_limits test.list test.st 10000
# bin/Sample_list_limits test.list test2.st 1000 Rmin=8.3 Rmax=8.7 zmin=-0.5 zmax=0.5 phimin=-0.025 phimax=0.025
bin/Coord_converter GCYfromGCA test.st test.st.Rzp
bin/Coord_converter HGPfromGCA test.st test.st.hgp
bin/Add_Uncert test.st.hgp test.st.hgp.obs 0 5 0

#end
