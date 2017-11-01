Code for creating and manipulating orbital tori.

The basic code for axisymmetric tori is described in Binney & McMillan MNRAS
456 1982 (2016)

The code for computing tori trapped around Omega_r=Omega_z is described in
Binney MNRAS 462 2792 (2016)

The code for computing non-axisymmetric tori including tori trapped at a CR
or Lindblad resonance is described in Binney MNRAS in press (2017)

To install the update, 

untar the tar file, go to directory utils and do make libutil.a 

Then go to the home directory and do make libtorus2.a 

Then to create tori traped at CR do make corot.exe This executable requires
as its argument the bar's pattern speed in Myr^{-1} times 1000. So for 0.04
Myr^{-1} do corot 40

To make tori trapped at OLR do make lindblad.exe and then lindblad 40, etc

The excutable nores.exe reads data files produced by corot.exe and
lindblad.exe and creates untrapped tori to produce a plot of the velocities
at which stars of known actions, trapped and untrapped, visit the Sun
