# Spiral Wave for TTP Model 

Another approach to spiral wave initialization is to 'rig' the simulation domain such that there are 4 sub-domains which are all initialized to different sets of ionic and gating variable states. A visualization of this simulation case is shown below

<p align="center">
   <img src="./SpiralWave.jpg" width="600">
</p>

Initial conditions for the Ten-tusscher-Panfilov (TTP) EP model need to be specified in the solver.xml to override the default variables that are hard-coded inside the solver. Each set of initial conditions parameters are defined in a seperate xml file that is included in the solver.xml via the Include_xml parameter. 
