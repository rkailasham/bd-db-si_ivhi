# bd-db-si_ivhi

si_ivhi.f uses a semi-implicit predictor-corrector method to integrate the stochastic differential equation (SDE) of a 
FENE dumbbell with internal viscosity (IV) and hydrodynamic interactions (HI), over an ensemble of trajectories.

The code needs an "inp.dat" file and a "tstep.dat" file to be present in the same folder it resides, for it to execute.
Given below is a sample of an "inp.dat" file. Copy and paste these contents to a new file titled "inp.dat" and 
make sure that the file is in the same directory as the FORTRAN code.


----------------- inp.dat ------------------------------------------

  0.0    1.0        100.0         10.00       1.00   0.30   2.00    1

------------------ End of inp.dat ------------------------------------

[Add spaces between numbers if necessary]

The code also needs a "tstep.dat" file in the same directory. Given below are the contents of a sample file.

------------------ tstep.dat ------------------------------------------

3 1000

200
500
1000

------------------- Endof tstep.dat -----------------------------------

[200,500,1000 must appear in consecutive lines]

For compiling the code : 
gfortran si_ivhi.f
This will generate an a.out file
./a.out will execute the code.

[Anatomy of inp.dat and tstep.dat files]

From left to right, the quantities mentioned in the inp.dat file are as follows : 
[Solvent quality] [d for EV] [FENE parameter] [Shear rate] [IV parameter] [HI parameter] [HI tensor] [Input style]  
[HI tensor]  : (1) corresponds to regularized ossen-burgers tensor
               (2) corresponds to RPY
               Any other option will set the HI parameter=0 and proceed without HI
 [Input Style] : (1) will read initial configs from equilibrated database
                 (2) will start from final configurations of a previous run.

In tstep.dat,

First number refers to number of timestep widths. Second number refers to number of trajectories.
The next lines indicate the number of BD steps.  

The code picks the initial configurations of the dumbbells from an equilibrated database. FPATH directs the code to the location of the database.
