# bd-db-si_ivhi
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

