This is a reduced version of the MATLAB code used to generate the March and November 2021 simulations of vaccinating 5-11 and 12-17 year olds in England.

The Matlab file `Vaccine_Compare.m` runs the simulation code and generates graphs of the impact of childhood vaccination. For compactness, this code only runs a single set of posterior values and uses surrogate vaccination and population size data.

`Vaccine_Compare.m` calls `Generate_Output.m` which in turn calls `Simulate_One_Region.m`; the underlying ODEs are within zLeakyVacc_ODEs.m`.

The code can be made to run far more rapidly by using `parfor` loops and by converting `LeakyVacc_ODEs.m` and `ODE_to_Observables.m` into [mex](https://uk.mathworks.com/help/matlab/call-mex-file-functions.html) files.
