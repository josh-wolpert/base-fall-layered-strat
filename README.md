# base-fall-layered-strat
Matlab functions to model 1-D river profile evolution with layered rock erodibilities governed by the stream power law and extract knickpoint data from model output. Parameter descriptions and examples are included in function headings. If you use this code in a publication, please cite: Wolpert, J.A. & Forte, A.M. (2021) Response of transient rock uplift and base level knickpoints to erosional efficiency contrasts in bedrock streams. Earth Surface Processes and Landforms, 1–18. Available from: https:// doi.org/10.1002/esp.5146
## SPIM_1D_UD
Function to solve the stream power incision model in 1-D using an explicit upwind differencing scheme. Users have full control over Hack's Law, SPIM parameters, and the geometry and number of planar contacts.
## extractknickdata
Companion function to SPIM_1D_UD. Takes in model output and finds knickpoints (defined as break points in channel steepness) at each timestep. Knickpoint elevations are plotted as a function of time, and users can isolate breakpoints to calculate vertical velocities, chi-space celerities, and knickpoint prominence (i.e., cross-knickpoint change in channel steepness) for individual knickpoints through time by drawing a line of at least two points nearest a single knickpoint elevation path. The function also outputs the time vs. elevation plot of each model run, and the user has the option to view the output of SPIM_1D_UD while extractknickdata runs.
