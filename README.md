# 234Th

## Optimization ##

This is the code used to generate results of the 234Th based POC and TOC export paper. 

To do parameter optimization, go to the CODE directory. You will find multiple .m files that are matlab files.

optimization.m is the main function, which sets up the environment and defines constant model parameters. 

The first function being called by optimization.m is neglogpost(2).m, which builds up objective function and computes first- and second-derivatives.

The neglogpost function calls eqPcycle.m and eqThcycle.m to get results from phosphorus and thorium-234 cycle models.
Function buildPFD.m is the function that is used to build particle flux divergence operator.

In summary, the structure of the modle is as follows,

                                       
optimization.m  --> neglogpost(2).m --> eqPcycle.m (-->buildPFD.m)
                                    --> eqThcycle.m (-->buildPFD.m)
                                      
                                      
                                      
## Error estimations ##                                      

You will find a couple of Matlab functions in the errorbar directory. These functions are used to generate model errorbars usnig Monte Carlo method.

The model structure is as follows

EB_main.m first calls neglogpost.m that produces phosphorus and Th-234 fields based on randomly drawn model parameters from a normal distribution with mean of optimal parameters and variance of second derivative (Hession) evaluated at optimal parameter values. EB_main.m then call POC_export_errorbar.m and totC_export_errorbar.m to produce new POC and TOC export fluxes.


## Supplementary Methods #

Latex files to generate Supplementary is located in NOTES directory.
