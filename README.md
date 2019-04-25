# 234Th

##Optimization

This is the code used to generate results of the 234Th based POC and TOC export paper. 

To do parameter optimization, go to the CODE directory. You will find multiple .m files that are matlab files.

optimization.m is the main function, which sets up the environment and defines constant model parameters. 

The first function being called by optimization.m is neglogpost(2).m, which builds up objective function and computes first- and second-derivatives.

The neglogpost function calls eqPcycle.m and eqThcycle.m to get results from phosphorus and thorium-234 cycle models.
Function buildPFD.m is the function that is used to build particle flux divergence operator.

In summary, the structure of the modle is as follows,

                                       
optimization.m  --> neglogpost(2).m --> eqPcycle.m (-->buildPFD.m)
                                    --> eqThcycle.m (-->buildPFD.m)
                                      
                                      
                                      
##Error estimations                                      
