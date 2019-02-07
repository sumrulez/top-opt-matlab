# top-opt-matlab
Topology optimization for plain stress problem in Matlab. Using Hyperworks to prepare the model and view results.

What it inside:

1) FE solver for plain stress 
2) Topology optimization algorythm usign SIMP method
3) import and export of results to HyperWorks format

P.S It was my diploma project for M.Sc. in Engineering. If you consider to use this code and need assistance, please contact me.

Here is the pictire explaining input and result programm gives:
 - on the left there is a bracket modeled with quad elements and constrained in the bottom. Forces is applied via the lug.
 - on the right side is the result of topology optimization (maximixation of stiffness with constrained volume fraction). Red color indicates that material should be here, blue - material to be removed.

![alt text](https://github.com/sumrulez/top-opt-matlab/blob/master/docs/matlab_topo_2.JPG)
