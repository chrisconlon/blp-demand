*** INSTRUCTIONS
This package contains the following files

*extract_params.m: 
this is a user-configurable file that MUST be edited and determines the 
correct specification passed to the optimization routine. It is where you specify the random
components of utility and the corresponding derivatives.

*solveRCBLPpar.m: 
this is the main routine that optimizes the BLP-GMM objective function,
computes the weighting matrix, and produces the standard errors. As an intermediary output
it produces first-step.mat for diagnostic purposes. Only the first few lines of the file are 
considered user configurable.

** These are not considered user-configurable files
*solveAllShares.m:
This is the main routine which concentrates out the mean utility parameters by solving the
share equations market by market. If the computer has multiple cores, it is capable of spreading
this task in parallel over several cores. It takes on two user-configurable options 'newton'
sovles for the shares using Newton's Method (fsolve) and 'fixed-point' which solves for the
shares using a modified fixed-point algorithm.

*solveNewton.m:
This solves the system of J share equations for a single market using Newton's method (fsolve)
and requires computing the Jacobian of the shares with respect to the mean utilities at each
iteration

*fp_squarem.m:
This is a port of the SQUAREM fixed-point algorithm in the R package "SQUAREM". It does not
have as many options and only handles the simplest cases. It is a generic implmentation of 
SQUAREM and is not BLP specific. In this context it is used to solve the system of share equations
for the mean utilities market by market.

*rc_share_safe.m:
This is an under/over-flow safe implementation of the random coefficients logit probability
for a SINGLE market. This should be re-written in C++/mex for optimal performance.

*RCBLP_Jacobian.m
This computes the derivative of the shares with respect to the mean utilities (deltas) and
nonlinear parameters (theta) for a single market. It has two calling modes. In the 'mpec'
mode it returns both the delta and theta derivatives. In the 'sterr' mode it returns the
d s / d theta = [d s/ d delta]^-1 d s/d theta. For more information consult the appendix 
to Nevo (2000).

Copyright (c) 2014, Christopher T. Conlon (cconlon@columbia.edu)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.
