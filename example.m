%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example using povided sample data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;
delete example.diary;
diary  example.diary;

%==============================================================================
% Step 1: Read in the data
%==============================================================================
load sample.mat

%==============================================================================
% Step 2: Set starting values
%==============================================================================
% set starting values:
% should be the number of random coefficients 
% (3 here for the 'correlated_normal' % example; 
% should be 2 for the 'indep_normal' example)
startval = rand(3,1); 

%==============================================================================
% Step 3: Produce draws (already included as a structure in sample.mat)
%==============================================================================
disp('Draws are already included as a structure in sample.mat');

%==============================================================================
% Step 4: Call solveRCBLPpar
%==============================================================================
dtable.delta = -10*rand(size(dtable.s0t)); %initialize delta as a placeholder
results = solveRCBLPpar(dtable,draws,startval,'correlated_normal');

%==============================================================================
% Step 5: Print output
%==============================================================================
% random coefficients (covariance on x2 parameters)
[results.theta results.SE(1:length(startval)) results.theta./results.SE(1:length(startval))]
% point estimates for x1 parameters
[results.beta results.SE(length(startval)+1:end) results.beta./results.SE(length(startval)+1:end)]

diary off