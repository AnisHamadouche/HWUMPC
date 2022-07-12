%% Copyright @ Heriot-Watt University UDRC WP 2.1
%% Authoer: Yun Wu

%% 
clear;
clc;
close all;

%% Parameter
N = 1; % horison length
ROW = 4 * N;
COL = 4 * N;
EXPONENT  = 5; % exponent of floating point
SIGNIFICAND = 4; % significand of floating point
BIT_WIDTH = 24; % total bit width of fixed point
INTE_WIDTH = 12; % integer bit width of fixed point
TOTALBITS = 10; % total bit width of posit 
MAX_ITER = 30; % gradient descent maximum iteration number

%% !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! 
% this is the switch to turn on/off the mex compilation, the mex file must
% be re-compiled once any above parameters are change
buildmexfile = 0; 

%% !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! 
% this is the switch to turn on/off delete of mex file at the end
deletemexfile = 0;

%% !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! 
% this is the option to choose different precision mex functions
precision_opt = 0; % 0: standard float
                   % 1: arbitrary float
                   % 2: fixed point
                   % 3: arbitrary posit

%% build mex files
if buildmexfile==1
    mexmain( ...
           ROW,COL, ...
           EXPONENT,SIGNIFICAND, ...
           BIT_WIDTH,INTE_WIDTH, ...
           TOTALBITS, ...
           MAX_ITER, ...
           matlabinc);
end

%% load the model
[H, F, x0, m, n, l, xmax, xmin, ymax, ymin] = model_spacecraft(N);

%% cvx
runcvx;

%% MEX LPGD example
runlpgdx;

%% Plot Example
plotfigs;

%% Delete mex files
if deletemexfile==1
    delete('*.mexa64');
    delete('*.mexw64');
end