% Play with the below: Coupled, TwoCarriers, RC (all set to 1 by default)

Coupled = 1;
TwoCarriers = 1;
RC = 1;

nx = 201; % nx is 201
l = 1e-6;

x =linspace(0,l,nx);
dx = x(2)-x(1);
xm = x(1:nx-1) + 0.5*dx;

Nd = 1e16 * 1e6; % Const. 1/cm3 (100 cm/m)^3

% ALTERATION 1 START: Assign linear gradient doping profile
%   doping profile is a row vector of size nx
%   Need Nd to be a vector of 201 points from 1e16 to 20e16 

Nd = (1e16 * 1e6):(0.095e16 * 1e6):(20e16 * 1e6); 

% the (0.095e16 * 1e6) was set to have the correct size (nx = 201)

                                 % Linear Gradient from: 
                                 % 1/cm3 (100 cm/m)^3     to 
                                 % 20/cm3 (2000 cm/m ?)^3
                                 
% probably would have been a good idea to just use linspace instead 
                                 
% ALTERATION 1 END:
                                 
% ALTERATION 2 START: Assign exponential gradient doping profile
%   doping profile (Nd) is a row vector of size nx
%   Need Nd to be a vector of 201 points from 1e16 to 20e16 

Nd = (1e16 * 1e6):(0.095e16 * 1e6):(20e16 * 1e6); 

% use logspace 

Nd = logspace(log10(1e16 * 1e6), log10(20e16 * 1e6), 201);

% does this use e^x or 10^x ?
                                 % Linear Gradient from: 
                                 % 1/cm3 (100 cm/m)^3     to 
                                 % 20/cm3 (2000 cm/m ?)^3

% ALTERATION 2 END:

NetDoping = ones(1,nx).*Nd; % doping profile is a row vector of size nx

x0 = l/2;
nw = l/20;
npDisturbance = 0 ; % 1e16*1e6*exp(-((x-x0)/nw).^2);

LVbc = 0;
RVbc = 0;

TStop = 14200000*1e-18;
PlDelt = 100000*1e-18;

% PlotYAxis = {[-1e-15 2e-15] [-2e-9 2e-9] [-1.5e-12 1.5e-12]...
%     [1e22 2e22] [0 1e22] [0 20e43]...
%     [-20e33 15e33] [-2.5e34 2e34] [-1.1e8 1.1e8] ...
%     [-1e8 1e8] [-10e-3 10e-3] [0 2e22]};

doPlotImage = 0;
PlotFile = 'Gau2CarRC.gif';
