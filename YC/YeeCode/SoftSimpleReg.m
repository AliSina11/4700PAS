winstyle = 'docked';
% winstyle = 'normal';

set(0,'DefaultFigureWindowStyle',winstyle)
set(0,'defaultaxesfontsize',18)
set(0,'defaultaxesfontname','Times New Roman')
% set(0,'defaultfigurecolor',[1 1 1])

% clear VARIABLES;
clear
global spatialFactor;
global c_eps_0 c_mu_0 c_c c_eta_0
global simulationStopTimes;
global AsymForcing
global dels
global SurfHxLeft SurfHyLeft SurfEzLeft SurfHxRight SurfHyRight SurfEzRight



dels = 0.75;
spatialFactor = 1;

c_c = 299792458;                  % speed of light
c_eps_0 = 8.8542149e-12;          % vacuum permittivity
c_mu_0 = 1.2566370614e-6;         % vacuum permeability
c_eta_0 = sqrt(c_mu_0/c_eps_0);


tSim = 200e-15  % set simulation time to 200 femtoseconds (fs) 
f = 230e12;     % set frequency to 230 TeraHertz (THz) (Altered)
lambda = c_c/f; % set wavelength to speed of light over frequency

xMax{1} = 20e-6;    % set xMax to 20 micrometers
nx{1} = 200;        % set nx to 200
ny{1} = 0.75*nx{1}; % set nx to 75% of nx


Reg.n = 1; % create Reg struct with Reg.n = 1

mu{1} = ones(nx{1},ny{1})*c_mu_0; % create vacuum permeability mapping

epi{1} = ones(nx{1},ny{1})*c_eps_0;  % create vacuum permitivvity mapping




% epi{1}(125:150,55:95)= c_eps_0*11.3; % adjust permitivitty within region
                                     % perhaps this is the "inclusion"?
                                     
% epi{1}(100:120,55:95)= c_eps_0*11.3; % we added another one

% epi{1}(70:90,55:95)= c_eps_0*11.3; % we added another one

% AMOGUS GRATING PATTERN

% FORMAT:
% epi{1}(xMin:xMax , yMin:yMax) = permeability

% epi{1}(60:140,45:120)= c_eps_0*11.3; % MAIN SQUARE

% epi{1}(50:60,70:90)= c_eps_0*11.3; % BACKPACK

% epi{1}(97:103,45:60)= c_eps_0; % Space between

% epi{1}(120:160,75:100)= c_eps_0*11.3*0.5; % VISOR

% VISOR up
% epi{1}(120:160,99:100)= c_eps_0*11.3*5;
% VISOR down
% epi{1}(120:160,75:76)= c_eps_0*11.3*5;
% VISOR left
% epi{1}(120:121,75:100)= c_eps_0*11.3*5;
% VISOR right
% epi{1}(159:160,75:100)= c_eps_0*11.3*5;

% LOSS

% Top Left Panel
%epi{1}(50:51,85:135)= c_eps_0*11.3;

% Top Right Panel
%epi{1}(135:136,85:135)= c_eps_0*11.3;
%epi{1}(165:166,85:130)= c_eps_0*11.3;

% Bottom Left Panel

%epi{1}(35:36,10:60)= c_eps_0*11.3;
%epi{1}(65:66,10:55)= c_eps_0*11.3;

% Bottom Right Panel

%epi{1}(135:136,10:60)= c_eps_0*11.3;
%epi{1}(140:180,10:11)= c_eps_0*11.3;

% PANEL BOUNDS
%{
epi{1}(100:101,2:148)= c_eps_0*11.3*10; % Y
 
epi{1}(1:198,75:76)= c_eps_0*11.3*10; % X
epi{1}(3:4,70:80)= c_eps_0*11.3*8; % X1
epi{1}(9:10,65:85)= c_eps_0*11.3*6; % X1
 
epi{1}(196:197,70:80)= c_eps_0*11.3*8; % X1
epi{1}(189:190,65:85)= c_eps_0*11.3*6; % X1
 
epi{1}(109:110,40:105)= c_eps_0*11.3*8; % X1
epi{1}(120:121,40:105)= c_eps_0*11.3*6; % X1
 
epi{1}(89:90,40:105)= c_eps_0*11.3*8; % X1
epi{1}(79:80,40:105)= c_eps_0*11.3*6; % X1
%}




% note that higher frequency waves will be more so reflected when incident
% on regions of higher permitivitty (Snell's Law Reflection).

sigma{1} = zeros(nx{1},ny{1}); % create conductivity mapping
sigmaH{1} = zeros(nx{1},ny{1}); % create conductivity mapping 2?

dx = xMax{1}/nx{1}; % dx increments = xMax over number of x increments
dt = 0.25*dx/c_c;   % dt increments = 0.25 times dx, over speed of light
nSteps = round(tSim/dt*2);  % create step per every 2*dt
yMax = ny{1}*dx;            % yMax = 0.75*xMax
nsteps_lamda = lambda/dx 

movie = 1;
Plot.off = 0;
Plot.pl = 0;
Plot.ori = '13';
Plot.N = 100;
Plot.MaxEz = 2;
Plot.MaxH = Plot.MaxEz/c_eta_0;
Plot.pv = [0 0 90];
Plot.reglim = [0 xMax{1} 0 yMax];

% The bc structure is used for background conditions such as
% sources, PMLs, etc. Basically anything other than the "sandbox"
% itself

% SOURCE 1
bc{1}.NumS = 1;
bc{1}.s(1).xpos = nx{1}/(4) + 0;
bc{1}.s(1).type = 'ss';
bc{1}.s(1).fct = @PlaneWaveBC;

% SOURCE 2

bc{1}.NumS = 2;
bc{1}.s(2).xpos = nx{1};
bc{1}.s(2).type = 'ss';
bc{1}.s(2).fct = @PlaneWaveBC;

% mag = -1/c_eta_0;
mag = 1;
phi = 0;
omega = (f*2*pi);
betap = 0;
t0 = 30e-15;
st = -0.05; % 15e-15 or -0.05
s = 0;
y0 = yMax/2;
sty = 1.5*lambda;

% source 1 params
bc{1}.s(1).paras = {1.2*mag,phi,omega,betap,t0,st,s,y0,sty,'s'};

% source 2 params (set mag as negative to make it go the other way)

bc{1}.s(2).paras = {-1.2*mag,phi,omega,betap,t0,st,s,y0,sty,'s'};

Plot.y0 = round(y0/dx);

bc{1}.xm.type = 'a';
bc{1}.xp.type = 'e'; %% 'a' or 'e'
bc{1}.ym.type = 'a';
bc{1}.yp.type = 'a';

pml.width = 20 * spatialFactor;
pml.m = 3.5;

Reg.n  = 1;
Reg.xoff{1} = 0;
Reg.yoff{1} = 0;

RunYeeReg
