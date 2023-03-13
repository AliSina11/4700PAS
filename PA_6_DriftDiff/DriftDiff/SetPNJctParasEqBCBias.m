Coupled = 1;
TwoCarriers = 1;
RC = 1;

nx = 201; % 201 position-x steps
l = 1e-6; % length of 1 um

x =linspace(0,l,nx); % create linspace of points along junction, called x
dx = x(2)-x(1); % dx is the space between two x(i) and x(i+1)
xm = x(1:nx-1) + 0.5*dx; % xm is the spaces halfway between x(i) and x(i+1)

ni = x < l/2; % ni is the left half of x
pi = x >= l/2; % pi is the right half of x

% Const. 1/cm3 (100 cm/m)^3
Nd = 4e16 * 1e6;        % donor doping
Na = 1e16 * 1e6;        % acceptor doping
NetDoping(ni) = Nd;     % NetDoping in ni is Nd 
NetDoping(pi) = -Na;    % NetDoping in pi is -Na (-ve because holes)

% ALTERATION 3 START - Assign Linear Doping Profile

NetDoping =  linspace((4e16 * 1e6), ((-1)*1e16 * 1e6), 201);

% ALTERATION 3 END



x0 = l/2;
nw = l/20;
npDisturbance = 0; % 0e16*1e6*exp(-((x-x0)/nw).^2);

JBC = 1;

RVbc = 0;

TStop = 80000000*1e-18;
PlDelt = 1000000*1e-18;

Phi =  C.Vt *log(Na*Nd/(niSi*niSi));
W  = sqrt(2*EpiSi*(Nd+Nd)*(Phi)/(C.q_0*Nd*Na));
Wn = W*Na/(Nd+Na);
Wp = (W - Wn);

LVbc = Phi;

PlotSS = 0;
% PlotYAxis = {[0 Phi+0.1] [0e5 40e5] [-20e2 40e2]...
%     [0e21 2.5e22] [0 1.1e22] [0 20e43]...
%     [-5e33 5e33] [-5e33 5e33] [-0e8 3e8] ...
%     [1e-3 1e8] [-3e6 1e6] [0 2.5e22]};
doPlotImage = 0;

SecondSim = 1;
LVbc2 = Phi-0.3;
TStop2 = TStop +  80000000*1e-18;

fprintf('Phi: %g W: %g Wn: %g Wp: %g \n',Phi,W,Wn,Wp) 
