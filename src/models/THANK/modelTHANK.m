function [G1,C,impact,eu,SDX,zmat,NY,NX] = modelTHANK(param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For given parameter values, this function
% 1) puts the model of Justiniano, Primiceri and Tambalotti (2009)
%    in Gensys' canonical form
% 2) solves the RE system of equations using Chris Sims' Gensys
% 3) augments the resulting state evolution equation with the observation
%    equation
% The solution takes the form:  x(t) = G1 * x(t-1) + impact * e(t)
%                               y(t) = zmat * x(t) + C
%                               e(t) ~ N(0,SDX*SDX')
%
% Alejandro Justiniano, Giorgio Primiceri and Andrea Tambalotti  4/30/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% INDEX for endogenous variables
% -------------------------------------------------------------------------
y           = 1;  % output
k           = 2;  % capital
L           = 3;  % hours
Rk          = 4;  % return on capital (rho)
w           = 5;  % real wages
p           = 6;  % inflation
mc          = 7;  % marginal cost
lambda      = 8;  % multiplier
c           = 9;  % consumption
R           = 10; % interst rate
u           = 11; % capital utilization
phi         = 12; % multiplier
i_s         = 13; % investment
kbar_s      = 14; % physical capital
wgap        = 15; % wage gap (w - mrs)
gdp         = 16; % gdp

% Shocks
z           = 17; % productivity shock
g           = 18; % government spending shock
miu         = 19; % investment shock
lambdap     = 20; % price markup shock
lambdaw     = 21; % wage markup shock
b           = 22; % intertemporal preference shock
mp          = 23; % monetary policy shock

% Expectational Variables
ep          = 24; 
ec          = 25;
elambda     = 26;
ephi        = 27;
eRk         = 28;
e_i_s       = 29;
ew          = 30;

% Lagged Variables 
gdp_1       = 31; 
c_1         = 32;
i_s_1       = 33;
w_1         = 34;
ARMAlambdap = 35; % auxiliary variables for parameterization of 
ARMAlambdaw = 36; % ... ARMA processes for markup shocks

% New variables:
lam_h = 37;
lam_s = 38;
c_h   = 39;
c_s   = 40;
t_h   = 41;
k_s   = 42;

ey     = 43; %% TODO: Checking here...
elam_s = 44;
elam_h = 45;

var_len = 45;

% variables for the flex prices and wages version of the model
% -------------------------------------------------------------------------
ystar       = var_len + 1;
kstar       = var_len + 2; 
Lstar       = var_len + 3;
Rkstar      = var_len + 4;
wstar       = var_len + 5;
mcstar      = var_len + 6;
lambdastar  = var_len + 7;
cstar       = var_len + 8;
Rstar       = var_len + 9;
ustar       = var_len + 10;
phistar     = var_len + 11;
i_s_star    = var_len + 12;
kbar_s_star    = var_len + 13;
wgapstar    = var_len + 14;
gdpstar     = var_len + 15;

ecstar      = var_len + 16;
elambdastar = var_len + 17;
ephistar    = var_len + 18;
eRkstar     = var_len + 19;
e_i_s_star  = var_len + 20;

% New variables:
lam_h_star = var_len + 21;
lam_s_star = var_len + 22;
c_h_star   = var_len + 23;
c_s_star   = var_len + 24;
t_h_star   = var_len + 25;

ey_star     = var_len + 26;
elam_h_star = var_len + 27;
elam_s_star = var_len + 28;

NY = var_len + 28;

% -------------------------------------------------------------------------
% INDEX for Shocks
% -------------------------------------------------------------------------
Rs       = 1;     % Monetary policy
zs       = 2;     % Technology
gs       = 3;     % Government spending
mius     = 4;     % Investment
lambdaps = 5;     % Price markup
lambdaws = 6;     % Wage markup
bs       = 7;     % Intertemporal preference

NX = 7;

% -------------------------------------------------------------------------
% INDEX for expectational errors
% -------------------------------------------------------------------------
pex      = 1;
cex      = 2;
lambdaex = 3;
phiex    = 4;
Rkex     = 5;
i_s_ex   = 6;
wex      = 7;

% ===
cexstar      = 8;
lambdaexstar = 9;
phiexstar    = 10;
Rkexstar     = 11;
i_s_ex_star  = 12;
% === * css / (expg * c_s_ss - h * css)

NETA = 7+5;

% -------------------------------------------------------------------------
% Index for the parameters
% -------------------------------------------------------------------------

% Calibrated parameters
gss       = 0.22;  % capital depreciation rate
delta     = 0.025; % steady state government spending to GDP ratio
t_h_0_Lss = 0.0;   % New parameter! Might estimate.

% Non SD parameters
% -------------------------------------------------------------------------
alpha     = param(1);  % share of capital in the prod. function
iotap     = param(2);  % price indexation (=0 static indexation, =1 dynamic)
iotaw     = param(3);  % wages indexation
gamma100  = param(4);  % steady state growth rate of technology (multiplied by 100)
h         = param(5);  % habit formation
lambdapss = param(6);  % steady state net price markup
lambdawss = param(7);  % steady state net wage markup
Lss       = param(8);  % steady state for log hours
pss100    = param(9);  % steady state quarterly net inflation (multiplied by 100)
Fbeta     = param(10); % 100 * ( inv(discount factor) - 1 )
niu       = param(11); % inverse Frisch elasticity
xip       = param(12); % price stickiness
xiw       = param(13); % wage stickiness
chi       = param(14); % elasticity of the capital utilization cost function
S         = param(15); % investment adjustment cost
fp        = param(16); % reaction to inflation
fy        = param(17); % reaction to output gap
fdy       = param(18); % reaction to output gap growth

rhoR = param(19); rhoz = param(20); rhog = param(21); rhomiu = param(22); 
rholambdap = param(23); rholambdaw = param(24); rhob = param(25); rhomp = param(26);
rhoARMAlambdap = param(27); rhoARMAlambdaw = param(28);

% Standard deviations
sdR = param(29); sdz = param(30); sdg = param(31); sdmiu = param(32); 
sdlambdap = param(33); sdlambdaw = param(34); sdb = param(35);
SDX = diag([sdR sdz sdg sdmiu sdlambdap sdlambdaw sdb]);

% New parameters
theta       = param(36);
sigma       = param(37);
sigma_prime = param(38);
tau_D       = param(39); % Where all the cases come in
tau_k       = param(40);

numpar = 35 + 5;  % Number of parameters
ncof   = 28 + 5;  % Number of coefficients not corresponding to standard deviations

% -------------------------------------------------------------------------
% Computation of the steady state
% -------------------------------------------------------------------------
gamma  = gamma100 / 100;
expg   = exp(gamma);
beta   = 100 / (Fbeta + 100);
rss    = expg / beta - 1; % rss100, pss100 pop into constants
rss100 = rss * 100; % REMOVED: pss    = pss100 / 100;
gss    = 1 / (1-gss);

expLss = exp(Lss);
Rkss   = (expg/beta-1+delta); % RHO
mcss   = 1 / (1 + lambdapss);
wss    = (mcss*((1-alpha)^(1-alpha))/((alpha^(-alpha))*Rkss^alpha))^(1/(1-alpha));
% Compute ratios wrt L_ss
kLss   = (wss/Rkss) * alpha/(1-alpha);
FLss   = (kLss^alpha - Rkss*kLss - wss);
yLss   = kLss^alpha - FLss;

i_s_Lss = (1 - (1-delta) * exp(-gamma)) * expg * kLss / (1-theta); %%%
cLss    = yLss/gss - (1-theta) * i_s_Lss; %%%

c_h_Lss = wss + t_h_0_Lss + tau_k * Rkss * kLss; %%%
c_s_Lss = (1/(1-theta))*cLss - (theta/(1-theta))*c_h_Lss; %%%

lam_h_Lss = expg / (expg*c_h_Lss - h*cLss);
lam_s_Lss = expg / (expg*c_s_Lss - h*cLss);

% Multiply by L_ss again
kss    = kLss * expLss;
F      = FLss * expLss;
yss    = yLss * expLss;
css    = cLss * expLss;

% New equations
c_h_ss = c_h_Lss * expLss;
c_s_ss = c_s_Lss * expLss;
i_s_ss = i_s_Lss * expLss;

lam_h_ss = lam_h_Lss / expLss;
lam_s_ss = lam_s_Lss / expLss;

R_ss = (pi*expg/beta) * lam_s_ss/(sigma*lam_s_ss + (1-sigma)*lam_h_ss);

% -------------------------------------------------------------------------
% System Matrices [stars = flexible price equilibrium]
% -------------------------------------------------------------------------
GAM0 = zeros(NY,NY);
GAM1 = zeros(NY,NY);
PSI  = zeros(NY,NX);
PPI  = zeros(NY,NETA);
C    = zeros(NY,1);

% eq 1 and 2, production function (y, ystar)
% -------------------------------------------------------------------------
GAM0(y, y) = 1;
GAM0(y, k) = -((yss + F) / yss) * alpha;     % (1 - lambda_p) ?
GAM0(y, L) = -((yss + F) / yss) * (1-alpha); % (1 - lambda_p) ?
% ===
GAM0(ystar, ystar) = 1;
GAM0(ystar, kstar) = -((yss + F) / yss) * alpha;
GAM0(ystar, Lstar) = -((yss + F) / yss) * (1-alpha);
% ===


% eq 3 and 4, cost minimization (L and Lstar)
% -------------------------------------------------------------------------
GAM0(L, L)  = -1;
GAM0(L, Rk) =  1;
GAM0(L, w)  = -1;
GAM0(L, k)  =  1;
% ===
GAM0(Lstar, Lstar)  = -1;
GAM0(Lstar, Rkstar) =  1;
GAM0(Lstar, wstar)  = -1;
GAM0(Lstar, kstar)  =  1;
% ===


% eq 5 and 6, marginal cost (s and mcstar )
% -------------------------------------------------------------------------
GAM0(mc, mc) = 1;
GAM0(mc, Rk) = -alpha;
GAM0(mc, w)  = -(1-alpha);
%===
GAM0(mcstar, mcstar) = 1;
GAM0(mcstar, Rkstar) = -alpha;
GAM0(mcstar, wstar)  = -(1-alpha);
%===


% eq 7 and 8, Phillips curve (p and Rstar)
% -------------------------------------------------------------------------
GAM0(p, p)       = 1;
GAM0(p, ep)      = -beta / (1 + iotap*beta);
GAM1(p, p)       = iotap / (1 + iotap*beta);
GAM0(p, mc)      = -((1 - beta*xip) * (1 - xip) / ((1 + iotap * beta) * xip)); % kappa_p
GAM0(p, lambdap) = -1; % normalized

%===
GAM0(Rstar, mcstar) = 1;
%===

%% TODO - CONSUMPTION
% eq 9 and 10, consumption foc (c and cstar)
% -------------------------------------------------------------------------
GAM0(c, lambda) = (expg - h*beta)*(expg - h);
GAM0(c, b)      = -(expg - h*beta*rhob)*(expg - h) / [(1-rhob)*(expg-h*beta*rhob)*(expg-h) / (expg*h+expg^2+beta*h^2)];
GAM0(c, z)      = -(beta*h*expg*rhoz - h*expg);
GAM0(c, c)      = (expg^2 + beta*h^2);
GAM0(c, ec)     = -beta * h * expg;
GAM1(c, c)      =         h * expg;

%===
GAM0(cstar, lambdastar) = (expg-h*beta)*(expg-h);
GAM0(cstar, b) = -(expg-h*beta*rhob)*(expg-h)/[(1-rhob)*(expg-h*beta*rhob)*(expg-h)/(expg*h+expg^2+beta*h^2)];
GAM0(cstar, z) = -(beta*h*expg*rhoz-h*expg);
GAM0(cstar, cstar) = (expg^2+beta*h^2);
GAM0(cstar, ecstar) = -beta*h*expg;
GAM1(cstar, cstar) = expg*h;
%===
%%

% eq 11 and 12, Euler equation (lambda and lambdastar)
% -------------------------------------------------------------------------
GAM0(lambda, lambda) = 1;
GAM0(lambda, lam_h)  = -theta;
GAM0(lambda, lam_s)  = -(1-theta);

GAM0(lam_s, lam_s) =  1;
GAM0(lam_s, b)     = -1;
GAM0(lam_s, z)     = h * css / (expg * c_s_ss - h * css);
GAM0(lam_s, c_s)   = expg * c_s_ss / (expg * c_s_ss - h * css);
GAM1(lam_s, c)     = h * css / (expg * c_s_ss - h * css);

GAM0(lam_h, lam_h) =  1;
GAM0(lam_h, b)     = -1;
GAM0(lam_h, z)     = h * css / (expg * c_h_ss - h * css);
GAM0(lam_h, c_h)   = expg * c_h_ss / (expg * c_h_ss - h * css);
GAM1(lam_h, c)     = h * css / (expg * c_h_ss - h * css);

GAM0(lam_s, R)      = -1;
GAM0(lam_s, z)      = rhoz;
GAM0(lam_s, ep)     = 1;
GAM0(lam_s, elam_s) = -(sigma*lam_s_ss)     / (sigma*lam_s_ss + (1-sigma)*lam_h_ss);
GAM0(lam_s, elam_h) = -((1-sigma)*lam_h_ss) / (sigma*lam_s_ss + (1-sigma)*lam_h_ss);
GAM0(lam_s, ey)     = -(sigma_prime)*sigma*(lam_s_ss - lam_h_ss) / ... 
                        (sigma*lam_s_ss + (1-sigma)*lam_h_ss); 
% ===
GAM0(lambdastar, lambdastar) = GAM0(lambda, lambda);
GAM0(lambdastar, lam_h_star) = GAM0(lambda, lam_h);
GAM0(lambdastar, lam_s_star) = GAM0(lambda, lam_s);

GAM0(lam_s_star, lam_s_star) = GAM0(lam_s, lam_s);
GAM0(lam_s_star, bstar)      = GAM0(lam_s, b);
GAM0(lam_s_star, zstar)      = GAM0(lam_s, z);
GAM0(lam_s_star, c_s_star)   = GAM0(lam_s, c_s);
GAM1(lam_s_star, cstar)      = GAM1(lam_s, c);

GAM0(lam_h_star, lam_h_star) = GAM0(lam_h, lam_h);
GAM0(lam_h_star, bstar)      = GAM0(lam_h, b);
GAM0(lam_h_star, zstar)      = GAM0(lam_h, z);
GAM0(lam_h_star, c_h_star)   = GAM0(lam_h, c_h);
GAM1(lam_h_star, cstar)      = GAM1(lam_h, c);

GAM0(lam_s_star, Rstar)       = GAM0(lam_s, R);
GAM0(lam_s_star, zstar)       = GAM0(lam_s, z);
GAM0(lam_s_star, epstar)      = GAM0(lam_s, ep);
GAM0(lam_s_star, elam_s_star) = GAM0(lam_s, elam_s);
GAM0(lam_s_star, elam_h_star) = GAM0(lam_s, elam_h);
GAM0(lam_s_star, ey_star)     = GAM0(lam_s, ey);     
% ===


% transfers
% -------------------------------------------------------------------------
GAM0(t_h, t_h) = theta;
GAM0(t_h, y)   = -tau_D;
GAM0(t_h, mc)  =  tau_D;
GAM0(t_h, L)   =  tau_D * (1-alpha);
GAM0(t_h, k)   =  tau_D * alpha - tau_k * (Rkss * kss / yss);
GAM0(t_h, Rk)  = -tau_k * (Rkss * kss / yss);
% ===
GAM0(t_h_star, t_h_star) = GAM0(t_h, t_h);
GAM0(t_h_star, ystar)    = GAM0(t_h, y);
GAM0(t_h_star, mcstar)   = GAM0(t_h, mc);
GAM0(t_h_star, Lstar)    = GAM0(t_h, L);
GAM0(t_h_star, kstar)    = GAM0(t_h, k);
GAM0(t_h_star, Rkstar)   = GAM0(t_h, Rk);
% ===


% eq 13 and 14, capital utilization foc (Rk and Rkstar)
% -------------------------------------------------------------------------
GAM0(Rk, Rk) = 1;
GAM0(Rk, u)  = -chi;
% ===
GAM0(Rkstar, Rkstar) = GAM0(Rk, Rk);
GAM0(Rkstar, ustar)  = GAM0(Rk, u);
% ===


% eq 15 and 16, capital foc (phi and phistar)
% -------------------------------------------------------------------------
GAM0(phi, phi)    = 1;
GAM0(phi, ephi)   = -beta * exp(-gamma) * (1-delta);
GAM0(phi, z)      = rhoz;
GAM0(phi, elam_s) = -(1 - beta*exp(-gamma) * (1-delta));
GAM0(phi, eRk)    = -(1 - beta*exp(-gamma) * (1-delta));
% ===
GAM0(phistar, phistar)     = GAM0(phi, phi);
GAM0(phistar, ephistar)    = GAM0(phi, ephi);
GAM0(phistar, z)           = GAM0(phi, z);
GAM0(phistar, elam_s_star) = GAM0(phi, elam_s);
GAM0(phistar, eRkstar)     = GAM0(phi, eRk);
% ===


% eq 17 and 18, investment foc (i and istar)
% -------------------------------------------------------------------------
GAM0(i_s, i_s)   = (1 + beta);
GAM0(i_s, lam_s) =  1 / (S * expg^2);
GAM0(i_s, phi)   = -1 / (S * expg^2);
GAM0(i_s, miu)   = -1 / (S * expg^2);
GAM0(i_s, z)     = (1 - beta*rhoz);
GAM0(i_s, e_i_s) = -beta;
GAM1(i_s, i_s)   = 1;
% ===
GAM0(i_s_star, i_s_star)   = GAM0(i_s, i_s);
GAM0(i_s_star, lam_s_star) = GAM0(i_s, lam_s);
GAM0(i_s_star, phistar)    = GAM0(i_s, phi);
GAM0(i_s_star, miu)        = GAM0(i_s, miu);
GAM0(i_s_star, z)          = GAM0(i_s, z);
GAM0(i_s_star, e_i_s_star) = GAM0(i_s, e_i_s);
GAM1(i_s_star, i_s_star)   = GAM1(i_s, i_s);
% ===


% eq 19 and 20, capital input (k and kstar)
% -------------------------------------------------------------------------
GAM0(k, k)      =  1;
GAM0(k, u)      = -1;
GAM0(k, z)      =  1;
GAM1(k, kbar_s) =  1;
% ===
GAM0(kstar, kstar)       = GAM0(k, k);
GAM0(kstar, ustar)       = GAM0(k, u);
GAM0(kstar, z)           = GAM0(k, z);
GAM1(kstar, kbar_s_star) = GAM1(k, kbar_s);
% ===


% eq 21 and 22, capital accumulation (kbar_s and kbar_s_star)
% -------------------------------------------------------------------------
GAM0(kbar_s, kbar_s) = 1;
GAM0(kbar_s, miu)    = -(1-(1-delta) * exp(-gamma));
GAM0(kbar_s, i_s)    = -(1-(1-delta) * exp(-gamma));
GAM1(kbar_s, kbar_s) =     (1-delta) * exp(-gamma);
GAM0(kbar_s, z)      =     (1-delta) * exp(-gamma);

% ===
GAM0(kbar_s_star, kbar_s_star) = 1;
GAM0(kbar_s_star, miu)         = -(1-(1-delta) * exp(-gamma));
GAM0(kbar_s_star, i_s_star)    = -(1-(1-delta) * exp(-gamma));
GAM1(kbar_s_star, kbar_s_star) =     (1-delta) * exp(-gamma);
GAM0(kbar_s_star, z)           =     (1-delta) * exp(-gamma);
% ===


% eq 23 and 24, wage Phillips curve (w and wstar)
% -------------------------------------------------------------------------
kappaw = (1-xiw*beta) * (1-xiw) / (xiw * (1+beta) * (1 + niu*(1 + 1/lambdawss)));
GAM0(w, w)       = 1;
GAM0(w, ew)      = -beta / (1 + beta);
GAM0(w, wgap)    = kappaw;
GAM0(w, p)       = (1 + beta*iotaw) / (1 + beta);
GAM0(w, ep)      = -beta / (1 + beta);
GAM0(w, z)       = (1 + beta*iotaw - beta*rhoz)/(1 + beta);
GAM0(w, lambdaw) = -1;
GAM1(w, w)       =     1 / (1+beta);
GAM1(w, p)       = iotaw / (1+beta);
GAM1(w, z)       = iotaw / (1+beta);
% ===
GAM0(wstar, wgapstar) = 1;
% ===


% eq 25 and 26, wage gap (wgap and wgapstar)
% -------------------------------------------------------------------------
GAM0(wgap,wgap) = 1;
GAM0(wgap,w) = -1;
GAM0(wgap,b) = 1/[(1-rhob)*(expg-h*beta*rhob)*(expg-h)/(expg*h+expg^2+beta*h^2)];
GAM0(wgap,L) = niu;
GAM0(wgap,lambda) = -1;

% ===
GAM0(wgapstar,wgapstar) = 1;
GAM0(wgapstar,wstar) = -1;
GAM0(wgapstar,b) = 1/[(1-rhob)*(expg-h*beta*rhob)*(expg-h)/(expg*h+expg^2+beta*h^2)];
GAM0(wgapstar,Lstar) = niu;
GAM0(wgapstar,lambdastar) = -1;
% ===


% eq 27, monetary policy rule (R)
% ----------------------------------------
GAM0(R,R) = 1;
GAM1(R,R) = rhoR;
GAM0(R,p) = -(1-rhoR)*fp;
GAM0(R,gdp) = -(1-rhoR)*fy-fdy;
GAM0(R,gdpstar) = (1-rhoR)*fy+fdy;
GAM1(R,gdp) = -fdy;
GAM1(R,gdpstar) = fdy;
GAM0(R,mp) = -1;


% eq 28 and 29,definition of gdp (gdp and gdpstar)
% -------------------------------------------------------------------------
GAM0(gdp,gdp) =  -1;
GAM0(gdp,y) = 1;
GAM0(gdp,u) =  -kss*Rkss/yss;

% ===
GAM0(gdpstar,gdpstar) =  -1;
GAM0(gdpstar,ystar) = 1;
GAM0(gdpstar,ustar) =  -kss*Rkss/yss;
% ===


% eq 30 and 31, market clearing (u and ustar)
% -------------------------------------------------------------------------
GAM0(u,c) = css/yss;
GAM0(u,i_s) = i_s_ss/yss;
GAM0(u,y) = -1/gss;
GAM0(u,g) = 1/gss;
GAM0(u,u) = kss*Rkss/yss;

% ===
GAM0(ustar,cstar) = css/yss;
GAM0(ustar,i_s_star) = i_s_ss/yss;
GAM0(ustar,ystar) = -1/gss;
GAM0(ustar,g) = 1/gss;
GAM0(ustar,ustar) = kss*Rkss/yss;
% ===


% eq 32 - 40, exogenous shocks
% -------------------------------------------------------------------------
GAM0(z,z) = 1; GAM1(z,z) = rhoz; PSI(z,zs) = 1;
GAM0(g,g) = 1; GAM1(g,g) = rhog; PSI(g,gs) = 1;

GAM0(lambdaw,lambdaw) = 1; GAM1(lambdaw,lambdaw) = rholambdaw; GAM0(lambdaw,ARMAlambdaw) = -1; GAM1(lambdaw,ARMAlambdaw) = -rhoARMAlambdaw;
GAM0(ARMAlambdaw,ARMAlambdaw) = 1; PSI(ARMAlambdaw,lambdaws) = 1;

GAM0(miu,miu) = 1; GAM1(miu,miu) = rhomiu; PSI(miu,mius) = 1;

GAM0(lambdap,lambdap) = 1; GAM1(lambdap,lambdap) = rholambdap; GAM0(lambdap,ARMAlambdap) = -1; GAM1(lambdap,ARMAlambdap) = -rhoARMAlambdap;
GAM0(ARMAlambdap,ARMAlambdap) = 1; PSI(ARMAlambdap,lambdaps) = 1;

GAM0(b,b) = 1; GAM1(b,b) = rhob; PSI(b,bs) = 1;
GAM0(mp,mp) = 1; GAM1(mp,mp) = rhomp; PSI(mp,Rs) = 1;


% eq 41 - 52, expectational terms
% -------------------------------------------------------------------------
GAM0(ep,p) = 1; GAM1(ep,ep) = 1; PPI(ep,pex) = 1;
GAM0(ec,c) = 1; GAM1(ec,ec) = 1; PPI(ec,cex) = 1;
GAM0(elambda,lambda) = 1; GAM1(elambda,elambda) = 1; PPI(elambda,lambdaex) = 1;
GAM0(ephi,phi) = 1; GAM1(ephi,ephi) = 1; PPI(ephi,phiex) = 1;
GAM0(eRk,Rk) = 1; GAM1(eRk,eRk) = 1; PPI(eRk,Rkex) = 1;
GAM0(e_i_s,i_s) = 1; GAM1(e_i_s,e_i_s) = 1; PPI(e_i_s,i_s_ex) = 1;
GAM0(ew,w) = 1; GAM1(ew,ew) = 1; PPI(ew,wex) = 1;

% ===
GAM0(ecstar,cstar) = 1; GAM1(ecstar,ecstar) = 1; PPI(ecstar,cexstar) = 1;
GAM0(elambdastar,lambdastar) = 1; GAM1(elambdastar,elambdastar) = 1; PPI(elambdastar,lambdaexstar) = 1;
GAM0(ephistar,phistar) = 1; GAM1(ephistar,ephistar) = 1; PPI(ephistar,phiexstar) = 1;
GAM0(eRkstar,Rkstar) = 1; GAM1(eRkstar,eRkstar) = 1; PPI(eRkstar,Rkexstar) = 1;
GAM0(e_i_s_star,i_s_star) = 1; GAM1(e_i_s_star,e_i_s_star) = 1; PPI(e_i_s_star,i_s_ex_star) = 1;
% ===

% eq 53 - 56, lagged variables
%--------------------------------------------------------------------------
GAM0(gdp_1,gdp_1) = 1; GAM1(gdp_1,gdp) = 1;
GAM0(c_1,c_1) = 1; GAM1(c_1,c) = 1;
GAM0(i_s_1,i_s_1) = 1; GAM1(i_s_1,i_s) = 1;
GAM0(w_1,w_1) = 1; GAM1(w_1,w) = 1;


% Solution of the RE system of equations using Chris Sims' Gensys
% -------------------------------------------------------------------------
[G1,C,impact,fmat,fwt,ywt,gev,eu] = GENSYS(GAM0,GAM1,C,PSI,PPI) ;


% Matrix that maps endogenous variables to observables
% -------------------------------------------------------------------------
zmat = zeros(7,size(G1,2));
zmat(1,gdp) = 1; zmat(1,gdp_1) = -1; zmat(1,z) = 1;
zmat(2,c) = 1; zmat(2,c_1) = -1; zmat(2,z) = 1;
zmat(3,i_s) = 1; zmat(3,i_s_1) = -1; zmat(3,z) = 1;
zmat(4,L) = 1;
zmat(5,w) = 1; zmat(5,w_1) = -1; zmat(5,z) = 1;
zmat(6,p) = 1;
zmat(7,R) = 1;


% constants in the observation equation
% -------------------------------------------------------------------------
C = zeros(7,1);
C(1) = gamma100;
C(2) = gamma100;
C(3) = gamma100;
C(4) = Lss;
C(5) = gamma100;
C(6) = pss100;
C(7) = pss100+rss100;
