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
s           = 7;  % marginal cost
lambda      = 8;  % multiplier
c           = 9;  % consumption
R           = 10; % interst rate
u           = 11; % capital utilization
phi         = 12; % multiplier
i_s         = 13; % investment
kbar        = 14; % physical capital
wgap        = 15; % wage gap (w - mrs)
gdp         = 16; % gdp

z           = 17; % productivity shock
g           = 18; % government spending shock
miu         = 19; % investment shock
lambdap     = 20; % price markup shock
lambdaw     = 21; % wage markup shock
b           = 22; % intertemporal preference shock
mp          = 23; % monetary policy shock

ep          = 24; % expectational variables
ec          = 25;
elambda     = 26;
ephi        = 27;
eRk         = 28;
e_i_s       = 29;
ew          = 30;
gdp_1       = 31;       % lagged variables
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
t_H   = 41;

var_len = 41;

% variables for the flex prices and wages version of the model
% -------------------------------------------------------------------------
ystar       = var_len + 1;
kstar       = var_len + 2; 
Lstar       = var_len + 3;
Rkstar      = var_len + 4;
wstar       = var_len + 5;
sstar       = var_len + 6;
lambdastar  = var_len + 7;
cstar       = var_len + 8;
Rstar       = var_len + 9;
ustar       = var_len + 10;
phistar     = var_len + 11;
i_s_star    = var_len + 12;
kbarstar    = var_len + 13;
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
t_H_star   = var_len + 25;


NY = var_len + 25;

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
% ===

NETA = 7+5;

% -------------------------------------------------------------------------
% Index for the parameters
% -------------------------------------------------------------------------

% calibrated parameters
gss   = 0.22;         % capital depreciation rate
delta = 0.025;        % steady state government spending to GDP ratio


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

numpar = 35;  % Number of parameters
ncof   = 28;  % Number of coefficients not corresponding to
              % standard deviations


% -------------------------------------------------------------------------
% Computation of the steady state
% -------------------------------------------------------------------------
gamma  = gamma100 / 100;
beta   = 100 / (Fbeta + 100);
rss    = exp(gamma) / beta - 1;
rss100 = rss * 100;
pss    = pss100 / 100;
gss    = 1 / (1-gss);

expLss = exp(Lss);
Rkss   = (exp(gamma)/beta-1+delta);
sss    = 1/(1+lambdapss);
wss    = (sss*((1-alpha)^(1-alpha))/((alpha^(-alpha))*Rkss^alpha))^(1/(1-alpha));
kLss   = (wss/Rkss)*alpha/(1-alpha);
FLss   = (kLss^alpha-Rkss*kLss-wss);
yLss   = kLss^alpha-FLss;
kss    = kLss*expLss;
i_s_ss = (1-(1-delta)*exp(-gamma))*kss*exp(gamma);
F      = FLss*expLss;
yss    = yLss*expLss;
css    = yss/gss-i_s_ss;


% -------------------------------------------------------------------------
% System Matrices
% -------------------------------------------------------------------------
GAM0 = zeros(NY,NY) ;
GAM1 = zeros(NY,NY) ;
PSI  = zeros(NY,NX) ;
PPI  = zeros(NY,NETA) ;
C    = zeros(NY,1) ;


% eq 1 and 2, production function (y, ystar)
% -------------------------------------------------------------------------
GAM0(y,y) = 1;
GAM0(y,k) = -((yss+F)/yss)*alpha;
GAM0(y,L) = -((yss+F)/yss)*(1-alpha);

% ===
GAM0(ystar,ystar) = 1;
GAM0(ystar,kstar) = -((yss+F)/yss)*alpha;
GAM0(ystar,Lstar) = -((yss+F)/yss)*(1-alpha);
% ===


% eq 3 and 4, cost minimization (L and Lstar)
% -------------------------------------------------------------------------
GAM0(L,Rk) = 1;
GAM0(L,w)  = -1;
GAM0(L,k)  = 1;
GAM0(L,L)  = -1;

% ===
GAM0(Lstar,Rkstar) = 1;
GAM0(Lstar,wstar)  = -1;
GAM0(Lstar,kstar)  = 1;
GAM0(Lstar,Lstar)  = -1;
% ===


% eq 5 and 6, marginal cost (s and sstar)
% -------------------------------------------------------------------------
GAM0(s,s)  = 1;
GAM0(s,Rk) = -alpha;
GAM0(s,w)  = -(1-alpha);

%===
GAM0(sstar,sstar)  = 1;
GAM0(sstar,Rkstar) = -alpha;
GAM0(sstar,wstar)  = -(1-alpha);
%===


% eq 7 and 8, Phillips curve (p and Rstar)
% -------------------------------------------------------------------------
GAM0(p,p) = 1;
GAM0(p,ep) = -beta/(1+iotap*beta);
GAM1(p,p) = iotap/(1+iotap*beta);
GAM0(p,s) = -((1-beta*xip)*(1-xip)/((1+iotap*beta)*xip));
GAM0(p,lambdap) = -1;

%===
GAM0(Rstar,sstar) = 1;
%===


% eq 9 and 10, consumption foc (c and cstar)
% -------------------------------------------------------------------------
expg = exp(gamma);
GAM0(c,lambda) = (expg-h*beta)*(expg-h);
GAM0(c,b) = -(expg-h*beta*rhob)*(expg-h)/[(1-rhob)*(expg-h*beta*rhob)*(expg-h)/(expg*h+expg^2+beta*h^2)];
GAM0(c,z) = -(beta*h*expg*rhoz-h*expg);
GAM0(c,c) = (expg^2+beta*h^2);
GAM0(c,ec) = -beta*h*expg;
GAM1(c,c) = expg*h;

%===
expg = exp(gamma);
GAM0(cstar,lambdastar) = (expg-h*beta)*(expg-h);
GAM0(cstar,b) = -(expg-h*beta*rhob)*(expg-h)/[(1-rhob)*(expg-h*beta*rhob)*(expg-h)/(expg*h+expg^2+beta*h^2)];
GAM0(cstar,z) = -(beta*h*expg*rhoz-h*expg);
GAM0(cstar,cstar) = (expg^2+beta*h^2);
GAM0(cstar,ecstar) = -beta*h*expg;
GAM1(cstar,cstar) = expg*h;
%===


% eq 11 and 12, Euler equation (lambda and lambdastar)
% -------------------------------------------------------------------------
GAM0(lambda,lambda) = 1;
GAM0(lambda,R) = -1;
GAM0(lambda,elambda) = -1;
GAM0(lambda,ep) = 1;
GAM0(lambda,z) = rhoz;

% ===
GAM0(lambdastar,lambdastar) = 1;
GAM0(lambdastar,Rstar) = -1;
GAM0(lambdastar,elambdastar) = -1;
GAM0(lambdastar,z) = rhoz;
% ===


% eq 13 and 14, capital utilization foc (Rk and Rkstar)
% -------------------------------------------------------------------------
GAM0(Rk,Rk) = 1;
GAM0(Rk,u) = -chi;

% ===
GAM0(Rkstar,Rkstar) = 1;
GAM0(Rkstar,ustar) = -chi;
% ===


% eq 15 and 16, capital foc (phi and phistar)
% -------------------------------------------------------------------------
GAM0(phi,phi) = 1;
GAM0(phi,ephi) = -beta*exp(-gamma)*(1-delta);
GAM0(phi,z) = rhoz;
GAM0(phi,elambda) = -(1-beta*exp(-gamma)*(1-delta));
GAM0(phi,eRk) = -(1-beta*exp(-gamma)*(1-delta));

% ===
GAM0(phistar,phistar) = 1;
GAM0(phistar,ephistar) = -beta*exp(-gamma)*(1-delta);
GAM0(phistar,z) = rhoz;
GAM0(phistar,elambdastar) = -(1-beta*exp(-gamma)*(1-delta));
GAM0(phistar,eRkstar) = -(1-beta*exp(-gamma)*(1-delta));
% ===


% eq 17 and 18, investment foc (i and istar)
% -------------------------------------------------------------------------
GAM0(i_s,lambda) = 1/(S*expg^2);
GAM0(i_s,phi) = -1/(S*expg^2);
GAM0(i_s,miu) = -1/(S*expg^2);
GAM0(i_s,i_s) = (1+beta);
GAM0(i_s,z) = (1-beta*rhoz);
GAM0(i_s,e_i_s) = -beta;
GAM1(i_s,i_s) = 1;

% ===
GAM0(i_s_star,lambdastar) = 1/(S*expg^2);
GAM0(i_s_star,phistar) = -1/(S*expg^2);
GAM0(i_s_star,miu) = -1/(S*expg^2);
GAM0(i_s_star,i_s_star) = (1+beta);
GAM0(i_s_star,z) = (1-beta*rhoz);
GAM0(i_s_star,e_i_s_star) = -beta;
GAM1(i_s_star,i_s_star) = 1;
% ===


% eq 19 and 20, capital input (k and kstar)
% -------------------------------------------------------------------------
GAM0(k,k) = 1;
GAM0(k,u) = -1;
GAM0(k,z) = 1;
GAM1(k,kbar) = 1;

% ===
GAM0(kstar,kstar)    = 1;
GAM0(kstar,ustar)    = -1;
GAM0(kstar,z)        = 1;
GAM1(kstar,kbarstar) = 1;
% ===


% eq 21 and 22, capital accumulation (kbar and kbarstar)
% -------------------------------------------------------------------------
GAM0(kbar,kbar) = 1;
GAM0(kbar,miu) = -(1-(1-delta)*exp(-gamma));
GAM0(kbar,i_s) = -(1-(1-delta)*exp(-gamma));
GAM1(kbar,kbar) = (1-delta)*exp(-gamma);
GAM0(kbar,z) = (1-delta)*exp(-gamma);

% ===
GAM0(kbarstar,kbarstar) = 1;
GAM0(kbarstar,miu) = -(1-(1-delta)*exp(-gamma));
GAM0(kbarstar,i_s_star) = -(1-(1-delta)*exp(-gamma));
GAM1(kbarstar,kbarstar) = (1-delta)*exp(-gamma);
GAM0(kbarstar,z) = (1-delta)*exp(-gamma);
% ===


% eq 23 and 24, wage Phillips curve (w and wstar)
% -------------------------------------------------------------------------
kappaw = [(1-xiw*beta)*(1-xiw)]/[xiw*(1+beta)*(1+niu*(1+1/lambdawss))];
GAM0(w,w) = 1;
GAM0(w,ew) = -beta/(1+beta);
GAM0(w,wgap) = kappaw;
GAM0(w,p) = (1+beta*iotaw)/(1+beta);
GAM0(w,ep) = -beta/(1+beta);
GAM0(w,z) = (1+beta*iotaw-beta*rhoz)/(1+beta);
GAM0(w,lambdaw) = -1;
GAM1(w,w) = 1/(1+beta);
GAM1(w,p) = iotaw/(1+beta);
GAM1(w,z) = iotaw/(1+beta);

% ===
GAM0(wstar,wgapstar) = 1;
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