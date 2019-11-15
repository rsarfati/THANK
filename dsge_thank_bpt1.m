clear all;

%DSGE-THANK 
%Number of Variables nx is the number of variables in the whole system, nz is the number of non-predetermined variables
nx = 35; nz = 7; nu = nz;
%/****************************************************/
%     @ Parameters of the model;@
%/****************************************************/
%Deep Parameters
beta=.9988;
delta=0.025;
rho=beta^(-1)-(1-delta);
alfa=0.17;
netmkup=0.29;
netmkupw=0.29;
h=0.84;
xip=0.77;
xiw=0.73;
iotap=0.05;
iotaw=0.08;
niu=1.96;
invelas=3.58;
chi=4.95;
rhor=0.78;
phipi=2.10;
phix=0.94;
phigrowth=0.0;

phia=0.9;
phiepsilon=0.8;
phimpshock=0.6;

loop=1; %loop over parameters; 
for loop =1:2
if loop==1
theta = 0.3;
s=1;
eta=0;
th0_L=0;
tauk=0;
taud=0;
w=(((alfa^alfa)*((1-alfa)^(1-alfa)))/((1+netmkup)*(rho^alfa)))^(1/(1-alfa));
k_L=alfa*w/((1-alfa*rho));
y_L=(1+netmkup)^(-1)*(k_L)^(alfa);
is_L=delta*k_L/(1-theta);
c_L=y_L-(1-theta)*is_L;
ch_L=w+th0_L+tauk*rho*k_L;
cs_L=((1-theta)^(-1))*c_L+((1-theta)^(-1))*theta*ch_L;
lamh_L=(ch_L-h*c_L)^(-1);
lams_L=(cs_L-h*c_L)^(-1);
lam_L=theta*lamh_L+(1-theta)*lams_L;
kappaw = (1-xiw*beta) * (1-xiw) / (xiw * (1+beta) * (1 + niu*(1 + 1/netmkupw)));

%STEADY STATE RATIOS
chc=ch_L/c_L;
csc=cs_L/c_L;
rhok_y=rho*k_L/y_L;

else
theta = 0.3;
s=1;
eta=0;
th0_L=0;
tauk=0;
taud=0;
w=(((alfa^alfa)*((1-alfa)^(1-alfa)))/((1+netmkup)*(rho^alfa)))^(1/(1-alfa));
k_L=alfa*w/((1-alfa*rho));
y_L=(1+netmkup)^(-1)*(k_L)^(alfa);
is_L=delta*k_L/(1-theta);
c_L=y_L-(1-theta)*is_L;
ch_L=w+th0_L+tauk*rho*k_L;
cs_L=((1-theta)^(-1))*c_L+((1-theta)^(-1))*theta*ch_L;
lamh_L=(ch_L-h*c_L)^(-1);
lams_L=(cs_L-h*c_L)^(-1);
lam_L=theta*lamh_L+(1-theta)*lams_L;
kappaw = (1-xiw*beta) * (1-xiw) / (xiw * (1+beta) * (1 + niu*(1 + 1/netmkupw)));

%STEADY STATE RATIOS
chc=ch_L/c_L;
csc=cs_L/c_L;
rhok_y=rho*k_L/y_L;

end
%THE Equilibrium Conditions

A = zeros(nx,nx);  B = zeros(nx,nx);  C = zeros(nx,nu);
phi = zeros(nu,nu); %exogenous shock processes

ypos = 1;  %
kpos = 2; %
lpos = 3; %
cpos = 4;  %
rhopos=5; % 
mcpos = 6;  %
pipos = 7; %
epi1pos=8;
lamspos=9;
cspos=10;
lamhpos=11;
chpos=12;
lampos=13; % 
rpos=14;
elams1pos=15;
thpos=16;
upos=17;
phipos=18;
ispos=19;
eis1pos=20;
wpos=21;
ew1pos=22;
gwpos=23;
xpos=24;
ephi1pos=25;
erho1pos=26;
elamh1pos=27;
ey1pos=28;

%Predetermined Variables
nk =7;

kspos = 29; %K(t)
wlagpos=30;
pilagpos=31;
clagpos=32;
islagpos=33;
rlagpos=34;
xlagpos=35;



%Shocks
apos = 1;
mpshockpos = 2;
epsilonpos=3; %investment-specific technlogical change

%Definiton of variables
%Equation 1: Prod function
B(1,ypos)=1;
B(1,kpos)=-(1+netmkup)*alfa;
B(1,lpos)=-(1+netmkup)*(1-alfa);

%@Equation 2: rental rate@
B(2,rhopos)=1;
B(2,wpos)=-1;
B(2,lpos)=-1;
B(2,kpos)=1;

%@ mg cost @
B(3,lampos)=1;
B(3,rhopos)=-alfa;
B(3,wpos)=alfa-1;

% philips curve
B(4, pipos)     = 1; 
B(4, epi1pos)     = -beta/(1+iotap*beta);
B(4, pilagpos)     = -iotap/(1+iotap*beta);
B(4, mcpos)       = -((1 - beta*xip) * (1 - xip) / ((1 + iotap * beta) * xip)); % kappa_p;

%@ Mg utility S  @
B(5,lamspos) = 1;
B(5,cspos)=csc/(csc-h);
B(5,clagpos)=-h/(csc-h);

%@ Mg utility H  @
B(6,lamhpos) = 1;
B(6,chpos)=chc/(chc-h);
B(6,clagpos)=-h/(chc-h);

% aggregate c
B(7, lampos)     = lam_L; 
B(7, lamhpos)  = -theta*lamh_L;
B(7, lamspos)  = -(1-theta)*lams_L;

%@ euler self-insurance  @
B(8,lamspos) = 1;
B(8,rpos)=-1;
B(8,epi1pos)=1;
B(8,elams1pos)=-s*lams_L/(s*lams_L+(1-s)*lamh_L);
B(8,elamh1pos)=-(1-s)*lamh_L/(s*lams_L+(1-s)*lamh_L);
B(8,ey1pos)=-eta*s*(lams_L-lamh_L)/(s*lams_L+(1-s)*lamh_L);

% transfer function
B(9, thpos)     = 1; 
B(9, ypos)  = -taud/theta;
B(9, mcpos)  = taud/theta;
B(9, kpos)  = (alfa*taud/theta)-tauk*rhok_y/theta;
B(9, lpos)  = (1-alfa)*taud/theta;
B(9, rhopos)  = -tauk*rhok_y/theta;

% c of H under ARBITRARY redistribution
B(10, chpos)     = 1; 
B(10, wpos)  = -w/ch_L;
B(10, lpos)  = -w/ch_L;
B(10, thpos)  = y_L/ch_L;

% aggregate c
B(11, cpos)     = 1; 
B(11, chpos)  = -theta*chc;
B(11, cspos)  = -(1-theta)*csc;

% utilization
B(12, rhopos)     = 1; 
B(12, upos)  = -chi;

%@ Euler Capital tobin's phi-q  @
B(13,phipos)=1;
B(13,ephi1pos)=-beta*(1-delta);
B(13,elams1pos)=-(1-beta*(1-delta));
B(13,erho1pos)=-(1-beta*(1-delta));

%@ Investment choice  @
B(14,lamspos)=1;
B(14,phipos)=-1;
B(14,ispos)=invelas*(1+beta);
B(14,islagpos)=-invelas;
B(14,eis1pos)=-invelas*beta;

% Kapital
B(15, kpos)     = 1; 
B(15, upos)  = -1;
B(15, kspos)  = -1;

%@capital accumulation  @
 A(16,kspos)=1;
 B(16,ispos)=delta;
 B(16,kspos)=(1-delta);

 %@wage PC  @
B(17,wpos)=1+beta;
B(17,wlagpos)=-1;
B(17,ew1pos)=-beta;
B(17,gwpos)=(1+beta)*kappaw;
B(17,pilagpos)=-iotaw;
B(17,pipos)=1+beta*iotaw;
B(17,epi1pos)=-beta;

%@MRS - wage gap  @
B(18,gwpos)=1;
B(18,wpos)=-1;
B(18,lpos)=niu;
B(18,lampos)=-1;

% Taylor
B(19, rpos)     = 1; 
B(19, rlagpos)     = -rhor; 
B(19, pipos)  = -(1-rhor)*phipi;
B(19, xpos)  = -(1-rhor)*phix - phigrowth;
B(19, xlagpos)  = phigrowth;
C(19, mpshockpos)  = 1; % RATE CUT

% gdp
B(20,xpos) = 1;
B(20,ypos) = -1;
B(20,upos) = rhok_y;

% market clearing
B(21,ypos) = 1;
B(21,cpos) = -c_L/y_L;
B(21,ispos) = -(1-theta)*is_L/y_L;
B(21,upos) = -rhok_y;

%Define E(t)pi(t+1)
A(22,pipos) = 1;
B(22,epi1pos) = 1;

%Define E(t)cs(t+1)
A(23,lamspos) = 1;
B(23,elams1pos) = 1;

%Define E(t)cs(t+1)
A(24,ispos) = 1;
B(24,eis1pos) = 1;

%Define E(t)cs(t+1)
A(25,wpos) = 1;
B(25,ew1pos) = 1;

%Define E(t)cs(t+1)
A(26,phipos) = 1;
B(26,ephi1pos) = 1;

%Define E(t)cs(t+1)
A(27,rhopos) = 1;
B(27,erho1pos) = 1;

%Define E(t)cs(t+1)
A(28,lamhpos) = 1;
B(28,elamh1pos) = 1;

%Define E(t)cs(t+1)
A(29,ypos) = 1;
B(29,ey1pos) = 1;

% lagged w
A(30,pilagpos) =1;
B(30,pipos)    =1;

% lagged w
A(31,clagpos) =1;
B(31,cpos)    =1;

% lagged w
A(32,islagpos) =1;
B(32,ispos)    =1;

% lagged w
A(33,wlagpos) =1;
B(33,wpos)    =1;

% lagged w
A(34,rlagpos) =1;
B(34,rpos)    =1;

% lagged w
A(35,xlagpos) =1;
B(35,xpos)    =1;



phi(apos,apos) = phia; %technology process introduced here
phi(mpshockpos,mpshockpos) = phimpshock;
phi(epsilonpos,epsilonpos) = phiepsilon;

[m,n,p,q,z22h,s,t,lambda] = solvek(A,B,C,phi,nk);
bigmn = [m n];
bigpq = [p q];
bigp = phi;
bigpsi = eye(nz,nu);

%%%%%%%%%%%%%%%%%% IRF analysis %%%%%%%%%%%%%%%%%%%%%%%%
% Have to specify ires and ishock, the index values for the 
% responding variable and the shock
% Using the solution of the model in state space form
% x(t+1) = Ax(t) + Bu(t+1)
% y(t) = Cx(t) + Du(t)
A = [p q;zeros(nz,nk) phi]; 
C = [m n];
D = zeros(nx-nk,nu);
B = [zeros(nk,nu);bigpsi];
% impo.m
ishock = 2; % .

npts = 40; % no of points plotted

% Loop over parameters, assign calculated solution for each case to one IRF
if loop ==1
A1 = [p q;zeros(nz,nk) phi]; 
C1 = [m n];
D1 = zeros(nx-nk,nu);
B1 = [zeros(nk,nu);bigpsi];
%[Y,X]=dimpulse(A,B,C,D,ishock,npts+1); 
%this does not work when have capital defined as here
[Y1,X1]=dimpulse(A1,B1,C1,D1,ishock,npts+2);
YP1=Y1(2:npts+2,:);
XP1=X1(2:npts+2,:);
else
A2 = [p q;zeros(nz,nk) phi]; 
C2 = [m n];
D2 = zeros(nx-nk,nu);
B2 = [zeros(nk,nu);bigpsi];
%[Y,X]=dimpulse(A,B,C,D,ishock,npts+1); 
%this does not work when have capital defined as here
[Y2,X2]=dimpulse(A2,B2,C2,D2,ishock,npts+2);
YP2=Y2(2:npts+2,:);
XP2=X2(2:npts+2,:);
end %ends if statement
%move on after first case:
loop = loop+1;
end %ends 'for' statement

jj=[0:npts];
%i1 = Y(:,ires);  % column index is the element of y you want

subplot(4,4,1)
plot(jj,YP1(:,ypos), jj,YP2(:,ypos),'-.r')
title('Output')
%axis([0 20 -.5 .5])
legend('stickyW-TANK (chi=1)','flex-W TANK (chi=1)')
text(0,YP1(1,ypos)+0.7,'FIGURE 1: Interest rate cut ') 
subplot(4,4,2)
plot(jj,YP1(:,cpos),jj,YP2(:,cpos),'-.r') 
title('Consumption')
%axis([0 20 -.25 .25]) 
subplot(4,4,3)
plot(jj,YP1(:,ispos),jj,YP2(:,ispos),'-.r')
title('Investment')
%axis([0 20 -.5 .5]) 
subplot(4,4,4)
plot(jj,YP1(:,lpos),jj,YP2(:,lpos),'-.r')
title('Hours worked')
%axis([0 20 -.5 .5]) 
subplot(4,4,5)
plot(jj,YP1(:,wpos),jj,YP2(:,wpos),'-.r') 
title('Real Wage')
%axis([0 20 -.25 .25]) 
subplot(4,4,6)
plot(jj,XP1(:,1),jj,XP2(:,1),'-.r')
title('Capital')
%axis([0 20 -.5 .5]) 
subplot(4,4,7)
plot(jj,YP1(:,rpos),jj,YP2(:,rpos),'-.r')
title('Nominal rate')
%axis([0 20 -.5 .5]) 
subplot(4,4,8)
plot(jj,YP1(:,pipos),jj,YP2(:,pipos),'-.r') 
title('Inflation')
%axis([0 20 -.25 .25]) 
subplot(4,4,9)
%plot(jj,XP1(:,3),jj,XP2(:,3),'-.r')
%title('Govt. Spending')
plot(jj,YP1(:,phipos),jj,YP2(:,phipos),'-.r')
title('Tobins q')
%axis([0 20 -.5 .5]) 
subplot(4,4,10)
plot(jj,YP1(:,cspos),jj,YP2(:,cspos),'-.r')
title('Consumption Savers')
%axis([0 20 -.5 .5]) 
subplot(4,4,11)
plot(jj,YP1(:,chpos),jj,YP2(:,chpos),'-.r') 
title('Consumption Hand-to-Mouth')
%axis([0 20 -.25 .25]) 
subplot(4,4,12)
%plot(jj,XP1(:,3),jj,XP2(:,3),'-.r')
%title('Govt. Spending')
plot(jj,YP1(:,rhopos),jj,YP2(:,rhopos),'-.r')
title('rental rate')
%axis([0 20 -.5 .5]) 
subplot(4,4,13)
plot(jj,YP1(:,mcpos),jj,YP2(:,mcpos),'-.r')
title('Marginal cost')
%axis([0 20 -.5 .5]) 
subplot(4,4,14)
plot(jj,YP1(:,xpos),jj,YP2(:,xpos),'-.r')
title('GDP x')
%axis([0 20 -.5 .5]) 
subplot(4,4,15)
plot(jj,YP1(:,upos),jj,YP2(:,upos),'-.r')
title('Utilization')
%axis([0 20 -.5 .5]) 
subplot(4,4,16)
plot(jj,YP1(:,kpos),jj,YP2(:,kpos),'-.r')
title('Installed K')
%axis([0 20 -.5 .5]) 

