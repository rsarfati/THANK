% chile4.m; this is buba2 withl net export balance
% It uses Fuhrer-Moore pricing. 

nx = 24; 
nz = 7;  nu = nz;
 
A = zeros(nx,nx); B = zeros(nx,nx); C = zeros(nx,nz);
 phi = zeros(nz,nz);
 
b1 = 1.; b2 = 0.333; b3 = 1.; b4 = .333;  %b4 is elas of sub [b1, b2 for exports
  
sh1 = .5; sh2 = .25; sh3 = 0.25; %sh2 is fraction of output exported; sh3 is govt. 
c2 =  sh2; gamma = b4*(-c2/(1-c2))*1.2;


 
theta = 6;%implies markup of 6/(6-5) = 1.2
beta = .99;
h = 0.8; sif = 2.5; %note: sif is inverse of sigma
B1 = h*(1 - sif);  B2 = beta*h*(1 + h*(1-sif)) - sif;
 
iy = 1;
iybar = 2;
iygap = 3;
iR = 4;
iq = 5;
is = 6;
ic = 7;
ilam = 8;
iexp = 9;
idx = 10;
ip = 11;
idp = 12;
iedp1 = 13; %     Etdp(t+1)
inetexp = 14;
iclag = 15;    
iRlag = 16;
iylag = 17;
ie1dx = 18;%Et-1dx
ie1ygap = 19; %Et-1ygap
ie1y = 20; 
idplag = 21;
iplag = 22;
ie1dp1 = 23; % Et-1dp(t+1)
ie1dp = 24;  %Et-1dp(t)


nk = 10;
 
iasc = 1;%scaled technology shock
iv = 2;
ier = 3;
ixi = 4;
iystar = 5;
iep = 6;
igov = 7;
 
%Nat income identity 
B(1,iy) = -1;
B(1,ic) = sh1;
B(1,iexp) = sh2;
C(1,igov) = sh3;
 
%defn  ygap = y - ybar
B(2,iygap) = -1;
B(2,iy) = 1;
B(2,iybar) = -1;
 
% pbar eqn  E(t)ygap(t+1) = alpha*ygap(t)
%A(3,iygap) = 1;
%B(3,iygap) = alpha;

%Fuhrer-Moore pricing eqn
B(3,idp) = -1;
B(3,idplag) = .5; 
B(3,iedp1) = .5;
B(3,iygap) = .02;
C(3,iep) = 1;


mu1 =.50;
mu2 = 0.4; 
mu3 = 0.8;
mu4 = .0;

% Policy rule  R(t) = 
B(4,iR) = -1;
%B(4,ie1dp1) = (1-mu3);  %*(1+mu1);
B(4,ie1dp) = (1-mu3)*(mu1+1);
%B(4,idp) = (1-mu3)*(mu1+1);
B(4,ie1ygap) = (1-mu3)*mu2/4; 
C(4,ier) = 1;
B(4,iRlag) = mu3;
B(4,ie1dx) = (1-mu3)*mu4;
 
% ybar eqn ybar(t) = asc(t) + gamma*q(t), yb is itself RW domestic part
B(5,iybar) = -1;
C(5,iasc) = 1;
B(5,iq) = gamma;

% defn p(t-1)
A(20,iplag) = 1;
B(20,ip) = 1;
 
%defn p(t) = dp(t) + p(t-1)
B(6,ip) = -1;
B(6,idp) = 1;
B(6,iplag) = 1;
 
%defn q(t) = s(t) - p(t) + pstar(t)
B(7,iq) = -1;
B(7,is) = 1;
B(7,ip) = -1;
%B(7,ipstar) = 1;
 
%UIP eqn  R(t) = Rstar(t) + E(t)s(t+1) - s(t) + xi(t)
A(8,is) = -1;
B(8,iR) = -1;
%B(8,Rstar) = 1;
B(8,is) = -1;
C(8,ixi) = 1;
 
%Cons. fn  beta*B1*Etc(t+1) = B2*c(t) - B1*c(t-1) - (1-h*beta)*lam(t)
A(9,ic) = beta*B1;
B(9,ic) = B2;
B(9,iclag) = -B1;
B(9,ilam) = -(1-h*beta);
C(9,iv) = 1;
 
%define clag
A(10,iclag) = 1;
B(10,ic) = 1;
 
% define Rlag
A(11,iRlag) = 1;
B(11,iR) = 1;
 
% lambda eqn  Etlam(t+1) + R(t) - Etdp(t+1) = lam(t)
A(12,ilam) = .9999;
A(12,idp) = -1;
B(12,ilam) = 1;
B(12,iR) = -1;
 
%log of  exports  net(t) = b1*ystar + b2*q  
%                            
B(13,iexp) = -1;
B(13,iq) = b2;
C(13,iystar) = b1;
 
% defn of dx(t): dx(t) = dp(t) + y(t) - y(t-1)
B(14,idx) = -1;
B(14,idp) = 1;
B(14,iy) = 1;
B(14,iylag) = -1;
 
%defn xlag
%A(15,ixlag) = 1;
%B(15,ix) = 1;
 
%defn ylag
A(15,iylag) = 1;
B(15,iy) = 1;
 
%defn Et-1dx(t)
A(16,ie1dx) = 1;
A(16,idx) = -1;

%defn Et-1ygap(t)
A(17,ie1ygap) = 1;
A(17,iygap) = -1;

%defn Et-1y(t)
A(18,ie1y) = 1;
A(18,iy) = -1;

%defn dp(t-1)
A(19,idplag) = 1;
B(19,idp) = 1;

%define Etdp(t+1)
B(21,iedp1) = 1;
A(21,idp) = 1;

%define Et-1dp(t+1)
A(22,ie1dp1) = 1;
A(22,iedp1) = -1;

%define Et-1dp(t)
A(23,ie1dp) = 1;
B(23,iedp1) = 1;

%define netexp
B(24,inetexp) = -1;
B(24,iexp) = 1;
B(24,iy) = -b3;
B(24,iq) = b4-1;


 
%AR parameters for exogenous variables yb, v, er, xi, inv, ep
phi(1,iasc) = .95;
phi(2,iv) = .0;
phi(3,ier) = .0;
phi(4,ixi) = 0.5;
phi(5,iystar) = .950;
phi(6,iep) = 0;
phi(7,igov) = .99;
 
%bigp = phi;
[m,n,p,q,z22h,s,t,lambda] = solvek(A,B,C,phi,nk);
%[m,n,p,q] = solvek(A,B,C,phi,nk);
bigmn = [m,n];
bigpq = [p,q];

bigpsi = zeros(nz,nu);
bigpsi(1,1) = 1;
bigpsi(2,2) = 1;
bigpsi(3,3) = 1;
bigpsi(4,4) = 1; 
bigpsi(5,5) = 1;
bigpsi(6,6) = 1;
bigpsi(7,7) = 1;
end; 
 
 bigr = [p q;zeros(nz,nk) phi];
bigpi = [m n];
 
 
