% sim33p.m   (Written by E. Nelson, B. McCallum, and C. Jensen)
% File for simulating models that have been solved with solvek. The
% output is standard deviations,covariances, and autocorrelations of
% selected endogenous variables (uses McGrattan's file autocor.)
% The user is assumed to have run solvek. The variables for which statistics
% are desired must be specified in the Y = [     ];    statement below. There
% are nvstat of them and their names should be entered. 
% Also, nlags specifies the number of lags in autocorrelation calculation.
% To simulate the model one puts the solution in state space form:
% s(t+1) = AS*s(t) + BS*u(t+1)
% y(t) = CS*s(t) + DS*u(t) 
% where s(t) is the state vector [k(t)' z(t)']', k(t)
% is the vector of predetermined endogenous variables,
% and z(t) is the vector of exogenous variables.
% Maximum value of nz is 6; minimum is 2.
% Autocorrelations are plotted; can print values also if desired
% Including more variable and more lags slows the calculation.


nlags = 2;  % nlags is the number of autocorrelations calculated
noplot = 0; % setting noplot = 0 will suppress plot of autocorrelations
auto = 1;   % setting auto = 1 will give printout of autocorr values
nvstat = 4; namev = ('      dp           y            ygap           R   ');
% nvstat is the number of variables with calculated statistics and
% namev lists their names. They must also be entered in Y defn below!!!

sig1 = 0.00;    % std dev of inn1
sig2 = 0.02;   % std dev of inn2
sig3 = 0.005;  % std dev of inn3
sig4 = 0.005;    % std dev of inn4
sig5 = 0.007;    % std dev of inn5
sig6 = 0.00;    % std dev of inn6

%...................................
nu = nz;
bigpsi = eye(nz,nu);
AS = [p q;zeros(nz,nk) phi];
CS = [m n];
BS = [zeros(nk,nu);bigpsi];  
%
[nrCS,ncCS] = size(CS);
DS = zeros(nrCS,nu);
%
randn('seed',sum(100*clock));
nn = 253;%length of each sample simulated
nsim = 100;% no. of simulations

inn1a = randn(nn,nsim);
inn1 = sig1*inn1a;
%
inn2a = randn(nn,nsim);
inn2 = sig2*inn2a;
%
inn3a = randn(nn,nsim);
inn3 = sig3*inn3a;
%
inn4a = randn(nn,nsim);
inn4 = sig4*inn4a;
%
inn5a = randn(nn,nsim);
inn5 = sig5*inn5a;
%
inn6a = randn(nn,nsim);
inn6 = sig6*inn6a;

AAA = zeros(nvstat*(nlags+1),nvstat);
BBB = zeros(nvstat,nvstat);

for i = 1:nsim;
if nu==2; inns = [inn1(:,i),inn2(:,i)]; else;
if nu==3; inns = [inn1(:,i),inn2(:,i),inn3(:,i)]; else;
if nu==4; inns = [inn1(:,i),inn2(:,i),inn3(:,i),inn4(:,i)]; else;
if nu==5; inns = [inn1(:,i),inn2(:,i),inn3(:,i),inn4(:,i),inn5(:,i)]; else;
 inns = [inn1(:,i),inn2(:,i),inn3(:,i),inn4(:,i),inn5(:,i),inn6(:,i)];
end; end; end; end; 
%
[ysim,stsim] = dlsim(AS,BS,CS,DS,inns);

%Y = ysim;%use this for all endog vars or next for selected

Y = [ysim(:,idp),ysim(:,iy),ysim(:,iygap),ysim(:,iR)];

y99 = [Y(51:250,:)];
Atemp = autocor(y99,nlags);
AAA = AAA + Atemp;
Btemp = cov(y99);
BBB = BBB + Btemp;
end;

AAA = AAA/nsim;
BBB = BBB/nsim;
CCC = sqrt(diag(BBB));
disp('standard deviations');
format short e
disp(CCC');
disp(namev);
disp(' ');
disp('covariances');
format short e;
disp(BBB);
format short;
[N,KK] = size(AAA);

if auto == 1 
disp('autocorrelations');   

for i = 1:nlags + 1;
disp(['       correlations at lag ' int2str(i-1)]);
disp(AAA(((i-1)*KK)+1:i*KK,:));
end;
else; end;

nl = N/KK;

if nlags == 0 break; end;
if noplot == 0 break; end;

B= zeros(nl,KK*KK);
jj = [0:nl-1];
for i = 1:KK;
for j = 1:KK;
for k = 1:nl;
B(k,(i-1)*KK+j) = AAA((k-1)*KK+i,j);
end
subplot(KK,KK,(i-1)*KK+j)
plot(jj,B(:,(i-1)*KK+j))
axis([0 nl-1 -1 1])
end
end



gtext('Autocorrelation functions for dp, y, ygap, R ')