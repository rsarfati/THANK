%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date: December 22, 2009
% This is the main code that controls the estimation of the baseline DSGE
% model in "Investment Shocks and Business Cycles," by
% Justiniano, Primiceri and Tambalotti
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% loading the relevant dataset
load DataJPT;
y=DataJPT;
T=length(y);


% initial guess for maximization algorithm
%guess=[0.25 0.21 0.15 0.53 0.85 0.25 0.12 0.5 0.75 0.1 4 0.85 0.75 5 2.5 2 0.1 0.25 0.8 0.25  0.98 0.7 0.95 0.98 0.7 0.15 0.75 0.95 0.2 0.9 0.35 5 0.15 0.2 0.05];

% ORIGINAL PAPER POSTMODE
guess=[0.1680    0.2087    0.1118    0.4836    0.7848    0.2361    0.1467    0.4159    0.7045    0.1147    3.6507    0.8315    0.6966    5.1221    2.5032    2.0552    0.0665    0.2377...
       0.8186    0.2327    0.9900    0.7270    0.9413    0.9837    0.6574    0.1368    0.7386    0.9380    0.2208    0.8740    0.3499    5.3592    0.1367    0.2127    0.0376    0.18];

%guess = [0.1675    0.2082    0.1011    0.4822    0.7464    0.2240    0.1442    0.4076    0.7006    0.1144    3.4775    0.8257    0.7062    5.0589    2.6513    1.9974    0.0635...
%    0.2470    0.8129    0.2868    0.9900    0.6909    0.9401    0.9835    0.6947    0.1367    0.7259    0.9433    0.2212    0.8737    0.3499    5.8661    0.1385    0.2152...
%    0.0358    0.1800];
   
x0=boundsINV(guess);

% posterior maximization
[fh,xh,gh,H,itct,fcount,retcodeh] = csminwel('logpostJPT',x0,.1*eye(length(x0)),[],10e-4,1000,T,y);


% processing the output of the maximization
postmode=bounds(xh);
JJ=jacobJPT(xh);
HH=JJ*H*JJ';


% preliminaries for the MCMC algorithm
const=.4;        % scaling constant for the inverse Hessian
M=2000;          % length of each of the two chains
N=200;          % number of discarded draws at the beginning of each chain


% Metropolis algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First chain
P1=zeros(M,length(xh));
logpostOLD=-1e+10;
while logpostOLD==-1e+10;
    P1(1,:)=mvnrnd(postmode,4*HH*const^2,1);
    logpostOLD=logpostJPT_MCMC(P1(1,:),T,y);
end

count=0;
for i=2:M
    if i==100*floor(.01*i);
        i
        ACCrate1=count/i
    end
    P1(i,:)=mvnrnd(P1(i-1,:),HH*const^2,1);
    logpostNEW=logpostJPT_MCMC(P1(i,:),T,y);
    if logpostNEW>logpostOLD
        logpostOLD=logpostNEW;
        count=count+1;
    else
        if rand(1)<exp(logpostNEW-logpostOLD);
            logpostOLD=logpostNEW;
            count=count+1;
        else
            P1(i,:)=P1(i-1,:);
        end
    end
end
ACCrate1=count/M;

% Second chain
P2=zeros(M,length(xh));
logpostOLD=-1e+10;
while logpostOLD==-1e+10;
    P2(1,:)=mvnrnd(postmode,4*HH*const^2,1);
    logpostOLD=logpostJPT_MCMC(P2(1,:),T,y);
end
% Metropolis algorithm
count=0;
for i=2:M
    if i==100*floor(.01*i);
        i
    end
    P2(i,:)=mvnrnd(P2(i-1,:),HH*const^2,1);
    logpostNEW=logpostJPT_MCMC(P2(i,:),T,y);
    if logpostNEW>logpostOLD
        logpostOLD=logpostNEW;
        count=count+1;
    else
        if rand(1)<exp(logpostNEW-logpostOLD);
            logpostOLD=logpostNEW;
            count=count+1;
        else
            P2(i,:)=P2(i-1,:);
        end
    end
end
ACCrate2=count/M;

% third chain
P3=zeros(M,length(xh));
logpostOLD=-1e+10;
while logpostOLD==-1e+10;
    P3(1,:)=mvnrnd(postmode,4*HH*const^2,1);
    logpostOLD=logpostJPT_MCMC(P3(1,:),T,y);
end

% Metropolis algorithm
count=0;
for i=2:M
    if i==100*floor(.01*i);
        i
    end
    P3(i,:)=mvnrnd(P3(i-1,:),HH*const^2,1);
    logpostNEW=logpostJPT_MCMC(P3(i,:),T,y);
    if logpostNEW>logpostOLD
        logpostOLD=logpostNEW;
        count=count+1;
    else
        if rand(1)<exp(logpostNEW-logpostOLD);
            logpostOLD=logpostNEW;
            count=count+1;
        else
            P3(i,:)=P3(i-1,:);
        end
    end
end
ACCrate3=count/M;

P=[P1(N+1:end,:);P2(N+1:end,:);P3(N+1:end,:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Processing output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MED=median(P)';
prctile5=prctile(P,5)';
prctile95=prctile(P,95)';
[MED prctile5 prctile95]           % parameter estimates (median and 90% CI)
