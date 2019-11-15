%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date: August 13, 2019
% This is the main code that controls the estimation of the baseline DSGE
% model in "..." by Bilbiie, Primiceri and Tambalotti
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

addpath('gensys')
addpath('csminwel')

% loading the relevant dataset
load DataTHANK;
y = DataTHANK;
T = length(y);

% initial guess for maximization algorithm
guess = [0.25 0.21 0.15 0.53 0.85 0.25 0.12 0.50 0.75 0.10 4.00 0.85 0.75 5.00 ... 
         2.50 2.00 0.10 0.25 0.80 0.25 0.98 0.70 0.95 0.98 0.70 0.15 0.75 0.95 ... 
         0.20 0.90 0.35 5.00 0.15 0.20 0.05 ...
         0.25 .98 0 0.1 0.1]; % New parameters
     
poscal=[38 40]';
     
x0 = boundsINV(guess);

% posterior maximization
[fh, xh, gh, H, tct, fcount, retcodeh] = csminwel('logpostTHANK', x0, ... 
                          0.1*10 * eye(length(x0)), [], 10e-5, 1000, poscal, T, y);

% processing the output of the maximization
postmode = bounds(xh);
logpost = logpostTHANK(xh',poscal,T,y);

% extracting the states
SmoothedStates=SmoothedStatesJPT(postmode',T,y);
disturbances=SmoothedStates(2:end,48:54);

% plotting inequality and risk
figure
subplot(2,2,1); plot([1954.5:.25:2004.75],SmoothedStates(2:end,45));
title('inequality')
subplot(2,2,2); plot([1954.5:.25:2004.75],bandpass(SmoothedStates(2:end,45),6,32)); hold on
subplot(2,2,2); plot([1954.5:.25:2004.75],bandpass(y(:,4),6,32));
legend('ineqiality','hours')
title('detrended inequality and hours')

subplot(2,2,3); plot([1954.5:.25:2004.75],SmoothedStates(2:end,46));
title('risk')
subplot(2,2,4); plot([1954.5:.25:2004.75],bandpass(SmoothedStates(2:end,46),6,32)); hold on
subplot(2,2,4); plot([1954.5:.25:2004.75],bandpass(y(:,4),6,32)/100);
legend('risk','hours/100')
title('detrended risk and hours/100')

% reconstructing the observables
[G1,C,impact,eu,SDX,zmat,NY,NX]=modelTHANKcycle(postmode');
States=SmoothedStates(1,:)';
STATES=zeros(T,NY);
detOBS=zeros(T,size(zmat,1));
for t=1:T
    States=G1*States+impact*disturbances(t,:)';
    STATES(t,:)=States;
    detOBS(t,:)=zmat*States;
end

% reconstructing the observables under counterfactual scenarios
% TANK
pv=postmode'; pv(37)=1;
[G1,C,impact,eu,SDX,zmat,NY,NX]=modelTHANKcycle(pv);
States=SmoothedStates(1,:)';
STATEStank=zeros(T,NY);
detOBStank=zeros(T,size(zmat,1));
for t=1:T
    States=G1*States+impact*disturbances(t,:)';
    STATEStank(t,:)=States;
    detOBStank(t,:)=zmat*States;  
end

% reconstructing the observables under counterfactual scenarios
% RANK
pv=postmode'; pv(37)=1; pv(36)=.01; pv(38)=0; pv(39)=0; pv(40)=0;
[G1,C,impact,eu,SDX,zmat,NY,NX]=modelTHANKcycle(pv);
States=SmoothedStates(1,:)';
STATESrank=zeros(T,NY);
detOBSrank=zeros(T,size(zmat,1));
for t=1:T
    States=G1*States+impact*disturbances(t,:)';
    STATESrank(t,:)=States;
    detOBSrank(t,:)=zmat*States;
end

figure
subplot(2,2,1); plot([1954.5:.25:2004.75],[detOBS(:,1) detOBStank(:,1) detOBSrank(:,1)]);
title('GDP growth')
legend('THANK','TANK','RANK')
subplot(2,2,2); plot([1954.5:.25:2004.75],[detOBS(:,2) detOBStank(:,2) detOBSrank(:,2)]);
title('C growth')
subplot(2,2,3); plot([1954.5:.25:2004.75],[detOBS(:,3) detOBStank(:,3) detOBSrank(:,3)]);
title('I growth')
subplot(2,2,4); plot([1954.5:.25:2004.75],[detOBS(:,4) detOBStank(:,4) detOBSrank(:,4)]);
title('log hours')

figure
subplot(1,2,1); plot([1954.5:.25:2004.75],[STATES(:,45) STATEStank(:,45)]);
title('inequality')
legend('THANK','TANK')
subplot(1,2,2); plot([1954.5:.25:2004.75],[bandpass(STATES(:,45),6,32) bandpass(STATEStank(:,45),6,32)]);
title('detrended inequality')
legend('THANK','TANK')



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% JJ       = jacobTHANK(xh);
% HH       = JJ * H * JJ';
% 
% % preliminaries for the MCMC algorithm
% const = .4;    % scaling constant for the inverse Hessian
% M     = 2000;  % length of each of the two chains
% N     = 200;   % number of discarded draws at the beginning of each chain
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Metropolis algorithm
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % First chain
% P1 = zeros(M,length(xh));
% logpostOLD = -1e+10;
% while logpostOLD == -1e+10
%     P1(1,:) = mvnrnd(postmode, 4 * HH * const^2, 1);
%     logpostOLD = logpostTHANK_MCMC(P1(1,:), T, y);
% end
% 
% count = 0;
% for i=2:M
%     if i==100*floor(.01*i)
%         i
%         ACCrate1 = count/i
%     end
%     P1(i,:)=mvnrnd(P1(i-1,:),HH*const^2,1);
%     logpostNEW=logpostTHANK_MCMC(P1(i,:),T,y);
%     
%     if logpostNEW > logpostOLD
%         logpostOLD=logpostNEW;
%         count=count+1;
%     else
%         if rand(1) < exp(logpostNEW - logpostOLD)
%             logpostOLD = logpostNEW;
%             count = count+1;
%         else
%             P1(i,:) = P1(i-1,:);
%         end
%     end
% end
% ACCrate1=count/M;
% 
% % Second chain
% P2=zeros(M,length(xh));
% logpostOLD=-1e+10;
% while logpostOLD == -1e+10
%     P2(1,:) = mvnrnd(postmode, 4 * HH * const^2, 1);
%     logpostOLD = logpostTHANK_MCMC(P2(1,:), T, y);
% end
% 
% % Metropolis algorithm
% count=0;
% for i=2:M
%     if i == 100*floor(.01*i)
%         i
%     end
%     P2(i,:) = mvnrnd(P2(i-1,:), HH * const^2, 1);
%     logpostNEW = logpostTHANK_MCMC(P2(i,:), T, y);
%     
%     if logpostNEW > logpostOLD
%         logpostOLD = logpostNEW;
%         count = count+1;
%     else
%         if rand(1) < exp(logpostNEW - logpostOLD)
%             logpostOLD = logpostNEW;
%             count = count+1;
%         else
%             P2(i,:) = P2(i-1,:);
%         end
%     end
% end
% ACCrate2 = count/M;
% 
% % Third chain
% P3 = zeros(M,length(xh));
% logpostOLD = -1e+10;
% while logpostOLD == -1e+10
%     P3(1,:) = mvnrnd(postmode,4*HH*const^2,1);
%     logpostOLD = logpostTHANK_MCMC(P3(1,:),T,y);
% end
% 
% % Metropolis algorithm
% count=0;
% for i=2:M
%     if i == 100*floor(.01*i)
%         i
%     end
%     P3(i,:) = mvnrnd(P3(i-1,:),HH*const^2,1);
%     logpostNEW = logpostTHANK_MCMC(P3(i,:),T,y);
%     
%     if logpostNEW > logpostOLD
%         logpostOLD = logpostNEW;
%         count = count + 1;
%     else
%         if rand(1) < exp(logpostNEW-logpostOLD) 
%             logpostOLD = logpostNEW;
%             count = count+1;
%         else
%             P3(i,:) = P3(i-1,:);
%         end
%     end
% end
% ACCrate3=count/M;
% 
% P = [P1(N+1:end,:); P2(N+1:end,:); P3(N+1:end,:)];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Processing output
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MED=median(P)';
% prctile5=prctile(P,5)';
% prctile95=prctile(P,95)';
% [MED prctile5 prctile95]           % parameter estimates (median and 90% CI)
