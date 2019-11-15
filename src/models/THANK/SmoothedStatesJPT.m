function SmoothedStates=SmoothedStatesJPT(param,T,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date: November 1, 2015
% this function computes the smoothed states for the model
% in Justiniano, Primiceri and Tambalotti (2010) for a given parameter
% value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solution of the model
[G1,C,M,eu,SDX,H,NY,NX]=modelTHANKcycle(param);

% initializing the Kalman filter
shat0=zeros(NY,1); shat=shat0;
Q=M*SDX*SDX'*M';
sig0=disclyap_fast(G1,Q); sig=sig0;
SHAT=zeros(NY,T);
SIG=zeros(NY,NY,T);

% Kalman filter recursion delivering log-likelihood
for t=1:T
    [shat,sig,loglh]=kfilter(y(t,:)',H,C,G1,shat,sig,0,Q);
    SHAT(:,t)=shat;
    SIG(:,:,t)=sig;
end

% Kalman smoother
SmoothedStates=zeros(NY,T); SmoothedStates(:,T)=shat;
bnT=shat';
SnT=sig;
for t=T-1:-1:1
    [bnT,SnT]=ksmooth_const_pseudo(SHAT(:,t)',squeeze(SIG(:,:,t)),bnT,SnT,G1,zeros(NY,1),Q);
    SmoothedStates(:,t)=bnT;
end
[bnT,junk]=ksmooth_const_pseudo(shat0',sig0,bnT,SnT,G1,zeros(length(G1),1),Q);
SmoothedStates=[bnT' SmoothedStates]';

% extracting the shock processes and the innovations
% processes=SmoothedStates(:,17:23);
% shocks=processes(2:end,:)-(diag(postmode(20:26))*processes(1:end-1,:)')';
% SmoothedStates(:,35:36);

    


function [btT,StT]=ksmooth_const_pseudo(btt,Stt,bnT,SnT,A,C,omega)
%[btT StT]=ksmooth(btt,Stt,bnT,SnT,A,omega)
% Smoothing recursion.  State evolution equation is
%    bn=A*bt+e,  Var(e)=omega
%    bt|t ~ N(btt,Stt) -- from Kalman Filter
%    bn|T ~ N(bnT,SnT) -- distribution of bn given full sample. From 
%                         KF if n=T, otherwise from this recursion
%    bt|T ~ N(btT,StT)
AS=A*Stt;
G=AS*A'+omega;
[u,d,v]=svd(G);
first0=min(find(diag(d)<1e-5));
if isempty(first0),first0=min(size(G))+1;end
u=u(:,1:first0-1);
v=v(:,1:first0-1);
d=diag(d);d=diag(1./d(1:first0-1)); 
pseudoinvG=v*d*u';
SAGI=AS'*pseudoinvG;
btT=(SAGI*(bnT'-A*btt'-C))'+btt;
StT=Stt-SAGI*AS+SAGI*SnT*SAGI';