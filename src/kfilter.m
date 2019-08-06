function [shatnew,signew,loglh]=kfilter(y,H,C,F,shat,sig,R,Q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL
% y(t)=C+H*s(t)+e(t)
% s(t)=F*s(t-1)+v(t)
% V(e(t))=R
% V(v(t))=Q
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(y);
omega=F*sig*F'+Q;
sigma=H*omega*H'+R;
k=omega*H'/sigma; 

shatnew=F*shat+k*(y-C-H*F*shat);
signew=omega-k*H*omega;

loglh=-.5*log(det(sigma))-.5*((y-C-H*F*shat)'/sigma)*(y-C-H*F*shat)-.5*n*log(2*pi);
