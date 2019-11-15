% impo.m
% Have to specify ires and ishock, the index values for the 
% responding variable and the shock
% Using the solution of the model
% in state space form
% x(t+1) = Ax(t) + Bu(t+1)
% y(t) = Cx(t) + Du(t)
A = [p q;zeros(nz,nk) phi]; %AS;%bigr;
C = [m n];%bigpi;
D = zeros(nx-nk,nu);
B = [zeros(nk,nu);bigpsi];

ishock = 4; % 2 is IS, 3 is MP, 4 is PC, 5 is tech.

npts = 20; % no of points plotted

%
[Y,X]=dimpulse(A,B,C,D,ishock,npts+1);
jj=[0:npts];
%i1 = Y(:,ires);  % column index is the element of y you want


subplot(2,2,1)
plot(jj,Y(:,iy))
title('y response')
%axis([0 20 -.5 .5]) 
 
subplot(2,2,2)
plot(jj,Y(:,idp)) 
title('dp response')
%axis([0 20 -.25 .25]) 
 
subplot(2,2,3)
plot(jj,Y(:,iygap))
title('ygap response')
%axis([0 20 -.5 .5]) 

subplot(2,2,4)
plot(jj,Y(:,iR))
title('R response')
%axis([0 20 -1.5 1.5]) 
 
   

gtext('FIGURE 1: RESPONSES TO UNIT SHOCK TO PC ') 
%end