% impoch2.m
% Have to specify ires and ishock, the index values for the 
% responding variable and the shock
% Using the solution of the model
% in state space form
% x(t+1) = Ax(t) + Bu(t+1)
% y(t) = Cx(t) + Du(t)
A = bigr;
C = bigpi;
D = zeros(nx-nk,nu);
B = [zeros(nk,nu);bigpsi];
%
ishock = 3;
%
[Y,X]=dimpulse(A,B,C,D,ishock,21);
jj=[0:20];
%i1 = Y(:,ires);  % column index is the element of y you want
%plot(jj,i1(1:21,:)),
%xlabel('Quarters after shock'),
%ylabel('Percent deviation'),
%title('Response of y to shock')


subplot(3,2,1)
plot(jj,Y(:,iy),'w')
ylabel('y')
%title('y response')
%axis([0 20 -.2 1.])
 
subplot(3,2,2)
plot(jj,Y(:,iq),'w') 
ylabel('q')
%title('q response')
%axis([0 20 -1 4])
 
subplot(3,2,6)
plot(jj,Y(:,inetexp),'w')
%plot(jj,X(:,1),'w')
ylabel ('net')
%title('   p response')
%axis([0 20 -2 .5]) 

subplot(3,2,4)
plot(jj,Y(:,is),'w')
ylabel('s')
%title('s response')
%axis([0 20 -1 4])
 
subplot(3,2,3)
plot(jj,Y(:,idp),'w')
ylabel('delta p')
%title('dp response')
%axis([0 20 -.2 .2]) 
 
subplot(3,2,5)
plot(jj,Y(:,iR),'w')
ylabel('R')
%title('R response'),
%axis([0 20 -.5 .5]), 
   
gtext('Responses to Unit Shock to MP, mu1 = .5, mu2 = .4, mu3 = .8') 
end