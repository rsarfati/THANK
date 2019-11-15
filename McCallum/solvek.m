function [m,n,p,q,z22h,s,t,lambda] = solvek(a,b,c,phi,nk)
% function [m,n,p,q,z22h,s,t,lambda] = solvek(a,b,c,phi,nk)
%
% A*E[x(t+1)|F(t)] = B*x(t) + C*z(t)
% where x(t) = [y(t)' k(t)'], and the nk-dimensional
% vector k(t) is predetermined, while
% E[z(t+1)|F(t)] = phi*z(t)
% 
% This function produces the Markov decision rule
% 
% y(t) = M*k(t) + N*z(t)
% k(t+1) = P*k(t) + Q*z(t)
%  (second line of above corrected EN May 27 1997)
% The supporting function qzdiv checks whether the matrix pencil
% b-za is singular.
%
% Written by Paul Klein April 1997, modified by B.McCallum and E. Nelson
% The algorithm follows Paul Klein: "Using the generalized Schur form
% to solve a system of linear expectational difference equations".
%
nx = size(a,1);
nz = size(c,2);
nd = nx-nk;
%
% Do the generalized Schur decomposition
%
[t,s,z,q,v] = qz(b',a');
[s,t,z,q] = reorder(s,t,z,q);%replacing qzdiv Apr 27 98
%[t,s,z,q] = qzdiv(1,t,s,z,q);%check whether l or 1
s = s';
t = t';
q = q';
zh = z';
%
%At this stage we have qaz' = s and qbz' = t where q and z are unitary,
% s and t are lower triangular, and the pairs (s(i,i),t(i,i)) are
% ordered in descending absolute value.
%
% Checking whether the no. of unstable eigenvalues is equal to nd 
global lambda alph bet
alpha = diag(s);
beta = diag(t);
alph = alpha;
bet = beta;
alpha = alpha + eps;
lambda = beta./alpha;
%
eigviol = 0;
unstable = lambda(abs(lambda)>1.000001);
stable = lambda(abs(lambda)<0.99999);
if (size(stable,1)>nk);
 disp('Multiple stable solutions');
end
if (size(unstable,1)>nd);
 disp('No stable solution');
end
%if ~(size(unstable,1)==nd);
% disp('Wrong number of stable eigenvalues.');
%  eigviol = 1;
%end
%
% Tranforming the system and the variables; partitioning
d = q*c;
d1 = d(1:nd,:);
d2 = d(nd+1:nx,:);
%
z11 = z(1:nd,1:nd);
z12 = z(1:nd,nd+1:nx);
z21 = z(nd+1:nx,1:nd);
z22 = z(nd+1:nx,nd+1:nx);
%
z11h = zh(1:nd,1:nd);
z12h = zh(1:nd,nd+1:nx);
z21h = zh(nd+1:nx,1:nd);
z22h = zh(nd+1:nx,nd+1:nx); 
%
s11 = s(1:nd,1:nd);
s21 = s(nd+1:nx,1:nd);
s22 = s(nd+1:nx,nd+1:nx);
%
t11 = t(1:nd,1:nd);
t21 = t(nd+1:nx,1:nd);
t22 = t(nd+1:nx,nd+1:nx);
%
% Calculating the decision rule for the transformed variables
%
vecl = (kron(phi',s11)-kron(eye(nz),t11))\d1(:);
l = reshape(vecl,nd,nz);
%
if rank(z22h) < nk;
   error('Rank condition not satisfied');
end
%
p = s22\t22;
r = s22\(t21*l-s21*l*phi+d2)+z22h\z21h*l*phi;
n = -z22h\z21h*l;
%
% Calculating the decision rule for the original variables
%
pre = inv(eye(nd)-z12h*z21);
mtilde = pre*z12h*z22;
ltilde = pre*z11h*l;
%
ptilde = z22h*p*(z21*mtilde+z22);
rtilde = z22h*(p*z21*ltilde+r);
m = real(mtilde);
n = real(ltilde);
%
p = real(ptilde);
q = real(rtilde);

%m = mtilde;
%n = ltilde;
%p = ptilde;
%q = rtilde;

