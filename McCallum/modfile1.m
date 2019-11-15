%This uses solvek to solve the "standard" macro model.

nx = 10; nz = 5; nu = nz;

b0 = 0; b1 = -.5; beta = .99;  alpha = .02;
mu0 = 0; mu1 = 1.99; mu2 = .0; mu3 = 0.8;

A = zeros(nx,nx);  B = zeros(nx,nx);  C = zeros(nx,nu);
phi = zeros(nu,nu);

iy = 1;
idp = 2;
iey1 = 3;
iR = 4;
iygap = 5;
ie1y = 6;
iRlag = 7;
ie1ygap = 8;
idplag = 9;
iygaplag = 10;

nk = 5;

icon = 1;
iv = 2;
ie = 3;
iu = 4;
iybar = 5;

%Define E(t)y(t+1)
A(1,iy) = 1;
B(1,iey1) = 1;

%Define E(t-1)y(t)
A(2,ie1y) = 1;
B(2,iey1) = 1;

%IS function
B(3,iy) = -1;
B(3,iR) = b1;
A(3,idp) = b1;
B(3,iey1) = 1;
C(3,iv) = 1;
C(3,icon) = b0;

%Price adjustment function
B(4,idp) = -1;
A(4,idp) = -0.5; %beta;
B(4,idplag) = 0.5;
B(4,iygap) = alpha;
C(4,iu) = 1;

%Policy rule
B(5,iR) = -1;
B(5,idp) = (1-mu3)*(mu1); %could enter mu1 differently
B(5,ie1ygap) = (1-mu3)*mu2; %should divide mu2 by 4
B(5,iRlag) = mu3;
C(5,ie) = 1;
C(5,icon) = mu0;

%Define output gap
B(6,iygap) = -1;
B(6,iy) = 1;
C(6,iybar) = -1;

%Define R(t-1)
A(7,iRlag) = 1;
B(7,iR) = 1;

%Define E(t-1)ygap
A(8,ie1ygap) = 1;
A(8,iygap) = -1;

%Define dp(t-1)
A(9,idplag) = 1;
B(9,idp) = 1;

%Define ygap(t-1)
A(10,iygaplag) = 1;
B(10,iygap) = 1;


phi(icon,icon) = 1;
phi(iv,iv) = 0.;
phi(iu,iu) = 0.;
phi(iybar,iybar) = .95;

   [m,n,p,q,z22h,s,t,lambda] = solvek(A,B,C,phi,nk);

bigmn = [m n];
bigpq = [p q];
%bigp = phi;
bigpsi = eye(nz,nu);

%end;