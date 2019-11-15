function  [yout,x] = dlsim(a, b, c, d, u, x0)
%DLSIM	Simulation of discrete-time linear systems.
%	DLSIM(A,B,C,D,U)  plots the time response of the discrete system:
%
%		x[n+1] = Ax[n] + Bu[n]
%		y[n]   = Cx[n] + Du[n]
%
%	to input sequence U.  Matrix U must have as many columns as there
%	are inputs, u.  Each row of U corresponds to a new time point.
%	DLSIM(A,B,C,D,U,X0) can be used if initial conditions exist.
%
%	DLSIM(NUM,DEN,U) plots the time response of the transfer function
%	description  G(z) = NUM(z)/DEN(z)  where NUM and DEN contain the 
%	polynomial coefficients in descending powers of z.  If 
%	LENGTH(NUM)=LENGTH(DEN) then DLSIM(NUM,DEN,U) is equivalent to 
%	FILTER(NUM,DEN,U).  When invoked with left hand arguments,
%		[Y,X] = DLSIM(A,B,C,D,U)
%		[Y,X] = DLSIM(NUM,DEN,U)
%	returns the output and state time history in the matrices Y and X.
%	No plot is drawn on the screen.  Y has as many columns as there 
%	are outputs and LENGTH(U) rows.  X has as many columns as there 
%	are states and LENGTH(U) rows. 
%
%	See also: DSTEP,DIMPULSE,DINITIAL and LSIM.

%	J.N. Little 4-21-85
%	Revised 7-18-88 JNL
%	Revised 7-31-90  Clay M. Thompson
%	Copyright (c) 1986-93 by the MathWorks, Inc.

error(nargchk(3,6,nargin));
if (nargin==4), error('Wrong number of input arguments.'); end

nargs = nargin;

if (nargs == 3)		% transfer function description 
  [num,den] = tfchk(a,b);
  u = c(:);
  [m,n] = size(num);
  if ((m == 1)&(nargout == 1))	% Use filter, it's more efficient
    y = filter(num,den,u);
    yout = y; 
    return
  else  			% Convert to state space
    [a,b,c,d] = tf2ss(num,den);
    nargs = 5;
  end
end

[ny,nu] = size(d);
if ny*ny==0, x = []; if nargout~=0, yout = []; end, return, end

[ns,nx] = size(a);
if (nargs == 5)
  x0 = zeros(1,ns);
end
error(abcdchk(a,b,c,d));
if min(size(u)) == 1
  u = u(:);	% Make sure u is a column vector
end

% Make sure u is has the right number of columns
[n,mu]=size(u);
if mu~=nu, error('U must have same number of columns as inputs.'); end

if isempty(a),
  x = [];
  y = u * d.';
else
  x = ltitr(a,b,u,x0);
  y = x * c.' + u * d.';
end

if nargout==0,	% With no output arguments, plot graph
  t = [0:length(y)-1];
  stairs(t,y)
  xlabel('No. of Samples'), ylabel('Amplitude')

  return  % Suppress output
end

yout = y; 

