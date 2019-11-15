function [yout,x,n] = dimpulse(a,b,c,d,iu,n)
%DIMPULSE Impulse response of discrete-time linear systems.
%	 DIMPULSE(A,B,C,D,IU)  plots the response of the discrete system:
%
%		x[n+1] = Ax[n] + Bu[n]
%		y[n]   = Cx[n] + Du[n]
%
%	to an unit sample applied to the single input IU.  The number of
%	points is determined automatically.
%
%	DIMPULSE(NUM,DEN)  plots the impulse response of the polynomial
%	transfer function  G(z) = NUM(z)/DEN(z)  where NUM and DEN contain
%	the polynomial coefficients in descending powers of z.
%
%	DIMPULSE(A,B,C,D,IU,N) or DIMPULSE(NUM,DEN,N) uses the user-
%	supplied number of points, N.  When invoked with left hand 
%	arguments,
%		[Y,X] = DIMPULSE(A,B,C,D,...)
%		[Y,X] = DIMPULSE(NUM,DEN,...)
%	returns the output and state time history in the matrices Y and X.
%	No plot is drawn on the screen.  Y has as many columns as there 
%	are outputs and X has as many columns as there are states.
%
%	See also: DSTEP,DINITIAL,DLSIM and IMPULSE.

%	J.N. Little 4-21-85
%	Revised CMT 7-31-90, ACWG 5-30-91 
%	Copyright (c) 1986-93 by the MathWorks, Inc.

nargin = 6;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ADDED BY BTM
if nargin==0, eval('exresp(''dimpulse'');'), return,end

error(nargchk(2,6,nargin));

if (nargin==2),		% Transfer function without number of points
	[num,den] = tfchk(a,b);
	iu = 1;
	[a,b,c,d] = tf2ss(num,den);

elseif (nargin==3),	% Transfer function with number of points
	[num,den] = tfchk(a,b);
	n = c;
	iu = 1;
	[a,b,c,d] = tf2ss(num,den);

elseif (nargin>=4)
	error(abcdchk(a,b,c,d));

end

[ny,nu] = size(d);
if (nu*ny==0), x = []; n = []; if nargout~=0, yout=[]; end, return, end

% Work out number of points if not supplied.
if nargin==4 | nargin==5 | nargin==2,
	if isempty(a),
		n = 3;
	else
		% The next line controls the number of samples in the plot if N not specified
		st=0.001; % Set settling time bound  = 0.1%
		[n,m]=size(b);
		x0=(eye(n,n)-a)\(b*ones(m,1));

		% Cater for pure integrator case
		infind=find(~finite(x0));
		x0(infind)=ones(length(infind),1);
		n=dtimvec(a,b,c,x0,st);
	end
	if nargin==4			% Multivariable systems
		[iu,nargin,y]=mulresp('dimpulse',a,b,c,d,n,nargout,0);
		if ~iu, if nargout, yout = y; end, return, end
	end
end

if (nargin <= 3)&(nargout <= 1),	% Transfer function description
	y = dlsim(num,den,[1;zeros(n-1,1)]);	% More efficient: uses FILTER
else
	if ~isempty(b), b=b(:,iu); end, d = d(:,iu);
	[y,x] = dlsim(a,b,c,d,[1;zeros(n-1,1)]);
end

if nargout==0,		% With no output arguments, plot graph
	status = ishold;

	stairs([0:n-1],y)
	hold on
	plot([0,n-1],[0;0],':')

	xlabel('No. of Samples')
	ylabel('Amplitude')

	if ~status, hold off, end	% Return hold to previous status
	return % Suppress output
end
yout = y;

