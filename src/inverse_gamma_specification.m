function [s,nu] = inverse_gamma_specification(mu,sigma,type)
%
% mu    --> expectation
% sigma --> standard deviation
% type = 1 --> Inverse Gamma 1
% type = 2 --> Inverse Gamma 2
%
% stephane.adjemian@cepremap.cnrs.fr [01/14/2004]


sigma2 = sigma^2;
mu2 = mu^2;

if type == 2;       % Inverse Gamma 2   
   nu   = 2*(2+mu2/sigma2);
   s    = 2*mu*(1+mu2/sigma2);
elseif type == 1;   % Inverse Gamma 1 
    if sigma2 < Inf;
      nu = sqrt(2*(2+mu2/sigma2));
      nu2 = 2*nu;
      nu1 = 2;
      err = 2*mu2*gamma(nu/2)^2-(sigma2+mu2)*(nu-2)*gamma((nu-1)/2)^2;
      while abs(nu2-nu1) > 1e-12
	if err > 0
	  nu1 = nu;
	  if nu < nu2
	    nu = nu2;
	  else
	    nu = 2*nu;
	    nu2 = nu;
	  end
	else
	  nu2 = nu;
	end
	nu =  (nu1+nu2)/2;
	err = 2*mu2*gamma(nu/2)^2-(sigma2+mu2)*(nu-2)*gamma((nu-1)/2)^2;
      end
      s = (sigma2+mu2)*(nu-2);
    else;
        nu  = 2;
        s   = 2*mu2/pi;
    end;   
else;
   s  = -1;
   nu = -1;
end;

% 01/18/2004 MJ replaced fsolve with secant
%               suppressed chck
%	        changed order of output parameters