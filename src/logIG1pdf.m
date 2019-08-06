function logIG1pdf = logIG1pdf(x,mu,sig);

% computing the deg and scale that parameterize the IG1 density as a
% function of mean and std
[scale,deg]=inverse_gamma_specification(mu,sig,1); 

% evaluating the desity
[junk logIG1pdf]=pdf_igone(x,deg,scale);