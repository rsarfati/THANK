function logGammapdf = logGammapdf(x,mu,sig);

% computing the alpha and beta that parameterize the Gamma density as a
% function of mean and std
beta=(sig^2)/mu; 
alpha=mu/beta;

% evaluating the desity
logGammapdf = log(gampdf(x,alpha,beta));