function logBetapdf = logBetapdf(x, mu, sig)

% computing the alpha and beta that parameterize the Beta density as a
% function of mean and std
alpha = mu     * [mu * (1-mu) / sig^2 - 1];
beta  = (1-mu) * [mu * (1-mu) / sig^2 - 1];

% evaluating the desity
logBetapdf = log(betapdf(x, alpha, beta));
