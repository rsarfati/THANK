function logpost = logpostTHANK(param,T,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Date: December 22, 2009
% this function computes the value of the posterior density for the model
% in Justiniano, Primiceri and Tambalotti (2010)
% To be used in the minimization algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% re-defining the parameter to transform the unconstrained into a
% constrained minimization
param = bounds(param);

% prior density
logprior = logpriorTHANK(param);

% posterior calculation
if logprior == -1e10
    logpost = 1e10;
    return    
else
    % Solution of the model
    [G1, C, M, eu, SDX, H, NY, NX] = modelTHANK(param);
    if eu ~= [1 1]'
        logpost = 1e10;
        return
    else
        % initializing the Kalman filter
        shat  = zeros(NY, 1);
        Q     = M * SDX * SDX' * M';
        sig   = disclyap_fast(G1, Q);
        LOGLH = 0;

        % Kalman filter recursion delivering log-likelihood
        for t=1:T
            [shat, sig, loglh] = kfilter(y(t,:)', H, C, G1, shat, sig, 0, Q);
            LOGLH = LOGLH + loglh;
        end
        logpost = -(LOGLH + logprior); % logpost = -logpost b/c using minimization algorithm (as opposed to maximization)
    end
end