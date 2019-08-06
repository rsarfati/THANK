function logpriorJPT = logpriorJPT_MCMC (param);

loprior = zeros(length(param));

prior(1)  = log(normpdf(param(1),0.30,0.05));       % alpha
prior(2)  = logBetapdf(param(2),0.5,0.15);          % iotap
prior(3)  = logBetapdf(param(3),0.5,0.15);          % iotaw
prior(4)  = log(normpdf(param(4),0.5,0.025));       % gamma100
prior(5)  = logBetapdf(param(5),0.5,0.1);           % h
prior(6)  = log(normpdf(param(6),0.15,0.05));       % lambdapss
prior(7)  = log(normpdf(param(7),0.15,0.05));       % lambdawss
prior(8)  = log(normpdf(param(8),0,0.5));           % Lss
prior(9)  = log(normpdf(param(9),0.5,0.1));         % pss
prior(10) = logGammapdf(param(10),0.25,0.1);        % Fbeta
prior(11) = logGammapdf(param(11),2,0.75);          % niu
prior(12) = logBetapdf(param(12),0.66,0.1);         % xip
prior(13) = logBetapdf(param(13),0.66,0.1);         % xiw
prior(14) = logGammapdf(param(14),5,1);             % chi 
prior(15) = logGammapdf(param(15),4,1);             % S
prior(16) = log(normpdf(param(16),1.7,0.3));        % fp
prior(17) = log(normpdf(param(17),0.125,0.05));     % fy
prior(18) = log(normpdf(param(18),0.125,0.05));     % fdy
 
prior(19) = logBetapdf(param(19),0.6,0.2);          % rhoR
prior(20) = logBetapdf(param(20),0.6,0.2);          % rhoz
prior(21) = logBetapdf(param(21),0.6,0.2);          % rhog
prior(22) = logBetapdf(param(22),0.6,0.2);          % rhomiu
prior(23) = logBetapdf(param(23),0.6,0.2);          % rholambdap
prior(24) = logBetapdf(param(24),0.6,0.2);          % rholambdaw
prior(25) = logBetapdf(param(25),0.6,0.2);          % rhob
prior(26) = logBetapdf(param(26),0.4,0.2);          % rhomp
prior(27) = logBetapdf(param(27),0.5,0.2);          % rhoARMAlambdap
prior(28) = logBetapdf(param(28),0.5,0.2);          % rhoARMAlambdaw
if any(param(19:28)>0.99); prior(19:28)=log(0); end

prior(29) = logIG1pdf(param(29),0.1,1);             % sdR
prior(30) = logIG1pdf(param(30),0.5,1);             % sdz
prior(31) = logIG1pdf(param(31),0.5,1);             % sdg
prior(32) = logIG1pdf(param(32),0.5,1);             % sdmiu
prior(33) = logIG1pdf(param(33),0.1,1);             % sdlambdap
prior(34) = logIG1pdf(param(34),0.1,1);             % sdlambdaw
prior(35) = logIG1pdf(param(35),0.1,1);             % sdb



if all(isfinite(prior))==0 | all(isreal(prior))==0;
    logpriorJPT=-1e10;
else
    logpriorJPT=sum(prior);
end