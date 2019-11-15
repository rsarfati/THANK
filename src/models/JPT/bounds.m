function theta = bounds (param);


alpha=bound01(param(1));     
iotap=bound01(param(2));     
iotaw=bound01(param(3));     
gamma100=param(4);  
h=bound01(param(5));  
lambdapss=param(6); 
lambdawss=param(7); 
Lss=param(8);       
pss100=param(9);   
Fbeta=bound0(param(10));
niu=bound0(param(11));
xip=bound01(param(12));
xiw=bound01(param(13));
chi=bound0(param(14));
S=bound0(param(15));
fp=param(16); 
fy=param(17);       
fdy=param(18);      
rhoR=bound0B(param(19),.99);
rhoz=bound0B(param(20),.99);
rhog=bound0B(param(21),.99);
rhomiu=bound0B(param(22),.99);
rholambdap=bound0B(param(23),.99);
rholambdaw=bound0B(param(24),.99);
rhob=bound0B(param(25),.99);
rhomp=bound0B(param(26),.99);
rhoARMAlambdap=bound0B(param(27),.99);
rhoARMAlambdaw=bound0B(param(28),.99);
sdR=bound0(param(29));
sdz=bound0(param(30)); 
sdg=bound0(param(31)); 
sdmiu=bound0(param(32));
sdlambdap=bound0(param(33)); 
sdlambdaw=bound0(param(34)); 
sdb=bound0(param(35));

eta=param(36);

theta = [alpha iotap iotaw gamma100 h lambdapss lambdawss Lss pss100 Fbeta...
    niu xip xiw chi S fp fy fdy rhoR rhoz rhog rhomiu rholambdap rholambdaw...
    rhob rhomp rhoARMAlambdap rhoARMAlambdaw sdR sdz sdg sdmiu sdlambdap sdlambdaw sdb eta];



function rho = bound01(param);
rho = 1-1/(1+exp(param));

function rho = bound0B(param,B);
rho = B*[1-1/(1+exp(param))];

function sigma = bound0(param);
sigma = exp(param);