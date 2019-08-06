function theta = boundsINV (param);


alpha=bound01INV(param(1));     
iotap=bound01INV(param(2));     
iotaw=bound01INV(param(3));     
gamma100=param(4);  
h=bound01INV(param(5));  
lambdapss=param(6); 
lambdawss=param(7); 
Lss=param(8);       
pss100=param(9);   
Fbeta=bound0INV(param(10));
niu=bound0INV(param(11));
xip=bound01INV(param(12));
xiw=bound01INV(param(13));
chi=bound0INV(param(14));
S=bound0INV(param(15));
fp=param(16); 
fy=param(17);       
fdy=param(18);      
rhoR=bound0BINV(param(19),.99);
rhoz=bound0BINV(param(20),.99);
rhog=bound0BINV(param(21),.99);
rhomiu=bound0BINV(param(22),.99);
rholambdap=bound0BINV(param(23),.99);
rholambdaw=bound0BINV(param(24),.99);
rhob=bound0BINV(param(25),.99);
rhomp=bound0BINV(param(26),.99);
rhoARMAlambdap=bound0BINV(param(27),.99);
rhoARMAlambdaw=bound0BINV(param(28),.99);
sdR=bound0INV(param(29));
sdz=bound0INV(param(30)); 
sdg=bound0INV(param(31)); 
sdmiu=bound0INV(param(32));
sdlambdap=bound0INV(param(33)); 
sdlambdaw=bound0INV(param(34)); 
sdb=bound0INV(param(35));


theta = [alpha iotap iotaw gamma100 h lambdapss lambdawss Lss pss100 Fbeta...
    niu xip xiw chi S fp fy fdy rhoR rhoz rhog rhomiu rholambdap rholambdaw...
    rhob rhomp rhoARMAlambdap rhoARMAlambdaw sdR sdz sdg sdmiu sdlambdap sdlambdaw sdb];



function rho = bound01INV(param);
rho = log(param/(1-param));

function rho = bound0BINV(param,B);
rho = log(param/B/(1-param/B));

function sigma = bound0INV(param);
sigma = log(param);