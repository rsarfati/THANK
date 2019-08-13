function J = jacobTHANK(param)

J = zeros(length(param),1);

J(1)  = bound01prime(param(1));     
J(2)  = bound01prime(param(2));     
J(3)  = bound01prime(param(3));     
J(4)  = 1;  
J(5)  = bound01prime(param(5));  
J(6)  = 1; 
J(7)  = 1; 
J(8)  = 1;       
J(9)  = 1;   
J(10) = bound0prime( param(10));
J(11) = bound0prime( param(11));
J(12) = bound01prime(param(12));
J(13) = bound01prime(param(13));
J(14) = bound0prime( param(14));
J(15) = bound0prime( param(15));
J(16) = 1; 
J(17) = 1;       
J(18) = 1;      
J(19) = bound0Bprime(param(19), .99);
J(20) = bound0Bprime(param(20), .99);
J(21) = bound0Bprime(param(21), .99);
J(22) = bound0Bprime(param(22), .99);
J(23) = bound0Bprime(param(23), .99);
J(24) = bound0Bprime(param(24), .99);
J(25) = bound0Bprime(param(25), .99);
J(26) = bound0Bprime(param(26), .99);
J(27) = bound0Bprime(param(27), .99);
J(28) = bound0Bprime(param(28), .99);
J(29) = bound0prime( param(29));
J(30) = bound0prime( param(30)); 
J(31) = bound0prime( param(31)); 
J(32) = bound0prime( param(32));
J(33) = bound0prime( param(33)); 
J(34) = bound0prime( param(34)); 
J(35) = bound0prime( param(35));

J(36) = bound0Bprime(param(36), .99);
J(37) = bound0Bprime(param(37), .99);
J(38) = bound01prime(param(38));
J(39) = bound0Bprime(param(39), .99);
J(40) = bound0Bprime(param(40), .99);

J = diag(J);

function y = bound01prime(x)
y = exp(x)/(1+exp(x))^2;

function y = bound0Bprime(x,B)
y = B*exp(x)/(1+exp(x))^2;

function y = bound0prime(x)
y = exp(x);
