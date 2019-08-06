  function X=disclyap_fast(G,V,ch)
% function X=disclyap_fast(G,V,ch)
% 
% Solve the discrete Lyapunov Equation 
% X=G*X*G'+V 
% Using the Doubling Algorithm 
%
% If ch is defined then the code will check if the resulting X 
% is positive definite and generate an error message if it is not 
% 
% Joe Pearlman and Alejandro Justiniano 
% 3/5/2005 
% =================================================================
if nargin == 2 | isempty( ch ) == 1 
    flag_ch = 0; 
else 
    flag_ch = 1; 
end 
s=size(G,1); 

tol = 1e-16; 

P0=V; 
A0=G; 

matd=1; 
while matd > tol 
    P1=P0+A0*P0*A0'; 
    A1=A0*A0;  
    matd=max( max( abs( P1 - P0 ) ) ); 
    P0=P1; 
    A0=A1; 
end 
clear A0 A1 P1; 

X=(P0+P0')/2; 

% Check that X is positive definite 
if flag_ch==1 
    [C,p]=chol(X); 
    if p ~= 0 
        error('X is not positive definite')
    end 
end 