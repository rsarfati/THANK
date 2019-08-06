  function [pd,logpdf]=pdf_igone(sd,v,s);
% function [pdf,logpdf]=pdf_igone(sd,v,s);
% Computes the pdf of an InverseGamma1 
% with degrees of freedom v and scale parameter s 
% This is for the standard deviation 
% Use invwishlub.m for the pdf of an IG2 
% AJ 4/29/04
num=log(2) - (v+1)*log(sd) -s/( 2*(sd^2) ); 
den=gammaln(v/2) + (v/2)*log(2/s); 
logpdf=num-den; 
pd=exp(logpdf); 
