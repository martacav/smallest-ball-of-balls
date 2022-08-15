function [x, flag] = infQ_2pts (p1, p2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       : Marta Cavaleiro
% Description  : Solves the InfQ for two points (columns of P)
% Input        : P ~ matrix with the points in its columns
% Output       : x ~ optimal solution
%                flag ~ 0 if the solution is "the middle"; 
%                1 if it is p2; 2 if it is p1 (this is to instruct what
%                point to be removed from the support set)                 
% Last revised : April 26, 2018

n = size(p1, 1); n = n-1;

x = zeros(n+1,1);
x0 = min([p1(n+1), p2(n+1), ( p1(n+1)+p2(n+1)-norm(p1(1:n)-p2(1:n)) )/2]);
x(n+1) = x0;
x(1:n) = ( (p1(n+1)-x0)*p2(1:n)+(p2(n+1)-x0)*p1(1:n) ) / ...
         ( (p1(n+1)-x0)+(p2(n+1)-x0) );
     
if x0 == p1(n+1)
	flag = 2;
elseif x0 == p2(n+1)
	flag = 1;
else 
	flag = 0;
end