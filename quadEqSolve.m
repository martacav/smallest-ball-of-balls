function x = quadEqSolve(a, b, c)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       : Marta Cavaleiro
% Description  : Solves the following quadratic equation 
%                ax^2 + bx + c = 0 
%                and it retunrs  -Inf if the solution is imaginary
% Input        : Constants a, b, c
% Output       : x ~ solution(s) (may be more than 1) or -Inf if imaginary
% Last revised : June 25, 2019


%SLOWER:
% x = roots([a,b,c]);
% if ~isreal(x)
%     x = -Inf;
% end


epsTol = 10^-12; %the precision here must be high


d = b^2 - 4*a*c;

if d < 0                
    x = -Inf;
    return
end

if abs(a) <= epsTol     
    x = -c/b;
    return
end

if abs(d) < epsTol      
    x = -b / (2*a);
    return
end


%there are two real solutions       
if b>0
    x(1) = (-b - sqrt(d)) / (2*a);      
else
	x(1) = (-b + sqrt(d)) / (2*a);    
end
x(2) = c/(a*x(1));     


end