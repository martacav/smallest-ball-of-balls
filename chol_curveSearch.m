function [xNew, facet, z] = chol_curveSearch (L, S, M, p, x, iterLog)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       : Marta Cavaleiro
% Description  : Performs one step of the curve search. Either stops at a
%                facet or at the point where p becomes feasible
% Input        : L ~ Cholesky factor of M'*M
%                x ~ current solution
%                S ~ n x r matrix of r affinely independent points in R^n
%                p ~ point in R^n s.t. {S, p} is affinely independent
%                iterLog ~ TRUE/FALSE whether we want a detailed output
%                report being printed on the console
% Output       : x_new ~ new solution
%                facet ~ 0 if intersects the bisectors related to p first;
%                  OR = i if intersects the facet opposite to point i
% Last revised : Jun 22, 2019



%SPECIAL CASE: 2 points
if size(S,2) == 1
    [xNew, facet] = infQ_2pts(S, p);
    z = 0;    
    return 
end


%Auxiliary variables from input
[n, s] = size(S); 
n = n-1;
pb = p(1:n);
pb1 = S(1:n,1); p10 = S(n+1,1);
normsSb2 = sum(S(1:n, :).^2, 1);
c = S(n+1,2:s)' - p10*ones(s-1,1);
b = (1/2) * ( normsSb2(2:s)' - normsSb2(1)*ones(s-1,1) - ...
            S(n+1,2:s)'.^2 + p10^2*ones(s-1,1) );
    

%I. GET VECTORS u, v, w, z        
opts1.LT = true;
aux1 = linsolve (L, [b-M'*pb1, c, -M'*(pb-pb1)], opts1);
opts2.UT = true;
aux2 = linsolve (L', aux1, opts2);
u = aux2(:,1);
v = aux2(:,2);
w = aux2(:,3);
z = M*w + (pb-pb1);


%II. CALCULATE STEPS
Mv = M*v;
Mu = M*u;
xFS = fullStep(x, p, S(:,1), z, Mu, Mv, iterLog);
[xPS, kPS] = partialStep (x, S(:,1), u, v, w, z, Mu, Mv, iterLog);


%III. TAKE A STEP
if xPS(n+1) >= xFS(n+1)
    xNew = xPS;
    facet = kPS;
else
    xNew = xFS;
    facet = 0;
end

end