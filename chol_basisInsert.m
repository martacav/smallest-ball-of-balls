function [Lnew, S, M] = chol_basisInsert (L, S, M, p, z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       : Marta Cavaleiro
% Description  : Updates the Cholesky factorization when a point is inserted in the basis
% Last revised : Jun 22, 2019


if isempty(S) %pstar was the solution of the curve search
    Lnew = []; 
    S = [p];
    M = [];
    return
end

[n, s] = size(S);
n = n-1;

mStar = p(1:n)-S(1:n, 1);      
if s==1    
    M = mStar;
    S = [S, p];
    Lnew = sqrt(M'*M); 
        
else          
    Lnew = zeros(s,s);
    Lnew(1:s-1,1:s-1) = L;
    opts.LT = true;
    aux = linsolve (L, M'*mStar, opts);    
    Lnew(s,:) = [aux', sqrt(mStar'*z)];
    M = [M, mStar];
    S = [S, p];   
   
end