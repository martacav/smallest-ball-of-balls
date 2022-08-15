function L2 = chol_lowRankUpdate (L, u, v)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       : Marta Cavaleiro
% Description  : Given the Cholesky factorization of B = LL', it returns 
%                the Cholesky factorization of
%                      B2 := (I + vu')B(I + uv')
%               (corresponding to a 2-rank update to matrix B)
%                Follows Goldfarb "Factorized Variable Metric Methods for
%                Unconstrained Optimization", 1976
% Input        : L ~ Cholesky factor (lower triangular) of B     
%                v, u ~ vectors that generate the B2
% Output       : L2 ~ Cholesky factor of B2
%                Y ~ matrix with the dual variables y_i on the columns
% Last revised : 22 Jun 2019

n = size(L,1);


%1. Get z, w s.t. Lz=v and w=L'u
opts.LT = true;
z = linsolve (L, v, opts);
w = L'*u;


%2. Get beta, gamm, lamb
%A. Recurrence 1 of Goldfarb's paper
s = zeros(n-1,1);
betabar = [zeros(n-1,1); 1/w(n)];
for j=n-1:-1:1
    r = betabar(j+1)*w(j);
    s(j) = (r^2+1)^(-1/2);
    c = r*s(j);
    betabar(j) = s(j)*betabar(j+1);
    betabar(j+1) = c*betabar(j+1);
end

%B. Recurrence 2 of Goldfarb's paper
beta = zeros(n-1,1);
gama = zeros(n-1,1);
lamb = zeros(n,1);
gamabar = zeros(n-1,1);
lambbar = zeros(n,1);

gamabar(1) = 1/betabar(1);
lambbar(1) = betabar(1)*w(1)+gamabar(1)*z(1);
for j=1:n-1
    lamb(j)=sqrt(lambbar(j)^2+s(j)^2);
    cbar=lambbar(j)/lamb(j);
    sbar=s(j)/lamb(j);
    beta(j)=cbar*betabar(j)-sbar*betabar(j+1);
    gama(j)=cbar*gamabar(j);
    gamabar(j+1)=sbar*gamabar(j);
    betabar(j+1)=sbar*betabar(j)+cbar*betabar(j+1);
    lambbar(j+1)=betabar(j+1)*w(j+1)+gamabar(j+1)*z(j+1);
end
lamb(n)=lambbar(n);


%3. Get Ltil and then L2
Ltil = diag(lamb)+[zeros(1,n); tril(w(2:n)*beta'+z(2:n)*gama',0), zeros(n-1,1)];
L2 = L*Ltil;


end