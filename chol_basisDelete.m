function [Lnew, S, M] = chol_basisDelete (L, S, M, flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       : Marta Cavaleiro
% Description  : Updates the Cholesky factorization when a point is deleted from the basis
% Last revised : Jun 22, 2019


[n, s] = size(S);
n = n-1;

if s<=2
   M=[];
   Lnew=[]; 
   S(:,flag)=[];
   return
end

switch flag
    case 1      
        a = [zeros(s-2,1); -1];
        b = ones(s-1,1);      
        Laux = chol_lowRankUpdate (L, a, b);
        Lnew = Laux(1:s-2, 1:s-2);
        M = M - M(:,s-1)*ones(1,s-1);
        M(:,s-1)=[];
        S(:,1)=S(:,s);
        S(:,s)=[];
            
    case s
        Lnew = L(1:s-2, 1:s-2);
        M(:,flag-1)=[];
        S(:,flag)=[];
            
    otherwise
        k=flag;
        Lnew = zeros(s-2, s-2);
        Lnew(1:k-2,1:k-2) = L(1:k-2,1:k-2);
        Lnew(k-1:s-2,1:k-2) = L(k:s-1,1:k-2);
        Lnew(k-1:s-2, k-1:s-2) = (cholupdate(L(k:s-1,k:s-1)', L(k:s-1, k-1)))';
        M(:,flag-1)=[];
        S(:,flag)=[];
end

end
            