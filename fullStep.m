function xFS = fullStep(xNow, pStar, p1, z, Mu, Mv, iterLog)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       : Marta Cavaleiro
% Description  : Calculates the full step
% Output       : xFS ~ the point corresponding the full step
% Last revised : Jul 9, 2019


epsTol = 10^-6;
n = length(p1) - 1; 

%I. CASE z = 0 (S U p* aff dep): NO FULL STEP
if norm(z) < epsTol
    xFS = -Inf*ones(n+1,1);
    return
end


%Auxiliary variables
pb1 = p1(1:n);
pStarb = pStar(1:n);

cStar = pStar(n+1) - p1(n+1);
bStar = (1/2) * ( pStarb'*pStarb - pb1'*pb1 - pStar(n+1)^2 + p1(n+1)^2);

ppp = pStarb - pb1; 
pptz = ppp'*z;

q = Mu + ( (bStar - ppp'*(Mu + pb1)) / pptz)*z;
r = Mv + ( (cStar - ppp'*Mv) / pptz)*z;


%II. SOLVE QUADRATIC
aa = r'*r - 1;
bb = 2*( q'*r + p1(n+1) );
cc = q'*q - p1(n+1)^2;
x0FS = quadEqSolve(aa, bb, cc);


%III. REMOVE SOLUTIONS THAT DO NOT SATISFY <=pstar0 AND <=x0
x0FS = x0FS ( x0FS <= pStar(n+1)+epsTol & x0FS <= xNow(n+1)+epsTol );
 

%IV. REMOVE SOLUTIONS THAT MAKE alphastar NEGATIVE
x0temp = x0FS;
for j = 1:length(x0temp)
    x0 = x0temp(j);

    alphaStar = (bStar + x0*cStar - ppp'*(Mu + x0*Mv + pb1)) / pptz;
    if alphaStar < -epsTol
        x0FS = x0FS(x0FS~=x0);
    end
end

if isempty(x0FS)
    x0FS = -Inf; %x0FS may be -Inf if no solution exists
end


%V. FIND MAXIMUM x0 AND GET xb
%keep the maximum (there may be two values)
x0FS = max(x0FS);
xbFS = q + x0FS*r + pb1;
xFS = [xbFS; x0FS];

end