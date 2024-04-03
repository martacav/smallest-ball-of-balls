function [xPS, kPS] = partialStep (xNow, p1, u, v, w, z, Mu, Mv, iterLog)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       : Marta Cavaleiro
% Description  : Calculates the partial step
% Output       : xPS ~ the point on the curve where a facet is first intersected
%                kPS ~ index of the point corresponding to the facet 
% Last revised : Jul 9, 2019


epsTol = 10^-6;
s = length(u)+1;
n = length(p1) - 1; 
x0 = xNow(n+1);

%I. CASE z = 0 (S U p* aff dep):
if norm(z)<epsTol
        
    numrtor = [-1+ones(s-1,1)'*(u+x0*v); -u-x0*v];
    denomntr = [-ones(s-1,1)'*w-1; w];
    ratios = numrtor ./ denomntr;
    ratios(denomntr>0) = Inf;   
    [~, kPS] = min(ratios);
    xPS = xNow;
    
    return 
end



%Auxiliary variables
pb1 = p1(1:n); p01 = p1(n+1);
MvtMv = Mv'*Mv;
MutMu = Mu'*Mu;
MutMv = Mu'*Mv;
ztz = z'*z;
sqrtzz = sqrt(ztz);

x0PS = [];
x0kFirst = -Inf;

%II. Calculation at the points where the each dual variable becomes 0 first
%coefficients of the quadratic for alpha_2:s  
aux = w.^2/ztz;
AA = v.^2 - (1-MvtMv)*aux;
BB = 2*(u.*v + (p01 + MutMv)*aux);
CC = u.^2 - (p01^2-MutMu)*aux;


for k=1:s  
      
    %-----> alpha_2:s    
    if k>1  
        i = k-1;

        %i. SOLVE QUADRATIC        
        x0temp = quadEqSolve(AA(i), BB(i), CC(i));
        x0k = [];
        if x0temp > -Inf
        %ii. REMOVE EXTRANEOUS SOLUTIONS      
            for j = 1:length(x0temp)
                x0 = x0temp(j);
                aux = Mu+x0*Mv;
                error = ( u(i) + x0*v(i) )  + w(i)*sqrt( ...
                    (p01 - x0)^2 - (aux)'*(aux)    )/sqrtzz;
                if abs(error) < epsTol
                    x0k = [x0k; x0temp(j)];
                end
            end                      
        end
    
    %----> alpha_1:   
    else
        %i. SOLVE QUADRATIC       
        v1 = ones(s-1,1)'*v;
        u1 = ones(s-1,1)'*u;
        w1 = ones(s-1,1)'*w;
        aux = (w1 + 1)^2 / ztz;
        aa = v1^2 -  aux*(1 - MvtMv);
        bb = 2*( -(1 - u1)*v1 + aux*(p01 + MutMv) );
        cc = (1-u1)^2 - aux*(p01^2 - MutMu);   
        x0temp = quadEqSolve(aa, bb, cc);        
        x0k = [];
        if x0temp > -Inf
        %ii. REMOVE EXTRANEOUS SOLUTIONS
            for j = 1:length(x0temp)
                x0 = x0temp(j);
                aux = Mu+x0*Mv;
                error = (1-u1-x0*v1) - (w1 + 1)*sqrt (...
                    (p01 - x0)^2 - (aux)'*(aux)    )/sqrtzz;
                if abs(error) < epsTol                    
                    x0k = [x0k; x0temp(j)]; 
                end
            end
        end
        
    end %if
    
     
    %iii. REMOVE SOLUTIONS THAT DO NOT SATISFY <=x0   
    x0k = x0k ( x0k <= xNow(n+1)+epsTol );
    
    
    %iv. KEEP THE MAXIMUM (first from x0k - they may be two values)
    x0k = max(x0k);
    
    %v. KEEP THE MAXIMUM OVER ALL k
    if x0k > x0kFirst 
        x0kFirst = x0k;
        kPS = k;        
        x0PS = x0k;
    end
end


%III. FIND xb
alphaStar = sqrt((p01 - x0PS)^2 - MutMu -2*x0PS*MutMv-MvtMv*x0PS^2)/sqrt(ztz);
xbPS = Mu + x0PS*Mv + alphaStar*z + pb1;
xPS = [xbPS; x0PS];

end