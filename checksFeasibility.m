function [k, maxGap, infeasCount] = checksFeasibility (Pb, P0, normsOfPbsq, x, option, iterLog)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       : Marta Cavaleiro
% Description  : Checks if (p_i-x)>=_Q 0 for all p_i\in P
% Note         : k=0 if x is feasible; If x is infeasible:
%                if option = 1, k is the index of the most infeasible 
%                point; if option = 2, k is the index of the first infeasible point found
% Last revised : Jun 24, 2019

global epsTol;
[n, m] = size(Pb);


xbtxb = x(1:n)'*x(1:n);

maxGap = NaN;
infeasCount = NaN;
gaps = P0 - sqrt(abs(normsOfPbsq - 2*x(1:n)'*Pb + xbtxb*ones(1, m)));


switch option 
    case 1 %the most infeasible one
        [mingap, k] = min(gaps); %finds the most infeasible one
                
        if mingap > x(n+1) - epsTol
            k = 0;
        end
        
        if iterLog 
            maxGap = x(n+1) - mingap;
            infeasCount = size(gaps(x(n+1)-gaps>=epsTol), 2); %gives the number of infeas constr.
        end
        
    case 2 %the first infeasible found
        k = find(gaps < x(n+1) - epsTol, 1); %the first infeasible
        if isempty(k)
            k=0;
        end
        
        if iterLog
            infeasCount = size(gaps(x(n+1)-gaps>=epsTol), 2); %gives the number of infeas constr.
            [mingap, ~] = min(gaps);
            maxGap = x(n+1) - mingap;
        end
end

end