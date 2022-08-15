function [x, S, output] = main(P, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       : Marta Cavaleiro
% Description  : Solves the InfQ for the points in P:
%                       max x_{n+1}  s.t.  (p_i-x) in SOC 
% Input        : P ~ matrix with the points in its columns
%                options ~ a structure with input options (See createinputstruct.m)
%			   When options is empty, the algorithm creates one.
% Output       : x ~ optimal solution
%                S ~ basis (indexes of the points)
%                output ~ a structure with output data (See createoutputstruct.m)
% Last revised : Jul 17, 2019
%
% If used, please cite: "Cavaleiro, M., Alizadeh, F. A dual simplex-type algorithm for 
%         the smallest enclosing ball of balls. Comput Optim Appl 79, 767–787 (2021)."



%NOTE:
%-The height of the point is considered to be the last coordinate. We consider 
% that any point is such that:
% x = (\bar{x}, x_0) with \bar{x}\in\RR^{n} and x_0\in\RR, so x\in\RR^{n+1}
% Pb = {\bar{p_i}, i=1,..,m} and P0 = {p_{i0}, i=1,..,m}


%options
if nargin == 1  %options not given
    options = createInputStruct (10^-8, 1, Inf, true);                        
end

global epsTol
epsTol = options.feasTol;

tic;
[n, m] = size(P); n = n-1;
Pb = P(1:n, :); P0 = P(n+1, :);
normsOfPbsq = sum(Pb.^2, 1);

if options.iterLog
    fprintf('dimension: %d\nconstraints: %d\n', n+1, m);   
    fprintf('----------------------------------------------------------------\n');
    fprintf('Iter  max Infeas Gap   Infeas Constr  norm(x_(j+1)-x_j)  size(S)\n' );
end

%----------------------- INITIAL CHECKS -----------------------
%i. solution when there's only one point
if m==1
    x = P;
    S = P;
    output = createOutputStruct (0, 0, 0, 1, 1, 'OPTIMAL');
    return
end

%ii. check if the solution is one of the given points
[~, kk] = min(P(n+1,:));
[ip, xjGap, infeasCount] = checksFeasibility (Pb, P0, normsOfPbsq, P(:,kk), 2, options.iterLog);
if ip == 0
    x = P(:,kk);
    S = x;
    tt = toc;
    output = createOutputStruct (1, 0, tt, 1, 1, 'OPTIMAL');
    if options.iterLog
        fprintf(' %d    %12.8f         %d              N/A              %d\n', output.iters, xjGap, infeasCount,  size(S,2));          
        fprintf('----------------------------------------------------------------\n');                
        fprintf('Status     Iter Count   Bs Upd Count     Time      Opt Basis Size\n');
        fprintf('%s       %d              %d          %.5f       %d\n', output.status, output.iters, output.basisUpdates, output.time, output.optBasisSize);
        fprintf('----------------------------------------------------------------\n');
    end    
    return
end
  

%------------------------- ALGORITHM ----------------------

%0. Initialization
x = P(:,1);
S = P(:,1);
M = [];
L = []; 
isOpt = 0;
output = createOutputStruct (0, 0, 0, 0, 0,'');

while isOpt == 0

    output.iters = output.iters + 1;    
    x_j = x;  %x at the beginning of the iter
        
%I. Optimality check and choice of an infeasible point
    [ip, xjGap, infeasCount] = checksFeasibility (Pb, P0, normsOfPbsq, x, options.feasOption, options.iterLog);   
    if ip == 0  
        if options.iterLog
            fprintf(' %d    %12.8f         %d           N/A             %d\n', output.iters, xjGap, infeasCount,  size(S,2));
        end
        break
    end   
    p = P(:, ip);

  
%II. INNER LOOP
    flag = 1;
    while flag ~= 0 && size(S, 2)>0
        
        [x, flag, z] = chol_curveSearch (L, S, M, p, x, options.iterLog);  
        if flag == 0
            %Leave the loop:
            break 
        else  
            %Update basis and Cholesky factor 
            [L, S, M] = chol_basisDelete (L, S, M, flag);
            
            output.basisUpdates = output.basisUpdates + 1;
            if toc > options.maxTime
               output.time = toc;
               output.status = 'MAXTIME';
               return 
            end
        end
    end
    
    
    %Update support set and Cholesky factor
    [L, S, M] = chol_basisInsert (L, S, M, p, z);    
    output.basisUpdates = output.basisUpdates + 1;
    output.maxBasisSize = max(size(S,2), output.maxBasisSize);
    
    if options.iterLog
        fprintf(' %d    %12.8f         %d         %12.8f      %d\n', output.iters, xjGap,infeasCount, norm(x-x_j), size(S,2));
    end
end


output.status = 'OPTIMAL';
output.time = toc;
output.optBasisSize = size(S,2);
if options.iterLog
    fprintf('----------------------------------------------------------------\n');                
    fprintf('Status     Iter Count   Bs Upd Count     Time      Opt Basis Size\n');
    fprintf('%s       %d           %d         %.5f       %d\n', output.status, output.iters, output.basisUpdates, output.time, output.optBasisSize);
    fprintf('----------------------------------------------------------------\n');
end

end