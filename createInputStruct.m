function in = createInputStruct (feasTol, feasOption, maxTime, iterLog)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       : Marta Cavaleiro
% Description  : Creates a struct with input options.
% Input        : feasTol ~ Primal feasibility tolerance 
%   		 feasOption ~ 1: it uses the most infeasible point; 2: it uses the first infeasible point found;
%		 maxTime ~ maximum running time allowed;
%		 iterLog ~ Whether to print a log of each iteration on the console during the algorithm's run;
% Output       : The struct data structure.
% Last revised : Jul 1, 2019



in = struct('feasTol', feasTol, 'feasOption', feasOption, 'maxTime', maxTime, ...
    'iterLog', iterLog);

end