function out = createOutputStruct (it, bu, t, mbs, obs, status)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author       : Marta Cavaleiro
% Description  : Creates a struct that collects statistical data on the algorithm's run
% Input        : it ~ number of iterations
%		 bu ~ number of basis updates
%		 t ~ running time
%		 mbs ~ maximum size observed of a basis
%		 obs ~ size of the optimal basis
%		 status ~ Either MAXTIME or OPTIMAL
% Output       : The struct data structure.
% Last revised : Jun 29, 2019

out = struct('iters', it, 'basisUpdates', bu, 'time', t, ...
    'maxBasisSize', mbs, 'optBasisSize', obs, 'status', status);

end