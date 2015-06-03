function RWGMH(TargetFun,ProposalFun,xparam1,vv,mh_bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_)
% Random walk Metropolis-Hastings algorithm. 
% 
% INPUTS 
%   o TargetFun  [char]     string specifying the name of the objective
%                           function (posterior kernel).
%   o xparam1    [double]   (p*1) vector of parameters to be estimated (initial values).
%   o mh_bounds  [double]   (p*2) matrix defining lower and upper bounds for the parameters. 
%   o dataset_              data structure
%   o options_              options structure
%   o M_                    model structure
%   o estim_params_         estimated parameters structure
%   o bayestopt_            estimation options structure
%   o oo_                   outputs structure
%
% ALGORITHM 
%   Metropolis-Hastings.       
%
% SPECIAL REQUIREMENTS
%   None.
%
% PARALLEL CONTEXT
% The most computationally intensive part of this function may be executed
% in parallel. The code sutable to be executed in
% parallel on multi core or cluster machine (in general a 'for' cycle),
% is removed from this function and placed in random_walk_metropolis_hastings_core.m funtion.
% Then the DYNARE parallel package contain a set of pairs matlab functions that can be executed in
% parallel and called name_function.m and name_function_core.m. 
% In addition in parallel package we have second set of functions used
% to manage the parallel computation.
%
% This function was the first function to be parallelized, later other
% functions have been parallelized using the same methodology.
% Then the comments write here can be used for all the other pairs of
% parallel functions and also for management funtions.

% Copyright (C) 2006-2013 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
%
% AUTHOR(S) frederic DOT karame AT univ DASH lemans DOT fr
%           stephane DOT adjemian AT univ DASH lemans DOT fr


% In Metropolis, we set penalty to Inf to as to reject all parameter sets triggering error in target density computation
global objective_function_penalty_base
objective_function_penalty_base = Inf;

options_.proposal_distribution = 'RWGMH_proposal';

% Initialization of the random walk metropolis-hastings chains.
[ ix2, ilogpo2, ModelName, MetropolisFolder, fblck, fline, npar, nblck, nruns, NewFile, MAX_nruns ] = ...
    RWGMH_initialization(TargetFun, xparam1, mh_bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_);

InitSizeArray = min([repmat(MAX_nruns,nblck,1) fline+nruns-1],[],2);

% Load last mh history file
load_last_mh_history_file(MetropolisFolder, ModelName);

% Only for test parallel results!!!

% To check the equivalence between parallel and serial computation!
% First run in serial mode, and then comment the follow line.
%   save('recordSerial.mat','-struct', 'record');

% For parllel runs after serial runs with the abobe line active.
%   TempRecord=load('recordSerial.mat');
%   record.Seeds=TempRecord.Seeds;

% Snapshot of the current state of computing. It necessary for the parallel
% execution (i.e. to execute in a corretct way portion of code remotely or
% on many core). The mandatory variables for local/remote parallel
% computing are stored in localVars struct.

vv = zeros(length(xparam1),length(xparam1)) ;
d = zeros(length(xparam1),length(xparam1)) ;

localVars =   struct('TargetFun', TargetFun, ...
                     'ProposalFun', ProposalFun, ...
                     'xparam1', xparam1, ...
                     'vv', vv, ...
                     'mh_bounds', mh_bounds, ...
                     'ix2', ix2, ...
                     'ilogpo2', ilogpo2, ...
                     'ModelName', ModelName, ...
                     'fline', fline, ...
                     'npar', npar, ...
                     'nruns', nruns, ...
                     'NewFile', NewFile, ...
                     'MAX_nruns', MAX_nruns, ...
                     'd', d, ...
                     'InitSizeArray',InitSizeArray, ...
                     'record', record, ...
                     'dataset_', dataset_, ...
                     'dataset_info', dataset_info, ...
                     'options_', options_, ...
                     'M_',M_, ...
                     'bayestopt_', bayestopt_, ...
                     'estim_params_', estim_params_, ...
                     'oo_', oo_,...
                     'varargin',[]);


% The user don't want to use parallel computing, or want to compute a
% single chain. In this cases Random walk Metropolis-Hastings algorithm is
% computed sequentially.

fout = RWGMH_core(localVars, fblck, nblck, 0);
record = fout.record;

irun = fout(1).irun;
NewFile = fout(1).NewFile;

update_last_mh_history_file(MetropolisFolder, ModelName, record);

skipline()
disp(['Estimation::RWGMH: Number of mh files: ' int2str(NewFile(1)) ' per block.'])
disp(['Estimation::RWGMH: Total number of generated files: ' int2str(NewFile(1)*nblck) '.'])
disp(['Estimation::RWGMH: Total number of iterations: ' int2str((NewFile(1)-1)*MAX_nruns+irun-1) '.'])
disp(['Estimation::RWGMH: Current acceptance ratio per chain: '])
for i=1:nblck
    if i<10
        disp(['                                                       Chain  ' num2str(i) ': ' num2str(100*record.AcceptanceRatio(i)) '%'])
    else
        disp(['                                                       Chain ' num2str(i) ': ' num2str(100*record.AcceptanceRatio(i)) '%'])
    end
end