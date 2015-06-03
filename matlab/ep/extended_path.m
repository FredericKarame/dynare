function ts = extended_path(initial_conditions,sample_size)
% Stochastic simulation of a non linear DSGE model using the Extended Path method (Fair and Taylor 1983). A time
% series of size T  is obtained by solving T perfect foresight models.
%
% INPUTS
%  o initial_conditions     [double]    m*nlags array, where m is the number of endogenous variables in the model and
%                                       nlags is the maximum number of lags.
%  o sample_size            [integer]   scalar, size of the sample to be simulated.
%
% OUTPUTS
%  o time_series            [double]    m*sample_size array, the simulations.
%
% ALGORITHM
%
% SPECIAL REQUIREMENTS

% Copyright (C) 2009-2015 Dynare Team
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
global M_ options_ oo_

options_.verbosity = options_.ep.verbosity;
verbosity = options_.ep.verbosity+options_.ep.debug;

% Set maximum number of iterations for the deterministic solver.
options_.simul.maxit = options_.ep.maxit;

% Prepare a structure needed by the matlab implementation of the perfect foresight model solver
pfm = setup_stochastic_perfect_foresight_model_solver(M_,options_,oo_);

exo_nbr = M_.exo_nbr;
ep = options_.ep;
steady_state = oo_.steady_state;
dynatol = options_.dynatol;

% Set default initial conditions.
if isempty(initial_conditions)
    if isempty(M_.endo_histval)
        initial_conditions = oo_.steady_state;
    else
        initial_conditions = M_.endo_histval;
    end
end


% Set the number of periods for the perfect foresight model
periods = options_.ep.periods;
pfm.periods = periods;
pfm.i_upd = pfm.ny+(1:pfm.periods*pfm.ny);

% keep a copy of pfm.i_upd
i_upd = pfm.i_upd;

% Set the algorithm for the perfect foresight solver
options_.stack_solve_algo = options_.ep.stack_solve_algo;

% Set check_stability flag
do_not_check_stability_flag = ~options_.ep.check_stability;

% Compute the first order reduced form if needed.
%
% REMARK. It is assumed that the user did run the same mod file with stoch_simul(order=1) and save
% all the globals in a mat file called linear_reduced_form.mat;

dr = struct();
if options_.ep.init
    options_.order = 1;
    [dr,Info,M_,options_,oo_] = resol(1,M_,options_,oo_);
end

% Do not use a minimal number of perdiods for the perfect foresight solver (with bytecode and blocks)
options_.minimal_solving_period = 100;%options_.ep.periods;

% Initialize the exogenous variables.
% !!!!!!!! Needs to fixed
options_.periods = periods;
make_ex_;

% Initialize the endogenous variables.
make_y_;

% Initialize the output array.
time_series = zeros(M_.endo_nbr,sample_size);

% Set the covariance matrix of the structural innovations.
variances = diag(M_.Sigma_e);
positive_var_indx = find(variances>0);
effective_number_of_shocks = length(positive_var_indx);
stdd = sqrt(variances(positive_var_indx));
covariance_matrix = M_.Sigma_e(positive_var_indx,positive_var_indx);
covariance_matrix_upper_cholesky = chol(covariance_matrix);

% (re)Set exo_nbr
%exo_nbr = effective_number_of_shocks;

% Set seed.
if options_.ep.set_dynare_seed_to_default
    set_dynare_seed('default');
end

% Set bytecode flag
bytecode_flag = options_.ep.use_bytecode;

% Simulate shocks.
switch options_.ep.innovation_distribution
  case 'gaussian'
    oo_.ep.shocks = transpose(transpose(covariance_matrix_upper_cholesky)*randn(effective_number_of_shocks,sample_size));
  otherwise
    error(['extended_path:: ' options_.ep.innovation_distribution ' distribution for the structural innovations is not (yet) implemented!'])
end


% Set waitbar (graphic or text  mode)
hh = dyn_waitbar(0,'Please wait. Extended Path simulations...');
set(hh,'Name','EP simulations.');

% hybrid correction
pfm.hybrid_order = options_.ep.stochastic.hybrid_order;
if pfm.hybrid_order
    oo_.dr = set_state_space(oo_.dr,M_,options_);
    options = options_;
    options.order = pfm.hybrid_order;
    pfm.dr = resol(0,M_,options,oo_);
else
    pfm.dr = [];
end

% number of nonzero derivatives
pfm.nnzA = M_.NNZDerivatives(1);

% setting up integration nodes if order > 0
if options_.ep.stochastic.order > 0
    [nodes,weights,nnodes] = setup_integration_nodes(options_.ep,pfm);
    pfm.nodes = nodes;
    pfm.weights = weights; 
    pfm.nnodes = nnodes;

    % compute number of blocks
    [block_nbr,pfm.world_nbr] = get_block_world_nbr(options_.ep.stochastic.algo,nnodes,options_.ep.stochastic.order,options_.ep.periods);
else
    block_nbr = options_.ep.periods
end


% set boundaries if mcp
[lb,ub,pfm.eq_index] = get_complementarity_conditions(M_, options_.ramsey_policy);
options_.lmmcp.lb = repmat(lb,block_nbr,1);
options_.lmmcp.ub = repmat(ub,block_nbr,1);
pfm.block_nbr = block_nbr;

% storage for failed draws
oo_.ep.failures.periods = [];
oo_.ep.failures.previous_period = cell(0);
oo_.ep.failures.shocks = cell(0);

% Initializes some variables.
t  = 0;
tsimul = 1;
% Main loop.
while (t<sample_size)
    if ~mod(t,10)
        dyn_waitbar(t/sample_size,hh,'Please wait. Extended Path simulations...');
    end
    % Set period index.
    t = t+1;
    % Put shocks in oo_.exo_simul (second line).
    exo_simul_1 = zeros(periods+2,exo_nbr);
    exo_simul_1(2,positive_var_indx) = oo_.exo_simul(2,positive_var_indx) + oo_.ep.shocks(t,:);
    if ep.init% Compute first order solution (Perturbation)...
        initial_path = simult_(initial_conditions,dr,exo_simul_1(2:end,:),1);
        endo_simul_1(:,1:end-1) = initial_path(:,1:end-1)*ep.init+endo_simul_1(:,1:end-1)*(1-ep.init);
    else
        if t==1
            endo_simul_1 = repmat(steady_state,1,periods+2);
        end
    end
    % Solve a perfect foresight model.
    % Keep a copy of endo_simul_1
    endo_simul = endo_simul_1;
    if verbosity
        save ep_test_1 endo_simul_1 exo_simul_1
    end
    if bytecode_flag && ~options_.ep.stochastic.order
        [flag,tmp] = bytecode('dynamic',endo_simul_1,exo_simul_1, M_.params, endo_simul_1, options_.ep.periods);
    else
        flag = 1;
    end
    if flag
        if options_.ep.stochastic.order == 0
            [flag,tmp,err] = solve_perfect_foresight_model(endo_simul_1,exo_simul_1,pfm);
        else
            switch(options_.ep.stochastic.algo)
              case 0
                [flag,tmp] = ...
                    solve_stochastic_perfect_foresight_model(endo_simul_1,exo_simul_1,pfm,options_.ep.stochastic.quadrature.nodes,options_.ep.stochastic.order);
              case 1
                [flag,tmp] = ...
                    solve_stochastic_perfect_foresight_model_1(endo_simul_1,exo_simul_1,options_,pfm,options_.ep.stochastic.order);
            end
        end
    end
    info_convergence = ~flag;
    if verbosity
        if info_convergence
                disp(['Time: ' int2str(t)  '. Convergence of the perfect foresight model solver!'])
        else
                disp(['Time: ' int2str(t)  '. No convergence of the perfect foresight model solver!'])
        end
    end
    endo_simul_1 = tmp;
    if info_convergence
        % Save results of the perfect foresight model solver.
        time_series(:,tsimul) = endo_simul_1(:,2);
        endo_simul_1(:,1:end-1) = endo_simul_1(:,2:end);
        endo_simul_1(:,1) = time_series(:,tsimul);
        endo_simul_1(:,end) = oo_.steady_state;
        tsimul = tsimul+1;
    else
        oo_.ep.failures.periods = [oo_.ep.failures.periods t];
        oo_.ep.failures.previous_period = [oo_.ep.failures.previous_period  endo_simul_1(:,1)];
        oo_.ep.failures.shocks = [oo_.ep.failures.shocks  shocks];
        endo_simul_1 = repmat(steady_state,1,periods+2);
        endo_simul_1(:,1) = time_series(:,tsimul-1);
    end
end% (while) loop over t

dyn_waitbar_close(hh);

if isnan(options_.initial_period)
    initial_period = dates(1,1);
else
    initial_period = optins_.initial_period;
end
if nargout
    ts = dseries(transpose([initial_conditions, time_series]),initial_period,cellstr(M_.endo_names));
else
    oo_.endo_simul = [initial_conditions, time_series];
    ts = dseries(transpose(oo_.endo_simul),initial_period,cellstr(M_.endo_names));
    dyn2vec;
end

 assignin('base', 'Simulated_time_series', ts);
