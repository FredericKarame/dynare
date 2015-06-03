function draw = RWGMH_proposal(Mean,b,options_,mh_bounds)
% Pseudo random draws from a multivariate normal distribution,
% \mathcal N_n(Mean,Sigma), with expectation Mean and variance Sigma.
%
% INPUTS 
%
%    Mean               [double]    number_of_chains*number_of_parameters matrix, expectation of the multivariate random variable.
%    
% OUTPUTS 
%    draw               [double]    number_of_chains*number_of_parameters matrix containing the proposal 
%        
% SPECIAL REQUIREMENTS

% Copyright (C) 2003-2009 Dynare Team
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

number_of_parameters = size(Mean,2) ;
number_of_chains = size(Mean,1) ;
draw = -1000*ones(1,number_of_parameters);
while ~(all( draw(:) > mh_bounds.lb ) && all( draw(:) < mh_bounds.ub ))
    indx1 = 1 ;
    indx2 = 1 ;
    while indx1==indx2 || indx1==b || indx2==b
      indx1 = randsample(1:number_of_chains,1) ;
      indx2 = randsample(1:number_of_chains,1) ;  
    end
    draw = Mean(b,:) ...
        + bsxfun(@times,options_.rwgmh.scale_chain,(2.38/(sqrt(2*(number_of_parameters))))*(Mean(indx1,:)-Mean(indx2,:)))  ...
        + bsxfun(@times,options_.rwgmh.scale_shock,(randn(1,number_of_parameters))) ;
end
