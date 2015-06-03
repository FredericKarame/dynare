function myoutput = RWGMH_core(myinputs,fblck,nblck,whoiam,ThisMatlab)
% PARALLEL CONTEXT
% This function contain the most computationally intensive portion of code in
% random_walk_metropolis_hastings (the 'for xxx = fblck:nblck' loop). The branches in 'for'
% cycle and are completely independent than suitable to be executed in parallel way.
%
% INPUTS
%   o myimput            [struc]     The mandatory variables for local/remote
%                                    parallel computing obtained from random_walk_metropolis_hastings.m
%                                    function.
%   o fblck and nblck    [integer]   The Metropolis-Hastings chains.
%   o whoiam             [integer]   In concurrent programming a modality to refer to the differents thread running in parallel is needed.
%                                    The integer whoaim is the integer that
%                                    allows us to distinguish between them. Then it is the index number of this CPU among all CPUs in the
%                                    cluster.
%   o ThisMatlab         [integer]   Allows us to distinguish between the
%                                    'main' matlab, the slave matlab worker, local matlab, remote matlab,
%                                     ... Then it is the index number of this slave machine in the cluster.
% OUTPUTS
%   o myoutput  [struc]
%               If executed without parallel is the original output of 'for b =
%               fblck:nblck' otherwise a portion of it computed on a specific core or
%               remote machine. In this case:
%                               record;
%                               irun;
%                               NewFile;
%                               OutputFileName
%
% ALGORITHM
%   Portion of Metropolis-Hastings.
%
% SPECIAL REQUIREMENTS.
%   None.

% PARALLEL CONTEXT
% The most computationally intensive part of this function may be executed
% in parallel. The code sutable to be executed in parallel on multi core or cluster machine,
% is removed from this function and placed in random_walk_metropolis_hastings_core.m funtion.
% Then the DYNARE parallel package contain a set of pairs matlab functios that can be executed in
% parallel and called name_function.m and name_function_core.m.
% In addition in the parallel package we have second set of functions used
% to manage the parallel computation.
%
% This function was the first function to be parallelized, later other
% functions have been parallelized using the same methodology.
% Then the comments write here can be used for all the other pairs of
% parallel functions and also for management funtions.


% Copyright (C) 2006-2014 Dynare Team
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

%if nargin<4,
%    whoiam=0;
%end

% reshape 'myinputs' for local computation.
% In order to avoid confusion in the name space, the instruction struct2local(myinputs) is replaced by:

TargetFun=myinputs.TargetFun;
ProposalFun=myinputs.ProposalFun;
xparam1=myinputs.xparam1;
vv=myinputs.vv;
mh_bounds=myinputs.mh_bounds;
ix2=myinputs.ix2;
ilogpo2=myinputs.ilogpo2;
ModelName=myinputs.ModelName;
fline=myinputs.fline;
npar=myinputs.npar;
nruns=myinputs.nruns;
NewFile=myinputs.NewFile;
MAX_nruns=myinputs.MAX_nruns;
d=myinputs.d;
InitSizeArray=myinputs.InitSizeArray;
record=myinputs.record;
dataset_ = myinputs.dataset_;
dataset_info = myinputs.dataset_info;
bayestopt_ = myinputs.bayestopt_;
estim_params_ = myinputs.estim_params_;
options_ = myinputs.options_;
M_ = myinputs.M_;
oo_ = myinputs.oo_;
varargin=myinputs.varargin;

% Necessary only for remote computing!
%if whoiam
%    Parallel=myinputs.Parallel;
    % initialize persistent variables in priordens()
%    priordens(xparam1,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7, bayestopt_.p3,bayestopt_.p4,1);
%end

MetropolisFolder = CheckPath('metropolis',M_.dname);
ModelName = M_.fname;
BaseName = [MetropolisFolder filesep ModelName];

options_.lik_algo = 1;
OpenOldFile = ones(nblck,1);
if strcmpi(ProposalFun,'rand_multivariate_normal')
    n = npar;
elseif strcmpi(ProposalFun,'rand_multivariate_student')
    n = options_.student_degrees_of_freedom;
end

%
% NOW i run the (nblck-fblck+1) metropolis-hastings chains
%
%???
if (options_.load_mh_file~=0) && (fline(b)>1) && OpenOldFile(b)
    load([BaseName '_mh' int2str(NewFile(b)) '_blck' int2str(b) '.mat'])
    x2 = [x2;zeros(InitSizeArray(b)-fline(b)+1,npar)];
    logpo2 = [logpo2;zeros(InitSizeArray(b)-fline(b)+1,1)];
    OpenOldFile(b) = 0;
else
    x2 = zeros(nruns(1),npar,nblck);
    logpo2 = [ilogpo2';zeros(nruns(1),nblck)];
    R_indicator = zeros(nruns(1),1) ;
end

logpost1 = -feval(TargetFun,xparam1(:),dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_);

set_dynare_seed('default');
search_for_convergence = 1 ;
ix3 = ix2 ;
isux = 0;
jsux = zeros(1,nblck);
irun = fline(1);
j = 1;
while j <= nruns(1)
    jloop=0;
    skipline()
    disp(['Draw # ' num2str(j)])
    skipline()
    %hh = dyn_waitbar(0,['RWGMH (' int2str(j) '/' int2str(options_.mh_nblck) ')...']);
    %set(hh,'Name','RWGMH');

    for b = fblck:nblck,
        jloop = jloop+1 ; 
        %loop sur les chaines
        %if whoiam
        %    prc0=(b-fblck)/(nblck-fblck+1)*(isoctave || options_.console_mode);
        %    hh = dyn_waitbar({prc0,whoiam,options_.parallel(ThisMatlab)},['RGWMH (' int2str(b) '/' int2str(options_.mh_nblck) ')...']);
        %else
        %    hh = dyn_waitbar(0,['RWGMH (' int2str(j) '/' int2str(options_.mh_nblck) ')...']);
        %    set(hh,'Name','RWGMH');
        %end
        
        par = RWGMH_proposal(ix2,b,options_,mh_bounds) ;

        if all( par(:) > mh_bounds.lb ) && all( par(:) < mh_bounds.ub )
            try
                logpost = -feval(TargetFun, par(:),dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_);
            catch
                logpost = -inf;
            end
        else
            logpost = -inf;
        end
        if (logpost > -inf) && ( log(rand) < logpost-logpo2(j,b) )
            x2(j,:,b) = par;
            logpo2(j+1,b) = logpost;
            ix3(b,:) = par ;
            isux = isux + 1;
            jsux(b) = jsux(b) + 1;
            if logpost>logpost1 
                 disp('New posterior mode found') 
                 disp('Parameters') 
                 format long 
                 disp(par')
                 disp(['Posterior log-density = ' num2str(logpost)])
                 skipline()
                 logpost1 = logpost ;
                 xparam1 = par ;
            end
        else
            x2(j,:,b) = ix2(b,:);
            logpo2(j+1,b) = logpo2(j,b);
            ix3(b,:) = ix2(b,:);
        end
    end
    %prtfrc = j/nruns(b);
    %if (mod(j, 3)==0 && ~whoiam) || (mod(j,50)==0 && whoiam)
    %    dyn_waitbar(prtfrc,hh,[ 'RWGMH (' int2str(b) '/' int2str(options_.mh_nblck) ') ' sprintf('Current acceptance ratio %4.3f', isux/j)]);
    %end
    if search_for_convergence && j>2 
        teta_bar_m = mean(x2(1:j,:,:),1) ;
        teta_bar = mean(teta_bar_m,3) ;
        teta_m = bsxfun(@minus,teta_bar_m,teta_bar) ;
        temp = bsxfun(@minus,x2(1:j,:,:),teta_bar_m) ;
        temp2 = permute(temp,[2 1 3]) ;
        W = zeros(npar,npar,nblck) ;
        B_N = zeros(npar,npar,nblck) ;
        for jj=1:nblck
            W(:,:,jj) = temp2(:,:,jj)*temp(:,:,jj) ;
            B_N(:,:,jj) = teta_m(:,:,jj)'*teta_m(:,:,jj) ;
        end    
        W = sum(W,3)/(nblck*(j-1)) ;
        B_N = sum(B_N,3)/(nblck-1) ;
        lambdamax = max(abs(eig(W\B_N))) ;
        R_indicator(j) = (j-1)/j + (1+1/nblck)*lambdamax ;
        if j>40 
            stock_acf = zeros(40,npar) ; 
            for lgg=1:nblck 
              sel = x2(1:j,:,lgg) ;
              temp = zeros(40,npar) ;
              for pp=1:npar
                  temp(:,pp) = acf(sel(:,pp),40) ;
              end
            stock_acf = stock_acf + temp ;
            end
            stock_acf = stock_acf/nblck ;
        end
    end     
    ix2 = ix3 ;
    format short ;
    if search_for_convergence && j>2
        disp(['Global R indicator: '  num2str(R_indicator(j))]) 
    end
    skipline() 
    if nblck<16
        disp('Current parameters') 
        disp(ix2') 
        disp('Associated posterior log-density')
        disp(logpo2(j+1,:)) 
        %disp('Variance') 
        %disp(std(logpo2(1:j+1,:),1)) 
        if search_for_convergence && j>40 
            disp('Averaged autocorrelation #40 in chains') 
            disp(stock_acf(40,:)) 
        end
        disp('Chain acceptance ratio') 
        disp(jsux*100/j)
    end
    disp(['Global acceptance ratio: ' num2str(100*sum(jsux)/(nblck*j)) '%'])
    if search_for_convergence && j>5 && R_indicator(j)<1.1 %&& all(abs(stock_acf(40,:)))<.5 && j>42 
        skipline() 
        disp('Convergence achevied.') 
        disp('Now begin the iterations to build the posterior distribution...') 
        search_for_convergence = 0 ;
        x2 = zeros(nruns(1),npar,nblck);
        logpo2 = [ilogpo2';zeros(nruns(1),nblck)];
        j = 0 ;
        irun = 0 ;
        jsux = zeros(1,nblck);
    end
    irun = irun + 1;
    j=j+1;
end
all_x2 = x2 ;
all_logpo2 = logpo2 ;
for b = fblck:nblck,
    %if (irun == InitSizeArray(b)) || (j == nruns(b)) % Now I save the simulations
        x2 = all_x2(:,:,b) ;
        logpo2 = all_logpo2(:,b) ;
        save([BaseName '_mh' int2str(NewFile(b)) '_blck' int2str(b) '.mat'],'x2','logpo2');
        fidlog = fopen([MetropolisFolder '/metropolis.log'],'a');
        fprintf(fidlog,['\n']);
        fprintf(fidlog,['%% Mh' int2str(NewFile(b)) 'Blck' int2str(b) ' (' datestr(now,0) ')\n']);
        fprintf(fidlog,' \n');
        fprintf(fidlog,['  Number of simulations.: ' int2str(nruns(b)) '\n']);
        fprintf(fidlog,['  Acceptance ratio......: ' num2str(jsux(b)/nruns(b)) '\n']);
        fprintf(fidlog,['  Posterior mean........:\n']);
        for i=1:npar 
            fprintf(fidlog,['    params:' int2str(i) ': ' num2str(mean(all_x2(:,i,b))) '\n']);
        end
        fprintf(fidlog,['    log2po:' num2str(mean(all_logpo2(:,b))) '\n']);
        fprintf(fidlog,['  Minimum value.........:\n']);
        for i=1:npar
            fprintf(fidlog,['    params:' int2str(i) ': ' num2str(min(all_x2(:,i,b))) '\n']);
        end
        fprintf(fidlog,['    log2po:' num2str(min(all_logpo2(:,b))) '\n']);
        fprintf(fidlog,['  Maximum value.........:\n']);
        for i=1:npar
            fprintf(fidlog,['    params:' int2str(i) ': ' num2str(max(all_x2(:,i,b))) '\n']);
        end
        fprintf(fidlog,['    log2po:' num2str(max(all_logpo2(:,b))) '\n']);
        fprintf(fidlog,' \n');
        fclose(fidlog);
 %       jsux = 0;
 %       if j == nruns(b)+1 % I record the last draw...
            record.LastParameters(b,:) = all_x2(end,:,b);
            record.LastLogPost(b) = all_logpo2(end,b);
 %       end
        % size of next file in chain b
%        InitSizeArray(b) = min(nruns(b)-j,MAX_nruns);
        % initialization of next file if necessary
%         if InitSizeArray(b)
%             x2 = zeros(InitSizeArray(b),npar);
%             logpo2 = zeros(InitSizeArray(b),1);
%             NewFile(b) = NewFile(b) + 1;
%             irun = 0;
%         end
    %end
end    
%        dyn_waitbar_close(hh);
record.AcceptanceRatio = jsux/j;
[record.LastSeeds(b).Unifor, record.LastSeeds(b).Normal] = get_dynare_random_generator_state();
OutputFileName(jloop,:) = {[MetropolisFolder,filesep], [ModelName '_mh*_blck' int2str(b) '.mat']};
myoutput.record = record;
myoutput.irun = irun;
myoutput.NewFile = NewFile;
myoutput.OutputFileName = OutputFileName;
