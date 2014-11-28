function [ ix2, ilogpo2, ModelName, MetropolisFolder, fblck, fline, npar, nblck, nruns, NewFile, MAX_nruns ] = ...
    RWGMH_initialization(TargetFun, xparam1, mh_bounds,dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_)
% Random Walk Genetic Metropolis-Hastings initialization.
% 
% INPUTS 
%   o TargetFun  [char]     string specifying the name of the objective
%                           function (posterior kernel).
%   o xparam1    [double]   (p*1) vector of parameters to be estimated (initial values).
%   o vv         [double]   (p*p) matrix, posterior covariance matrix (at the mode).
%   o mh_bounds  [double]   (p*2) matrix defining lower and upper bounds for the parameters. 
%   o dataset_              data structure
%   o options_              options structure
%   o M_                    model structure
%   o estim_params_         estimated parameters structure
%   o bayestopt_            estimation options structure
%   o oo_                   outputs structure
%  
% OUTPUTS 
%   None  
%
% SPECIAL REQUIREMENTS
%   None.

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

ix2 = [];
ilogpo2 = [];
ModelName = []; 
MetropolisFolder = [];
fblck = [];
fline = [];
npar = [];
nblck = [];
nruns = [];
NewFile = [];
MAX_nruns = [];

ModelName = M_.fname;
if ~isempty(M_.bvar)
    ModelName = [ModelName '_bvar'];
end

MetropolisFolder = CheckPath('metropolis',M_.dname);
BaseName = [MetropolisFolder filesep ModelName];

nblck = options_.mh_nblck;
nruns = ones(nblck,1)*options_.mh_replic;
npar  = length(xparam1);
MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);

if ~options_.load_mh_file && ~options_.mh_recover
    % Here we start a new metropolis-hastings, previous draws are discarded.
    disp('Estimation::RWGMH: Multiple chains mode.')
    % Delete old mh files if any...
    files = dir([BaseName '_rwgmh*_blck*.mat']);
    if length(files)
        delete([BaseName '_rwgmh*_blck*.mat']);
        disp('Estimation::mcmc: Old rwgmh-files successfully erased!')
    end
    % Delete old metropolis log file.
    file = dir([ MetropolisFolder '/metropolis.log']);
    if length(file)
        delete([ MetropolisFolder '/metropolis.log']);
        disp('Estimation::RWGMH: Old metropolis.log file successfully erased!')
        disp('Estimation::RWGMH: Creation of a new metropolis.log file.')
    end
    fidlog = fopen([MetropolisFolder '/metropolis.log'],'w');
    fprintf(fidlog,'%% RWGMH log file (Dynare).\n');
    fprintf(fidlog,['%% ' datestr(now,0) '.\n']);
    fprintf(fidlog,' \n\n');
    fprintf(fidlog,'%% Session 1.\n');
    fprintf(fidlog,' \n');
    fprintf(fidlog,['  Number of blocks...............: ' int2str(nblck) '\n']);
    fprintf(fidlog,['  Number of simulations per block: ' int2str(nruns(1)) '\n']);
    fprintf(fidlog,' \n');
    % Find initial values for the nblck chains...
    fprintf(fidlog,['  Initial values of the parameters:\n']);
    disp('Estimation::RWGMH: Searching for initial values...')
    ix2 = zeros(nblck,npar);
    ilogpo2 = zeros(nblck,1);
    for j=1:nblck
        validate    = 0;
        init_iter   = 0;
        trial = 1;
        while validate == 0 && trial <= 10
            candidate = rand_multivariate_normal( transpose(xparam1), options_.rwgmh_init_scale, npar);
            if all(candidate(:) > mh_bounds(:,1)) && all(candidate(:) < mh_bounds(:,2)) 
                ix2(j,:) = candidate;
                ilogpo2(j) = - feval(TargetFun,ix2(j,:)',dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,oo_);
                if ~isfinite(ilogpo2(j)) % if returned log-density is
                                         % Inf or Nan (penalized value)
                    validate = 0;
                else
                    fprintf(fidlog,['    Blck ' int2str(j) ':\n']);
                    for i=1:length(ix2(1,:))
                        fprintf(fidlog,['      params:' int2str(i) ': ' num2str(ix2(j,i)) '\n']);
                    end
                    fprintf(fidlog,['      logpo2: ' num2str(ilogpo2(j)) '\n']);
                    j = j+1;
                    validate = 1;
                end
            end
            init_iter = init_iter + 1;
            if init_iter > 100 && validate == 0
                disp(['Estimation::RWGMH: I couldn''t get a valid initial value in 100 trials.'])
                if options_.nointeractive
                    disp(['Estimation::RWGMH: I reduce rwgmh_init_scale by ten percent:'])
                    options_.mh_init_scale = .9*options_.rwgmh_init_scale;
                    disp(sprintf('Estimation::RWGMH: Parameter rwgmh_init_scale is now equal to %f.',options_.mh_init_scale))
                else
                    disp(['Estimation::RWGMH: You should Reduce rwgmh_init_scale...'])
                    disp(sprintf('Estimation::RWGMH: Parameter rwgmh_init_scale is equal to %f.',options_.rwgmh_init_scale))
                    options_.rwgmh_init_scale = input('Estimation::mcmc: Enter a new value...  ');
                end
                trial = trial+1;
            end
        end
        if trial > 10 && ~validate
            disp(['Estimation::RWGMH: I''m unable to find a starting value for block ' int2str(j)])
            return
        end
    end
    fprintf(fidlog,' \n');
    disp('Estimation::RWGMH: Initial values found!')
    skipline()
    fprintf(fidlog,' \n');
    fblck = 1;
    fline = ones(nblck,1);
    NewFile = ones(nblck,1);
    % Delete the mh-history files
    delete_mh_history_files(MetropolisFolder, ModelName);
    %  Create a new record structure
    fprintf(['Estimation::RWGMH: Write details about the RWGMH... ']);
    AnticipatedNumberOfFiles = ceil(nruns(1)/MAX_nruns);
    AnticipatedNumberOfLinesInTheLastFile = nruns(1) - (AnticipatedNumberOfFiles-1)*MAX_nruns;
    record.Nblck = nblck;
    record.MhDraws = zeros(1,3);
    record.MhDraws(1,1) = nruns(1);
    record.MhDraws(1,2) = AnticipatedNumberOfFiles;
    record.MhDraws(1,3) = AnticipatedNumberOfLinesInTheLastFile;
    record.AcceptanceRatio = zeros(1,nblck);
    for j=1:nblck
        % we set a different seed for the random generator for each block then we record the corresponding random generator state (vector)
        set_dynare_seed(options_.DynareRandomStreams.seed+j);
        % record.Seeds keeps a vector of the random generator state and not the scalar seed despite its name
        [record.InitialSeeds(j).Unifor,record.InitialSeeds(j).Normal] = get_dynare_random_generator_state();
    end
    record.InitialParameters = ix2;
    record.InitialLogPost = ilogpo2;
    record.LastParameters = zeros(nblck,npar);
    record.LastLogPost = zeros(nblck,1);
    record.LastFileNumber = AnticipatedNumberOfFiles ;
    record.LastLineNumber = AnticipatedNumberOfLinesInTheLastFile;
    fprintf('Ok!\n');
    id = write_mh_history_file(MetropolisFolder, ModelName, record);
    disp(['Estimation::RWGMH: Details about the MCMC are available in ' BaseName '_mh_history_' num2str(id) '.mat'])
    skipline()
    fprintf(fidlog,['  CREATION OF THE MH HISTORY FILE!\n\n']);
    fprintf(fidlog,['    Expected number of files per block.......: ' int2str(AnticipatedNumberOfFiles) '.\n']);
    fprintf(fidlog,['    Expected number of lines in the last file: ' int2str(AnticipatedNumberOfLinesInTheLastFile) '.\n']);
    fprintf(fidlog,['\n']);
    for j = 1:nblck,
        fprintf(fidlog,['    Initial state of the Gaussian random number generator for chain number ',int2str(j),':\n']);
        for i=1:length(record.InitialSeeds(j).Normal)
            fprintf(fidlog,['      ' num2str(record.InitialSeeds(j).Normal(i)') '\n']);
        end
        fprintf(fidlog,['    Initial state of the Uniform random number generator for chain number ',int2str(j),':\n']);
        for i=1:length(record.InitialSeeds(j).Unifor)
            fprintf(fidlog,['      ' num2str(record.InitialSeeds(j).Unifor(i)') '\n']);
        end
    end,
    fprintf(fidlog,' \n');
    fclose(fidlog);
elseif options_.load_mh_file && ~options_.mh_recover
    % Here we consider previous mh files (previous mh did not crash).
    disp('Estimation::RWGMH: I am loading past RWGMH simulations...')
    load_last_mh_history_file(MetropolisFolder, ModelName);
    mh_files = dir([ MetropolisFolder filesep ModelName '_rwgmh*.mat']);
    if ~length(mh_files)
        error('Estimation::RWGMH: I cannot find any MH file to load here!')
    end
    fidlog = fopen([MetropolisFolder '/metropolis.log'],'a');
    fprintf(fidlog,'\n');
    fprintf(fidlog,['%% Session ' int2str(length(record.MhDraws(:,1))+1) '.\n']);
    fprintf(fidlog,' \n');
    fprintf(fidlog,['  Number of blocks...............: ' int2str(nblck) '\n']);
    fprintf(fidlog,['  Number of simulations per block: ' int2str(nruns(1)) '\n']);
    fprintf(fidlog,' \n');
    past_number_of_blocks = record.Nblck;
    if past_number_of_blocks ~= nblck
        disp('Estimation::RWGMH: The specified number of blocks doesn''t match with the previous number of blocks!')
        disp(['Estimation::RWGMH: You declared ' int2str(nblck) ' blocks, but the previous number of blocks was ' int2str(past_number_of_blocks) '.'])
        disp(['Estimation::RWGMH: I will run the Metropolis-Hastings with ' int2str(past_number_of_blocks) ' blocks.' ])
        nblck = past_number_of_blocks;
        options_.mh_nblck = nblck;
    end
    % I read the last line of the last mh-file for initialization of the new metropolis-hastings simulations:
    LastFileNumber = record.LastFileNumber;
    LastLineNumber = record.LastLineNumber;
    if LastLineNumber < MAX_nruns
        NewFile = ones(nblck,1)*LastFileNumber;
        fline = ones(nblck,1)*(LastLineNumber+1);
    else
        NewFile = ones(nblck,1)*(LastFileNumber+1);
        fline = ones(nblck,1);
    end
    ilogpo2 = record.LastLogPost;
    ix2 = record.LastParameters;
    fblck = 1;
    NumberOfPreviousSimulations = sum(record.MhDraws(:,1),1);
    fprintf('Estimation::RWGMH: I am writing a new rwgmh-history file... ');
    record.MhDraws = [record.MhDraws;zeros(1,3)];
    NumberOfDrawsWrittenInThePastLastFile = MAX_nruns - LastLineNumber;
    NumberOfDrawsToBeSaved = nruns(1) - NumberOfDrawsWrittenInThePastLastFile;
    AnticipatedNumberOfFiles = ceil(NumberOfDrawsToBeSaved/MAX_nruns);
    AnticipatedNumberOfLinesInTheLastFile = NumberOfDrawsToBeSaved - (AnticipatedNumberOfFiles-1)*MAX_nruns;  
    record.LastFileNumber = LastFileNumber + AnticipatedNumberOfFiles;
    record.LastLineNumber = AnticipatedNumberOfLinesInTheLastFile;
    record.MhDraws(end,1) = nruns(1);
    record.MhDraws(end,2) = AnticipatedNumberOfFiles;
    record.MhDraws(end,3) = AnticipatedNumberOfLinesInTheLastFile;
    record.InitialSeeds = record.LastSeeds;
    write_mh_history_file(MetropolisFolder, ModelName, record);
    fprintf('Done.\n')
    disp(['Estimation::RWGMH: Ok. I have loaded ' int2str(NumberOfPreviousSimulations) ' simulations.'])
    skipline()
    fclose(fidlog);
elseif options_.mh_recover
    % The previous metropolis-hastings crashed before the end! I try to recover the saved draws...
    disp('Estimation::RWGMH: Recover mode!')
    load_last_mh_history_file(MetropolisFolder, ModelName);
    nblck = record.Nblck;% Number of "parallel" mcmc chains.
    options_.mh_nblck = nblck;
    if size(record.MhDraws,1) == 1
        OldMh = 0;% The crashed metropolis was the first session.
    else
        OldMh = 1;% The crashed metropolis wasn't the first session.
    end
    % Default initialization:
    if OldMh
        ilogpo2 = record.LastLogPost;
        ix2 = record.LastParameters;
    else
        ilogpo2 = record.InitialLogPost;
        ix2 = record.InitialParameters;
    end
    % Set NewFile, a nblck*1 vector of integers, and fline (first line), a nblck*1 vector of integers.
    if OldMh
        LastLineNumberInThePreviousMh = record.MhDraws(end-1,3);% Number of lines in the last mh files of the previous session.
        LastFileNumberInThePreviousMh = sum(record.MhDraws(1:end-1,2),1);% Number of mh files in the the previous sessions.
        if LastLineNumberInThePreviousMh < MAX_nruns% Test if the last mh files of the previous session are not complete (yes)
            NewFile = ones(nblck,1)*LastFileNumberInThePreviousMh;
            fline = ones(nblck,1)*(LastLineNumberInThePreviousMh+1);
        else% The last mh files of the previous session are complete.
            NewFile = ones(nblck,1)*(LastFileNumberInThePreviousMh+1);
            fline = ones(nblck,1);
        end
    else
        LastLineNumberInThePreviousMh = 0;
        LastFileNumberInThePreviousMh = 0;
        NewFile = ones(nblck,1);
        fline = ones(nblck,1);
    end
    % Set fblck (First block), an integer targeting the crashed mcmc chain.
    fblck = 1;
    % How many mh files should we have ?
    ExpectedNumberOfMhFilesPerBlock = sum(record.MhDraws(:,2),1);
    ExpectedNumberOfMhFiles = ExpectedNumberOfMhFilesPerBlock*nblck;
    % I count the total number of saved mh files...
    AllMhFiles = dir([BaseName '_rwgmh*_blck*.mat']);
    TotalNumberOfMhFiles = length(AllMhFiles);
    % And I quit if I can't find a crashed mcmc chain. 
    if (TotalNumberOfMhFiles==ExpectedNumberOfMhFiles)
        disp('Estimation::RWGMH: It appears that you don''t need to use the mh_recover option!')
        disp('                   You have to edit the mod file and remove the mh_recover option') 
        disp('                   in the estimation command')
        error()
    end
    % I count the number of saved mh files per block.
    NumberOfMhFilesPerBlock = zeros(nblck,1);
    for b = 1:nblck
        BlckMhFiles = dir([BaseName '_rwgmh*_blck' int2str(b) '.mat']);
        NumberOfMhFilesPerBlock(b) = length(BlckMhFiles);
    end
    % Is there a chain with less mh files than expected ? 
    while fblck <= nblck
        if  NumberOfMhFilesPerBlock(fblck) < ExpectedNumberOfMhFilesPerBlock
            disp(['Estimation::rwgmh: Chain ' int2str(fblck) ' is not complete!'])
            break
            % The mh_recover session will start from chain fblck.
        else
            disp(['Estimation::rwgmh: Chain ' int2str(fblck) ' is complete!'])
        end
        fblck = fblck+1;
    end
    % How many mh-files are saved in this block?
    NumberOfSavedMhFilesInTheCrashedBlck = NumberOfMhFilesPerBlock(fblck);
    % Correct the number of saved mh files if the crashed metropolis was not the first session (so
    % that NumberOfSavedMhFilesInTheCrashedBlck is the number of saved mh files in the crashed chain 
    % of the current session).  
    if OldMh
        NumberOfSavedMhFilesInTheCrashedBlck = NumberOfSavedMhFilesInTheCrashedBlck - LastFileNumberInThePreviousMh;
    end
    NumberOfSavedMhFiles = NumberOfSavedMhFilesInTheCrashedBlck+LastFileNumberInThePreviousMh;
    % Correct initial conditions.
    if NumberOfSavedMhFiles
        load([BaseName '_rwgmh' int2str(NumberOfSavedMhFiles) '_blck' int2str(fblck) '.mat']);
        ilogpo2(fblck) = logpo2(end);
        ix2(fblck,:) = x2(end,:);
    end
end