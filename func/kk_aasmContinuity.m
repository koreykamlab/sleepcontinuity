function [stats,fh] = kk_aasmContinuity(states,doPlot,doSave,fname,saveDir,epoch)
%function [stats,fh] = kk_aasmContinuity(states,doPlot,doSave,fname,saveDir,epoch)
% description: compute continuity of states from hypnogram
%
% default: <nStates> number of possible states (5 = AASM), not an input
%**huAASM states: 0-W,1-N1,2-N2,3-N3,5-REM
% Stage W (Wakefulness)
% Stage N1 (NREM 1)
% Stage N2 (NREM 2)
% Stage N3 (NREM 3)
% Stage R (REM)
% Stage 'S' (all sleep)
%
%
% INPUT:
%   <states>: vector of hypnogram called 'states'. for consistency to
%   mouse, keep as row.
%   <doPlot>: true/false
%   <doSave>: true/false
%   <fname>: if doSave true, give filename
%   <saveDir>: path of directory to save to
%   <epoch>: epoch length, default is 30sec for hu, 1sec for ms
%
% OUTPUT:
%   <stats>: strucutre containing continuity metrics/table of models
%   <fh>: store the computed figures
%
% USAGE:
%   >> [stats] = kk_aasmContinuity(states,doPlot,doSave,fname,saveDir,epoch)


%% TO DO
%{
%-add pareto curve fits
%-add timestamps for transitions (N-R, R-W, N-W)
% handle states when missing








%}

%% test case
%{
%% numeric example
fname = 'huHypnogram_AASMnum.mat';
load(fname)
saveDir = 'd:\GitHub\kk_codeRepo\_exampleData\hypnos';
tic
[stats] = kk_aasmContinuity(states,1,1,fname,saveDir)
toc
stats.stateStatsT.survmdls(6,1).exp.mdl.theta
stats.stateStatsT.first3min_epoch


%% example w/o REM
fname = '140122VA-d_hyp.mat';
saveDir = '/Users/koreykam/Google Drive/_DATA/hu/masraiHypnograms';
load(fname)
[stats] = kk_aasmContinuity(states,0,1,fname,saveDir)

%% string example



%}

%% handle defaults
disp('calling kk_aasmContinuity')

%dbstop if warning %enable for troubleshooting
if ~exist('states','var') || isempty(states)
    error('ERROR: states hypnogram vector does NOT exist')
end
if ~exist('doPlot','var') || isempty(doPlot)
    doPlot = 1;
end
if ~exist('doSave','var') || isempty(doSave)
    doSave = 1;
end
if doSave && ~exist('fname','var') || isempty(fname)
    error('ERROR: fname is missing');
end
if doSave && ~exist('saveDir','var')
    %[~,~,saveDir] = kk_determineComp;
    error('ERROR: no saveDir selected!')
end
if ~exist('epoch','var') || isempty(epoch)
    epoch = 30; %epoch length in seconds for AASM
end
if ~exist('thresh','var') || isempty(thresh)
    thresh = [180 300]; %in sec = 3 and 5 min
end
%sanity thresh
if ~le(thresh(1),thresh(2))
    thresh = fliplr(thresh); %if false, flip
end

%initialize fig count
fhcount = 0;
fh = [];

%% force to column vector
if isrow(states)
    states = states(:); %force to row vector
    %force row>> reshape(x, 1, prod(size(x)))
    disp('Converted input to COLUMN vector')
elseif iscolumn(states)
    disp('''states'' is a COLUMN vector')
else
    error('incorrect input size, expected ROW vector')
end

%% get state ID/names
stateID = unique(states);
disp('stateIDs'),celldisp({num2str(stateID)})
state_possible = [0;1;2;3;5];
disp('possible states'),celldisp({num2str(state_possible)})
stateNames = {'W','N1','N2','N3','REM','SLEEPall','NREMall'};

%% sanity check existing states
if ~all(ismember(states,[0,1,2,3,5])) %'states' can also NOT have all classes
    error('ERROR: unexpected state selection')
elseif ~all(ismember([0,1,2,3,5],states))
    warning('NOT all classes exist in "states"')
    stats.missingState = true;
end

%% force 9s to 1s 
try
    %find idx of zeros
    [idx] = find(states == 9);
    if ~isempty(idx)
        %assign as wake
        states(idx) = 0;
        fprintf('epoch# %i now wake\n',idx)
    end
catch ME
    disp('no unscored')
end

%% calc state logicals
stateLogic = zeros(length(states),length(stateNames));
stateDur_min = zeros(length(stateNames),1);

%logicalindex loop through hypnogram
for ii = 1:length(state_possible)
    stateLogic(:,ii) = states == state_possible(ii); %see kk_msContinuity if not all stateIDs exist in hypnogram, assign to length of stateNames
    
    %calc duration in min
    try
        stateDur_min(ii,1) = nnz(stateLogic(:,ii));
    catch
        stateDur_min(ii,1) = 0;
        fprintf('NOT done @stateLogical %s\n',state_possible(ii))
    end
    %tell me
    fprintf('done @stateLogical %i\n',state_possible(ii))
    if isequal(stateDur_min(ii,1),0)
        warning('...state %i has 0 duration...',state_possible(ii))
    end
end

%handle all SLEEP case
stateLogic(:,end-1) = states ~= state_possible(1); %find where NOT wake -> SLEEP
stateDur_min(end-1,1) = nnz(stateLogic(:,end-1)); %sum duration of SLEEP
fprintf('done @stateLogical %s\n',stateNames{end-1})

%handle all NREM case
stateLogic(:,end) = states ~= state_possible(1) & states ~= state_possible(5); %find where NOT wake AND NOT REM
stateDur_min(end,1) = nnz(stateLogic(:,end)); %sum duration of SLEEP
fprintf('done @stateLogical %s\n',stateNames{end})

%assign logical matrix as table
tmpT = table(stateLogic);
stateLogicT = splitvars(tmpT,'stateLogic','NewVariableNames',stateNames);
clearvars tmpT stateLogic

%assign state stats as table
stateStatsT = table(stateDur_min,'VariableNames',{'dur_min'},'RowNames',stateNames);

%% sanity check total duration of sleep
if ~isequal(sum(stateStatsT.dur_min(2:5)),stateStatsT.dur_min(6))
    error('all sleep stateLogical does not match sum of individual sleepLogicals')
else
    disp('duration sum PASSED...')
end

%% calc proportion of states
stateStatsT.proportion = stateStatsT.dur_min/sum(stateStatsT.dur_min(1:end-2));

if doPlot
    fhcount = fhcount + 1;
    fh(fhcount) = figure;
    pie(stateStatsT.proportion(1:end-2)) %except SLEEPall & NREMall
    legend(stateNames(1:end-2),'Location','southoutside','Orientation','horizontal')
    title('proportion of each state')
end

%% calculate duration of each state
%ref: https://www.mathworks.com/help/stats/nonparametric-and-empirical-probability-distributions.html?searchHighlight=ecdf&s_tid=doc_srchtitle

%temp cell containers for ecdf
ecdfF_cell = cell(1,width(stateLogicT)); %F out of ecdf
ecdfT_cell = cell(1,width(stateLogicT)); %T out of ecdf

%find length/indices of all detected regions
for ii = 1:width(stateLogicT)
    %find transitions, pre/post-pend 0 to COLUMN vector
    dtctsig = diff([0; stateLogicT{:,ii}; 0]);
    %detect positive diff
    startIdx = find(dtctsig > 0);
    %detect negative diff -1 sample
    endIdx = find(dtctsig < 0)-1;
    %sanity that length of indices match
    if ~isequal(length(startIdx),length(endIdx))
        %FIX THIS ERROR if possible: drop last endIdx?
        error('transitions of hypnogram do not match')
    end
    
    %calc duration b/w start-end epoch, multiply difference by epoch length
    %to get dur in min
    %must store in cell to input into table
    disp('make dur_vec in min for survival modeling')
    stateStatsT.dur_vec(ii,1) = {(endIdx-startIdx+1) * epoch ./ 60}; %convert to min
    stateStatsT.startIdx(ii,1) = {startIdx}; %store start idx
    stateStatsT.endIdx(ii,1) = {endIdx}; %store end idx
    
    %sanity check if dur_vec empty prior to ecdf
    if isempty(stateStatsT.dur_vec{ii,1})
        stateStatsT.dur_vec{ii,1} = 0;
        warning('dur_vec IS empty @ %s',stateNames{ii})
    end
    if isempty(stateStatsT.startIdx{ii,1})
        stateStatsT.startIdx{ii,1} = 0;
        warning('startIdx IS empty @ %s',stateNames{ii})
    end
    if isempty(stateStatsT.endIdx{ii,1})
        stateStatsT.endIdx{ii,1} = 0;
        warning('endIdx IS empty @ %s',stateNames{ii})
    end
    
    %compute ecdf and store in cell
    [ecdfF_cell{1,ii},ecdfT_cell{1,ii}] = ecdf(stateStatsT.dur_vec{ii,1},...
        'function','survivor','alpha',0.05);
end

if doPlot
    %do CDF survival
    fhcount = fhcount + 1;
    fh(fhcount) = figure;
    s = zeros(1,width(stateLogicT));
    for ii = 1:width(stateLogicT)
        s(ii) = subplot(1,width(stateLogicT),ii);
        ecdf(stateStatsT.dur_vec{ii,1},'function','survivor','alpha',0.05,'bounds','on');
        hold on,
        title(sprintf('%s',stateNames{ii})),
        xlabel(''); ylabel('');
        %add reflines
        %refline([0 0.5])
        Med = median(stateStatsT.dur_vec{ii,1});
        Mea = mean(stateStatsT.dur_vec{ii,1});
        plot([Med Med],[0 0.5],'Color',[0.2 0.2 0.2],'LineStyle','--');
        text(Med,0.51,sprintf('median=%.1f',Med),'vert','bottom','horiz','center','Color',[0.4 0.4 0.4]);
        text(Mea,0.55,sprintf('mean=%.1f',Mea),'vert','bottom','horiz','center','Color',[0.4 0.4 0.4]);
        hold off,
        %delete(subplot(1,6,5)) %figure out how to redraw with deleted
        %subplot
    end
    legend({'CDF','CI','CI','med'})
    sgtitle('continuity of states') %add superior title
    ylabel(s(1),'survival prob.'),
    xlabel(s(1),'duration (min)'), %want this in minutes
end

%% calculate percentiles in units of duration
tau = [25, 50, 75]; %desired percentiles
for ii = 1:width(stateLogicT)
    stateStatsT.prctle(ii,1) = {prctile(stateStatsT.dur_vec{ii,1},tau)}; %store in cell
    stateStatsT.median(ii,1) = {median(stateStatsT.dur_vec{ii,1})}; %store in cell
    stateStatsT.iqr(ii,1)    = {iqr(stateStatsT.dur_vec{ii,1})}; %store in cell
    stateStatsT.mean(ii,1)   = {mean(stateStatsT.dur_vec{ii,1})}; %store in cell
    stateStatsT.std(ii,1)    = {std(stateStatsT.dur_vec{ii,1})}; %store in cell
end

%% calculate transitions
%by loop
for ii = 1:width(stateLogicT)
    stateStatsT.transitionsto(ii,1)   = sum(diff(stateLogicT{:,ii}) == 1);
    stateStatsT.transitionsfrom(ii,1) = sum(diff(stateLogicT{:,ii}) == -1);
    firstEpoch = find(stateLogicT{:,ii}, 1, 'first');
    %handle empty case
    if isempty(firstEpoch)
        firstEpoch = 0;
        warning('firstEpoch IS empty @ %s',stateNames{ii})
    end
    
    %latencies
    stateStatsT.latencyto_min(ii,1)   = firstEpoch * epoch / 60; %calc latency to state in minutes
    try
        stateStatsT.first3min_epoch(ii,1)   = stateStatsT.startIdx{ii,1}(find(stateStatsT.dur_vec{ii,1}>=3,1,'first'));
    catch
        stateStatsT.first3min_epoch(ii,1) = 0;
    end
    try
        stateStatsT.first5min_epoch(ii,1)   = stateStatsT.startIdx{ii,1}(find(stateStatsT.dur_vec{ii,1}>=5,1,'first'));
    catch
        stateStatsT.first5min_epoch(ii,1) = 0;
    end
    try
        stateStatsT.first10min_epoch(ii,1)   = stateStatsT.startIdx{ii,1}(find(stateStatsT.dur_vec{ii,1}>=10,1,'first'));
    catch
        stateStatsT.first10min_epoch(ii,1) = 0;
    end
    
end

%% quant of transition pairs: eg N3-R, N2-R etc (NOT DONE)

%% calc proportion > or < threshold in units of duration
for ii = 1:width(stateLogicT)
    numEvents = length(stateStatsT.dur_vec{ii,1});
    idxBelow = find(stateStatsT.dur_vec{ii,1} < thresh(1)/60); %less than 'low' thresh, convert thresh to min
    idxAbove = find(stateStatsT.dur_vec{ii,1} > thresh(2)/60); %greater than 'high' thresh, convert thresh to min
    numBelow = length(idxBelow);
    numAbove = length(idxAbove);
    stateStatsT.propBelowThresh(ii,1) = numBelow/numEvents;
    stateStatsT.propAboveThresh(ii,1) = numAbove/numEvents;
end

if doPlot
    fhcount = fhcount + 1;
    fh(fhcount) = figure;
    for ii=1:width(stateLogicT)
        s(ii) = subplot(1,width(stateLogicT),ii);
        c = categorical({'<lowThres','>highThresh','w/inThresh'});
        Y = [stateStatsT.propBelowThresh(ii,1),stateStatsT.propAboveThresh(ii,1),1-stateStatsT.propBelowThresh(ii,1)-stateStatsT.propAboveThresh(ii,1)];
        bar(c,Y);
        text(1:length(Y),Y,num2str(round(Y',2)),'vert','bottom','horiz','center');
        xlabel(''); ylabel('');
        title(sprintf('%s',stateNames{ii})),
        ylim([0 1])
        box off
    end
    ylabel(s(1),'proportion @criterion'),
    sgtitle(sprintf('<lowThres=%i min, \n>highThresh=%i min',thresh(1)/60,thresh(2)/60))
end

%% arousal detection
%bouts b/w 1 and 10 sec with sleep preceeding and succeeding
%find > 1 & < 10
threshMinutes = [1,180]/60;

idx = find(stateStatsT.dur_vec{1,1} > threshMinutes(1) & stateStatsT.dur_vec{1,1} < threshMinutes(2)); %less than 'low' thresh, convert thresh to min
arousalCount = length(idx);

%%
%%
%%
%%
%% survival modeling (DONE)
%loop through all states
for ii = 1:width(stateLogicT)
    try
        doSanity = doPlot;
        %determine if 0
        if isequal(stateStatsT.dur_vec{ii,1},0)
            %fill as NaN
            stateStatsT.survmdls(ii,1).exp.mdl.theta = NaN;
            warning('unable to model states w/ 0 epochs...')
            
        else
            stateStatsT.survmdls(ii,1) = kk_survModels(minutes(stateStatsT.dur_vec{ii,1}),stateNames{ii},doSanity);
            if any(structfun(@isempty, stateStatsT.survmdls(ii,1)))
                warning('empty struct of survmodel object for state %i',ii)
                %find empties in struct and replace with NaN
                empties = structfun(@isempty,stateStatsT.survmdls(ii,1));
                %get fieldnames
                fnames=fieldnames(stateStatsT.survmdls(ii,1));
                %assign empties as NaN
                stateStatsT.survmdls(ii,1) = structfun(@(x) NaN, stateStatsT.survmdls(ii,1), 'UniformOutput', false);
            end
            %grab theta parameter of expfit
            fprintf('%s expfit theta = %.3f\n',stateNames{ii},stateStatsT.survmdls(ii,1).exp.mdl.theta)
            %grab gof errors for expfit
            stateStatsT.expgof(ii,1) = stateStatsT.survmdls(ii,1).exp.gof;
            if any(structfun(@isempty, stateStatsT.expgof(ii,1)))
                stateStatsT.expgof(ii,1) = structfun(@(x) NaN, stateStatsT.expgof(ii,1), 'UniformOutput', false);
            end
            %grab gof errors for weifit
            stateStatsT.weigof(ii,1) = stateStatsT.survmdls(ii,1).wei.gof;
            if any(structfun(@isempty, stateStatsT.weigof(ii,1)))
                stateStatsT.weigof(ii,1) = structfun(@(x) NaN, stateStatsT.weigof(ii,1), 'UniformOutput', false);
            end
        end
        
    catch ME
        %rethrow(ME)
        warning('%s',ME.identifier)
    end
end

%%
%%
%%
%% organize output to struct
stats.executiontime = strrep(datestr(now),':','-');
stats.stateStatsT = stateStatsT;
stats.stateLogicT = stateLogicT;
stats.stateID = stateID;
stats.state_possible = state_possible;
stats.stateNames = stateNames;
stats.thresh = thresh;
stats.arousalCount = arousalCount;

%% save all to disk
if doSave
    [~,basename,~] = fileparts(fname);
    save(fullfile(saveDir,sprintf('%s_continuity',basename)),'stats','-v7.3');
    warning('SAVED continuity struct for file: %s ',fname)
    
    %get all figs
    %save them as one .fig file
    savefig(fh,fullfile(saveDir,sprintf('%s_continuity',basename)),'compact')
    warning('SAVED continuity figs for file: %s ',fname)
end

%% end
end
