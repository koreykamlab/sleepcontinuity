
%% single sub example
addpath(genpath('~\sleepcontinuity'))
fname = 'huHypnogram_AASMnum.mat';
load(fname)
doPlot=1;
doSave=1;
saveDir = fullfile(getenv('USERPROFILE'), 'Downloads');
epoch=30;
tic
[stats,fh] = kk_aasmContinuity(states,doPlot,doSave,fname,saveDir,epoch);
toc

%% checks
display(stats.stateStatsT)
figure, stairs(states)

%ex: get continuity of all sleep states
stats.stateStatsT.median{6,1}
stats.stateStatsT.survmdls(6,1).exp.mdl.theta

%ex: states and latency to first 3min bout of each state
stats.stateNames
stats.stateStatsT.first3min_epoch

%% consider grp lvl example
load('huHypnogram_AASMnum_grp.mat')

%define Ns per group (varies for each group)
nA=5;
nB=5;

%define median bout# of desired state (varies for each group)
nbootA = 17;
nbootB = 16;

%do boot
bootA = randsample(grptoyDat.A.N3,nbootA*nA,1);
bootB = randsample(grptoyDat.B.N3,nbootB*nB,1);

[h,p] = kstest2(bootA,bootB)

%% 