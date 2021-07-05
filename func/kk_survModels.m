function mdls = kk_survModels(dur,stateName,doSanity)
%function mdls = kk_survModels(dur,doSanity)
%description: survival models from duration data
%
% INPUT:
%   <dur> duration data
%   <stateName> label for state
%   <doSanity> plot sanity checks
%
% OUTPUTS:
%   <mdls> stats structure
%
% USAGE:
%   >> mdls = kk_survModels(dur,doSanity)

%% REFS:
%{
https://www.mathworks.com/help/stats/survival-analysis.html
https://www.mathworks.com/help/stats/model-data-using-the-distribution-fitting-tool.html
%expdist: https://www.mathworks.com/help/predmaint/ref/reliabilitysurvivalmodel.html#d119e20422
theta: a*exp(-theta*x)
%weibll:
a	Scale parameter	a>0
b	Shape parameter	b>0

%}

%% TO DO
%{
%-add pareto curve fits




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bivariate survival models: needs X and Y, not just dur...
%-add cox regression (coxphfit):
%exp degrade: https://www.mathworks.com/help/predmaint/ref/exponentialdegradationmodel.html#d119e14762
%https://www.mathworks.com/help/stats/lifetime-of-light-bulbs.html
%covar: https://www.mathworks.com/help/predmaint/ref/covariatesurvivalmodel.html#mw_7b1b3e0b-1e9f-4f8e-a19d-54e5c9443c13
%}


%% test case
%{
dur = exprnd(10,[50 1]);
dur = minutes(dur);
mdls = kk_survModels(dur);
%}

%% params
disp('calling kk_survModels')
%sanity check duration
if ~isduration(dur)
    error('ERROR: missing duration data')
end
if ~exist('stateName','var') || isempty(stateName)
    stateName = 'default state';
end
if ~exist('doSanity','var') || isempty(doSanity)
    doSanity = 1;
end

%force to minutes
dur = duration(dur,'Format','m');
mdls.classofduration = class(dur);

%%
%%
%% force to minutes
%convert to double
X = minutes(dur);

%% fit distribution: exponential
%generate empirical xy
[y,x] = ecdf(X,'function','survivor');
%define explicit fitting func
mdls.exp.fitfunc = fittype('a*exp(-theta*x)');
%check finite
%https://www.mathworks.com/matlabcentral/answers/169861-help-with-error-in-fit-re-nan-and-inf-values
idx = isfinite(x) & isfinite(y);

%fit to exp function
[mdls.exp.mdl,mdls.exp.gof,mdls.exp.fitout] = fit(x(idx),y(idx),mdls.exp.fitfunc,'StartPoint',[0 0]);

%% get proportion at quantiles
mdls.tau = [25,50,75];
mdls.prctle_dur = prctile(x,mdls.tau);
mdls.iqr = iqr(x);

%interp the intercept
intercept = [5,10]; %in min
%mdls.intercept5_10 = interp1(x(2:end),y(2:end),intercept,'spline'); % Spline interpolation using not-a-knot end conditions interpolation
mdls.intercept5_10 = interp1(x(2:end),y(2:end),intercept,'pchip'); % Spline interpolation using not-a-knot end conditions interpolation

%% sanity plot: exponential
if doSanity
    figure,
    subplot(121),cdfplot(X); title(sprintf('ecdf: %s',stateName))
    subplot(122),plot(x,y); title(sprintf('exp fit: %s',stateName))
    hold on,
    plot(mdls.exp.mdl),
    text(mean(x)+std(x),0.50,sprintf('exp%s=%.3f','\theta',mdls.exp.mdl.theta),'vert','bottom','horiz','center','Color',[0.4 0.4 0.4]);
    %add intercept
    plot(intercept,mdls.intercept5_10,'ko');   % Plot interpolated points on ECDF
    plot([intercept(1) intercept(1)],[0 mdls.intercept5_10(1)],'Color',[0.4 0.4 0.4])
    plot([intercept(2) intercept(2)],[0 mdls.intercept5_10(2)],'Color',[0.4 0.4 0.4])
    legend({'data','expfit','5&10min intercept'})
    hold off,
    ylabel('Survival probability');
    xlabel('X (min.)')
end

%% weibull
%define explicit fitting func
mdls.wei.fitfunc = fittype('weibull');
%check finite
%https://www.mathworks.com/matlabcentral/answers/169861-help-with-error-in-fit-re-nan-and-inf-values
idx = isfinite(x) & isfinite(y);

%fit it
[mdls.wei.mdl,mdls.wei.gof,mdls.wei.fitout] = fit(x(idx),y(idx),mdls.wei.fitfunc,'StartPoint',[1 1]);

if doSanity
    figure,
    plot(x,y)
    hold on,
    plot(mdls.wei.mdl)
    text(mean(x)+std(x),0.50,sprintf('b|%s=%.3f, a|%s=%.3f','\kappa',mdls.wei.mdl.b,'\lambda',mdls.wei.mdl.a),'vert','bottom','horiz','center','Color',[0.4 0.4 0.4]);
    %add intercept
    plot(intercept,mdls.intercept5_10,'ko');   % Plot interpolated points on ECDF
    plot([intercept(1) intercept(1)],[0 mdls.intercept5_10(1)],'Color',[0.4 0.4 0.4])
    plot([intercept(2) intercept(2)],[0 mdls.intercept5_10(2)],'Color',[0.4 0.4 0.4])
    legend({'data','weibullfit'})
    hold off,
    ylabel('Survival probability');
    xlabel('X (min.)')
    title(sprintf('weibull fit: %s',stateName))
end

%% distributionFitter into table
%assemble table of models

mdls.T = table([fitdist(X,'Exponential');...
    fitdist(X,'Weibull');...
    fitdist(X,'Gamma');...
    fitdist(X,'Lognormal')],...
    'VariableNames',{'mdl'},'RowNames',{'exp','wei','gam','logn'});

%% compute errors on distfit
%negative log likelihood
mdls.T.nll = arrayfun(@negloglik,mdls.T.mdl);

%%
%%
%% fit stats
%https://www.mathworks.com/help/curvefit/evaluating-goodness-of-fit.html
%for exponential
mdls.exp.formula = formula(mdls.exp.mdl);
mdls.exp.ci = confint(mdls.exp.mdl);
%get fitoptions
mdls.exp.opts = fitoptions(mdls.exp.fitfunc);

%for weibull
mdls.wei.formula = formula(mdls.wei.mdl);
mdls.wei.ci = confint(mdls.wei.mdl);
%get fitoptions
mdls.wei.opts = fitoptions(mdls.wei.fitfunc);

%%
%%
%% create surival model/remaining useful life
mdls.RUL.mdl = reliabilitySurvivalModel("Exponential","LifeTimeUnit","minutes");
fit(mdls.RUL.mdl,dur,"minutes");

%% plot distribution
[mdls.RUL.est,mdls.RUL.ci,mdls.RUL.pdf] = predictRUL(mdls.RUL.mdl);

if doSanity
    figure,
    bar(mdls.RUL.pdf.RUL,mdls.RUL.pdf.ProbabilityDensity)
    ylabel('Survival probability');
    xlabel('X (min.)')
    title(sprintf('surival RUL: %s',stateName))
end

%%
%% end
end


