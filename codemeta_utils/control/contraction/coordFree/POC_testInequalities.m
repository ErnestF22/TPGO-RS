% Generate the inequality plots for data_contractionOpt_20190307_001958.mat
% for some kd

clear all; close all;

load('data_contractionOpt_20190307_001958.mat');

%% Beta = 0
% Generate plots for kd = 100 where opt problem is infeasible near kv=61.65
opts = {'minm2',0.013,'maxm2',0.017,'minm3',0,'maxm3',0.005,'maxsteps',200};
POC_EvalTermsContraction(100,[60,61.8,62,70],maxD,magW,opts{:});

% Generate plots for kd=61.65 where opt problem is infeasible near kv=42.49
opts = {'minm2',0.018,'maxm2',0.023,'minm3',0,'maxm3',0.005,'maxsteps',200};
POC_EvalTermsContraction(kd_list(62),[40,kv_list(43),44,48],maxD,magW,opts{:});

% Generate plots for kd=4.136 where opt problem is infeasible near kv=22.3
opts = {'minm2',.04,'maxm2',.09,'minm3',.002,'maxm3',.0025,'maxsteps',200};
POC_EvalTermsContraction(kd_list(5),[20,22,23,27],maxD,magW,opts{:});

% Generate plots for kd = 0.5 where opt problem is infeasible near kv=65.41
% opts = {'minm2',0.01425,'maxm2',0.01435,'minm3',0.00019,'maxm3',0.00021,'maxsteps',200};
opts = {'minm2',0.01510,'maxm2',0.01520,'minm3',0.0002,'maxm3',0.00024,'maxsteps',200};
POC_EvalTermsContraction(0.5,[65.4;66;70],maxD,magW,opts{:});

%% Beta small
opts = {'minm2',0.013,'maxm2',0.017,'minm3',0,'maxm3',0.005,...
    'maxsteps',200,'beta',1e-2};
POC_EvalTermsContraction(100,[70],maxD,magW,opts{:});