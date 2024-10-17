function iccv15_plots_testOutliers(idxFig)
figDir='../../../papers/vision/ICCV15-RotationAveragingSurvey/figures';
figDimWide=[580 130];
figDimSingle=[200 130];
flagSaveFig=2;
figExt='epsc';
axXLim=80;

if ~exist('idxFig','var') || isempty(idxFig)
    idxFig=1:30;
end

dataset='fountain_herzjesu_castleentry_herzjesularge_castle_castlelarge';
%dataset='fountain_herzjesu_castleentry_herzjesularge_castlelarge';
%dataset='fountain_herzjesu_castleentry_herzjesularge';
%dataset='castle';
iFig=0;

fs=setFigFontSize(8);
fn=setFigFont('Times');

%%
iFig=iFig+1;
if ismember(iFig,idxFig)
    ax=[0 axXLim 0 100];
    figure(iFig); 
    subplot(1,3,1)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalLaplacian','optsGlobalLaplacianAdjust','optsGlobalLaplaciannormalize','optsGlobalLaplaciannormalizeAdjust'},...
        'legendNames')
    axis(ax)
    subplot(1,3,2)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalSpectral','optsGlobalSpectralAdjust','optsGlobalSpectralNormalize'},...
        'legendNames')
    axis(ax)
    subplot(1,3,3)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalSpectral','optsGlobalLaplacian'},...
        'legendNames')
    axis(ax)
    disp('# Lesson (linear methods)')
    disp('Laplacian works better with normalization, Spectral without.')
    disp('Linear adjustement does not make a significant difference.')
    disp('Spectral works better than Laplacian')
    savefigure(fullfile(figDir,'linear'),figExt,figDimWide,flagSaveFig);
end


%%
iFig=iFig+1;
if ismember(iFig,idxFig)
    ax=[0 axXLim 0 50];
    figure(iFig); 
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalSDP','optsGlobalsdpprojectBeforeFactorization','optsGlobalSDPoptsLowRankspectrahedral'},...
        'legendNames')
    axis(ax)
    disp('# Lesson (SDP trace + pre-projection)')
    disp('Spectrahedral constraints make an insignificant difference.')
    disp('Projecting the recovered measurement matrix before factorizaiton does not help.')
    savefigure(fullfile(figDir,'sdpTrace'),figExt,figDimSingle,flagSaveFig);
end

%%
iFig=iFig+1;
if ismember(iFig,idxFig)
    figure(iFig); 
    subplot(1,3,1)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalALM','optsGlobalALMoptsLowRankinnerProj0'},...
        'legendNames')
    subplot(1,3,2)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalALMoptsLowRankcriterionl1','optsGlobalALMoptsLowRankcriterionl1innerProj0'},...
        'legendNames')
    subplot(1,3,3)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalALMoptsLowRankcriterionl12','optsGlobalALMoptsLowRankcriterionl12innerProj0'},...
        'legendNames')
    disp('# Lesson (ALM)')
    disp('Adding the intermediate projection step to ALM is significatively better.')
    disp('However, there might be convergence problems at the origin.')
    disp('This happens for sparse datasets, a single lambda does not work.')
end

%%
iFig=iFig+1;
if ismember(iFig,idxFig)
    figure(iFig); 
    subplot(1,3,1)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalalmoptsLowRankcriterionl1innerProj0optsOptimizerlowrankpriorSDP','optsGlobalalmoptsLowRankcriterionl1optsOptimizerlowrankpriorSDP'},...
        'legendNames')
    subplot(1,3,2)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalalmoptsLowRankcriterionl12innerProj0optsOptimizerlowrankpriorSDP','optsGlobalalmoptsLowRankcriterionl12optsOptimizerlowrankpriorSDP'},...
        'legendNames')
    subplot(1,3,3)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalSDP','optsGlobalalmoptsLowRankcriterionl1optsOptimizerlowrankpriorSDP','optsGlobalalmoptsLowRankcriterionl12optsOptimizerlowrankpriorSDP'},...
        'legendNames')
    disp('# Lesson (SDP, ALM)')
    disp('Also with SDP, adding the projection step in the middle of ALM improves performances')
end


%%
iFig=iFig+1;
if ismember(iFig,idxFig)
    ax=[0 axXLim 0 90];
    figure(iFig); 
    subplot(1,3,1)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalALM','optsGlobalALMoptsLowRankinnerProj0','optsGlobalSDP'},...
        'legendNames')
    axis(ax)
    subplot(1,3,2)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalALMoptsLowRankcriterionl1','optsGlobalALMoptsLowRankcriterionl1innerProj0',...
        'optsGlobalalmoptsLowRankcriterionl1optsOptimizerlowrankpriorSDP','optsGlobalalmoptsLowRankcriterionl1innerProj0optsOptimizerlowrankpriorSDP'},...
        'legendNames')
    axis(ax)
    subplot(1,3,3)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalALMoptsLowRankcriterionl12','optsGlobalALMoptsLowRankcriterionl12innerProj0',...
        'optsGlobalalmoptsLowRankcriterionl12optsOptimizerlowrankpriorSDP','optsGlobalalmoptsLowRankcriterionl12innerProj0optsOptimizerlowrankpriorSDP'},...
        'legendNames')
    axis(ax)
    savefigure(fullfile(figDir,'sdpVsLowrank'),figExt,figDimWide,flagSaveFig);

    disp('# Lesson (SDP+ALM, ALM)')
    disp('ALM for L2 with projection gives the same performances as trace SDP.')
    disp('Nuclear norm prior gives better results.')
end

%%
iFig=iFig+1;
if ismember(iFig,idxFig)
    ax=[0 axXLim 0 60];
    figure(iFig); 
%     sfm_poseCombine_testOutliers_plot(dataset,...
%         'methods',{'optsGlobalALM','optsGlobalALMoptsLowRankcriterionl1','optsGlobalALMoptsLowRankcriterionl12',...
%         'optsGlobalSDP','optsGlobalalmoptsLowRankcriterionl1optsOptimizerlowrankpriorSDP','optsGlobalalmoptsLowRankcriterionl12optsOptimizerlowrankpriorSDP'}, ...
%         'legendNames')
    sfm_poseCombine_testOutliers_plot(dataset,...
         'methods',{'optsGlobalALM','optsGlobalALMoptsLowRankcriterionl1','optsGlobalALMoptsLowRankcriterionl12'},...
        'legendNames')
    axis(ax)
    savefigure(fullfile(figDir,'lowrank'),figExt,figDimSingle,flagSaveFig);
    disp('# Lesson (ALM and SDP)')
    disp('SDP and ALM (with intermediate projection) formulations have virtually the same performances.')
    disp('Among the ALM methods, the L12 is better than the L1 which is better than the L2.')
end

%%
iFig=iFig+1;
if ismember(iFig,idxFig)
    figure(iFig);
    subplot(2,3,1)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalALM','localRefineoptsGlobalALM','localRefineoptsLocalmethodWeiszfeldoptsGlobalALM'},...
        'legendNames')
    subplot(2,3,2)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalALMoptsLowRankcriterionl1','localRefineoptsGlobalALMoptsLowRankcriterionl1','localRefineoptsLocalmethodWeiszfeldoptsGlobalALMoptsLowRankcriterionl1'},...
        'legendNames')
    subplot(2,3,3)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalALMoptsLowRankcriterionl12','localRefineoptsGlobalALMoptsLowRankcriterionl12','localRefineoptsLocalmethodWeiszfeldoptsGlobalALMoptsLowRankcriterionl12'},...
        'legendNames')
    subplot(2,3,4)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalSDP','localRefineoptsGlobalSDP','localRefineoptsLocalmethodWeiszfeldoptsGlobalSDP'},...
        'legendNames')
    subplot(2,3,5)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalSpectral','localRefineoptsGlobalSpectral','localRefineoptsLocalmethodWeiszfeldoptsGlobalSpectral'},...
        'legendNames')
    disp('# Lesson (local refinement)')
    disp('Local refinement significatively improves solution')
    disp('Reshaping function and Weiszfeld cost produce similar results (')
end

%%
iFig=iFig+1;
if ismember(iFig,idxFig)
    ax=[0 axXLim 0 70];
    figure(iFig);
    subplot(1,3,1)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalSpectral','localRefineoptsGlobalSpectral','localRefineoptsLocalmethodWeiszfeldoptsGlobalSpectral'},...
        'legendNames')
    axis(ax)
    subplot(1,3,2)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalSDP','localRefineoptsGlobalSDP','localRefineoptsLocalmethodWeiszfeldoptsGlobalSDP'},...
        'legendNames')
    axis(ax)
    subplot(1,3,3)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalALMoptsLowRankcriterionl12','localRefineoptsGlobalALMoptsLowRankcriterionl12','localRefineoptsLocalmethodWeiszfeldoptsGlobalALMoptsLowRankcriterionl12'},...
        'legendNames')
    axis(ax)
    savefigure(fullfile(figDir,'iterative'),figExt,figDimWide,flagSaveFig);
    disp('# Lesson (local refinement)')
    disp('Local refinement significatively improves solution')
    disp('Reshaping function and Weiszfeld cost produce similar results (')
end

%%
iFig=iFig+1;
if ismember(iFig,idxFig)
    figure(iFig); 
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'localRefineoptsGlobalALM','localRefineoptsGlobalALMoptsLowRankcriterionl1',...
        'localRefineoptsGlobalALMoptsLowRankcriterionl12','localRefineoptsGlobalSDP',...
        'localRefineoptsGlobalSpectral'},...
        'legendNames')
    disp('# Lesson (local refinement)')
    disp('Also with local refinement, the relation between different methods is the same. ALM L1 and ALM L12 are the best')
end

%%
iFig=iFig+1;
if ismember(iFig,idxFig)
    figure(iFig); 
    subplot(2,3,1)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalALM','inferenceOutliersoptsGlobalALM'},...
        'legendNames')
    subplot(2,3,2)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalALMoptsLowRankcriterionl1','inferenceOutliersoptsGlobalALMoptsLowRankcriterionl1'},...
        'legendNames')
    subplot(2,3,3)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalALMoptsLowRankcriterionl12','inferenceOutliersoptsGlobalALMoptsLowRankcriterionl12'},...
        'legendNames')
    subplot(2,3,4)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalSDP','inferenceOutliersoptsGlobalSDP'},...
        'legendNames')
    subplot(2,3,5)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalSpectral','inferenceOutliersoptsGlobalSpectral'},...
        'legendNames')
    disp('# Lesson (outlier inference)')
    disp('Outlier inference significatively improves solution')
end

%%
iFig=iFig+1;
if ismember(iFig,idxFig)
    figure(iFig); 
    ax=[0 axXLim 0 60];
    subplot(1,3,1)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalSpectral','inferenceOutliersoptsGlobalSpectral'},...
        'legendNames')
    axis(ax)
    subplot(1,3,2)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalSDP','inferenceOutliersoptsGlobalSDP'},...
        'legendNames')
    axis(ax)
    subplot(1,3,3)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'optsGlobalALMoptsLowRankcriterionl12','inferenceOutliersoptsGlobalALMoptsLowRankcriterionl12'},...
        'legendNames')
    axis(ax)
    savefigure(fullfile(figDir,'inference'),figExt,figDimWide,flagSaveFig);
    disp('# Lesson (outlier inference)')
    disp('Outlier inference significatively improves solution')
end

%%
iFig=iFig+1;
if ismember(iFig,idxFig)
    figure(iFig); 
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'inferenceOutliersoptsGlobalALM','inferenceOutliersoptsGlobalALMoptsLowRankcriterionl1',...
        'inferenceOutliersoptsGlobalALMoptsLowRankcriterionl12','inferenceOutliersoptsGlobalSDP',...
        'inferenceOutliersoptsGlobalSpectral'},...
        'legendNames')
    disp('# Lesson (outlier inference)')
    disp('With outlier inference, spectral performs relatively better. ALM L1 and ALM L12 are the best, but now ALM L1 is significatively better')

end
%%

iFig=iFig+1;
if ismember(iFig,idxFig)
    figure(iFig);
    subplot(2,3,1)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'inferenceOutliersoptsGlobalALM','localRefineoptsGlobalALM','inferenceOutlierslocalRefineoptsGlobalALM'},...
        'legendNames')
    subplot(2,3,2)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'inferenceOutliersoptsGlobalALMoptsLowRankcriterionl1','localRefineoptsGlobalALMoptsLowRankcriterionl1','inferenceOutlierslocalRefineoptsGlobalALMoptsLowRankcriterionl1'},...
        'legendNames')
    subplot(2,3,3)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'inferenceOutliersoptsGlobalALMoptsLowRankcriterionl12','localRefineoptsGlobalALMoptsLowRankcriterionl12','inferenceOutlierslocalRefineoptsGlobalALMoptsLowRankcriterionl12'},...
        'legendNames')
    subplot(2,3,4)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'inferenceOutliersoptsGlobalSDP','localRefineoptsGlobalSDP','inferenceOutlierslocalRefineoptsGlobalSDP'},...
        'legendNames')
    subplot(2,3,5)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'inferenceOutliersoptsGlobalSpectral','localRefineoptsGlobalSpectral','inferenceOutlierslocalRefineoptsGlobalSpectral'}...
    )
    disp('# Lesson (outlier inference + local refinement)')
    disp('Outlier inference together with local refinement strictly improves the solution over either one alone.')
    disp('The improvement is less marked with ALM L1 and ALM L12')
end
%%

iFig=iFig+1;
if ismember(iFig,idxFig)
    figure(iFig);
    ax=[0 axXLim 0 75];
    subplot(1,3,1)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'inferenceOutliersoptsGlobalSpectral','localRefineoptsGlobalSpectral','inferenceOutlierslocalRefineoptsGlobalSpectral'},...
        'legendNames')
    axis(ax)
    subplot(1,3,2)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'inferenceOutliersoptsGlobalSDP','localRefineoptsGlobalSDP','inferenceOutlierslocalRefineoptsGlobalSDP'},...
        'legendNames')
    axis(ax)
    subplot(1,3,3)
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'inferenceOutliersoptsGlobalALMoptsLowRankcriterionl12','localRefineoptsGlobalALMoptsLowRankcriterionl12','inferenceOutlierslocalRefineoptsGlobalALMoptsLowRankcriterionl12'},...
        'legendNames')
    axis(ax)
    savefigure(fullfile(figDir,'inferenceIterative'),figExt,figDimWide,flagSaveFig);
    disp('# Lesson (outlier inference + local refinement)')
    disp('Outlier inference together with local refinement strictly improves the solution over either one alone.')
    disp('The improvement is less marked with ALM L1 and ALM L12')
end

%%
iFig=iFig+1;
if ismember(iFig,idxFig)
    figure(iFig);
    ax=[0 axXLim 0 75];
    sfm_poseCombine_testOutliers_plot(dataset,...
        'methods',{'inferenceOutlierslocalRefineoptsGlobalALM','inferenceOutlierslocalRefineoptsGlobalALMoptsLowRankcriterionl1',...
        'inferenceOutlierslocalRefineoptsGlobalALMoptsLowRankcriterionl12','inferenceOutlierslocalRefineoptsGlobalSDP',...
        'inferenceOutlierslocalRefineoptsGlobalSpectral'},...
        'legendNames')
    axis(ax)
    savefigure(fullfile(figDir,'inferenceIterativeAll'),figExt,figDimSingle.*[1 1.2],flagSaveFig);
    
    disp('# Lesson (outlier inference + local refinement)')
    disp('The best solution is ALM L1, followed by ALM L12.')
    disp('Spectral and SDP perform surprisingly well.')
end
%%

setFigFontSize(fs);
setFigFont(fn);
