function cvpr13_plots_aVsR
%dest='paper';
dest='poster';
flagSaveFile=2;
flagPlots=true;

switch dest
    case 'paper'
        figDimR=[7 6]*savefigureResolution('cm');
        figDimT=[7 3.7]*savefigureResolution('cm');
        figDir='~/Documents/UPenn/papers/vision/CVPR13-EssentialManifold/figures/';
        flagMatchesImg=false;
        fs=setFigFontSize(7);
        fn=setFigFont('Times');
    case 'poster'
        figDimR=2*[8 5.5]*savefigureResolution('cm');
        figDimT=figDimR;
        %figDimMatches=[7 6]*savefigureResolution('cm');
        figDir='~/Documents/UPenn/presentations/vision/CVPR14-EssentialManifold/figures/';
        flagMatchesImg=true;
        fs=setFigFontSize(18);
        fn=setFigFont('Helvetica');
end

close all

load('cvpr13_fountain_data_aVsR')
idxImgAll=[data.match.idxImg];
idxImgValid=[data.matchFiltered.idxImg];
[~,idxIdxIntersect]=intersect(idxImgAll',idxImgValid','rows');

if flagPlots
    load('averagingVsRansac_trials_29-Oct-2013_12_59_04')
    disp(['Number of trials: ' num2str(iTrial)])
    errors=errors(:,1:iTrial);
    errRCell=cellfun(@(x) x.errR, errors(idxIdxIntersect),'uniformoutput',false);
    errRMean=mean(cat(3,errRCell{:}),3);
    errRMedian=median(cat(3,errRCell{:}),3);
    legendTextR=errors{1,1}.legendTextR;
    idxSamples=cellfun(@(x) strcmpi('Samples',x), legendTextR);

    errRMean(:,idxSamples)=[];%mean(errR(:,idxSamples));
    errRMedian(:,idxSamples)=[];%median(errR(:,idxSamples));
    legendTextR(idxSamples)=[];

    figure(1)
    plot(errRMean*180/pi)
    hold on
    plotInterp(errRMedian*180/pi)
    hold off
    legend(legendTextR)
    xlabel('# of samples')
    ylabel('Angle error [deg]')
    savefigure(fullfile(figDir,'averagingVsRansac_errR_mean'),'epsc',figDimR,flagSaveFile)

    colorOrder=get(gca,'ColorOrder');

    errTCell=cellfun(@(x) x.errT, errors,'uniformoutput',false);
    errTMean=mean(cat(3,errTCell{:}),3);
    errTMedian=median(cat(3,errTCell{:}),3);
    legendTextT=errors{1,1}.legendTextT;
    idxSamples=cellfun(@(x) strcmpi('Samples',x), legendTextT);
    errTMean(:,idxSamples)=[];%mean(errTMean(:,idxSamples));
    errTMedian(:,idxSamples)=[];%mean(errTMean(:,idxSamples));
    legendTextT(idxSamples)=[];

    figure(2)
    [~,idxRT]=intersect(char(legendTextR),char(legendTextT),'rows');
    plot(errTMean*180/pi)
    set(gca,'ColorOrder',colorOrder(sort(idxRT),:));
    set(gca,'NextPlot','replacechildren');
    plot(errTMean*180/pi)
    hold on
    plotInterp(errTMedian*180/pi)
    hold off
    legend(legendTextT)
    xlabel('# of samples')
    ylabel('Angle error [deg]')
    savefigure(fullfile(figDir,'averagingVsRansac_errT_mean'),'epsc',figDimT,flagSaveFile)
end

if flagMatchesImg
    matchIdxImg=sfm_getMatchIdxImg(data);
    iiMatch=3;
    idxImg=matchIdxImg(:,iiMatch);
%     copyfile(data.imgFileName{idxImg(1)},fullfile(figDir,'img1.png'))
%     copyfile(data.imgFileName{idxImg(2)},fullfile(figDir,'img2.png'))
    img1=sfm_getImageById(data,idxImg(1));
    img2=sfm_getImageById(data,idxImg(2));
    img1=imresize(img1,0.66);
    img2=imresize(img2,0.66);
    imwrite(img1,fullfile(figDir,'img1.png'),'png')
    imwrite(img2,fullfile(figDir,'img2.png'),'png')
    [x1,x2]=sfm_getFeatureLocationFromMatchId(data,iiMatch);
%     x1=data.imgSize(:,idxImg(1));
%     x1=[x1 [x1(1);0] [0;x1(2)]];
    x1(2,:)=-x1(2,:);%data.imgSize(2,idxImg(1))-x1(2,:);
    x2(2,:)=-x2(2,:);
    xscale=0.01;
    x1=x1*xscale;
    x2=x2*xscale;
    fid=fopen(fullfile(figDir,'img1data.txt'),'w');
    fprintf(fid,'%.10f %.10f\n',x1);
    fclose(fid);
    fid=fopen(fullfile(figDir,'img2data.txt'),'w');
    fprintf(fid,'%.10f %.10f\n',x2);
    fclose(fid);
    fid=fopen(fullfile(figDir,'img1fivePoints.tex'),'w');
    for iPoint=1:5
        fprintf(fid,'\\coordinate (img1-%d) at (%.10f,%.10f);\n',iPoint,x1(:,iPoint));
    end
    fclose(fid);
    fid=fopen(fullfile(figDir,'img2fivePoints.tex'),'w');
    for iPoint=1:5
        fprintf(fid,'\\coordinate (img2-%d) at (%.10f,%.10f);\n',iPoint,x2(:,iPoint));
    end
    fclose(fid);
%     figure(3)
%     sfm_rawDisplayFeature(img1,x1)
%     savefigure(fullfile(figDir,'imgFeatures1'),'epsc',figDimMatches,flagSaveFile)
%     figure(3)
%     sfm_rawDisplayFeature(img2,x2)
%     savefigure(fullfile(figDir,'imgFeatures2'),'epsc',figDimMatches,flagSaveFile)
end

setFigFontSize(fs)
setFigFont(fn)

function plotInterp(y)
[xInterp,yInterp]=interpHalf(y);
plot(xInterp,yInterp,'.','MarkerSize',3)

function [xInterp,yInterp]=interpHalf(y)
sz=size(y);
xInterp=1:0.5:size(y,1);
yInterp=cat(3,y,([y(2:end,:);zeros(1,sz(2))]+y)/2);
yInterp=reshape(permute(yInterp,[3 1 2]),sz(1)*2,sz(2));
yInterp(end,:)=[];
