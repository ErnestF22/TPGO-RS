function eccv14_plots_cameraDepthCalibration
oldFS=setFigFontSize(7);
pathTex='~/Documents/UPenn/papers/network/ECCV14-AnisotropicPoseAveraging/';
pathFigures=fullfile(pathTex,'figures');
flagSaveFigures=2;
figSize=[200,110]*0.9;
figSizeNetwork=[200,130]*0.9;

close all
load 'cameraDepthSensorPlaneCalibration_testTrials_Trials50_It2000.mat'

methodAbsolutePoses='references';

switch 1
    case 1
        NPlanesVis=3;
        E=[1 2; 2 3; 3 1];
        E=[E;E+3;E+6];
        E=[E; [11*ones(9,1) (1:9)']];
        EComp=[10*ones(9,1) (1:9)'];
        d=allDatasets{4};
        T=G2T(cat(3,d.GPlaneTagsNoise(:,:,1:3*NPlanesVis),d.GCamera,d.GDepth));

        draw3dcameraFromG(d.GPlaneTagsNoise(:,:,1:3*NPlanesVis),methodAbsolutePoses,'shape','april','flagAxesLabel',false)
        hold on
        draw3dcameraFromG(d.GCamera,methodAbsolutePoses,'color','g','flagAxesLabel',false)
        draw3dcameraFromG(d.GDepth,methodAbsolutePoses,'color','b','flagAxesLabel',false)
        draw3dPlane(rigidTransformG(d.GDepth,d.NVecPlanesDepth(:,1:NPlanesVis),'planes',methodAbsolutePoses,'cw'),'side',3,'patchopts',{'EdgeColor','k','FaceColor','none'},'flagAxesLabel',false)
        plot3([T(1,E(:,1));T(1,E(:,2))],[T(2,E(:,1));T(2,E(:,2))],[T(3,E(:,1));T(3,E(:,2))],'r-')
        plot3([T(1,EComp(:,1));T(1,EComp(:,2))],[T(2,EComp(:,1));T(2,EComp(:,2))],[T(3,EComp(:,1));T(3,EComp(:,2))],'b-')
        hold off
        axis equal
        view(72,20)
        set(gcf,'Renderer','Painters')
        savefigure(fullfile(pathFigures,'incompletePoseDataset'),'epsc',figSizeNetwork,flagSaveFigures)
        
    case 2

        eAvg=[allErrorsAvg{:}];
        eLinear=[allErrorsLinear{:}];

        eRot=[eAvg.poseRot;eLinear.poseRot]'*180/pi;
        %eTransl=[eAvg.poseTransl;eLinear.poseTransl]'*100;

        %compute translation angular error
        r=[allResults{:}];
        TDepthNormAvg=cnormalize(G2T(cat(3,r.GDepthAvg)));
        eTranslAvg=acos(TDepthNormAvg(2,:))';
        TDepthNormLinear=cnormalize(G2T(cat(3,r.GDepthLinear)));
        eTranslLinear=acos(TDepthNormLinear(2,:))';
        eTransl=[eTranslAvg eTranslLinear]*180/pi;

        figure(1)
        cumDistBoxPerc(eRot)
        ylabel('Percentage of trials')
        xlabel('Angular Error [deg]')
        %legend('Averaged','Linear')
        disp('Angular Error [deg]')
        disp(mean(eRot))
        disp(std(eRot))
        savefigure(fullfile(pathFigures,'errorIncompleteDistributionRotations'),'epsc',figSize,flagSaveFigures)

        figure(2)
        cumDistBoxPerc(eTransl)
        ylabel('Percentage of trials')
        xlabel('Translation Error [cm]')
        legend('Averaged','Linear')
        disp('Translation angular Error [deg]')
        savefigure(fullfile(pathFigures,'errorIncompleteDistributionTranslations'),'epsc',figSize,flagSaveFigures)
        disp(mean(eTransl))
        disp(std(eTransl))
    case 3
        t_node=buildTestTagNetwork('seed',11);
        draw3dcameraFromG(t_node.gitruth(:,:,3:end),'shape','april','flagAxesLabel',false);
        hold on
        draw3dcameraFromG(t_node.gitruth(:,:,1:2),'flagAxesLabel',false);
        testNetworkDisplay(t_node,'edgesonly')
        axis equal
        view(-8,20)
        set(gcf,'Renderer','Painters')
        savefigure(fullfile(pathFigures,'completePoseDataset'),'epsc',figSizeNetwork,flagSaveFigures)
end
setFigFontSize(oldFS)
