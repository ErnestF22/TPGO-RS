function testnet = testNetwork_params(testNum,num_nodes,mode,min_node_deg,sigmaR,sigmaT)
    randn('state',0)
    rand('state',0)
    close all
    
    if ~exist('testNum','var')
        testNum=4;
    end

    if ~exist('sigmaR','var')
        sigmaR=0.00;
    end

    if ~exist('sigmaT','var')
        sigmaT=0.00;
    end
    
    switch testNum
        case 1
            methodAbsolutePoses='references';
            structType='array';
        case 2
            methodAbsolutePoses='poses';
            structType='array';
        case 3
            methodAbsolutePoses='references';
            structType='single';
        case 4
            methodAbsolutePoses='poses';
            structType='single';
    end
    
    N=num_nodes;
%     A=adjgallery(N,'banded',3);
%     A=adjgallery(N,'banded',2); % 14 edges usual
%     A=adjgallery(N,'full'); %uncomment this for complete graph
    if (strcmp(mode, 'full'))
        A = adjgallery(N, mode);
    elseif (strcmp(mode, 'banded'))
        A = adjgallery(N, mode, min_node_deg);
    else
        error('Unknown testNetwork mode');
    end
    [G,X]=testNetworkCreateAbsolutePoses(N);
    
    t_node=testNetworkCreateStruct(A,'Type',structType);
    t_node=testNetworkAddGroundTruth(t_node,G,methodAbsolutePoses,'flagInvertG',true);
    
    t_node=testNetworkProjectImages(t_node,X,methodAbsolutePoses);
    
    %figure
    %set(gcf,'Name','Images in each camera')
    %testNetworkDisplayImages(t_node)
    
    t_node=testNetworkAddMeasurements(t_node,'Method','Essential','SigmaR',sigmaR,'SigmaT',sigmaT); %this adds gij
    
    t_node=testNetworkInitializeStates(t_node,'MethodR','NoisyTruth',...
        'MethodT','NoisyTruth','MethodScale','NoisyTruth',...
        'SigmaR',sigmaR,'SigmaT',sigmaT); %this adds gi
    
    %Note:
    % - measurements errors should be all zero
    % - errors should be in the order of 2 times the variances given during
    %   state initialization
    % figure(3)
    % subplot(2,1,1)
    % hist([rotErr rotErrMs])
    % subplot(2,1,2)
    % hist([translErr translErrMs])
    % disp(['Geometric std of ratios (est/mes)', num2str(geostd(scaleRatios)), '/', num2str(geostd(scaleRatiosMs))])
    %figure
    %testNetworkDisplayErrors(t_node,'rtsn','optsComputeErrors',{methodAbsolutePoses},'boxplots')
    
    %see inside testNetworkDisplayErrors for examples of use of
    %   testNetworkComputeErrors
    %   testNetworkDisplay

    % edge_1_i = t_node.E(1,1);
    % edge_1_j = t_node.E(1,2);
    % R_tmp_truth = (t_node.gitruth);
    % R_tmp_truth_i = R_tmp_truth(:,:,edge_1_i);
    % R_tmp_truth_j = R_tmp_truth(:,:,edge_1_j);
    % R_tmp_truth_ij = (t_node.gijtruth);
    % disp("[inv(R_tmp_truth_i) * R_tmp_truth_j, R_tmp_truth_ij(:,:,1)]")
    % disp([inv(R_tmp_truth_i) * R_tmp_truth_j, R_tmp_truth_ij(:,:,1)])
    
    % num_edges = size(t_node.E, 1);
    % for ee = 1:num_edges
    %     e_i = t_node.E(ee,1);
    %     e_j = t_node.E(ee,2);
    %     t_node.gij(:,:,ee) = inv(t_node.gi(:,:,e_i)) * t_node.gi(:,:,e_j);
    % end

    t_node.lambdaij = t_node.lambdaijtruth;

    testnet = t_node;



end %file function