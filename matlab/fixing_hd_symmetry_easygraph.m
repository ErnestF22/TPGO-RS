function fixing_hd_symmetry_easygraph

N = 5;
mindeg = 2;

% 1a) PW TRANSLATION DATA INPUT: R, T are the gt, Tijs_nois are the input data
% testdata = testNetwork_params(3, N, 'banded', mindeg); %4 would be the default
testdata = testNetwork_params_translpert(3, N, 'banded', mindeg); 

T = [0, 0, 1; -5,-2, 1; -1, -1, -1; 5, -2, 1; 3, 3, 4]';

testdata = testNetworkInitializeStates(testdata, 'T', T);

figure(702)
% transf_ssom = RT2G(R_out_noiseless, T_out_noiseless);
% testdata = testNetwork_params(3, N, 'banded', 2); %4 would be the default
% testdata.E = edges;
% testdata.NNodes = N;
% testdata.gi = transf_out2;
% % testdata.gijtruth = 
% testdata.lambdaij = lambdas_ssom_out2;
% testdata.gitruth = RT2G_stiefel(X_gt.R, X_gt.T);
% testdata.lambdaijtruth = X_gt.lambda;
% testdata = testNetworkCompensate(testdata);
hold on;
red=[65535	8567	0]/65535;
opts_draw_camera={'Color1',red,'Color2',red};
testNetworkDisplay(testdata, 'member','gi','optionsDrawCamera', opts_draw_camera, 'flagDisplayEdges', false, 'flagDisplayCameraNumber', true)
% green=[15934	35723	14392]/65535/0.6;           %camera color
% opts_draw_camera={'Color1',green,'Color2',green};  %options to pass to drawCamera
% testNetworkDisplay(testdata, 'member','gitruth', 'optionsDrawCamera', opts_draw_camera, 'flagDisplayEdges', false, 'flagDisplayCameraNumber', true)

% edgesStyle='g';
edges = testdata.E;
num_edges = size(edges, 1);
orange_vecs = zeros(3, num_edges);
purple_vecs = zeros(3, num_edges);

R = G2R(testdata.gi);
T = G2T(testdata.gi);
tijs = G2T(testdata.gij);
lambdaijs = testdata.lambdaij;
for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    R_i = R(:,:,ii);
    T_i = T(:,ii);
    T_j = T(:,jj);
    tij = tijs(:,ee);
    orange_vecs(:,ii) = R_i * tij;
    lambdaij = lambdaijs(ee);
    purple_vecs(:,ii) = (T_i - T_j) / lambdaij;
end

edgesStyle='r';
plotArrows(T(:,edges(:,1)), orange_vecs, edgesStyle) % T(:,edges(:,2))
edgesStyle='g';
plotArrows(T(:,edges(:,1)), purple_vecs, edgesStyle) % T(:,edges(:,2))
 
hold off;


end %file function