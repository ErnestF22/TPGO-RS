clc;
clear;
close all;

% 1a) PW TRANSLATION DATA INPUT: R, T are the gt, Tijs_nois are the input data
testdatas = [];

d = 3;
mu = 0.0;

sigmas = readmatrix("data/sigmas.txt");

for ii = 1:size(sigmas,1)

    s = sigmas(ii);

    N = 5;
    mindeg = 2;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdatas = [testdatas, testdata];
    
    N = 5;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdatas = [testdatas, testdata];
    
    N = 6;
    mindeg = 2;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdatas = [testdatas, testdata];
    
    N = 6;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdatas = [testdatas, testdata];
    
    N = 7;
    mindeg = 2;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdatas = [testdatas, testdata];
    
    N = 7;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdatas = [testdatas, testdata];
    
    N = 8;
    mindeg = 2;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdatas = [testdatas, testdata];
    
    N = 8;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdatas = [testdatas, testdata];
    
    N = 9;
    mindeg = 2;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdatas = [testdatas, testdata];
    
    N = 9;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdatas = [testdatas, testdata];
    
    % N = 10;
    % mindeg = 2;
    % testdata = [];
    % testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    % testdata.mindeg = mindeg;
    % testdata.sigma = s;
    % testdatas = [testdatas, testdata];
    
    N = 10;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdatas = [testdatas, testdata];
    
    N = 25;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, s, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdatas = [testdatas, testdata];
end



for t = 1:size(testdatas,2)
    tdata = testdatas(t);
    num_edges = size(tdata.E, 1);
    sigma = tdata.sigma;
    folder_name = strcat("data/ssom_testdata_noisy/harder/", ...
         "tdata_n", string(tdata.NNodes), ...
         "_mindeg", string(tdata.mindeg), ...
         "_sigma", sprintf( '%03d', sigma*100 ));
    [status, msg, msgID] = mkdir(folder_name);
    %edges
    % save('poc2degree_data/R_gt.mat', "R_globalframe")
    varEdges = tdata.E;
    writematrix(varEdges, ...
        convertStringsToChars(strcat(folder_name, "/edges.csv")), 'Delimiter', ',')
    %tijs
    Tijs = G2T(tdata.gijtruth);
    for ee = 1:num_edges
        Tijs(:,ee) = Tijs(:,ee) / tdata.lambdaijtruth(1,ee);
    end

    Tijs_nois = G2T(tdata.gij);
    Tijs_nois = Tijs_nois + sigma*(rand(size(Tijs_nois)));
    for ee = 1:num_edges
        Tijs_nois(:,ee) = Tijs_nois(:,ee) / tdata.lambdaij(1,ee);
    end
    
    writematrix(Tijs_nois, ...
        convertStringsToChars(strcat(folder_name, "/tijs.csv")), 'Delimiter', ',')
    writematrix(Tijs, ...
        convertStringsToChars(strcat(folder_name, "/tijs_truth.csv")), 'Delimiter', ',')
    %gt
    gt_truth_vec = [vec(G2R(tdata.gitruth)); vec(G2T(tdata.gitruth)); vec(tdata.lambdaijtruth)];
    writematrix(gt_truth_vec, convertStringsToChars(strcat(folder_name, "/Xgt_truth.csv")))
    gt_vec = [vec(G2R(tdata.gi)); vec(G2T(tdata.gi)); vec(tdata.lambdaij)];
    writematrix(gt_vec, convertStringsToChars(strcat(folder_name, "/Xgt.csv")))
    %n
    n = tdata.NNodes;
    writematrix(n, convertStringsToChars(strcat(folder_name, "/n.csv")))
    %num edges
    e = tdata.NEdges;
    writematrix(e, convertStringsToChars(strcat(folder_name, "/e.csv")))
    %
    R_initguess = randrot_som(d, tdata.NNodes);
    transl_initguess = 10 * rand(d, tdata.NNodes);
    lambdas_initguess = ones(e, 1);
    % transf_initguess = RT2G(R_initguess, transl_initguess);
    transf_initguess_vec = [R_initguess(:); transl_initguess(:); lambdas_initguess(:)];
    writematrix(transf_initguess_vec, convertStringsToChars(strcat(folder_name, "/ssom_x_start.csv")))
end

