clc;
clear;
close all;

% 1a) PW TRANSLATION DATA INPUT: R, T are the gt, Tijs_nois are the input data
testdatas = [];

d = 3;
mu = 0.0;

N = 5;
mindeg = 2;
testdata = testNetwork_params(3, N, 'banded', mindeg);
testdata.mindeg = mindeg;
testdatas = [testdatas, testdata];

N = 5;
mindeg = 3;
testdata = testNetwork_params(3, N, 'banded', mindeg);
testdata.mindeg = mindeg;
testdatas = [testdatas, testdata];

N = 6;
mindeg = 2;
testdata = testNetwork_params(3, N, 'banded', mindeg);
testdata.mindeg = mindeg;
testdatas = [testdatas, testdata];

N = 6;
mindeg = 3;
testdata = testNetwork_params(3, N, 'banded', mindeg);
testdata.mindeg = mindeg;
testdatas = [testdatas, testdata];

N = 7;
mindeg = 2;
testdata = testNetwork_params(3, N, 'banded', mindeg);
testdata.mindeg = mindeg;
testdatas = [testdatas, testdata];

N = 7;
mindeg = 3;
testdata = testNetwork_params(3, N, 'banded', mindeg);
testdata.mindeg = mindeg;
testdatas = [testdatas, testdata];

N = 8;
mindeg = 2;
testdata = testNetwork_params(3, N, 'banded', mindeg);
testdata.mindeg = mindeg;
testdatas = [testdatas, testdata];

N = 8;
mindeg = 3;
testdata = testNetwork_params(3, N, 'banded', mindeg);
testdata.mindeg = mindeg;
testdatas = [testdatas, testdata];

N = 9;
mindeg = 2;
testdata = testNetwork_params(3, N, 'banded', mindeg);
testdata.mindeg = mindeg;
testdatas = [testdatas, testdata];

N = 9;
mindeg = 3;
testdata = testNetwork_params(3, N, 'banded', mindeg);
testdata.mindeg = mindeg;
testdatas = [testdatas, testdata];

% N = 10;
% mindeg = 2;
% testdata = testNetwork_params(3, N, 'banded', mindeg);
% testdata.mindeg = mindeg;
% testdatas = [testdatas, testdata];

N = 10;
mindeg = 3;
testdata = testNetwork_params(3, N, 'banded', mindeg);
testdata.mindeg = mindeg;
testdatas = [testdatas, testdata];

N = 25;
mindeg = 3;
testdata = testNetwork_params(3, N, 'banded', mindeg);
testdata.mindeg = mindeg;
testdatas = [testdatas, testdata];

sigmas = readmatrix("data/sigmas.txt");
for tdata = testdatas
    for s = 1:length(sigmas)
        sigma = sigmas(s);
        folder_name = strcat("data/ssom_testdata_noisy/harder/", ...
             "tdata_n", string(tdata.NNodes), ...
             "_mindeg", string(tdata.mindeg), ...
             "_sigma", sprintf( '%02d', sigma*10 ));
        [status, msg, msgID] = mkdir(folder_name);
        %edges
        % save('poc2degree_data/R_gt.mat', "R_globalframe")
        varEdges = tdata.E;
        writematrix(varEdges, ...
            convertStringsToChars(strcat(folder_name, "/edges.csv")), 'Delimiter', ',')
        %tijs
        Tijs = G2T(tdata.gijtruth);
        Tijs_nois = Tijs + sigma*randn(size(Tijs)) + ...
                mu * ones(size(Tijs));
        writematrix(Tijs_nois, ...
            convertStringsToChars(strcat(folder_name, "/tijs.csv")), 'Delimiter', ',')
        writematrix(Tijs, ...
            convertStringsToChars(strcat(folder_name, "/tijs_truth.csv")), 'Delimiter', ',')
        %gt
        gt_vec = [vec(G2R(tdata.gitruth)); vec(G2T(tdata.gitruth)); vec(tdata.lambdaijtruth)];
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
        writematrix(transf_initguess_vec, convertStringsToChars(strcat(folder_name, "/startX.csv")))
    end
end

