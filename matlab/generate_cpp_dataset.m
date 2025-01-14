clc;
clear;
close all;

% 1a) PW TRANSLATION DATA INPUT: R, T are the gt, Tijs_nois are the input data
testdatas = [];

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


for tdata = testdatas
    folder_name = strcat("data/cpp_testdata/", ...
         "tdata_n", string(tdata.NNodes), "_mindeg", string(tdata.mindeg));
    [status, msg, msgID] = mkdir(folder_name);
    %edges
    % save('poc2degree_data/R_gt.mat', "R_globalframe")
    varEdges = tdata.E;
    writematrix(varEdges, ...
        convertStringsToChars(strcat(folder_name, "/edges.csv")), 'Delimiter', ',')
    %tijs
    varTijs = G2T(tdata.gijtruth);
    writematrix(varTijs, ...
        convertStringsToChars(strcat(folder_name, "/tijs.csv")), 'Delimiter', ',')
    %gt
    gt_vec = [vec(G2R(tdata.gitruth)); vec(G2T(tdata.gitruth))];
    writematrix(gt_vec, convertStringsToChars(strcat(folder_name, "/X_gt.csv")))
    %n
    n = tdata.NNodes;
    writematrix(n, convertStringsToChars(strcat(folder_name, "/n.csv")))
    %num edges
    e = tdata.NEdges;
    writematrix(e, convertStringsToChars(strcat(folder_name, "/e.csv")))

end

