clc;
clear;
close all;

% 1a) PW TRANSLATION DATA INPUT: R, T are the gt, Tijs_nois are the input data
testdatas = [];

% 0b) Noise PARAMS
%NOTE: sigmas, mus can be seen as couples for each test
sigmas = readmatrix("data/sigmas.txt"); %sigma = stdev, sigma.^2 = variance
mus = readmatrix("data/mus.txt"); %OBS. generally, mus can be d-dimensional; here, we just assume them as scalar (i.e. a d-dimensional vector with all coordinates equal)

% sigmas = [0.0; 0.1];

for ss = 1:size(sigmas,1)
    s = sigmas(ss);

    N = 5;
    mindeg = 2;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, 0.0, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 0; % compensation part
    testdata.a = 2.0; % log scale-compensation param    
    testdatas = [testdatas, testdata];
end

run_repeated_test(testdatas, sigmas, mus);
testdatas = [];

for ss = 1:size(sigmas,1)
    s = sigmas(ss);

    N = 5;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, 0.0, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 0;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];
end

run_repeated_test(testdatas, sigmas, mus);
testdatas = [];

for ss = 1:size(sigmas,1)
    s = sigmas(ss);

    N = 6;
    mindeg = 2;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, 0.0, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 0;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];
end

run_repeated_test(testdatas, sigmas, mus);
testdatas = [];

for ss = 1:size(sigmas,1)
    s = sigmas(ss);

    N = 6;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, 0.0, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 0;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];
end

run_repeated_test(testdatas, sigmas, mus);
testdatas = [];

for ss = 1:size(sigmas,1)
    s = sigmas(ss);

    N = 7;
    mindeg = 2;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, 0.0, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 0;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];
end

run_repeated_test(testdatas, sigmas, mus);
testdatas = [];

for ss = 1:size(sigmas,1)
    s = sigmas(ss);

    N = 7;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, 0.0, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 0;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];
end

run_repeated_test(testdatas, sigmas, mus);
testdatas = [];

for ss = 1:size(sigmas,1)
    s = sigmas(ss);

    N = 8;
    mindeg = 2;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, 0.0, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 0;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];
end

run_repeated_test(testdatas, sigmas, mus);
testdatas = [];

for ss = 1:size(sigmas,1)
    s = sigmas(ss);

    N = 8;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, 0.0, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 0;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];
end

run_repeated_test(testdatas, sigmas, mus);
testdatas = [];

for ss = 1:size(sigmas,1)
    s = sigmas(ss);

    N = 9;
    mindeg = 2;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, 0.0, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 0;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];
end

run_repeated_test(testdatas, sigmas, mus);
testdatas = [];

for ss = 1:size(sigmas,1)
    s = sigmas(ss);

    N = 9;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, 0.0, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 0;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];
end

run_repeated_test(testdatas, sigmas, mus);
testdatas = [];

for ss = 1:size(sigmas,1)
    s = sigmas(ss);

    N = 10;
    mindeg = 2;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, 0.0, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 0;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];
end

run_repeated_test(testdatas, sigmas, mus);
testdatas = [];

for ss = 1:size(sigmas,1)
    s = sigmas(ss);

    N = 10;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, 0.0, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 0;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];
end

run_repeated_test(testdatas, sigmas, mus);
testdatas = [];

for ss = 1:size(sigmas,1)      
    s = sigmas(ss);

    N = 25;
    mindeg = 3;
    testdata = [];
    testdata = testNetwork_params(3, N, 'banded', mindeg, 0.0, s);
    testdata.mindeg = mindeg;
    testdata.sigma = s;
    testdata.rho = 0;
    testdata.a = 2.0;
    testdatas = [testdatas, testdata];
end

run_repeated_test(testdatas, sigmas, mus);
testdatas = [];