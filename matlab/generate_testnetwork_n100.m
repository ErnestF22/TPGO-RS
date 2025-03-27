function generate_testnetwork_n100

N = 50;
mindeg = 2;
d = 3;

tdata = testNetwork_params(3, N, 'banded', mindeg); %3 is the mode
tdata.mindeg = mindeg;

sigmas = readmatrix("data/sigmas.txt");
mu = 0.0;

for s = 1:length(sigmas)
    sigma = sigmas(s);
    folder_name = strcat("data/big_testnetworks/", ...
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
    gt_vec = [vec(G2R(tdata.gitruth)); vec(G2T(tdata.gitruth))];
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
    % transf_initguess = RT2G(R_initguess, transl_initguess);
    transf_initguess_vec = [R_initguess(:); transl_initguess(:)];
    writematrix(transf_initguess_vec, convertStringsToChars(strcat(folder_name, "/startX.csv")))
end


end % file function
