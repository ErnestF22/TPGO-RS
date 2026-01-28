function quick_save_tdata(sigma,tdata,X_initguess,Tijs, Tijs_nois)

d = 3;
% sigma
folder_name = strcat("data/check_same_output/", ...
     "tdata_n", string(tdata.NNodes), ...     % "_mindeg", string(tdata.mindeg),
     "_sigma", sprintf( '%02d', sigma*10 ));
[status, msg, msgID] = mkdir(folder_name);
%edges
% save('poc2degree_data/R_gt.mat', "R_globalframe")
varEdges = tdata.E;
writematrix(varEdges, ...
    convertStringsToChars(strcat(folder_name, "/edges.csv")), 'Delimiter', ',')
%tijs
% Tijs = G2T(tdata.gijtruth);
% Tijs_nois = Tijs + sigma*randn(size(Tijs));
if ~exist("Tijs_nois", "var")
    Tijs_nois = Tijs;
end
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
% transf_initguess = RT2G(R_initguess, transl_initguess);
transf_initguess_vec = [X_initguess.R(:); X_initguess.T(:); X_initguess.lambda(:)];
writematrix(transf_initguess_vec, convertStringsToChars(strcat(folder_name, "/ssom_x_start.csv")))