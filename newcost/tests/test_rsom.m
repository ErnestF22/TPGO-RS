function problem=test_rsom()

testdata = testNetwork_som(3); %4 is the default option

nrs = 4;
d = 3;
N = testdata.NNodes;
edges = (testdata.E);
% num_edges = testdata.NEdges;

sz=[nrs,d,N];

Tijs = G2T(testdata.gijtruth);
R_gf = G2R(testdata.gitruth);
T_gf = G2T(testdata.gitruth);

T_gf_stief = cat_zero_row(T_gf);

[P, frct] = make_step1_p_fct(T_gf_stief, Tijs, edges);
[LR, PR, BR] = make_LR_PR_BR_noloops(R_gf, Tijs, edges);


problem=struct("sz",sz, ...
    'P', P, 'frct', frct, ...
    'PR', PR, 'LR', LR, 'BR', BR);



problem.cost=@(x) rsom_cost_rot_stiefel(x,problem);
problem.egrad=@(x) rsom_egrad_rot_stiefel(x,problem);
problem.grad=@(x) rsom_rgrad_rot_stiefel(x,problem);
problem.ehess=@(x,u) rsom_ehess_rot_stiefel(x,u,problem);
problem.hess=@(x,u) rsom_rhess_rot_stiefel(x,u,problem);


end %file function


