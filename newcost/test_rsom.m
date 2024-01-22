function problem=test_rsom()

testdata = testNetwork_som(3); %4 is the default option

nrs = 3;
d = 3;
N = testdata.NNodes;
edges = (testdata.E);
num_edges = testdata.NEdges;

sz=[nrs,d,N];

Tijs = G2T(testdata.gijtruth);
R_gf = G2R(testdata.gitruth);
T_gf = G2T(testdata.gitruth);

[P, frct] = compute_step1_p_fct(R_gf, T_gf, Tijs, edges);


A = 0;
B = 0;
problem=struct("sz",sz, ...
    'P',P, 'fixed_cost_term', frct, ...
    'A',A,'B',B);



problem.cost=@(x) rsom_cost_rot_stiefel(x,problem);
problem.egrad=@(x) rsom_egrad_rot_stiefel(x,problem);
problem.rgrad=@(x) rsom_rgrad_rot_stiefel(x,problem);
problem.ehess=@(x,u) rsom_ehess_rot_stiefel(x,u,problem);
problem.rhess=@(x,u) rsom_rhess_rot_stiefel(x,u,problem);


end %file function


