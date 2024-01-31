function h=rsom_rhess_rot_stiefel(x,u,problem)
% From Boumal's Intro book, formula (5.34) page 111:
% Hess f(x)[u] =
% Proj_x \bigl(\mP_u (grad \overline{f}(x)) \bigr) + <- part1
% Proj_x(Hess \overline{f}(x)[u] \bigr) <- part2
% Note that part2 = 0 in our case

ehess = rsom_ehess_rot_stiefel(x, u, problem);
ehess_proj = stiefel_tangentProj(x, ehess);


% G = matStack(rsom_egrad_rot_stiefel(x,problem));
% X = matStack(x);
% U = matStack(u);
G = rsom_egrad_rot_stiefel(x, problem);
% DGf = - U * 0.5 * (X' * G + G' * X) - X * 0.5 * (U' * G + G' * U);
term_1 = multiprod(u, ...
    0.5*multiprod(multitransp(x), G) + 0.5*multiprod(multitransp(G), x)); 
term_2 = multiprod(x, ...
    0.5*multiprod(multitransp(u), G) + 0.5*multiprod(multitransp(G), u));
DGf = - term_1 - term_2;
% h = stiefel_tangentProj(x, matUnstack(DGf, problem.sz(1)))+ ehess_proj;
h = stiefel_tangentProj(x, DGf) + ehess_proj;
end
