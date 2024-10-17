function h=rsom_rhess_rot_stiefel(x,u,problem)
% From Boumal's Intro book, formula (5.34) page 111:
% Hess f(x)[u] =
% Proj_x \bigl(\mP_u (grad \overline{f}(x)) \bigr) + <- part1
% Proj_x(Hess \overline{f}(x)[u] \bigr) <- part2
% Note that part2 = 0 in our case

% ehess = rsom_ehess_rot_stiefel(x, u, problem);
% ehess_proj = stiefel_tangentProj(x, ehess);

G = rsom_egrad_rot_stiefel(x, problem);
term_1 = multiprod(u, ...
    0.5*multiprod(multitransp(x), G) + 0.5*multiprod(multitransp(G), x)); 
term_2 = multiprod(x, ...
    0.5*multiprod(multitransp(u), G) + 0.5*multiprod(multitransp(G), u));
DGf = - term_1 - term_2;
h = stiefel_tangentProj(x, DGf); %ehess_proj = zeros(nrs,d,N)

end %file function
