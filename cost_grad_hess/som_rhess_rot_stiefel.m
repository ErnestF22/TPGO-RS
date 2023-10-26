function h=som_rhess_rot_stiefel(x,u,problem)
eh = ehess(x,u,problem);
X = matStack(x);
X_dot = matStack(u);
U = matStack(egrad(x,problem));
stief_proj_differential = X_dot * (X' * U + U' * X) + ...
    X * (X_dot' * U + U' * X_dot);
h = stp_boumal(x, matUnstack(stief_proj_differential, problem.sz(1))) + ...
    stp_boumal(x, eh);
end
