function Tijs = make_tijs_scaled(lambdas, tijs)
    num_edges = length(lambdas);

    Tijs = zeros(size(tijs));
    for ee = 1:num_edges
        Tijs(:,ee) = lambdas(ee) * tijs(:,ee);
    end
end