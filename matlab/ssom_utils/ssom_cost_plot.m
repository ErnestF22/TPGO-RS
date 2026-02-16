function cost_out = ssom_cost_plot(X, problem_data)

lambdas = X.lambda;
T = X.T;
R = X.R;


edges = problem_data.edges;
tijs = problem_data.tijs;
rho = problem_data.rho;

num_edges = size(edges, 1);

figure(5)

cost_out = 0.0;
for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    lambda_e = lambdas(ee);
    tij_e = tijs(:, ee);
    T_i = T(:, ii);
    T_j = T(:, jj);
    R_i = R(:, :, ii);
    a = T_i - T_j;
    b = R_i * tij_e;
    cost_ee = trace(a' * a + 2 * lambda_e * (a' * b) + lambda_e^2 * (b' * b)); 
    scale_compensation_ee = relu_som(ssom_relu_argument(lambda_e));
    cost_out = cost_out + cost_ee + rho * scale_compensation_ee * scale_compensation_ee;

    point1 = R_i * lambda_e * tij_e;
    point2 = T_i - T_j;

    clf;
    origin = [0,0,0];
    hold on;

    plot3([origin(1) point1(1)],[origin(2) point1(2)],[origin(3) point1(3)],'r-^', 'LineWidth',3);
    plot3([origin(1) point2(1)],[origin(2) point2(2)],[origin(3) point2(3)],'g-^', 'LineWidth',3);
    
    grid on;
    xlabel('X axis'), ylabel('Y axis'), zlabel('Z axis')
    set(gca,'CameraPosition',[1 2 3]);
    hold off;
end

end %file function