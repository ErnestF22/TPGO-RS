function [A, B] = make_a_b_T_lambda_qp(R, problem_data, ii, jj, ee)

num_edges = size(problem_data.edges, 1);
N = size(R, 3);
A = zeros(3, 3*N);
B = zeros(3, num_edges);

tijs = problem_data.tijs;

A(:, (ii-1)*3+1:ii*3) = eye(3);
A(:, (jj-1)*3+1:jj*3) = -eye(3);

B(:, ee) = R(:,:,ii) * tijs(:,ee);

end