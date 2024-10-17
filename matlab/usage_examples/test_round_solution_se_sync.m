nrs = 4;
d = 3;
N = 5;


problem_struct_round_solution.d = d;
problem_struct_round_solution.n = N; %we do them one at a time

Y_opt = repmat(cat_zero_row(rotz(30)), 1, N);

% Y_opt_rounded = zeros(d,d,N);

% for ii = 1:N
%     Y_opt_rounded_ii = round_solution_se_sync(Y_opt(:,:,ii), problem_struct_round_solution);
% 
%     Y_opt_rounded(:,:,ii) = Y_opt_rounded_ii;
% 
% end

Y_opt_rounded = round_solution_se_sync(Y_opt, problem_struct_round_solution);

disp("Y_opt")
disp(Y_opt)
disp("Y_opt_rounded")
disp(Y_opt_rounded)

