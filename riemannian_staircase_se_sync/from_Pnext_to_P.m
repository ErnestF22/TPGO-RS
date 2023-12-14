function P_prev = from_Pnext_to_P(P_next, problem_struct_prev)
%FROM_PNEXT_TO_P


nrs_prev = problem_struct_prev.sz(1);
d = problem_struct_prev.sz(2);
N = problem_struct_prev.sz(3);


P_prev = zeros(nrs_prev*N, d);

ids_P = reshape(1:nrs_prev*N, nrs_prev, []);
ids_to_sum = 0:size(ids_P,2)-1;
% ids_L_next = reshape(1:nrs_next, nrs_next, []);
ids_P_next = ids_P + ids_to_sum;

for ii = 1:N
    P_prev(ids_P(:,ii), :) = P_next(ids_P_next(:,ii), :);
end


end

