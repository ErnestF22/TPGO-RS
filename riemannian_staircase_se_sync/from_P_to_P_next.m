function P_next = from_P_to_P_next(P, problem_struct_next)
%FROM_P_TO_P_NEXT 


nrs_next = problem_struct_next.sz(1);
d = problem_struct_next.sz(2);
N = problem_struct_next.sz(3);

P_next = zeros(nrs_next*N, d);

nrs = nrs_next - 1;

ids_P = reshape(1:nrs*N, nrs, []);
ids_to_sum = 0:size(ids_P,2)-1;
% ids_L_next = reshape(1:nrs_next, nrs_next, []);
ids_P_next = ids_P + ids_to_sum;

for ii = 1:N
    P_next(ids_P_next(:,ii), :) = P(ids_P(:,ii), :);
end


end

