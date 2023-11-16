function L_next = from_L_to_L_next(L, problem_struct_next)
%FROM_L_TO_L_NEXT

nrs_next = problem_struct_next.sz(1);
d = problem_struct_next.sz(2);
N = problem_struct_next.sz(3);

L_next = zeros(nrs_next*N, nrs_next*N);

nrs = nrs_next - 1;

ids_L = reshape(1:nrs*N, nrs, []);
ids_to_sum = 0:size(ids_L,2)-1;
% ids_L_next = reshape(1:nrs_next, nrs_next, []);
ids_L_next = ids_L + ids_to_sum;

for ii = 1:N
    L_next(ids_L_next, ids_L_next) = L(ids_L, ids_L);
end


end

