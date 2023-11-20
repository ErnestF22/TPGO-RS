function L_prev = from_Lnext_to_L(L, problem_struct_prev)
%FROM_LNEXT_TO_L

nrs_prev = problem_struct_prev.sz(1);
d = problem_struct_prev.sz(2);
N = problem_struct_prev.sz(3);

L_prev = zeros(nrs_prev*N, nrs_prev*N);

ids_L = reshape(1:nrs_prev*N, nrs_prev, []);
ids_to_sum = 0:size(ids_L,2)-1;
% ids_L_next = reshape(1:nrs_next, nrs_next, []);
ids_L_next = ids_L + ids_to_sum;

for ii = 1:N
    L_prev(ids_L(:,ii), ids_L(:,ii)) = L(ids_L_next(:,ii), ids_L_next(:,ii));
end


end

