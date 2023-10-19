function c = compute_fixed_cost_term(tijs_vec, d)
%COMPUTE_FIXED_COST_TERM function that takes in input the costs on the
%edges that are given as data of the problem, and for each one of them
%computes the outer product, then summing them all and finally returning
%the trace of the matrix resulting from all the sums

c_mat = zeros(d,d); % \trace \sum_{e_{ij}} T_{ij}T_{ij}'
for ii = 1:size(tijs_vec, 2)
    c_mat_ij = tijs_vec(:, ii) * tijs_vec(:, ii)';
    c_mat = c_mat + c_mat_ij;
end
% disp("trace(c):");
% disp(trace(c));
c = trace(c_mat);

end %function