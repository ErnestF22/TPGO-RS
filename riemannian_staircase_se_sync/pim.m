function v_max = pim(A, thresh)
%PIM (acronym for Power Iteration Method) Iterative method that returns 
%an eigenvector associated to che maximum eigenvalue of matrix A

if ~exist('thresh','var')
  thresh=1e-5;
end

iterative_change = 1e+6;
d = size(A, 1);
v = rand(d,1);
iteration_num = 0;
while (iterative_change > thresh) && (iteration_num < 1000)
    iteration_num = iteration_num + 1;
    v_prev = v;
    v = v ./ norm(v);
    v = - A * v;
    iterative_change = min(norm(v_prev - v), norm(v_prev + v));
end

v_max = v;

end

