function [A,At,At_full] = reducedIncidenceMatrix(edges_id, n, m)
% For a directed graph with n+1 nodes and m edges, the function computes 
% the n x m reduced incidence matrix A, i.e., the incidence matrix of the 
% graph without the first row. It also returns At = A' and the complete
% incidence matrix At_full (n+1 x m). 
% - "edges_id" is an m x 2 matrix, such that each row [i j] encodes an edge
% starting from i and ending at j (assume the indices start from 1)
% - n = number of nodes - 1 (optional)
% - m = number of edges (optional).
% 
% Luca Carlone
% Politecnico di Torino, Italy
% Jan 27, 2011

add_ind = 0; % assume 1-based indices for nodes

if nargin < 2 % we have to compute m, n
  m = size(edges_id,1);
  min_node_id = min(min(edges_id(:,1)),min(edges_id(:,2))); % to check if a 0-based or 1-based nodes numeration is used
  if min_node_id==0
    error('reduced_incidence_matrix: 0-based indices detected. This is not supported in the current version')
    add_ind = 1;
  elseif min_node_id==1
    add_ind = 0;
  else
    error('reduced_incidence_matrix: minimum node id is neither 0 nor 1.')
  end
  
  n = max(  max(edges_id(:,1)), max(edges_id(:,2))  ) - 1; % assumes 1-based indices in edges
  % NB: n is (number of nodes - 1)
end

row_ids = [ [1:m] [1:m] ]; % first elements, second elements, on the same row
col_ids = [ edges_id(:,1); edges_id(:,2) ];
values = [ -ones(1,m) ones(1,m) ]; 
At_full = sparse(row_ids,col_ids,values,m,n+1);

At = At_full(:,2:end); %reduced incidence matrix
A = At';