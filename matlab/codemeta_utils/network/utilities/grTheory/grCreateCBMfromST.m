function C = createCBMfromST(edges_id,l)
% For a directed graph with n+1 nodes and m edges, the function computes 
% an l x m cycle basis matrix, with l = m-n. 
% - "edges_id" is an m x 2 matrix, such that each row [i j] encodes an edge
% starting from i and ending at j (assume the indices start from 1)
% - l = cyclomatic number
% 
% Luca Carlone
% Georgia Institute of Technology, USA
% May 20, 2013
%

m = size(edges_id,1);
n = m-l;
C_st_lc = spalloc(l,m,l*n+1); % we first build the part corresponding to the spanning tree

%% sparse adjacency matrix
Ad = sparse(zeros(n+1,n+1)); 
for k=1:m
  id1 = edges_id(k,1);
  id2 = edges_id(k,2);
  Ad(id1,id2) = 1;
  Ad(id2,id1) = 1;
end
%% Compute spanning tree
[tree pred] = graphminspantree(Ad);

%% Distinguish edges in the spanning tree from chords
edges_st_id = [];
for i=1:n+1
  id1 = pred(i);
  id2 = i;
  k = find(edges_id(:,1)==id1 & edges_id(:,2)==id2);
  if length(k)>0
    if (length(k)>1) error('createCBMfromST: graph is not simple');  end
    edges_st_id(end+1) = k;
  else
    k = find(edges_id(:,1)==id2 & edges_id(:,2)==id1);
    if length(k)>0
      if (length(k)>1) error('createCBMfromST: graph is not simple');  end
      edges_st_id(end+1) = k;
    end
  end
end
edges_lc_id = setdiff([1:m],edges_st_id);

id_st = edges_id(edges_st_id,:); % edges in the spanning tree
id_lc = edges_id(edges_lc_id,:); % chords

for i=1:l % we build the cycle basis matrix by rows
  id1 = id_lc(i,1);
  id2 = id_lc(i,2);
  C_st_lc(i,1:n) = findOrientedPathOnSTFast(Ad , id_st, id2, id1, pred);
  C_st_lc(i,n+i) = 1;
end

C = spalloc(l,m,l*n+1);
C(:,edges_st_id) = C_st_lc(:,1:n);
C(:,edges_lc_id) = C_st_lc(:,n+1:end);