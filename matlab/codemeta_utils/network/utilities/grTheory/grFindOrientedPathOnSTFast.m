function [orientedPath] = findOrientedPathOnSTFast( A , id_st, id_start, id_end, prec)
% find an oriented path on a spanning tree
% id_st is the collection of oriented edges forming a spanning tree with n+1 nodes (each row is id1->id2)
% id_start is the id of the starting node
% id_end is the id of the node to reach
% orientedPath is a vector of n elements in which the i-th element is zero
% if the node does not belong to the path, +1 if the edge is traversed 
% according to its orientation when travelling from id_start to id_end, and
% it is -1 if it is traversed in the opposite direction
%
% Luca Carlone - GeorgiaTech - 2013/10/6
n = size(A,1)-1; % size of the spanning tree
orientedPath = sparse(zeros(1,n));

%% we start at the last node
id2 = id_end;
while 1
  id1 = prec(id2);
  if id1 == 0
    break; % we reached the root, we have to start from id_start and reach the root from there
  end
  
  %% fill in the corresponding entry in orientedPath (with sign) 
  edge = find(id_st(:,1)==id1 & id_st(:,2)==id2);
  if isempty(edge)
    edge = find( id_st(:,1)==id2 & id_st(:,2)==id1);
    orientedPath(edge) = -1;
    if isempty(edge)
      error('findOrientedPathOnSTFast: incorrect edge')
    end
  else
    orientedPath(edge) = +1;
  end
  
  if id1 == id_start
    return; % id_start is along the path to the root
  end
  id2 = id1; 
end

%% we start from id_start
id2 = id_start;
while 1
  id1 = prec(id2);
  if id1 == 0
    return; % we reached the root
  end
  
  %% fill in the corresponding entry in orientedPath (with sign) 
  edge = find(id_st(:,1)==id1 & id_st(:,2)==id2); 
  if isempty(edge)
    edge = find( id_st(:,1)==id2 & id_st(:,2)==id1);
    if orientedPath(edge) == 0  % not visited in the first while loop
      orientedPath(edge) = +1; % note, the sign is inverted wrt before, to follow the cycle
    else
      orientedPath(edge) = 0;
    end
    if isempty(edge)
      error('findOrientedPathOnSTFast: incorrect edge')
    end
  else
    if orientedPath(edge) == 0  % not visited in the first while loop
      orientedPath(edge) = -1; % note, the sign is inverted wrt before, to follow the cycle
    else
      orientedPath(edge) = 0;
    end
  end
  id2 = id1; 
end

