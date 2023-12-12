function E=bnFranchi_graphEdges(NNodes)
E=[1 2; 2 1; [repmat([1;2],NNodes-2,1) reshape(repmat(3:NNodes,2,1),[],1)]];
