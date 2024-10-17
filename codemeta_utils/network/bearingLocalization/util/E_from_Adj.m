function E = E_from_Adj(Adj)


symAdj = (.5*(Adj + Adj') + 10^(-5)) > 1;   
   

[ii,jj] = find( triu(symAdj,1) );
E = [ii jj];
 

end

