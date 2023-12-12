%function [V,E]=graphSquashNonReferenced(E)
%Find the set V of vertices that appear in the edge list E and renumber
%them from one to length(V)
function [V,E]=graphSquashNonReferenced(E)
V=unique(E(:));
E=mapValues(E,[V (1:length(V))']);
