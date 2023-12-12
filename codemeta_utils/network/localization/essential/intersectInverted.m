%Find subset of a set of integers contained in another set
%function idxs2=intersectInverted(s1,s2)
%Returns the indeces of the elements in s2 that are contained in s1.
%Requires that s1 and s2 contains only positive integers
%The function works by making an inverted index for the elements in s1 and
%then looking up the elements in s2
function idxs2=intersectInverted(s1,s2)
N=max(max(s1),max(s2));
flags1=sparse(1,N);
flags1(s1)=1;
idxs2=find(flags1(s2));
