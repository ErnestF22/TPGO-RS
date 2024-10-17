%Get the second index of a match for a given first index
%function id2=sfm_rawGetMatchFromId(m,id1)
function id2=sfm_rawGetMatchFromId(m,id1)
k1=1;
k2=mod(k1,2)+1;

id2=m(k2,m(k1,:)==id1);
