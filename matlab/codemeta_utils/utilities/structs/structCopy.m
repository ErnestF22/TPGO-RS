%Copy fields from a structure to another
%function s1=structCopy(s1,s2)
%Makes a copy of each field in s2 into a field in s1 with the same name
%This function is useful when we want to copy values between two structures
%that have different types. The output will have the same type as s1.
function s1=structCopy(s1,s2)
fields=fieldnames(s2);
for f=fields'
    s1.(f{1})=s2.(f{1});
end
