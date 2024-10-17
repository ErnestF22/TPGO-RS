%function root=depPackageFindCommonRoot(nameList)
%Given a cell array of strings, find the common root among all of them
function root=depPackageFindCommonRoot(nameList)
%convert to char
charList=char(nameList);
%use diff to find how many columns have all the same character
rootLength=find(max(diff(double(charList))~=0),1,'first')-1;

root=nameList{1}(1:rootLength);