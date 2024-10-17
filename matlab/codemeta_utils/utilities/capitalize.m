%function s=capitalize(s)
%Turns into capitals the letters in s at the beginning of a word.
function s=capitalize(s)
if iscell(s)
    N=length(s);
    for iN=1:N
        s{iN}=capitalize(s{iN});
    end
else
    offset='A'-'a';
    for iChar=1:length(s)
        if (iChar==1 || s(iChar-1)==' ') && s(iChar)>='a' && s(iChar)<='z'
            s(iChar)=s(iChar)+offset;
        end
    end
end
