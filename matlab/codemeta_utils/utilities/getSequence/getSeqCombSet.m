%function h=getSeqCombSet(s,n)
%Returns a handle to a function h() which, when called, sequentially
%generates all the length(s)^n possible words of length n where each digit is taken
%from the set s
function h=getSeqCombSet(s,n)
nextState=getSeqComb(length(s)*ones(1,n));
h=@nextComb;

    function w=nextComb()
        state=nextState();
        w=s(state);
    end
end
