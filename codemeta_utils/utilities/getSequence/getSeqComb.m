%function h=getSeqComb(v)
%Returns a handle to a function h() which, when called, sequentially
%generates all combinations in {1:v(1) x 1:v(2) x ... x 1:v(length(v))}
function h=getSeqComb(v)
cnt=0;
h=@counter;

    function out=counter()
        cnt=cnt+1;
        out=myind2sub(v,cnt);
    end

end
