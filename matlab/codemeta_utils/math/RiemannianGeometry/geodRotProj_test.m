for(itrial=1:1000)
    R=rot(randn(3,1));
    R0=rot(randn(3,1));
    v=randn(3,1);
    
    [S,tfinalIter]=geodRotProj(R,R0,v,'iter');
    [S,tfinalClosed]=geodRotProj(R,R0,v);
    
    err(itrial)=abs(tfinalClosed-tfinalIter);
end

plot(err)