function edges = make_edges(num_edges, N)
    %for now, just a dummy function: TODO later -> improve it
    %right now, it just initializes 13 "random" (but always the same)
    %correspondences between 8 cameras
    edges = [];

    if num_edges~=13
        return;
    end
    if N~=8
        return;
    end
    
    edges = [1 2; 3 4; 5 6; 7 8; 1 8; 2 7; 3 6; 4 5; 1 3; 2 4; 3 5; 4 6; 5 7]'; %i-th column means that Tijs(i) is the translation between first row of the column and the second one
end