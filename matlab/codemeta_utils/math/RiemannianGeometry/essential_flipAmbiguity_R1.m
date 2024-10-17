function Q1=essential_flipAmbiguity_R1(Q1,k)
Rxpi=diag([1;-1;-1]);
switch k
    case {1,3}
        %nothing to do
    case {2,4}
        Q1=Rxpi*Q1;
    otherwise
        error('Value of k invalid')
end
