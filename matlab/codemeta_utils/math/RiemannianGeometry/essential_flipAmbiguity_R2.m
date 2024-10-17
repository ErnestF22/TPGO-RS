function Q2=essential_flipAmbiguity_R2(Q2,k)
Rxpi=diag([1;-1;-1]);
Rzpi=diag([-1;-1;1]);

switch k
    case 1
        %nothing to do
    case 2
        Q2=Rxpi*Rzpi*Q2;
    case 3
        Q2=Rzpi*Q2;
    case 4
        Q2=Rxpi*Q2;
    otherwise
        error('Value of k invalid')
end
