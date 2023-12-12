function [ output_args ] = contraction3D_plotGershgorinDisc( A )
%Given a square A matrix, plot its gershgorin disc and real eigenvalues
A
figure;
%plot the eigenvalues
eG = eig(A);
plot(real(eG), imag(eG), 'rx');
hold on
%plot the discs
x = 0:1e-3:2*pi;
for i = 1:length(A)
    r = sum(abs(A(i,:))) - abs(A(i,i)); %radius is sum of all off-diagonal terms
    center = [real(A(i,i)), imag(A(i,i))];
    plot(r*sin(x)+center(1), r*cos(x)+center(2));    
end

end

