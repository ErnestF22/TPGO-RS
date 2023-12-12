% POC that the Gershgorin circle is always moving left as ||omega||
% decreases for the Omega matrix if m creates pos-def matrix
close all;
A = rand(2,2);
B = (A+A')/2+2*eye(2);
m = [B(1,1);B(1,2);B(2,2)]
kd = -5;
kv = 10;
normW = linspace(0,10);
normW_max = max(normW)

alpha = normW*(m(3)*kv-m(2))/4;
alpha_max = normW_max*(m(3)*kv-m(2))/4;
magDiff_rad = ( abs(m(3)/8*normW_max^2 + alpha_max*i) - abs(m(3)/8*normW.^2 + alpha*i) ).^2;
magDiff_cent = ( m(2)/4*(normW_max^2-normW.^2) ).^2;

centroid = m(2)/4*normW.^2;
radius = abs(m(3)/8*normW.^2 + alpha*i);

% Compare magnitude of centroid change vs radius change
A = diff(centroid).^2;
B = diff(radius).^2;
any(B<A) % This should return false==0

figure
plot(normW, -centroid + radius)
hold on
plot(normW, centroid+ radius)

