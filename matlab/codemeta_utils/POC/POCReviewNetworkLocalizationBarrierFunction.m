function POCReviewNetworkLocalizationBarrierFunction

e_sigma= @(w1,w2,w3) rot_dist(eye(3),rot(w1)*rot(w2)*rot(w3));

V=@(w1,w2,w3) 1/3*(log(pi^2)-log(pi^2-e_sigma(w1,w2,w3)^2));

z=zeros(3,1);
u1=pi*sphere_randn();
u2=pi*sphere_randn();
u3=pi*sphere_randn();

function test(w1,w2,w3)
    disp('Angles')
    disp([norm(w1) norm(w2) norm(w3)])
    disp('V(w1,w2,w3)')
    disp(V(w1,w2,w3))
end

test(z,z,z)

test(u1,z,z)

test(u1,u2,z)

end