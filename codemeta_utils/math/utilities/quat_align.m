%function q=quat_align(q)
%Given a quaternion q, q and -q represent the same rotation. This function
%flips the sign of each column of q such that the scalar product with the
%first one is maximized

function q=quat_align(q)

N=size(q,2);
qref=q(:,1);

for(iq=1:N)
    if(q(:,iq)'*qref<0)
        q(:,iq)=-q(:,iq);
    end
end
