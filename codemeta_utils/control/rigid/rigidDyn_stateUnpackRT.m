
function [R,w,T,v]=rigidDyn_stateUnpackRT(x)
[R,w]=rotDyn_stateUnpack(x);
[T,v]=realDyn_stateUnpack(x);
