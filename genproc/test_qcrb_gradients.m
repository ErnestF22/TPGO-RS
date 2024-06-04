function test_qcrb_gradients

% node_deg = 2;
% 
% qc0=make_rand_stiefel_3d_array(4,4,1);
% if det(qc0) < 0
%     qc0(:,1) = - qc0(:,1);
% end
% rb0=make_rand_stiefel_3d_array(2,2,1);
% if det(rb0) < 0
%     rb0(:,1) = - rb0(:,1);
% end
% 
% disp("check_is_rotation(qc0)")
% disp(check_is_rotation(qc0))
% disp("check_is_rotation(rb0)")
% disp(check_is_rotation(rb0))
% 
% 
% % [qcRt, dqcRt, ~, ~, ~, ~] = rot_geodFun(qc0, []);
% % [rbRt, drbRt, ~, ~, ~, ~] = rot_geodFun(rb0, []);
% M.qc = qc0;
% M.rb = rb0;
% 
% [st,dst,s0,ds0,vVec,ddst,dvVec] = qcrb_randGeodFun(M);
% 
% 
% 
% Qa = [-1.64302828263923e-16	0.106320705792361	-0.889261908456345	-0.444869830049638;
%     0.993981726033715	-0.108890962053455	-0.0116398298358148	-0.00275700117918123;
%     2.81044311504080e-17	0.0248892039991507	-0.444885487710104	0.895241548605309;
%     0.109546010018785	0.988038052621039	0.105615696536622	0.0250160529834866];
% qcd_i = eye(4);
% Ri = [0.0864141652363848	0.0420429561776974	-0.993397284823312;
%     0.256372166990434	0.964519148326362	0.0626995624731963;
%     -0.962094412651845	0.260563503917528	-0.0701748073529320;
%     0.0343546966691522	-0.00647014480909629	0.0656208487016165];
% 
% % Qa_i_1 = [0.00282804903025583, 0.0399962827864656];
% % Qa_i_2 = [0.0626139221207881, -0.997232067403872];
% 
% 
% figure(1)
% f=@(t) mycost_qcrb(st(t), node_deg, Qa, qcd_i, Ri);
% egradf=@(t) myegrad_qcrb(st(t), node_deg, Qa, qcd_i, Ri);
% df=@(t) stiefel_metric([],egradf(t),dst(t));
% funCheckDer(f.qc,df.qc)

end %file function


