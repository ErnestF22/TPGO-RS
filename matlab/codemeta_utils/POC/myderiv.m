% This code was generated using ADiGator version 0.4.2
% ©2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function y = myderiv(x,k)
global ADiGator_myderiv
if isempty(ADiGator_myderiv); ADiGator_LoadData(); end
Gator1Indices = ADiGator_myderiv.myderiv.Gator1Indices;
Gator1Data = ADiGator_myderiv.myderiv.Gator1Data;
% ADiGator Start Derivative Computations
cada1f1dx = k(:).*x.dx;
cada1f1 = k.*x.f;
cada1td1 = zeros(3,1);
cada1td1(Gator1Indices.Index1) = cada1f1dx;
y.dx = cada1td1;
y.f = [cada1f1;1];
%User Line: y=[k.*x;1];
y.dx_size = [4,3];
y.dx_location = Gator1Indices.Index2;
end


function ADiGator_LoadData()
global ADiGator_myderiv
ADiGator_myderiv = load('myderiv.mat');
return
end