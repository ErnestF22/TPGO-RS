%Returns true if the style string contains a line style
%function flag=styleContainsLine(s)
function flag=styleContainsLine(s)
flag=isempty(regexp(s,'[\*\+sodxvph\-:\.]','once'));
