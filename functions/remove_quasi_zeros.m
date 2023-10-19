function m_out = remove_quasi_zeros(m_in, thr)
%REMOVE_QUASI_ZEROS sets to zero all elements below a certain threshold
%to improve output readability; if not given, default threshold is 0.001
if nargin == 1
    thr = 0.001;
end

m_in(abs(m_in)<=thr) = 0;

m_out = m_in;
end %function
