function cf = cf_uniform2( t, a, b )
%CF_UNIFORM returns the value of characteristic function of uniform distribution 
%
% Purpose:
%   This function returns the value of characteristic function of uniform
%   distribution.
%
% Arguments:
%   t ... a scalar or vector
%   a ... the lower boundary of the uniform distribution. Must be a real scalar.
%   b ... the upper limit of the uniform distribution. Must be a real
%   scalar.
%
% Usage:
%
%   cf = cf_uniform (t, a, b);
%
% History:
%   Created TD, 1-feb-13
%   Modified TD, 17-feb-13: Added code to deal with case when t == 0.
%
% ---------------------------------------------------------------------

    %% Check the arguments
%     if (nargin() ~= 3)
%         help cf_uniform
%         error ('For usage see above')
%         
%     elseif (a >= b)
%         help cf_uniform
%         error (char('The lower boundary of the uniform distribution must be ', ...
%             'less than the upper boundary'))
%         
%     elseif (isscalar(a) == false || isreal(a) == false)
%         error ('The lower boundary of uniform distribution must be a real scalar')
%                 
%     elseif (isscalar(b) == false || isreal(b) == false)
%         error ('The upper boundary of uniform distribution must be a real scalar')
%         
%     elseif (isvector(t) == false)
%         help cf_uniform
%         error ('The argument t can be a vector or a scalar')
%     end
    
    
    %% Do the calcualtion
    %  Modified TD, 17-feb-13: Added code to deal with the situation when 
    %  t == 0
%     temp = t == 0;
%     if (isempty(temp) == false)
%         cf(temp) = 1;
%     end
%     
%     temp = ~temp;
%     if (isempty(temp) == false)
%         cf(temp) = (exp(1i*t(temp)*b) - exp(1i*t(temp)*a)) ...
%             ./ (1i * t(temp) * (b - a));
%     end
    
cf = (exp(1i*t*b) - exp(1i*t*a)) ./ (1i * t * (b - a));

end

