function [pjt,pijt]=rc_share_safe(delta, params,draws,mkt)
% This gives the RC choice probabilities after integration (pjt)
% Also gives the individual choice probabilities (pijt) for derivatives
%
% This function should be re-written in mex/C++ for improved speed
    [ns k] = size(draws.v);
    u1 = repmat(delta,[1 ns]) + mkt.x2*params.betai';
    utils = exp(u1);
    pijt = bsxfun(@rdivide,utils,1+sum(utils));
    pjt=pijt*draws.w;
    
% under/overflow safety
% should write an under/overflow safe version of this function
%    m = max(u1); 
%    utils = exp(bsxfun(@minus,u1,m));
%    pijt = bsxfun(@rdivide,utils,exp(-m)+sum(utils));

end
