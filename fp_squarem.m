function [xnew]=fp_squarem(g,x,varargin)
% USAGE: [delta]=fp_squarem(f,x,optional arguments)
% Applies the SQUAREM acceleration method for fixed point iterations port
% of SQUAREM in R by Ravi Varadhan
% C. Conlon (November 2014)
% 
% f is a function handle that performs a single fixed point iteration
% x is an initial guess of the value such that f(x) = x
% 
% Optional Arguments
%
% algorithm: "squarem" (default), "contraction"
% alphaversion: (1, 2, 3) different step-length-schemes (default=3)
% noisy : 0 for no output, 1 for final output only, 2 for iteration level output (default=0)
% con_tol : tolerance for convergence (default 1e-13)
% max_iter: maximum number of iterations before failure (default 1e5)
% stepmin0/stepmax0: initial step min/max (default=1)
% mstep: accepted step scaling factor (default 4).

p = inputParser;
%addRequired(p,'g');
%addRequired(p,'x',@isvector);
defaultAlg = 'squarem';
validAlg = {'squarem','contraction'};
checkAlg = @(x) any(validatestring(x,validAlg));

addParameter(p,'algorithm',3);
addParameter(p,'alphaversion',3);
addParameter(p,'mstep',4);
addParameter(p,'con_tol',1e-14);
addParameter(p,'noisy',0);
addParameter(p,'max_iter',1e5);
addParameter(p,'stepmin0',1);
addParameter(p,'stepmax0',1);
p.KeepUnmatched = true;
parse(p,varargin{:});

iter=0;
fpevals=0;
xchng = Inf;

%display(p.Results);
stepmin = p.Results.stepmin0;
stepmax = p.Results.stepmax0;

while((iter < p.Results.max_iter)  & (xchng > p.Results.con_tol)  )
    iter=iter+1;
    x1=g(x);
    fpevals=fpevals+1;
    
    if(strcmp('contraction',p.Results.algorithm))
        xchng=mean(abs(x1-x));
        xnew=x1;
        display([' Iteration ' num2str(iter) ' and ' num2str(xchng) ])
        x=xnew;
        continue
    end
    q1=x1-x;
    x2 =g(x1);
    fpevals=fpevals+1;
    q2=x2-x1;
    
    % Form quadratic terms
    sr2 = q1'*q1;
    sq2 = sqrt(q2'*q2);
    sv2 = (q2-q1)'*(q2-q1);
    srv = q1'*(q2-q1);
    
    % Get the step-size
    [alpha]=compute_alpha(sv2,sr2,srv,stepmin,stepmax,p.Results.alphaversion);
    xtmp = x + 2 * alpha * q1+ alpha.^2 * (q2-q1);
    % Fixed point iteration beyond the quadratic step
    xnew = g(xtmp);
    fpevals=fpevals+1;
    
    if(any(isnan(xnew)))
        display('Error')
        xnew=x2;
    end
    if (alpha == stepmax)
        stepmax = p.Results.mstep*stepmax;
    end
    if ((alpha == stepmin) & alpha < 0)
        stepmin = p.Results.mstep*stepmin;
    end
    if(p.Results.noisy==2)
        display([' Iteration ' num2str(iter) ' delta change: ' num2str(norm(xnew-xtmp,Inf)) ...
            ' and alpha=' num2str(alpha) ' ( ' num2str(stepmin) ' , ' num2str(stepmax) ')'])
    end
    xchng=mean(abs(xnew-xtmp));
    x=xnew;
end

if(p.Results.noisy)
    display(['Fixed Point Evaluations: ' num2str(fpevals) ' Optimization Error ' num2str(xchng) ])
end

    function [alpha]=compute_alpha(sv2,sr2,srv,stepmin,stepmax,alphaversion)
        switch alphaversion
            case 1
                alpha = -srv/sv2;
            case 2
                alpha = -sr2/srv;
            case 3
                alpha = sqrt(sr2/sv2);
        end
        alpha = max(stepmin,min(stepmax,alpha));
    end
end