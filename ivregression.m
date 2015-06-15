 function [beta,resid]=ivregression(Y,X,Z,W)
        if nargin <4
            W = (Z'*Z) \ eye(size(Z,2));
        end
        beta=(X'*Z * W * Z'*X)\(X'*Z * W * Z'*Y);
        if nargout >1
            resid=Y-X * beta;
        end
 end
