function [param]=extract_params(pp,draws)
    % note since this is user customizable this does NOT do any input
    % checking

     % Each entry in betai needs to extract coefficients from pp and the
     % correct components from draws so that there is a corresponding value
     % of betai for each element in dtable.x2.
     
     %param.betai = [pp(2)*draws.v(:,1)  pp(3)*draws.v(:,1)+pp(3)*draws.v(:,2)];
     
     % The betamask is a matrix where each row corresponds to an element of
     % betai and consists of two columns.
     % column 1: the index of the corresponding column of x2
     % column 2: the index of the corresponding column of dbeta
     % the first column should never be greater than dim(x2)
     % the second column should never be greater than dim(dbeta)
     
     %param.betamask = [1 1; 2 1; 2 2];
     
     % Each entry in dbeta corresponds to the derivative of beta_i.
     % It is possible to have fewer entries in dbeta than there are in
     % betai.
     
     %param.dbeta = [draws.v(:,1:2)];
     
    % Consider the following case of correlated normal random coefficients
     
    %param.betai = [pp(1)*draws.v(:,1)  pp(2)*draws.v(:,1)+pp(3)*draws.v(:,2)];
    %param.dbeta = [draws.v(:,1:2)];
    %param.betamask = [1 1; 2 1; 2 2;];

    % Consider the following case where beta_i = pi * y_i
     
    param.betai = [pp(1)*draws.v(:,1) pp(2)*draws.v(:,1)+pp(3)*draws.v(:,2)];
    param.dbeta = [draws.v(:,1:2)];
    
    % 1-1 implies first X first draws
    % 2-2 implies second X second draw
    % 2-1 implies second X first draw ( this is the correlation term of the
    % cholesky root)
    param.betamask = [1 1; 2 1; 2 2;];
    
    % These are upper and lower bounds on parameter values. Some parameters
    % should only take on positive value (such as standard deviations).
    param.lb = [0 -Inf 0];
    param.ub = [Inf Inf Inf];
     
end