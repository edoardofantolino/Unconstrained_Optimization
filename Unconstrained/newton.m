function [xk, fk, gradfk_norm, k, xseq] = newton(x0, f, gradf, ...
    Hessf, alpha, kmax, tollgrad)
%
% [xk, fk, gradfk_norm, k, xseq] = newton(x0, f, Hessf, gradf, alpha, kmax,
% tollgrad)
%
% Function that performs the newton optimization method, for a 
% given function for the choice of the step length alpha.
%
% INPUTS:
% x0 = n-dimensional column vector;
% f = function handle that describes a function R^n->R;
% gradf = function handle that describes the gradient of f;
% Hessf = function handle that describes the Hessian of f;
% alpha = the initial factor that multiplies the descent direction at each
% iteration;
% kmax = maximum number of iterations permitted;
% tolgrad = value used as stopping criterion w.r.t. the norm of the
% gradient;
%
% OUTPUTS:
% xk = the last x computed by the function;
% fk = the value f(xk);
% gradfk_norm = value of the norm of gradf(xk)
% k = index of the last iteration performed
% xseq = n-by-k matrix where the columns are the xk computed during the 
% iterations
%


% Initializations
xseq = zeros(length(x0), kmax);

xk = x0;
k = 0;
gradfk_norm = norm(gradf(xk));

while k < kmax && gradfk_norm >= tollgrad
    % Compute the descent direction as solution of
    % Hessf(xk) p = - graf(xk)
    pk = -Hessf(xk)\gradf(xk);
    % Compute the new value for xk
    xk = xk + alpha * pk;
    
    % Compute the gradient of f in xk
    gradfk_norm = norm(gradf(xk));
    
    % Increase the step by one
    k = k + 1;
    
    % Store current xk in xseq
    xseq(:, k) = xk;
end

% Compute f(xk)
fk = f(xk);

% "Cut" xseq to the correct size
xseq = xseq(:, 1:k);

end