function [x] = gauss_seidel_solver(type, alpha, beta, gamma, x, b, params, iters)
% GAUSS_SEIDEL_SOLVER: This function uses a Gauss Seidel solver to solve
% the linear system A x = b given by the coefficients alpha, beta, gamma,
% and b.
%
% INPUT:
%    type  - A type indicator telling what type of field we are working on.
%    alpha - The diagonal values of the A-matrix
%    beta  - The off-diagonal values of the A-matrix
%    gamma - The off-diagonal values of the A-matrix
%        x - The field that we wish to solve for.
%        b - The right hand side field of the matrix equation.
%   params - The parameter values.
%OUTPUT:
%        x - The solution field.
%
% Copyright 2012, Kenny Erleben, DIKU.

x = set_boundary_conditions(type, x, params);

i = 2:params.I-1;
j = 2:params.J-1;

for iter=1:iters
    
  x(i,j) = (...
    b(i,j)...
    - beta.*( x(i+1,j) + x(i-1,j))...
    - gamma*(x(i,j+1) + x(i,j-1) )...
    )./alpha;
  
  x = set_boundary_conditions(type, x, params);
      
end

end
