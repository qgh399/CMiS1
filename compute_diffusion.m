function [x err] = compute_diffusion(type, x, x0, params)
% COMPUTE_DIFFUSION: This function solves a implicity time integration of
% the diffusion term of the velocity field.
%
% INPUT:
%
%    type  -- A type indicator telling us whether we are working on the u
%             or v components of the velocity field.
%    x     -- A warmstart value for the resulting field. 
%    x0    -- The current value of the velocity field component.
%  params  -- Parameter values.
%
% OUTPUT:
%
%   x  -- The diffused velocity field component.
%
% Copyright 2012, Kenny Erleben, DIKU.



%-------------------------------------------------------------------------
% We use a first order implicit time integration
%
%  \frac{d \vecu}{dt} = nu nabla^2 \vec u
%
%  \vec u^{t+1} = \vec u_t + dt*nu*(\frac{d^2}{d x} u^{t+1} + \frac{d^2}{dy} u^{t+1})
%
% Using 2. order central finite difference approximations
%
%   f''(x) approx  (f(x+h) - 2f(x) + f(x-h)) / h^2
%
% Then the u-compnent (v-component is similar) we have
%
%    u(i,j) = u0(i,j) + dt*nu ( ( u(i+1,j) - 2 u(i,j) + u(i-1,j) )/dx*dx 
%                             + ( u(i,j+1) - 2 u(i,j) + u(i,j-1) )/dy*dy )
%
% Now define
%
%    beta   = - \frac{ dt \nu }{ dx^2 }
%    gamma  = - \frac{ dt \nu }{ dy^2 }
%    alpha  = 1 - 2 (beta+gamma) 
%
% Then we have
%
%  alpha*u(i,j) + beta ( u(i+1,j) + u(i-1,j))  
%                            + gamma (u(i,j+1) + u(i,j-1)  ) = u0(i,j)
%
% If we concatenate all u(i,j) into one vector x then we have a linear
% system
%
%     A x = b
%
% The A-matrix  will be a banded matrix with alpha on its diagonal
% and beta and gamma on off diagonals. The b vector will correspond to u0.
%-------------------------------------------------------------------------

gamma = - (params.dt*params.viscosity)/(params.dy*params.dy);
beta  = - (params.dt*params.viscosity)/(params.dx*params.dx);
alpha = 1 - 2*(beta+gamma);

[x, err] = gauss_seidel_solver(type, alpha, beta, gamma, x, x0, params);

end