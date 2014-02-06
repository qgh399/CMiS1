function [u, v] = compute_pressure_projection(u, v, params)
% COMPUTE_PRESSURE_PROJECTION: This function compute the pressure
% projection that will make the input velocity field divergence free
% upon return.
%
% INPUT:
%
%     u   - The x-component of the velocity field.
%     v   - The y-component of the velocity field.
% params  - The parameter values.
%
% OUTPUT:
%
%    u   -- The updated x-component of the velocity field.
%    v   -- The updated y-component of the velocity field.
%
% Copyright 2012, Kenny Erleben, DIKU.


%--- The outline of the projection step ----------------------------------
% The main ideas behind the pressure projection is that we are the left
% with the PDE term
%
%   \vec u^{t+1} = \vec u^t -  dt/rho*( nabla p ) 
%
% Here \vec u = (u,v)^ is the velocity field. Also we wish that
%
%   \nabla \cdot \vec u^{t+1} = 0
%
% Now computing the divergence of the first equation gives us
%
%   \nabla \cdot\vec u^{t+1} = \nabla \cdot\vec u^t -  dt*( \nabla \cdot nabla p ) 
% 
% Which we simpplify to the Poisson equation
%
%  \nabla^2 p = (rho/dt)  \nabla \cdot\vec u^t
% 
% As rho/dt is a postive constant for the entire fluid domain we can simply
% let the pressure field arbsorp this term. That is we redefine 
%
%    p <- p*dt/rho
%
% That way we do not need to care about this constant value.

%--- Calculate the divergence of current velocity field -------------------
% From calculus we know 
%
%    div( (u,v) ) = nabla \cdot (u,v)^ = d/dx u + d/dy v
%
% Using a CD approximation we have that
%
%    d/dx u_{i,j} \approx (u_{i+1,j}-u_{i-1,j})/ (2 dx) 
%
% Similar for d/dy approximations.

b = zeros(params.I, params.J);
i = 2:params.I-1;
j = 2:params.J-1;

b(i,j) = ( u(i+1,j) - u(i-1,j) ) / (2*params.dx)...
       + ( v(i,j+1) - v(i,j-1) ) / (2*params.dy);
             
b      = set_boundary_conditions(0,b, params);

%--- Solve for pressure field ---------------------------------------------
% Now we actually solve the Poisson equation, one should look into the
% solver for details on how this is done.
%
% We have
%
%   d^2/dx^2 p  + d^2/dy^2 p = nabla \cdot \vec u^t  
%
% The right hand side is the minus b-field that we computed above.
%
% Using 2nd order central difference approximations we get
%
%  (p(i+1,j) - 2 p(i,j) + p(i-1,j))/dx^2 
%          + (p(i,j+1) - 2 p(i,j) + p(i,j-1))/dy^2 = b(i,j)
%
% Cleaning up
%
%  -2(1/dx*dx + 1/dy*dy) p(i,j) + 1/dx*dx ( p(i+1,j) + p(i-1,j) )
%                               + 1/dy*dy ( p(i,j+1) + p(i,j-1) ) = b(i,j)
%
% Defining
% 
%  beta  = 1/(dx*dx)
%  gamma = 1/(dy*dy)
%  alpha = - 2 (beta+gamma)
%
% We have
%
%  alpha p(i,j) + beta ( p(i+1,j)+p(i-1,j)) 
%               + gamma (p(i,j+1)+p(i,j-1)) = b(i,j)
%
% Concatenating p(i,j) into a vecto x we find a linear system
%
%      A x = b
%
% where the A-matrix is a banded matrix having alpha on its
% diagonal and beta and gamma on off diagonals.
% -------------------------------------------------------------------------
beta   = 1/(params.dx*params.dx);
gamma  = 1/(params.dy*params.dy);
alpha  = -2*(beta + gamma);

p = zeros(params.I,params.J);
p = gauss_seidel_solver(0, alpha, beta, gamma, p, b, params);

%--- Make velocity field divegence free by adding pressure forces ---------
% Basically we are using Helmholtz-Hodge decomposition. To compute the
% gradient of the pressure field 
%
% grad p = nabla p = (dp/dx, dp/dy)^
%
% We apply a CD scheme as we did above for the velocity field components.
%

nabla_px = (p(i+1,j)-p(i-1,j))/(2*params.dx);
nabla_py = (p(i,j+1)-p(i,j-1))/(2*params.dy);

u(i,j) = u(i,j) - nabla_px;
v(i,j) = v(i,j) - nabla_py;

%--- Finally we make sure that boundary conditions are correct ------------
u = set_boundary_conditions(1,u,params);
v = set_boundary_conditions(2,v,params);

end