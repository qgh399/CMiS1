function phi = compute_advection(type, phi, u, v, params)
% COMPUTE_ADVECTION: Advects the velocity field using a semi-Lagrangian
% time integration technique. (a.k.a.method of characteristics)
%
% INPUT:
%    type - The type of field we are working on.
%    phi  - The current value of the field.
%      u  - The current value of the x-component of the velocity field.
%      v  - The current value of the y-component of the velocity field.
%  params - The parameter values.
% OUTPUT:
%    phi  - The updated field value.
%
% Copyright 2012, Kenny Erleben, DIKU.


%-------------------------------------------------------------------------
% For details read
%
%    Jos Stam, "Stable Fluids", In SIGGRAPH 99 Conference Proceedings,
%    Annual Conference Series, August 1999, 121-128.
%
%    Jos Stam, "Real-Time Fluid Dynamics for Games". Proceedings of the
%    Game Developer Conference, March 2003.
%
%    Mark J. Harris , William V. Baxter , Thorsten Scheuermann , Anselmo
%    Lastra, Simulation of cloud dynamics on graphics hardware, Proceedings of the ACM SIGGRAPH/EUROGRAPHICS conference on Graphics hardware, July 26-27, 2003, San Diego, California
%
%   Mark Harris, Fast Fluid Dynamics Simulation on the GPU, GPU Gems, Chapter 38
%   Addison-Wesley, 2004 (http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html)
%
%-------------------------------------------------------------------------

%--- Convert grid nodes into particles -----------------------------------
[y_cur, x_cur] = meshgrid(1:params.I,1:params.J);

%--- Back track particles to the position they came from -----------------
x_old = x_cur - (u*(params.dt/params.dx));
y_old = y_cur - (v*(params.dt/params.dy));

%--- Particles that is traced outside the domain is projected back onto
%--- the domain 
x_old(x_old>(params.I)) = params.I;
x_old(x_old<1)          = 1;

y_old(y_old>(params.J)) = params.J;
y_old(y_old<1)          = 1;

%--- Update the value of phi to the value at the old positions ------------
phi = interp2(y_cur, x_cur, phi, y_old, x_old,'*linear');

%--- Apply proper boundary conditions for the updated field ---------------
phi = set_boundary_conditions(type,phi,params);
end
