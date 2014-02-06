function [u v smoke diff_err_u diff_err_v press_err1 press_err2] =  compute_fractional_step(u, v, smoke, params)
% COMPUTE_FRACTIONAL_STEP: Solves an Navier Stokes equation iteration
% for velocity and smoke field.
% INPUT:  
%    u  -- The current value of the x-component of the velocity field.
%    v  -- The current value of the y-component of the velocity field.
% smoke -- The current value of the smoke density field.
% Output:
%    u  -- The updated value of the x-component of the velocity field.
%    v  -- The updated value of the y-component of the velocity field.
% smoke -- The updated value of the smoke density field.
%
% Copyright 2012, Kenny Erleben, DIKU.

%--- First phase we time integrate the velocity field of the fluid --------

%--- First step: integrate buoyancy force ---------------------------------
fv = compute_buoyancy_force(smoke, params);

v = v + fv*params.dt; 

%--- Second step integrate vorticity force --------------------------------
[fu fv] = compute_vorticity_confinement(u,v,params);

v = v + fv*params.dt; 
u = u + fu*params.dt; 

%--- Third step: Integrate diffusion term ---------------------------------
u0 = u;
v0 = v;
[u, diff_err_u]  = compute_diffusion(1, u, u0, params);

[v, diff_err_v]  = compute_diffusion(2, v, v0, params);

%--- Intermediate Step: extra projection for better accuracy --------------
[u, v, press_err1] = compute_pressure_projection(u,v,params);

%--- Fourth Step: Integrate advection term --------------------------------
u0 = u;
v0 = v;
u  = compute_advection(1,u,u0,v0,params);
v  = compute_advection(2,v,u0,v0,params);

%--- Fifth Step: Finaly projection step to make fluid divergence free -----
[u, v, press_err2] = compute_pressure_projection(u,v,params);

%--- Second phase, we time integrate our trace (smoke) density field ------

%--- First step: We add some more new smoke -------------------------------

source = zeros(params.I,params.J);

source( params.smoke_source_i, params.smoke_source_j ) = params.smoke_rate; 

smoke = smoke + source*params.dt; 

%--- Second step: We advect the smoke by the fluid velocity ---------------
smoke = compute_advection(0,smoke,u,v,params);

end