function [u v smoke elapsed] =  compute_fractional_step(u, v, smoke, params)
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
elapsed = zeros(1,8);
i_elapsed = 1;
%--- First step: integrate buoyancy force ---------------------------------

t = tic;
fv = compute_buoyancy_force(smoke, params);

v = v + fv*params.dt; 
elapsed(i_elapsed) = toc(t);
i_elapsed = i_elapsed+1;

%--- Second step integrate vorticity force --------------------------------
t = tic;
[fu fv] = compute_vorticity_confinement(u,v,params);

v = v + fv*params.dt; 
u = u + fu*params.dt; 
elapsed(i_elapsed) = toc(t);
i_elapsed = i_elapsed+1;

%--- Third step: Integrate diffusion term ---------------------------------
t = tic;
u0 = u;
v0 = v;
u  = compute_diffusion(1, u, u0, params);
v  = compute_diffusion(2, v, v0, params);
elapsed(i_elapsed) = toc(t);
i_elapsed = i_elapsed+1;

%--- Intermediate Step: extra projection for better accuracy --------------
t = tic;
[u v] = compute_pressure_projection(u,v,params);
elapsed(i_elapsed) = toc(t);
i_elapsed = i_elapsed+1;

%--- Fourth Step: Integrate advection term --------------------------------
t = tic;
u0 = u;
v0 = v;
u  = compute_advection(1,u,u0,v0,params);
v  = compute_advection(2,v,u0,v0,params);
elapsed(i_elapsed) = toc(t);
i_elapsed = i_elapsed+1;

%--- Fifth Step: Finaly projection step to make fluid divergence free -----
t = tic;
[u v] = compute_pressure_projection(u,v,params);
elapsed(i_elapsed) = toc(t);
i_elapsed = i_elapsed+1;

%--- Second phase, we time integrate our trace (smoke) density field ------

%--- First step: We add some more new smoke -------------------------------
t = tic;
source = zeros(params.I,params.J);

source( params.smoke_source_i, params.smoke_source_j ) = params.smoke_rate; 

smoke = smoke + source*params.dt; 
elapsed(i_elapsed) = toc(t);
i_elapsed = i_elapsed+1;

%--- Second step: We advect the smoke by the fluid velocity ---------------
t = tic;
smoke = compute_advection(0,smoke,u,v,params);
elapsed(i_elapsed) = toc(t);
i_elapsed = i_elapsed+1;

end