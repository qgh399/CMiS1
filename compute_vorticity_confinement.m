function [fu fv] = compute_vorticity_confinement(u,v,params)
% COMPUTE_VORTICITY_CONFINEMENT -- Computes a vorticity force that will try
% and induce more 'rotational' motion into the fluid. Note: there is NO
% gridsize dependent scaling implemented! 
%
% INPUT:
%
%   u    - The current value of the x-component of the velocity field.
%   v    - The current value of the y-component of the velocity field.
% params - The parameter values.
%
% OUTPUT:
%
%   fu   - The current value of the x-component of vorticity confinement force.
%   fv   - The current value of the y-component of vorticity confinement force.
%
% Copyright 2012, Kenny Erleben, DIKU.



%-------------------------------------------------------------------------
% For details read 
%
%  @inproceedings{Fedkiw:2001:VSS:383259.383260,
%   author = {Fedkiw, Ronald and Stam, Jos and Jensen, Henrik Wann},
%   title = {Visual simulation of smoke},
%   booktitle = {Proceedings of the 28th annual conference on Computer graphics and interactive techniques},
%   series = {SIGGRAPH '01},
%   year = {2001},
%   isbn = {1-58113-374-X},
%   pages = {15--22},
%   numpages = {8},
%   url = {http://doi.acm.org/10.1145/383259.383260},
%   doi = {http://doi.acm.org/10.1145/383259.383260},
%   acmid = {383260},
%   publisher = {ACM},
%   address = {New York, NY, USA},
%   keywords = {Euler equations, Navier-Stokes equations, computational fluid dynamics, participating media, semi-Lagrangian methods, smoke, stable fluids, vorticity confinement},
%  }
%

i = 2:params.I-1;
j = 2:params.J-1;

%--- w =  nabla times (u,v)^T --------------------------------------------
w      = zeros(params.I,params.J);
w(i,j) = (v(i+1,j)-v(i-1,j))/(2*params.dx)...
       - (u(i,j+1)-u(i,j-1))/(2*params.dy);
w = set_boundary_conditions(0,w,params);
%--- n = nabla |w| -------------------------------------------------------
nu = zeros(params.I,params.J);
nv = zeros(params.I,params.J);

aw      = abs(w);
nu(i,j) = (aw(i+1,j)-aw(i-1,j))/(2*params.dx);
nv(i,j) = (aw(i,j+1)-aw(i,j-1))/(2*params.dy);

nu = set_boundary_conditions(0,nu,params);
nv = set_boundary_conditions(0,nv,params);
%--- N = n / | n | -------------------------------------------------------
norm_n = sqrt(nu.*nu + nv.*nv);
Nu     = nu ./ (norm_n+eps);
Nv     = nv ./ (norm_n+eps);

%--- f = N x w -----------------------------------------------------------
h  = min(params.dx, params.dy);
fu = (Nv .*  w)*(params.vorticity*h);
fv = (Nu .* -w)*(params.vorticity*h);
end