function [ fv ] = compute_buoyancy_force(rho, params)
% COMPUTE_BUOYANCY_FORCE: Imperical buoyancy force generation works to
% some extent but is not physically founded.
%
% INPUT:
%
%    rho  -   The density field for which the bouncy force is wanted.
%  params -   Parameter values.
%
% OUTPUT:
%
%       fv -  The y-component of the resulting boouncy force.
%
% Copyright 2012, Kenny Erleben, DIKU.

T_ambient = sum( rho(:) );
T_ambient = T_ambient ./ (params.I*params.J);
fv = - params.drop .* rho + params.lift .* (rho - T_ambient);
end