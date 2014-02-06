function draw_fluid(u, v, smoke, params)
% DRAW_FLUID: Draws the current fluid solution
% INPUT:
%
%        u - The x component of the velocity field.
%        v - The y component of the velocity field.
%    smoke - Smoke density field.
%   params - The parameter values.
%
% Copyright 2012, Kenny Erleben, DIKU.

color =  smoke;
color(color<0) = 0;
c = color(params.idx);
patch_color = [c, c, c, c];
patch( params.patch_x', params.patch_y', patch_color' );
shading interp;
end
