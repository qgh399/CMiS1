function [phi] = set_boundary_conditions(type, phi, params)
% SET_BOUNDARY_CONDITIONS: Sets the value of the ghost nodes such that
% boundarcy conditions are satisfied.
%
% We always have a boxed domain where a static solid wall is placed between
% the outer one ring ghost nodes and inner domain nodes. 
%
% INPUT
%
%    type  -  Field type indicator, can be 0 (pressure field), 
%             1 (u-component of velocity field), 2(v-component
%             of velocity field)
%
%    phi   - The field where boundaryc ocnditions should be applied to. 
%   params - The parameter values.
%
% OUTPUT
%
%      phi - The updated field wih ghost values set appropriately.
%
% Copyright 2012, Kenny Erleben, DIKU

if nargin~=3
  error('incorrect number of arguments given')
end
if (type<0) || (type>2)
  error('illegal field type was given')
end

I = params.I;
J = params.J;

i   = 2:I-1;
j   = 2:J-1;

%-------------------------------------------------------------------------
% With out loss of generality let us take a look at a left solid wall. The
% ghost node to the left is given by index i-1 and the fluid cell just
% right of the wall is given by index i.
%
% A pressure value at the wall would be estimated using a CD scheme
%
%  p_wall =  (p_{i-1} + p_i)/2
%
%  D(p_wall,\vec n) = \nabla p_wall^T \vec n \approx   (p_i - p_{i-1})/ dx
%
% Similar formulas hold for u_wall and v_wall.

if type>0
    
  %--- No slip condition means that velocity is zero at walls, (v_wall,u_wall)^T = 0
  phi(i,1) = - phi(i,2);      
  phi(i,J) = - phi(i,J-1);

  phi(1,j) = - phi(  2,j);
  phi(I,j) = - phi(I-1,j);
  
else
  
  %--- Pressure difference across a solid wall is zero, dp_wall/dn = 0
  phi(i,1) =  phi(i,2);
  phi(i,J) =  phi(i,J-1);
  
  phi(1,j) = phi(  2,j);
  phi(I,j) = phi(I-1,j);
  
end

%--- Finally we fix ghost corners -----------------------------------------
phi(1,1) = 0.5 * (phi(2,   1) + phi(1,     2));
phi(1,J) = 0.5 * (phi(2,   J) + phi(1,   J-1));
phi(I,1) = 0.5 * (phi(I,   2) + phi(I-1,   1));
phi(I,J) = 0.5 * (phi(I-1, J) + phi(I,   J-1));

end
