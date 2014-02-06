function params = create_params(width, height, I, J)
% CREATE_PARAMS -- This function sets up a structure containg all essential
% parameter values for the fluid simulation. This makes it easier to pass
% parameter values around to the different functions that are called during
% simulation.
%
% OUTPUT:
%
%  width  - The physical width of the domain
%  height - The physical height of the domain
%  I      - The number of grid nodes along the x-axis.
%  J      - The number of grid nodes along the y-axis.
%  params - A struct of parameter values.
%
% Copyright 2012, Kenny Erlbeen, DIKU.


fps        = 30;           % Frames per second
T          = 1;            % Total simulated time
dt         = 0.01;         % Time step size
viscosity  = 0.0000001;    % Viscosity coefficient
rel_tol    = 0.01;         % Relative tolerance for linear system solver
abs_tol    = 0.01;         % Absolute tolerance for linear system solver
max_iter   = 250;           % Maximum number of iterations of linear system solver
drop       = 0.000625;     % Buoyancy density drop coefficient
lift       = 0.25;         % Buoyancy temparature lift coefficient
dx         = width/(I-1);  % Grid spaing between two grid nodes.
dy         = height/(J-1); % Grid spaing between two grid nodes.
I          = I+2;          % Add ghost nodes to domain.
J          = J+2;          % Add ghost nodes to domain.
smoke_rate = 1000;         % The rate of smoke added to the simulation
vorticity  = 75;           % The amount of small scale that should be added by the vorticity forces.

smoke_source_i = floor(I/2)-floor(I/10):floor(I/2)+floor(I/10);
smoke_source_j = floor(J/3)-floor(J/10):floor(J/3)+floor(J/10);


%--- Create some auxiliary precomputed stuff to make rendering more fast --
[subJ, subI] = meshgrid(2:I-1, 2:J-1);
idx = sub2ind( [I, J], subI(:), subJ(:));

span_x = -dx:dx:width+dx;
span_y = -dy:dy:height+dy;
[full_yc, full_xc] = meshgrid(span_x, span_y);

xc = full_xc(idx);
yc = full_yc(idx);

x1 = xc - dx/2;
x2 = xc + dx/2;
x3 = xc + dx/2;
x4 = xc - dx/2;

y1 = yc - dy/2;
y2 = yc - dy/2;
y3 = yc + dy/2;
y4 = yc + dy/2;

patch_x = [x1, x2, x3, x4];
patch_y = [y1, y2, y3, y4];
%--------------------------------------------------------------------------

params = struct(...
  'T',T,...
  'fps',fps,...
  'dt',dt,...
  'viscosity',viscosity,...
  'rel_tol',rel_tol,...
  'abs_tol',abs_tol,...
  'max_iter',max_iter,...
  'drop',drop,...
  'lift',lift,...
  'I',I,...
  'J',J,...
  'dx',dx,...
  'dy',dy,...
  'smoke_rate',smoke_rate,...
  'smoke_source_i',smoke_source_i,...
  'smoke_source_j',smoke_source_j,...
  'vorticity', vorticity,...
  'width',width,...
  'height',height,...
  'idx',idx,...
  'patch_x',patch_x,...
  'patch_y',patch_y...
  );
end