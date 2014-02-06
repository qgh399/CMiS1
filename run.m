% Solves the incompressible Navier Stokes equations in 2D 
% Copyright 2012, Kenny Erleben, DIKU.
clear all;
close all;
clc;

params = create_params(1,1,100,100);

u      = zeros(params.I, params.J);
v      = zeros(params.I, params.J);
smoke  = zeros(params.I, params.J);

T_wanted = params.T;
frame    = 1;

iter = 1;
diff_err = zeros(1, params.max_iter - 1);
press_err = zeros(1, params.max_iter - 1);

while T_wanted>0
  
  dt_wanted = 1/ params.fps;
  
  min_dist = min(params.dx, params.dy );
  
  grid_speed = min_dist / dt_wanted;  
    
  while dt_wanted > 0
  
    %--- Apply a crude CFL condition as a rule of thumb to avoid
    %--- stability issues. For details read: 
    %
    % @book{osher2003level,
    %  title={Level set methods and dynamic implicit surfaces},
    %  author={Osher, S. and Fedkiw, R.P.},
    %  isbn={9780387954820},
    %  lccn={70022436},
    %  series={Applied mathematical sciences},
    %  url={http://books.google.dk/books?id=SQQI2vqWR7gC},
    %  year={2003},
    %  publisher={Springer}
    %  }
    %
    %  Around equation 3.10 and 3.11
    %
    speed = sqrt( u.^2 + v.^2 );
    max_speed = max( [ speed(speed>0);  grid_speed ]  );
    cfl_dt = min( dt_wanted, 0.9*(min_dist/max_speed) );
    
    params.dt = cfl_dt;
    
    [u, v, smoke, diff_err_u, diff_err_v, press_err1, press_err2] = ...
        compute_fractional_step(u,v,smoke,params);
    
    diff_err = diff_err + (diff_err_u + diff_err_v)/2;
    press_err = press_err + (press_err1 + press_err2)/2;
    iter = iter + 1;

    dt_wanted = dt_wanted - cfl_dt;
  end
  
  figure(100);
  clf;
  draw_fluid(u,v,smoke,params);
  axis equal tight;
  colormap(gray);
  title(['fame no.'  num2str(frame) ]);

  T_wanted = T_wanted - 1/params.fps;
  frame = frame + 1;
end
