function [ phi ] = semi_lagrangian( iters )

[X,Y] = meshgrid(1:101, 1:101);
[N,M] = size(X);
phi = peaks(N);
dt = 0.01;

for i=1:iters
    u_x = Y - 51;
    u_y = -X + 51;
    x_old = X - dt*u_x;
    y_old = Y - dt*u_y;
    
    x_old(x_old>N) = N; 
    x_old(x_old<1) = 1;

    y_old(y_old>M) = M;
    y_old(y_old<1) = 1;
    
    phi = interp2(phi, x_old, y_old);
    
    surf(X,Y,phi)
    drawnow
end

end

