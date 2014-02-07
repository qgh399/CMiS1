function [ phi ] = semi_lagrangian( iters )

[X,Y] = meshgrid(1:101, 1:101);
[N,M] = size(X);
phi = peaks(N);
phi_orig = phi;
err = zeros(1,iters);
dt = 0.01;
    
u_x = Y - 51;
u_y = -X + 51;

for i=1:iters
    x_old = X - dt*u_x;
    y_old = Y - dt*u_y;
    
    x_old(x_old>N) = N; 
    x_old(x_old<1) = 1;

    y_old(y_old>M) = M;
    y_old(y_old<1) = 1;
    
    phi = interp2(phi, x_old, y_old);
    
    RMS_err = norm(phi_orig - phi);
    err(i) = RMS_err;
    
    figure(1);
    subplot(1,2,1);
    surf(X,Y,phi)
    drawnow
    
    subplot(1,2,2);
    plot(1:i, err(1:i))
    grid on;
    drawnow
end

figure(2);
subplot(1,2,1);
surf(X,Y,phi_orig)
subplot(1,2,2);
surf(X,Y,phi)


end

