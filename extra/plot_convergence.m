% The implicit diffusion convergence plot
figure(1)
diff_data = load('../../data/implicit_diffusion_convergence_iter983.mat');
diff_val = diff_data.avg_diff_err;
semilogy(1:8, diff_val(1:8));
xlabel('# Iterations', 'FontSize', 18);
ylabel('\epsilon', 'FontSize', 18);
grid on

figure(2)
press_data = load('../../data/pressure_projection_convergence_iter983.mat');
press_val = press_data.avg_press_err;
semilogy(1:length(press_val), press_val);
xlabel('# Iterations', 'FontSize', 18);
ylabel('\epsilon', 'FontSize', 18);
grid on