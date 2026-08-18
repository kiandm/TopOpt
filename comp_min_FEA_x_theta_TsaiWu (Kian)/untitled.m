% h_fd = 1e-6;
%% dc/dx
dcdx = importdata('finitediff_dcdx.txt', ' ');
elements1 = dcdx.data(:,1);
analytic1 = dcdx.data(:,2);
numeric1 = dcdx.data(:,3);
abs_err1 = dcdx.data(:,4);
rel_err1 = dcdx.data(:,5) * 100;

% Plot relative error per element
figure(1);
tiledlayout(2,1);
nexttile
plot(elements1, rel_err1, '-x');
xlabel('Elements');
ylabel('Relative Error (%)');
xlim([0,1600]);
title('Relative Error for Each Element (dc/dx)');
grid on;

% Plot absolute error per element
nexttile;
plot(elements1, abs_err1, '-x');
xlabel('Elements');
ylabel('Absolute Error');
xlim([0,1600]);
title('Absolute Error for Each Element (dc/dx)');
grid on;
drawnow;

%% dc/dtheta
dcdtheta = importdata('finitediff_dcdtheta.txt', ' ');
elements2 = dcdtheta.data(:,1);
analytic2 = dcdtheta.data(:,2);
numeric2 = dcdtheta.data(:,3);
abs_err2 = dcdtheta.data(:,4);
rel_err2 = dcdtheta.data(:,5) * 100;

% Plot relative error per element
figure(2);
tiledlayout(2,1);
nexttile;
plot(elements2, rel_err2, '-x');
xlabel('Elements');
ylabel('Relative Error (%)');
xlim([0,1600]);
title('Relative Error for Each Element (dc/dtheta)');
grid on;

% Plot absolute error per element
nexttile;
plot(elements2, abs_err2, '-x');
xlabel('Elements');
ylabel('Absolute Error');
xlim([0,1600]);
title('Absolute Error for Each Element (dc/dtheta)');
grid on;
drawnow;

%% dg/dx
dgdx = importdata('finitediff_dgdx.txt', ' ');
elements3 = dgdx.data(:,1);
analytic3 = dgdx.data(:,2);
numeric3 = dgdx.data(:,3);
abs_err3 = dgdx.data(:,4);
rel_err3 = dgdx.data(:,5) * 100;

% Plot relative error per element
figure(3);
tiledlayout(2,1);
nexttile;
plot(elements3, rel_err3, '-x');
xlabel('Elements');
ylabel('Relative Error (%)');
xlim([0,1600]);
title('Relative Error for Each Element (dg/dx)');
grid on;

% Plot absolute error per element
nexttile;
plot(elements3, abs_err3, '-x');
xlabel('Elements');
ylabel('Absolute Error');
xlim([0,1600]);
title('Absolute Error for Each Element (dg/dx)');
grid on;
drawnow;

%% dg/dtheta
dgdtheta = importdata('finitediff_dgdtheta.txt', ' ');
elements4 = dgdtheta.data(:,1);
analytic4 = dgdtheta.data(:,2);
numeric4 = dgdtheta.data(:,3);
abs_err4 = dgdtheta.data(:,4);
rel_err4 = dgdtheta.data(:,5) * 100;

% Plot relative error per element
figure(4);
tiledlayout(2,1);
nexttile;
plot(elements4, rel_err4, '-x');
xlabel('Elements');
ylabel('Relative Error (%)');
xlim([0,1600]);
title('Relative Error for Each Element (dg/dtheta)');
grid on;

% Plot absolute error per element
nexttile;
plot(elements4, abs_err4, '-x');
xlabel('Elements');
ylabel('Absolute Error');
xlim([0,1600]);
title('Absolute Error for Each Element (dg/dtheta)');
grid on;
drawnow;

%%
figure(5); clf;
patch('Faces', conn', ...
      'Vertices', coords', ...
      'FaceVertexCData', rel_err1, ...
      'FaceColor', 'flat', ...
      'EdgeColor', 'none');
set(gcf, 'Color', 'white')
axis equal off;
colorbar;
%clim([0 100]);   % 1 = failure limit
title('Relative error for dc/dx');
drawnow;

figure(6); clf;
patch('Faces', conn', ...
      'Vertices', coords', ...
      'FaceVertexCData', rel_err2, ...
      'FaceColor', 'flat', ...
      'EdgeColor', 'none');
set(gcf, 'Color', 'white')
axis equal off;
colorbar;
%clim([0 100]);   % 1 = failure limit
title('Relative error for dc/dtheta');
drawnow;

figure(7); clf;
patch('Faces', conn', ...
      'Vertices', coords', ...
      'FaceVertexCData', rel_err3, ...
      'FaceColor', 'flat', ...
      'EdgeColor', 'none');
set(gcf, 'Color', 'white')
axis equal off;
colorbar;
%clim([0 100]);   % 1 = failure limit
title('Relative error for dg/dx');
drawnow;

figure(8); clf;
patch('Faces', conn', ...
      'Vertices', coords', ...
      'FaceVertexCData', rel_err4, ...
      'FaceColor', 'flat', ...
      'EdgeColor', 'none');
set(gcf, 'Color', 'white')
axis equal off;
colorbar;
%clim([0 100]);   % 1 = failure limit
title('Relative error for dg/dtheta');
drawnow;
%%
figure(9); clf;
patch('Faces', conn', ...
      'Vertices', coords', ...
      'FaceVertexCData', abs_err1, ...
      'FaceColor', 'flat', ...
      'EdgeColor', 'none');
set(gcf, 'Color', 'white')
axis equal off;
colorbar;
%clim([0 100]);   % 1 = failure limit
title('Absolute error for dc/dx');
drawnow;

figure(10); clf;
patch('Faces', conn', ...
      'Vertices', coords', ...
      'FaceVertexCData', abs_err2, ...
      'FaceColor', 'flat', ...
      'EdgeColor', 'none');
set(gcf, 'Color', 'white')
axis equal off;
colorbar;
%clim([0 100]);   % 1 = failure limit
title('Absolute error for dc/dtheta');
drawnow;

figure(11); clf;
patch('Faces', conn', ...
      'Vertices', coords', ...
      'FaceVertexCData', abs_err3, ...
      'FaceColor', 'flat', ...
      'EdgeColor', 'none');
set(gcf, 'Color', 'white')
axis equal off;
colorbar;
%clim([0 100]);   % 1 = failure limit
title('Absolute error for dg/dx');
drawnow;

figure(12); clf;
patch('Faces', conn', ...
      'Vertices', coords', ...
      'FaceVertexCData', abs_err4, ...
      'FaceColor', 'flat', ...
      'EdgeColor', 'none');
set(gcf, 'Color', 'white')
axis equal off;
colorbar;
%clim([0 100]);   % 1 = failure limit
title('Absolute error for dg/dtheta');
drawnow;