%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HOMEWORK #7
% Joshua Julian Damanik (20194701)
% AE551 - Introduction to Optimal Control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc, close all;
addpath('lib');

eps = 0.1;
k = 15;

rho_list = [10, 100, 1000];
delta_rho = 2;

iter = 0;

X_init = [0; 10];

f = @(X) (0.5*(X(1)-1).^2+10*(X(2)-1).^2);
c = @(X) (X(1)-2)^2+2-X(2);

for jj = 1:3
    X = X_init;
    rho = rho_list(jj);
    iter = iter + 1;
    
    La = @(X, lambda) (f(X) + lambda'*c(X)+rho.*c(X)'*c(X));

    lambda = zeros(length(1),1);
    X_data{jj} = X;
    lambda_data{jj} = lambda;
    
    for ii = 2:10000
        LaX = @(X) La(X,lambda);
        X_bar = minimize(X, eps, k, LaX);

        lambda = lambda + 2*rho*c(X_bar);
        fprintf('jj=%d,ii=%d\n',jj,ii);
        
        X_data{jj}(:,ii) = X_bar;
        lambda_data{jj}(:,ii) = lambda;
        
        X = X_bar;
        
        if norm(X_data{jj}(:,ii)-X_data{jj}(:,ii-1))/norm(X_data{jj}(:,ii)) < 1e-10
            break;
        end
    end
end

%% Function Graph

N_grid = 50;
axis_xy = [-2 6 0 15];
x_cont = linspace(axis_xy(1), axis_xy(2), N_grid);
y_cont = linspace(axis_xy(3), axis_xy(4), N_grid);
[X_cont, Y_cont] = meshgrid(x_cont, y_cont);

F_cont = zeros(N_grid);
Cy_cont = zeros(N_grid);

for i=1:N_grid
    for j=1:N_grid
        F_cont(i,j) = f([X_cont(i,j), Y_cont(i,j)]');
    end
end

figure(1);
s = contour(X_cont, Y_cont, F_cont);
colorbar;
hold on;

%% Constraint Graph

C_data = (x_cont-2).^2+2;
plot(x_cont, C_data, 'r--');
axis(axis_xy);

%% Search Path Graph
color = [1, 0.3, 0.3;
         0.3, 0.5, 0.3;
         0.3, 0.3, 1];
for j=1:3
    points = X_data{j};
    for i=1:size(points,2)-1
        F_data = f(points(:,i));
        qlen = [points(:,i+1) - points(:,i)];% f(X_data(:,i+1))-F_data];
        quiver(points(1,i), points(2,i), ...% F_data, ...
                    qlen(1), qlen(2), ... % qlen(3), ...
                    'r', 'AutoScale', 'off', 'LineWidth', 1, ...
                    'MaxHeadSize', min(1 / norm(qlen),1), ...
                    'color', color(j,:));
    end
end
xlabel('x');
ylabel('y');
zlabel('z');
legend('Function', 'Constraint', 'Search path', 'Location', 'SouthEast');

%% Lambda Graph
figure(2); hold on;
for j=1:3
    p=plot(1:length(lambda_data{j}), lambda_data{j}, 'Color', color(j,:));
end
grid on;
xlabel('Iteration');
ylabel('Lambda');
legend('rho=10', 'rho=100', 'rho=1000');
