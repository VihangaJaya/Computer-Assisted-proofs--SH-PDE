 

tic
L = 15;
%alpha = 10+(15)*rand(1,1);
alpha = 15.8; 
m = 19;
a0 = ones(m+1,1);

tspan  = [0,10];
%m = length(a0)-1;
%initialValuesMatrix(i, :) = [L, alpha, a0'];

[t,a_sol] = ode45(@(t, a) function_handle(a,m, L, alpha), tspan, a0');
disp('Solution matrix:');

%solution_vector = a_sol(m+1,:)


figure;
plot(t, a_sol(:, 1:m+1));  % Change the index if needed based on your specific case
title('ODE45 Solution');
xlabel('Time');
ylabel('Solution Values');
%legend('a_0', 'a_1', 'a_2', ..., 'a_m');  % Adjust legends based on the number of elements in a_sol



disp('Solution matrix:');
solution_vector = a_sol(end, :);
disp(solution_vector);

main;


toc
%%
%x = [-1.0532,1.4593,0.7171];
%alpha  = 9.2735;



%function_handle(x, 2, 100, alpha)
%%
p = @(r) mid(Z3).*r.^3 + (Z2).*r.^2 + ((Z1)+(Z0)-1).*r + (Y0);
x = linspace(rmin,rmax,1000);
plot(p(x))


%%
% Data from the LaTeX table
m_values = [20, 29, 40];
time_values = [6.897597, 116.828626,  60.620375];

% Plotting
figure;
plot(m_values, time_values, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
grid on;

% Adding labels and title
xlabel('Value of m');
ylabel('Time (s)');
title('Time vs. Value of m');

% Adding a legend if needed
legend('Data Points');

% Display the plot
disp('Time vs. Value of m Plot:');
disp(' ');
disp('Press any key to continue...');
pause;
