Gb = 100;       
Ib = 10;        
p1 = 0.03;      % Glucose decay rate (1/min)
p2 = 0.02;      % Insulin action decay rate (1/min)
p3 = 0.01;      % Insulin effectiveness rate (1/min)
n = 0.1;        % Insulin clearance rate (1/min)
variability = 0.15; % variability range (±15%)
r=100;
G_init = Gb;       % Initial glucose concentration (mg/dL)
I_init = Ib;        % Initial plasma insulin concentration (mU/L)
X_init = 0;         % Initial insulin action
D_max_big = 300; % Maximum glucose influx for big meals (mg/dL/min)
D_max_small = 120; % Maximum glucose influx for small meals (mg/dL/min)
sigma = 30; % Spread of the spike (minutes)
t_meals_big = [4 * 60, 12 * 60, 18 * 60]; % Big meals: breakfast, lunch, dinner
t_meals_small = [1 * 60, 7 * 60, 10 * 60, 14 * 60, 20 * 60]; % Smaller meals: snacks
D_t = zeros(size(time));
% Larger Meals
for i = 1:length(t_meals_big)
    D_t = D_t + D_max_big * exp(-(time - t_meals_big(i)).^2 / (2 * sigma^2));
end
% Smaller Meals 
for i = 1:length(t_meals_small)
    D_t = D_t + D_max_small * exp(-(time - t_meals_small(i)).^2 / (2 * sigma^2));
end
D_t_interp = @(t) interp1(time, D_t, t, 'linear', 0);
%% Linearized Model 
A = [-p1 -Gb   0;
     0   -p2   p3;
     0    0,   -n]; % State matrix 
C = [1, 0, 0]; % Output matrix 
B = [0; 0; 1]; % Input matrix 
E = [1; 0; 0]; % Disturbance matrix 

% Determining L values
desired_poles = [-2, -2.1, -2.2];
L = design_observer_gain(A, B,C, E, p1, p2, p3, Ib, Gb, 3, desired_poles);

% Determining Kt and Ki values
A_aug = [A, zeros(3,1); -C, 0]; % Add integral state
B_aug = [B; 0];                % Add input for integral state
C_aug = [C, 0];
E_aug = [E;0];
desired_eigs= [-1,-1.1,-1.2,-1.3];
Ka=place(A_aug,B_aug,desired_eigs);
Kt = Ka(1:3); % State feedback gains
Ki = Ka(4);   % Integral feedback gain
disp('State Feedback Gains (Kt):');
disp(Kt);
disp('Integral Gain (Ki):');
disp(Ki);
% Define parameters
params.p1 = p1; params.p2 = p2; params.p3 = p3; 
params.n = n; params.Gb = Gb; params.Gref = Gref;
params.Ib=Ib;

% Define gains
gains.Kt = Kt; 
gains.Ki = Ki;

% Define initial conditions
initial_conditions.G = G_init;
initial_conditions.I = I_init;
initial_conditions.X = X_init;

tspan = [0, 1000]; % Time in minutes
x0 = [120; 0; 11;120;0;11;0]; % Initial conditions: [G, X, I]

[t,x_lin] = ode45(@(t,x_lin) linodefcn(t,x_lin,A,B,C,E,Kt,Ki,r,L,D_t_interp(t)),tspan,x0);

plot_results(t,x_lin)

function dxdt_lin=linodefcn(t,x,A,B,C,E,Kt,Ki,r,L,D_t)
    x_hat=x(4:6);
    z=x(7);
    x=x(1:3);
    y=C*x;
    y_hat=C*x_hat;
    %D=D_t;
    D=0;

    %u=-Kt*x_hat+Ki*z;
    u=0;
    
    dxdt=A*x+B*u+E*D;
    dx_hatdt=A*x_hat+B*u+L*(y-y_hat);
    dzdt=r-y;
    
    dxdt_lin=[dxdt;dx_hatdt;dzdt];
end
%% Nonlinear Dynamics


%% Functions
function L = design_observer_gain(A, B, C, E, p1, p2, p3, Ib, Gb, n, desired_poles)
    % Sanity Checker for the mentally deranged 
    if length(desired_poles) ~= n
        error('Number of desired eigenvalues must match the size of A.');
    end
    % Compute the observer gain L
    L = place(A', C', desired_poles)'; 
    % Verify eigenvalues
    observer_eigenvalues = eig(A - L * C);
    disp('Eigenvalues of A - LC:');
    disp(observer_eigenvalues);
    % Display Results 
    disp('Values of L:');
    disp(L);
end


function plot_results(t, x_lin)
    % Function to plot glucose, insulin action, and insulin concentration results.
    % Parameters:
    % t - Time vector
    % x - State matrix where each column represents a different state variable

    % Plot the first set of states (G, X, I)
    figure;
    subplot(3, 1, 1);
    plot(t/60, x_lin(:,1), 'r', 'LineWidth', 1.5); 
    grid on;
    xlabel('Time (minutes)');
    ylabel('Glucose (G)');
    legend('True G', 'Estimated G');
    title('Glucose Levels');

    subplot(3, 1, 2);
    plot(t/60, x_lin(:,2), 'g', 'LineWidth', 1.5); 
    grid on;
    xlabel('Time (minutes)');
    ylabel('Insulin Action (X)');
    legend('True X', 'Estimated X');
    title('Insulin Action');

    subplot(3, 1, 3);
    plot(t/60, x_lin(:,3), 'b', 'LineWidth', 1.5); 
    grid on;
    xlabel('Time (minutes)');
    ylabel('Insulin (I)');
    legend('True I', 'Estimated I');
    title('Insulin Concentration');

    % Plot the second set of states (G, X, I)
    figure;
    subplot(3, 1, 1);
    plot(t/60, x_lin(:,4), 'r', 'LineWidth', 1.5);
    grid on;
    xlabel('Time (minutes)');
    ylabel('Glucose (G)');
    legend('True G', 'Estimated G');
    title('Glucose Levels');

    subplot(3, 1, 2);
    plot(t/60, x_lin(:,5), 'g', 'LineWidth', 1.5); 
    grid on;
    xlabel('Time (minutes)');
    ylabel('Insulin Action (X)');
    legend('True X', 'Estimated X');
    title('Insulin Action');

    subplot(3, 1, 3);
    plot(t/60, x_lin(:,6), 'b', 'LineWidth', 1.5); 
    grid on;
    xlabel('Time (minutes)');
    ylabel('Insulin (I)');
    legend('True I', 'Estimated I');
    title('Insulin Concentration');
end


