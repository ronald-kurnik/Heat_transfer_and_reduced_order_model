%% 04_heat_transfer_and_reduced_order.m
% 1-D unsteady heat conduction + lumped thermal network
% Relevant to both roles: heat transfer analysis and generating
% lower-order models for system-level simulations.
%
% Demonstrates:
%   - Method of lines for 1-D heat equation
%   - Reduction to a few-node thermal network
%   - Validation of reduced model against full PDE

clear; close all; clc;

%% Physical parameters (steel-like wall)
k_th   = 40;          % W/m-K
rho    = 7800;        % kg/m3
cp     = 500;         % J/kg-K
alpha  = k_th/(rho*cp); % thermal diffusivity

L      = 0.01;        % m  wall thickness
T_init = 300;         % K
T_hot  = 800;         % K  (gas side)
T_cool = 350;         % K  (coolant side)
h_hot  = 500;         % W/m2-K
h_cool = 2000;        % W/m2-K

%% 1-D heat equation via Method of Lines
Nx = 41;
x  = linspace(0,L,Nx)';
dx = x(2)-x(1);
T  = T_init * ones(Nx,1);

% Time
dt = 0.002;
Tfinal = 30;          % s
Nt = round(Tfinal/dt);
t  = (0:Nt)*dt;

T_hist = zeros(Nx,Nt+1);
T_hist(:,1) = T;

for n = 1:Nt
    % Interior
    d2T = (T(1:end-2) - 2*T(2:end-1) + T(3:end))/dx^2;
    T(2:end-1) = T(2:end-1) + dt * alpha * d2T;
    
    % Robin BCs (explicit)
    % left (x=0, hot gas): -k dT/dx = h_hot (T_hot - T(1))
    % right (x=L, coolant): -k dT/dx = h_cool (T(end) - T_cool)
    T(1)   = (k_th/dx * T(2)   + h_hot  * T_hot)  / (k_th/dx + h_hot);
    T(end) = (k_th/dx * T(end-1) + h_cool * T_cool) / (k_th/dx + h_cool);
    
    T_hist(:,n+1) = T;
end

%% Reduced-order thermal network (3 nodes)
% Node 1: gas-side surface, Node 2: bulk wall, Node 3: coolant-side surface
% Capacities and conductances derived from geometry

C1 = rho*cp*(L/6)*1;     % approximate
C2 = rho*cp*(L/3)*1;
C3 = rho*cp*(L/6)*1;
G12 = k_th / (L/2);      % W/K  (per unit area)
G23 = k_th / (L/2);
Gh  = h_hot;
Gc  = h_cool;

% State: [T1; T2; T3]
A_rom = [ -(Gh+G12)/C1 ,  G12/C1 ,  0 ;
           G12/C2 , -(G12+G23)/C2 , G23/C2 ;
           0 , G23/C3 , -(G23+Gc)/C3 ];
B_rom = [ Gh/C1 , 0 ; 0 , 0 ; 0 , Gc/C3 ];
u_rom = [T_hot; T_cool];

% Integrate ROM with same dt
T_rom = T_init * ones(3,1);
T_rom_hist = zeros(3,Nt+1);
T_rom_hist(:,1) = T_rom;

for n = 1:Nt
    T_rom = T_rom + dt * (A_rom*T_rom + B_rom*u_rom);
    T_rom_hist(:,n+1) = T_rom;
end

%% Comparison plots
figure('Name','Heat Transfer PDE vs ROM','Position',[50 50 1000 600]);

subplot(2,2,1);
imagesc(t, x*1000, T_hist);
colorbar; set(gca,'YDir','normal');
xlabel('t [s]'); ylabel('x [mm]'); title('Full 1-D temperature field [K]');

subplot(2,2,2);
plot(t, T_hist(1,:), 'b', t, T_hist(round(Nx/2),:), 'r', t, T_hist(end,:), 'g', 'LineWidth',1.3);
hold on;
plot(t, T_rom_hist(1,:), 'b--', t, T_rom_hist(2,:), 'r--', t, T_rom_hist(3,:), 'g--', 'LineWidth',1.5);
legend('PDE surface hot','PDE mid','PDE surface cool', ...
       'ROM T1','ROM T2','ROM T3','Location','best');
xlabel('t [s]'); ylabel('T [K]'); title('Surface & mid-plane temperatures');
grid on;

subplot(2,2,3);
plot(x*1000, T_hist(:,end), 'k-', 'LineWidth',1.5); hold on;
plot([0 L/2 L]*1000, T_rom_hist(:,end), 'ro-', 'LineWidth',1.5, 'MarkerSize',8);
xlabel('x [mm]'); ylabel('T [K]'); title('Final temperature profile');
legend('PDE','ROM nodes'); grid on;

subplot(2,2,4);
% Heat flux comparison
q_hot_pde = h_hot * (T_hot - T_hist(1,:));
q_hot_rom = h_hot * (T_hot - T_rom_hist(1,:));
plot(t, q_hot_pde, 'b', t, q_hot_rom, 'r--', 'LineWidth',1.3);
xlabel('t [s]'); ylabel('q [W/m^2]'); title('Gas-side heat flux');
legend('PDE','ROM'); grid on;

fprintf('\n=== Heat transfer & model reduction demo complete ===\n');
fprintf('Talking points for interview:\n');
fprintf(' - Full CFD / FEA gives spatial detail; system models need 0-D/1-D ROMs\n');
fprintf(' - Validate reduced models against high-fidelity and against test data\n');
fprintf(' - Thermal networks are standard for real-time and controls co-simulation\n');
