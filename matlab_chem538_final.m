% Hydrogen atom energies and radial wavefunctions
clear;
clc;
tic
% Constants
Z = 1; % Atomic number of hydrogen
a0 = 0.0529; % Bohr radius in nm
eV_to_cm1 = 8065.54; % Conversion factor from eV to cm^-1
n_max = 4; % Maximum principal quantum number
m_e = 9.10938356e-31; % Electron mass (kg)
e = 1.602176634e-19; % Elementary charge (C)
epsilon_0 = 8.8541878128e-12; % Vacuum permittivity (F/m)
h = 6.63e-34; % Planck's constant (Js)
hbar = h/(2*pi); % (Js)
c = 2.9979e10; % (cm/s)
JWN=1/(h*c);

% Parameters
R = 10; % Box size (nm) 
dr = 0.001; % step size for integrals 
r_num = (0:dr:R)'; % Radial grid
rho = 2 * Z * r_num / a0;
N = length(r_num); % Number of radial grid points
n_basis = 500; % Increased basis set size
n_num = 1:n_basis;
num_states = 10; % Total states from 1s to 4f
[rr,nn]=meshgrid(r_num,n_num);

% Initialize storage for energies and wavefunctions
Energy_analytic = zeros(num_states, 3);
WF_anal = zeros(num_states, length(r_num));
Expectation_r = zeros(1,10);
state_index = 1;
% Calculate energies and wavefunctions for all states from 1s to 4f
for n = 1:n_max
    for l = 0:(n-1)
        % Energy in eV
        E_n = -13.6 / n^2*eV_to_cm1;
        % Store energy and quantum numbers
        Energy_analytic(state_index, :) = [E_n, n, l];
        % Radial wavefunction
        rho = 2 * Z * r_num / (n * a0);
        L = laguerreL(n-l-1, 2*l+1, rho); % Associated Laguerre polynomial
        R_r = sqrt((2*Z/(n*a0))^3 * factorial(n-l-1) / (2*n * factorial(n+l))) ...
        .* exp(-rho/2) .* rho.^(l) .* L;
        % Store wavefunction
        WF_anal(state_index, :) = R_r;
        % Calculate expectation value <r> using analytical formula
        Expectation_r(state_index) = (a0 / 2) * (3 * n^2 - l * (l + 1)); % In angstroms
        % Increment state index
        state_index = state_index + 1;
    end
end

% Display energies
fprintf('State\t n\t l\t Energy (cm-1)\n');
for i = 1:num_states
    fprintf('%d\t %d\t %d\t %.4f\n', i, Energy_analytic(i, 2), Energy_analytic(i, 3), Energy_analytic(i, 1));
end

% Plot all wavefunctions
figure;
hold on;
for i = 1:num_states
    plot(r_num, WF_anal(i, :), 'LineWidth', 1.5, 'DisplayName', ...
    sprintf('n = %d, l = %d', Energy_analytic(i, 2), Energy_analytic(i, 3)));
end
xlabel('r (a_0)');
ylabel('R_{n,l}(r)');
title('Radial Wavefunctions for States 1s to 4f');
legend show;
grid on;
hold off;

basis=sqrt(2/R)*sin((pi/R)*nn.*rr);
ddbasis=-((pi*nn)/R).^2.*basis;
l_max=3;
total_xwf=cell(1,l_max+1);
total_ener=cell(1,l_max+1);
for l=0:l_max
    K=-(((hbar^2)/2*m_e).*ddbasis)*basis'*dr*JWN; 
    V_coulomb = (-(e^2)./(4*pi*epsilon_0*rr*1e-9) + ((hbar^2*(l*(l+1)))./(2*m_e*(rr*1e-9).^2)));
    % Avoid the singularity at x = 0
    V_coulomb(rr == 0) = 0; % Set potential at x = 0 to 0
    V=(V_coulomb.*basis)*basis'*dr*JWN;
    H = K + V;
    [xvec, xeval] = eig(H, 'vector');
    % Sorting eigenvalues (outputting sorted eigenvalues and index of it)
    [xevalsort, ind] = sort(xeval);
    total_ener{l+1}=xevalsort;
    % Coefficients
    xvecsort = xvec(:, ind);
    xwf = xvecsort.' * basis;
    total_xwf{l+1}=xwf;
end

% Extract 4 rows from the first cell (l = 0)
rows_l0 = total_xwf{1}(1:4, :);
S=total_ener{1}(1:4,:);
% Extract 3 rows from the second cell (l = 1)
rows_l1 = total_xwf{2}(1:3, :);
P=total_ener{2}(1:3,:);
% Extract 2 rows from the third cell (l = 2)
rows_l2 = total_xwf{3}(1:2, :);
D=total_ener{3}(1:2,:);
% Extract 1 row from the fourth cell (l = 3)
rows_l3 = total_xwf{4}(1, :);
F=total_ener{4}(1,:);
ENERGY=[S;P;D;F];
WF_num=[rows_l0;rows_l1;rows_l2;rows_l3];

% Initialize the energy matrix
ENERGY_MATRIX = zeros(10, 3);
% Assign energies (1st column)
ENERGY_MATRIX(:, 1) = ENERGY;

% Define principal quantum number (n) and angular quantum number (l)
% Rows correspond to states: 1s, 2s, 2p, 3s, 3p, 3d, 4s, 4p, 4d, 4f
n_values = [1; 2; 3; 4; 2; 3; 4; 3; 4; 4];
l_values = [0; 0; 0; 0; 1; 1; 1; 2; 2; 3];

% Create a matrix that combines energy, n, and l values
energy_and_quantum = [ENERGY, n_values, l_values];

% Sort the matrix by n and l (first by n, then by l)
sorted_energy_and_quantum = sortrows(energy_and_quantum, [2, 3]);

% Now extract the sorted values back into ENERGY_MATRIX
ENERGY_MATRIX(:, 1) = sorted_energy_and_quantum(:, 1); % Sorted energies
ENERGY_MATRIX(:, 2) = sorted_energy_and_quantum(:, 2); % n values
ENERGY_MATRIX(:, 3) = sorted_energy_and_quantum(:, 3); % l values

desired_order = [1, 2, 5, 3, 6, 8, 4, 7,9,10];
WF_num = WF_num(desired_order, :);

% Numeric expectation value 
exp_r = (WF_num.*rr(1:10,:)) * WF_num' *dr ;
diag_exp_r = diag(exp_r); 
 
% Figure 1 

figure(50)
plot(n_num(1:6),Energy_analytic(1:6,1),'r-O','MarkerFaceColor','r'); 
hold on 
plot(n_num(1:6),ENERGY_MATRIX(1:6,1),'b-S','MarkerFaceColor','b');
hold off
xlabel('State index'); 
ylabel('Energy (cm^-1)'); 
title('Orbital Energies');
legend('Analytic Energies','Numeric Energies','Location','best')

% State labels for legend
state_labels = {'1s', '2s', '2p', '3s', '3p', '3d'};

% Figure 2
figure(12);
hold on;

% Plot analytical wavefunctions
for i = 1:6
    plot(r_num, WF_anal(i, :), 'LineWidth', 1.5, 'DisplayName', sprintf('Analytical %s', state_labels{i}));
end

% Plot numerical wavefunctions
for i = 1:6
    plot(r_num, WF_num(i, :), '-.', 'LineWidth', 1.5, 'DisplayName', sprintf('Numerical %s', state_labels{i}));
end

hold off;

% Add labels, title, and legend
xlabel('r (nm)'); xlim([0 1]);
ylabel('\psi_{nl}(r)');
title('Radial Wavefunctions for 1s to 3d States');
legend('show', 'Location', 'best');

% Figure 3 

figure(51)
plot(n_num(1:6),Expectation_r(1,1:6),'r-O','MarkerFaceColor','r'); 
hold on 
plot(n_num(1:6),diag_exp_r(1:6,1),'b-S','MarkerFaceColor','b');
hold off
xlabel('State index'); 
ylabel('<r> (nm)'); 
title('Expectation Values of the radius operator');
legend('Analytic Energies','Numeric Energies','Location','northwest');

theta=linspace(0,pi,50);
phi=linspace(0,2*pi,99);
dx=theta(2)-theta(1);
l=2;
m=2;

[Theta,Phi]=meshgrid(theta,phi);

x=cos(Theta);

LEGj=legendre(l,cos(theta));

Plm=LEGj(m+1,:);

Plm_mat=repmat(Plm,length(phi),1);


num=factorial(l-abs(m))*(2*l+1);
denom=factorial(l+abs(m))*4*pi;

normal=sqrt(num/denom);

Ylm=(1i)^(m+abs(m))*normal*Plm_mat.*exp(1i*m*Phi);

Norm=sum(abs(Ylm).^2.*sin(Theta),"all").*dx^2;

[X,Y,Z]=sph2cart(Phi,theta-pi/2,abs(Ylm).^2);
maxR=max(abs(Ylm).^2,[],"all")*1.1;



CO=abs(Ylm).^2;

figure(1)
surf(X,Y,Z,CO);
title("Probability density in any l and m state ")
xlim([-maxR,maxR]);
ylim([-maxR,maxR]);
zlim([-maxR,maxR]);
xlabel('x');
ylabel('y');
zlabel('z');
axis equal; colorbar;

if(abs(m)>0)
    orbital=1/sqrt(2)*(Ylm+conj(Ylm));
else
    orbital=Ylm;
end

Norm2=sum(abs(orbital).^2.*sin(Theta),"all").*dx^2;

maxorb=max(abs(orbital).^2,[],"all")*1.1;

[X1,Y1,Z1]=sph2cart(Phi,theta-pi/2,abs(orbital).^2);
CO2=abs(orbital).^2;

figure(2)
surf(X1,Y1,Z1,CO2);
title("Probability density of orbitals obtained after summation")
xlim([-maxorb,maxorb]);
ylim([-maxorb,maxorb]);
zlim([-maxorb,maxorb]);
xlabel('x');
ylabel('y');
zlabel('z');
axis equal; colorbar;

toc
