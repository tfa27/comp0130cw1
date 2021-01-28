clear all;
clc;

format long;

%read csv files
filename1 = 'Workshop2_Pseudo_ranges.csv';
Pseudo_ranges = csvread(filename1);


filename2 = 'Workshop2_Pseudo_range_rates.csv';
Pseudo_range_rates = csvread(filename2);

filename3 = 'Workshop2_GNSS_Pos_ECEF.csv';
GNSS_pos = csvread(filename3);

%initialize the constants
Define_Constants;

%Task 2A 

%a -----------------------------------------------
%Initialise the Kalman filter state vector estimate
[x_est,P_matrix] = Initialise_GNSS_KF

x_k_1_plus = x_est;
P_k_1_plus = P_matrix;

%b ----------------------------------------------
%Compute the transition matrix 

tau_s = 1; %s
phi_k_1 = eye(8);


for i = 1:3
    phi_k_1(i,3+i) = tau_s *1;
end
phi_k_1(7,8) = tau_s *1;

%c ----------------------------------------
%Compute the system noise covariance matrix

S_a_e = 5 %m^2s^-3
S_c_phi_a = 0.01;
S_cf_a = 0.04;
Q_k_1 = eye(8);

for i = 1:3
    Q_k_1(i,i) = 1/3 * S_a_e * tau_s^3;
    Q_k_1(i,3+i) = 1/2 * S_a_e * tau_s^2;
    Q_k_1(3+i,i) = 1/2 * S_a_e * tau_s^2;
    Q_k_1(3+i,3+i) = S_a_e * tau_s;
end

Q_k_1(7,7) = S_c_phi_a * tau_s + 1/3 * S_cf_a * tau_s^3;
Q_k_1(7,8) = 1/2 * S_cf_a * tau_s^2;
Q_k_1(8,7) = 1/2 * S_cf_a * tau_s^2;
Q_k_1(8,8) = S_cf_a * tau_s;


%d ---------------------------------------------------------
%Use the transition matrix to propagate the state estimates

x_k_minus = phi_k_1 * x_k_1_plus

%e -----------------------------------------------------------
% propagate the error covariance matrix

P_k_minus = phi_k_1 * P_k_1_plus * transpose(phi_k_1) + Q_k_1

%f ----------------------------------------------------------
% Predict the ranges from the approximate user position
%Earth rotation rate
omega_ie = 7.292115 * 10^-5; %rad/s
%speed of light
c = 299792458; %m/s

%initilaize a vector for user - sat ranges
r_aj_minus = zeros(10,1);

%Compute the Cartesian ECEF positions of the satellites at time 0
time = 0;
num_sat = Pseudo_ranges(1,2:11);
%calculate the statelite position and velocity wrt Earth center
% initialize a matrix for psoitions and velocities
r_ej = zeros(3,10);
v_ej = zeros(3,10);

for r = 1:10
    [r_ej(:,r),v_ej(:,r)] = Satellite_position_and_velocity(time,num_sat(r));  
end

%predicted user pos
r_ea_minus = x_k_minus(1:3);
v_ea_minus = x_k_minus(4:6);

% line-of-sight unit vector
u_aj = zeros(3,10);

%Predict the range rates from the approximate user position to each satellite 
r_aj_dot_minus = zeros(10,1);

for r = 1:10
    Ce = eye(3);
    r_aj_minus(r) = sqrt(transpose(Ce * r_ej(:,r) - r_ea_minus)*(Ce * r_ej(:,r) - r_ea_minus)); 
    temp = omega_ie * r_aj_minus(r)/c;
    Ce = [1 temp 0; -temp 1 0; 0 0 1];
    r_aj_minus(r) = sqrt(transpose(Ce * r_ej(:,r) - r_ea_minus)*(Ce * r_ej(:,r) - r_ea_minus)); 
    %g ------------------------------------------------------------------------------------
    u_aj(:,r) = (Ce * r_ej(:,r) - r_ea_minus)/ r_aj_minus(r)
    %h ------------------------------------------------------------------
    r_aj_dot_minus(r) = transpose(u_aj(:,r))* (Ce *(v_ej(:,r) + Omega_ie * r_ej(:,r)) - (v_ea_minus + Omega_ie * r_ea_minus))
end

%i ----------------------------------------------------------------------
% Compute the measurement matrix
H_k = zeros(20,8);

for r = 1:10
    H_k(r,1:3) = - transpose(u_aj(:,r));
    H_k(r,7) = 1;
    H_k(10 + r,8) = 1;
    H_k(10 + r,4:6) = - transpose(u_aj(:,r));
   
end

%j -----------------------------------------------------------------------
%  Compute the measurement noise covariance matrix 
R_k = zeros(20,20);
sigma_p = 10;
sigma_r = 0.05;

for r = 1:10
    R_k(r,r) = sigma_p^2;
    R_k(10 + r,10 + r) = sigma_r^2;
end

%k -----------------------------------------------------------------------
%Compute the Kalman gain matrix 
K_k = P_k_minus * transpose(H_k)* inv(H_k * P_k_minus * transpose(H_k) + R_k)


%i ------------------------------------------------------------------------
%Formulate the measurement innovation vector
delta_z_minus = zeros(20,1);

%reciever clock offset estimate
delta_rho_c = x_k_minus(7)
%receiver clock drift estimate
delta_rho_dot_c = x_k_minus(8)

%pseudo range
rho_a = Pseudo_ranges(2,2:11)

rho_dot_a = Pseudo_range_rates(2,2:11)

for r = 1:10
    delta_z_minus(r) = rho_a(r) - r_aj_minus(r) - delta_rho_c;
    delta_z_minus(10+r) = rho_dot_a(r) - r_aj_dot_minus(r) - delta_rho_dot_c;
end
delta_z_minus
%m -----------------------------------------------------------------------
% Update the state estimates using

x_k_plus = x_k_minus + K_k * delta_z_minus

%n ------------------------------------------------------------- 
%Update the error covariance matrix
I = eye(8);
P_k_plus = (I - K_k * H_k)*P_k_minus

%o ---------------------------------------------------------------
%Convert this Cartesian ECEF position solution to latitude, longitude and height

[L_b_f,lambda_b_f,h_b_f,v_eb_n_f] = pv_ECEF_to_NED(x_k_plus(1:3),x_k_plus(4:6))

latitude_f = rad_to_deg * L_b_f
longitude_f = rad_to_deg * lambda_b_f




%TASK 2B -----------------------------------------------------------

%Initialise the Kalman filter state vector estimate
[x_est,P_matrix] = Initialise_GNSS_KF

x_k_1_plus = x_est;
P_k_1_plus = P_matrix;

for j = 2:182
    disp('time step:');
    time = j - 2
    tau_s = 1; %s
    phi_k_1 = eye(8);


    for i = 1:3
        phi_k_1(i,3+i) = tau_s *1;
    end
    phi_k_1(7,8) = tau_s *1;

    %c ----------------------------------------
    %Compute the system noise covariance matrix

    S_a_e = 5; %m^2s^-3
    S_c_phi_a = 0.01;
    S_cf_a = 0.04;
    Q_k_1 = eye(8);

    for i = 1:3
        Q_k_1(i,i) = 1/3 * S_a_e * tau_s^3;
        Q_k_1(i,3+i) = 1/2 * S_a_e * tau_s^2;
        Q_k_1(3+i,i) = 1/2 * S_a_e * tau_s^2;
        Q_k_1(3+i,3+i) = S_a_e * tau_s;
    end

    Q_k_1(7,7) = S_c_phi_a * tau_s + 1/3 * S_cf_a * tau_s^3;
    Q_k_1(7,8) = 1/2 * S_cf_a * tau_s^2;
    Q_k_1(8,7) = 1/2 * S_cf_a * tau_s^2;
    Q_k_1(8,8) = S_cf_a * tau_s;


    %d ---------------------------------------------------------
    %Use the transition matrix to propagate the state estimates

    x_k_minus = phi_k_1 * x_k_1_plus;

    %e -----------------------------------------------------------
    % propagate the error covariance matrix

    P_k_minus = phi_k_1 * P_k_1_plus * transpose(phi_k_1) + Q_k_1;

    %f ----------------------------------------------------------
    % Predict the ranges from the approximate user position
    %Earth rotation rate
    omega_ie = 7.292115 * 10^-5; %rad/s
    %speed of light
    c = 299792458; %m/s

    %initilaize a vector for user - sat ranges
    r_aj_minus = zeros(10,1);

    %Compute the Cartesian ECEF positions of the satellites at time 0
    
    num_sat = Pseudo_ranges(1,2:11);
    %calculate the statelite position and velocity wrt Earth center
    % initialize a matrix for psoitions and velocities
    r_ej = zeros(3,10);
    v_ej = zeros(3,10);

    for r = 1:10
        [r_ej(:,r),v_ej(:,r)] = Satellite_position_and_velocity(time,num_sat(r));  
    end

    %predicted user pos
    r_ea_minus = x_k_minus(1:3);
    v_ea_minus = x_k_minus(4:6);

    % line-of-sight unit vector
    u_aj = zeros(3,10);

    %Predict the range rates from the approximate user position to each satellite 
    r_aj_dot_minus = zeros(10,1);

    for r = 1:10
        Ce = eye(3);
        r_aj_minus(r) = sqrt(transpose(Ce * r_ej(:,r) - r_ea_minus)*(Ce * r_ej(:,r) - r_ea_minus)); 
        temp = omega_ie * r_aj_minus(r)/c;
        Ce = [1 temp 0; -temp 1 0; 0 0 1];
        r_aj_minus(r) = sqrt(transpose(Ce * r_ej(:,r) - r_ea_minus)*(Ce * r_ej(:,r) - r_ea_minus)); 
        %g ------------------------------------------------------------------------------------
        u_aj(:,r) = (Ce * r_ej(:,r) - r_ea_minus)/ r_aj_minus(r);
        %h ------------------------------------------------------------------
        r_aj_dot_minus(r) = transpose(u_aj(:,r))* (Ce *(v_ej(:,r) + Omega_ie * r_ej(:,r)) - (v_ea_minus + Omega_ie * r_ea_minus));
    end

    %i ----------------------------------------------------------------------
    % Compute the measurement matrix
    H_k = zeros(20,8);

    for r = 1:10
        H_k(r,1:3) = - transpose(u_aj(:,r));
        H_k(r,7) = 1;
        H_k(10 + r,8) = 1;
        H_k(10 + r,4:6) = - transpose(u_aj(:,r));

    end

    %j -----------------------------------------------------------------------
    %  Compute the measurement noise covariance matrix 
    R_k = zeros(20,20);
    sigma_p = 10;
    sigma_r = 0.05;

    for r = 1:10
        R_k(r,r) = sigma_p^2;
        R_k(10 + r,10 + r) = sigma_r^2;
    end

    %k -----------------------------------------------------------------------
    %Compute the Kalman gain matrix 
    K_k = P_k_minus * transpose(H_k)* inv(H_k * P_k_minus * transpose(H_k) + R_k);


    %i ------------------------------------------------------------------------
    %Formulate the measurement innovation vector
    delta_z_minus = zeros(20,1);

    %reciever clock offset estimate
    delta_rho_c = x_k_minus(7);
    %receiver clock drift estimate
    delta_rho_dot_c = x_k_minus(8);

    %pseudo range
    rho_a = Pseudo_ranges(j,2:11);

    rho_dot_a = Pseudo_range_rates(j,2:11);

    for r = 1:10
        delta_z_minus(r) = rho_a(r) - r_aj_minus(r) - delta_rho_c;
        delta_z_minus(10+r) = rho_dot_a(r) - r_aj_dot_minus(r) - delta_rho_dot_c;
    end
    delta_z_minus;
    %m -----------------------------------------------------------------------
    % Update the state estimates using

    x_k_plus = x_k_minus + K_k * delta_z_minus;

    %n ------------------------------------------------------------- 
    %Update the error covariance matrix
    I = eye(8);
    P_k_plus = (I - K_k * H_k)*P_k_minus;

    %o ---------------------------------------------------------------
    %Convert this Cartesian ECEF position solution to latitude, longitude and height

    [L_b_f,lambda_b_f,h_b_f,v_eb_n_f] = pv_ECEF_to_NED(x_k_plus(1:3),x_k_plus(4:6));
    
    latitude = rad_to_deg * L_b_f
    longitude = rad_to_deg * lambda_b_f
    height = h_b_f
    
    v_North = v_eb_n_f(1)
    v_East = v_eb_n_f(2)
    v_Down = v_eb_n_f(3)

    
    
    %update
    x_k_1_plus = x_k_plus;
    P_k_1_plus = P_k_plus;
    
end