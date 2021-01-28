clear all;
clc;

format long;

filename1 = 'Workshop3_Speed_heading.csv';
speed_heading = csvread(filename1);

filename2 = 'Workshop3_GNSS_Pos_Vel_NED.csv';
GNSS_pos_vel = csvread(filename2);

%initialize the constants
Define_Constants;

%initalize a state
states = zeros(351,5);

%initial values
%latitude
L_k_1 = deg_to_rad * 50.4249580; %deg
%longitude
lambda_k_1 = deg_to_rad * -3.5957974; %deg

%initial velocities
psi_0 = deg_to_rad * speed_heading(1,3);
V_k_1 = speed_heading(1,2);
V_N_k_1 = V_k_1 * cos(psi_0);
V_E_k_1 = V_k_1 * sin(psi_0);


 states(1,1) = speed_heading(1,1);
 states(1,2) = 50.4249580;
 states(1,3) = -3.5957974;
 states(1,4) = V_N_k_1;
 states(1,5) = V_E_k_1;

for r = 1:350
    %average velocity between epochs k−1 and k
    psi_k_1 = deg_to_rad * speed_heading(r,3);
    psi_k = deg_to_rad * speed_heading(r+1,3);
    v_k = speed_heading(r+1,2);
    
    %initialize a matrix for the velocities
    v_avg_k = zeros(2,350);
    
    v_avg_k(:,r) = 1/2 * [cos(psi_k) + cos(psi_k_1); sin(psi_k) + sin(psi_k_1)] * v_k;

    V_avg_N_k = v_avg_k(1,r);
    V_avg_E_k = v_avg_k(2,r);

    %latitude, lonbgitude calculation
    %height 
    h = 37.4;
    
    %time
    t_k_1 = speed_heading(r,1);
    t_k = speed_heading(r+1,1);
    
    [R_N,R_E]= Radii_of_curvature(L_k_1);
    %latitude
    L_k = L_k_1 + V_avg_N_k *(t_k - t_k_1)/(R_N + h);
    lambda_k = lambda_k_1 + V_avg_E_k *(t_k - t_k_1)/((R_E + h)*cos(L_k));
    
    %compute the damped instantaneous DR velocity
    
    V_N_k = 1.7 * V_avg_N_k - 0.7 * V_N_k_1;
    V_E_k = 1.7 * V_avg_E_k - 0.7 * V_E_k_1;
    
    %update
    L_k_1 = L_k;
    lambda_k_1 = lambda_k;
    V_N_k_1 = V_N_k;
    V_E_k_1 = V_E_k;
    states(r+1,1) = speed_heading(r+1,1);
    states(r+1,2) = rad_to_deg * L_k;
    states(r+1,3) = rad_to_deg * lambda_k;
    states(r+1,4) = V_N_k;
    states(r+1,5) = V_E_k;
end

%clearvars -except states
