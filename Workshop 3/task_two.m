alternate_script; % <----- TASK 1
GNSS_stuff = csvread('Workshop3_GNSS_Pos_Vel_NED.csv');
deg2rad = pi / 180;
rad2deg = 180 / pi;

% col 1 time in seconds
% col 2 geodetic latitude in degrees
% col 3 geodetic longitude in degrees
% col 4 geodetic height in metres
% col 5 to 7 earth referenced velocity in metres per second, resolved along north east and down respectively 

% final step of coursework!!!

%Part 1
delta_v_n = 0;
delta_v_e = 0;
delta_L = 0;
delta_lambda = 0;
x_init = [delta_v_n; delta_v_e; delta_L; delta_lambda];

% Part 2
sigma_sub_v = 0.1;
sigma_sub_r = 10;

h_k_min_1 = GNSS_stuff(1, 4);
L_k_min_1 = GNSS_stuff(1,2) * deg2rad;
[R_N,R_E]= Radii_of_curvature(GNSS_stuff(1,2) * deg2rad);
P_sub_zero_super_plus = eye(4);
P_sub_zero_super_plus(1:2, 1:2) = P_sub_zero_super_plus(1:2, 1:2) *  sigma_sub_v^2;
P_sub_zero_super_plus(3,3) = sigma_sub_r^2 / ((R_N + h_k_min_1)^2);
P_sub_zero_super_plus(4,4) = sigma_sub_r^2 / ((R_E + h_k_min_1)^2 * (cos(L_k_min_1))^2);
S_dr = 0.2;
x_k_min_1_carat_plus = x_init;
res = [];
time = 0;
res = [res; time, transpose(x_init)];


for i = 2:351
    
    [R_N,R_E]= Radii_of_curvature(GNSS_stuff(i,2) * deg2rad); % maybe minus 1 ???
    
    
    % EQUATION 6
    tau_s = 0.5;
    phi_k_min_1 = eye(4);

    phi_k_min_1(3,1) = tau_s / (R_N + h_k_min_1);
    phi_k_min_1(4,2) = tau_s / ((R_E + h_k_min_1) * cos(L_k_min_1));
    
    
    % EQUATION 7
    qew_k_min_1 = eye(4);
    qew_k_min_1(1:2, 1:2) = qew_k_min_1(1:2, 1:2) * (S_dr * tau_s);
    qew_k_min_1(3,3) = (1/3)  * ( (S_dr * tau_s^3) / (R_N + h_k_min_1) );
    qew_k_min_1(4,4) = (1/3) * ( (S_dr * tau_s^3) / (((R_E + h_k_min_1)^2) * (cos(L_k_min_1))) );
    qew_k_min_1(1,3) = (1/2) * ( (S_dr * tau_s^2) / (R_N + h_k_min_1) );
    qew_k_min_1(2,4) = (1/2) * ( (S_dr * tau_s^2) / ((R_E + h_k_min_1) * cos(L_k_min_1)) );
    qew_k_min_1(3,1) = (1/2) * ( (S_dr * tau_s^2) / (R_N + h_k_min_1) );
    qew_k_min_1(4,2) = (1/2) * ( (S_dr * tau_s^2) / ((R_E + h_k_min_1) * cos(L_k_min_1)) );
    
    
    % EQUATION 8
    x_k_carat_minus = phi_k_min_1 * x_k_min_1_carat_plus;
    
    
    % EQUATION 9
    P_k_minus = phi_k_min_1 * P_sub_zero_super_plus * transpose(phi_k_min_1) + qew_k_min_1;
    
    
    % EQUATION 10
    H_k = zeros(4,4);
    H_k(1,3) = -1;
    H_k(3,1) = -1;
    H_k(4,2) = -1;
    H_k(2,4) = -1;
    
    
    sigma_gr = 0.02;
    % EQUATION 11
    R_k = eye(4);
    R_k(1,1) = sigma_gr^2 / ((R_N + GNSS_stuff(i, 4))^2);
    R_k(2,2) = sigma_gr^2 / ((R_E + GNSS_stuff(i, 4))^2 * (cos(GNSS_stuff(i,2) * deg2rad)^2));
    R_k(3:4, 3:4) = R_k(3:4, 3:4) * sigma_gr^2;
    
    
    % EQUATION 12
    K_k = P_k_minus * transpose(H_k) * (H_k * P_k_minus * transpose(H_k) + R_k);
    
    
    % EQUATION 13 - HAVE TO USE TASK 1
    L_sub_k_super_g = GNSS_stuff(i, 2) * deg2rad;
    L_sub_k_super_d = states(i, 2) * deg2rad;
    lambda_sub_k_super_g = GNSS_stuff(i,3) * deg2rad;
    lambda_sub_k_super_d = states(i,3) * deg2rad;
    
    v_N_G = GNSS_stuff(i,5);
    v_E_G = GNSS_stuff(i,6);
    v_N_D = states(i,4);
    v_E_D = states(i,5);
    
    vec = [L_sub_k_super_g - L_sub_k_super_d; lambda_sub_k_super_g - lambda_sub_k_super_d; v_N_G - v_N_D; v_E_G - v_E_D];
    
    delta_zee_k = vec - (H_k * x_k_carat_minus);
    
    
    % EQUATION 14
    x_carat_plus_subk = x_k_carat_minus + K_k * delta_zee_k;
    
    x_k_min_1_carat_plus = x_carat_plus_subk;
    
    vector2 = zeros(4,1);
    vector2(1,1) = (states(i,2)*deg2rad - x_carat_plus_subk(1)) * rad2deg;
    vector2(2,1) = (states(i,3)*deg2rad - x_carat_plus_subk(2)) * rad2deg;
    vector2(3,1) = v_N_D - x_carat_plus_subk(3);
    vector2(4,1) = v_E_D - x_carat_plus_subk(4);
    
    
    % EQUATION 15
    P_sub_zero_super_plus = (eye(4) - K_k * H_k) * P_k_minus;
    
    L_k_min_1 = GNSS_stuff(i,2) * deg2rad;
    h_k_min_1 = GNSS_stuff(i,4) * deg2rad;
    
    time = time + 0.5;
    
    
    % EQUATION 16
    res = [res; time, transpose(vector2)];
end

clearvars -except states res GNSS_stuff
