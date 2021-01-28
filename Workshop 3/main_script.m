% The file Workshop3_Speed_Heading.csv contains 175s of car odometry and magnetic
% compass data. Column 1 of this comma-separated variable (CSV) format file contains time
% in seconds, column 2 contains the forward speed in metres per second and column 3
% contains heading in degrees. The compass provides a measurement of the instantaneous
% heading. However, the odometer counts wheel rotations from which it calculates the
% distance travelled. The odometer speed output is therefore the average speed since the
% previous speed measurement, not the instantaneous speed.
% The initial geodetic latitude and longitude (at time 0) are 50.4249580 and −3.5957974,
% respectively. Use the speed and heading measurements to compute a position solution for
% the rest of the car’s trajectory.

% column 3 is psi
clear all

% PART 1

speed_heading = csvread('Workshop3_Speed_Heading.csv');
init_lat = 50.4249580;
init_long = -3.5957974;

velocities = [];
psi_k_min_1 = 0;
deg2rad = pi / 180;
for i = 1:351   % CONVERT TO RADIANS
   if (i == 1)
       vel = [(cos(speed_heading(i, 3)*deg2rad)); (sin(speed_heading(i,3)*deg2rad))];
   else
       vel = [(cos(speed_heading(i, 3)*deg2rad)+cos(psi_k_min_1*deg2rad)); (sin(speed_heading(i,3)*deg2rad) + sin(psi_k_min_1*deg2rad))]; 
   end
   vel = vel * 0.5 * speed_heading(i,2);
   velocities = [velocities; transpose(vel)];
   psi_k_min_1 = speed_heading(i,3);
end


% PART 2
lat_long = [];

h = 37.4;
last_lat = init_lat * deg2rad;
last_long = init_long * deg2rad;
rad2deg = 180 / pi;

time = 1;
for i = 1:351
    
    [R_N,R_E]= Radii_of_curvature(last_lat);
    
    latitude = last_lat + ((velocities(i, 1) * (0.5)) / (R_N + h));
    longitude = last_long + ((velocities(i,2)*(0.5)) / ((R_E + h) * cos(latitude)));
    
    lat_long = [lat_long; latitude * rad2deg, longitude * rad2deg];
    
    last_lat = latitude;
    last_long = longitude;
    time = time + 0.5;
end

clear last_lat last_long latitude longitude R_E R_N i h time vel psi_k_min_1

% PART 3

damped_dr_vel = [];
last_vel_1 = 0;
last_vel_2 = 0;
for i = 1:351.
    if (i == 1)
        damped_dr_vel = [damped_dr_vel; velocities(i,1), velocities(i,2)];
    else
        damped_dr_vel = [damped_dr_vel; (1.7 * velocities(i, 1) - 0.7 * last_vel_1), (1.7 * velocities(i,2) - 0.7 * last_vel_2)];
    end
    last_vel_1 = velocities(i,1);
    last_vel_2 = velocities(i,2);
end
clear i last_vel_1 last_vel_2 


results = [lat_long, damped_dr_vel];
