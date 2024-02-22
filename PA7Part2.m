% Amarachi Akunna Onyekachi
% akunna1@live.unc.edu
% 06/13/2020
% PA7Part2.m
%
% Records video of bouncing ball and calculates the energy and coefficient
% of restitution

clc
clear
close all

%% Declarations
% Constant
G = -9.81; % acceleration due to gravity (m/s^2)

% Ball properties
ballMass = 0.40; % ball mass (kg)
ballDIn = 8.63;  % ball diameter (in)
ballDm = ballDIn*0.0254; % ball diameter (m)
ballDPx = 191; % ball diameter (px)
px2m = ballDm/ballDPx; % pixel to meter conversion factor
bounceCount = 12;

% Video information
vidFile = 'AkunnaBounce.mp4'; % with extension

% Load video
vid = VideoReader(vidFile);    % reads in .mp4 file
frameRate = vid.FrameRate;
vid_duration = vid.Duration;
total_frames = frameRate * vid_duration;
vid_width = vid.Width;
vid_height = vid.Height;
vid_bits = vid.BitsPerPixel;
frameStart = floor(0.402*frameRate);
frameStop = floor(5.654*frameRate); % the ball rolls off after 5.654ls

% Cropping rectangle dimension
crop_x_start = 7; % x-coordinate to begin cropping image
crop_y_start = 460; % y-coordinate to begin cropping image
crop_x_size = 733; % shifts crop_x_start by crop_x_size along x-axis
crop_y_size = 1460; % shifts crop_y_start by crop_y_size along y-axis

% Equations
PE =  @(h) ballMass * G * -h; % Potential energy formula
KE = @(v) 0.5 * ballMass * v.^2; % Kinteic Energy formula


%% Threshold Video
n = 1; % place holder for centroid index
% Step through each frame 

for k = frameStart:frameStop
    frameSlice = read(vid,k); % loads current frame into frameSlice
    % frameSlice = imrotate(frameSlice,-90); % uncomment this if needed
    
    % Crop image
    frameSlice_crop = imcrop(frameSlice, [crop_x_start,crop_y_start,crop_x_size,crop_y_size]); 
    
    % Changes everything execpt ball to black and changes ball to white
    % by filtering out unwanted colors like greyish-green, brown, white &
    % black
    img_arr_bw = ~((frameSlice_crop(:,:,1) > 120) & (frameSlice_crop(:,:,1) < 200) &...
        (frameSlice_crop(:,:,2) > 120) & (frameSlice_crop(:,:,2) < 200) &...
        (frameSlice_crop(:,:,3) > 120) & (frameSlice_crop(:,:,3) < 200)) &...
        ~((frameSlice_crop(:,:,1) > 55) & (frameSlice_crop(:,:,1) < 80) &...
        (frameSlice_crop(:,:,2) > 30) & (frameSlice_crop(:,:,2) < 55) &...
        (frameSlice_crop(:,:,3) > 30) & (frameSlice_crop(:,:,3) < 60)) &...
        ~((frameSlice_crop(:,:,1) > 140) & (frameSlice_crop(:,:,1) < 250) &...
        (frameSlice_crop(:,:,2) > 140) & (frameSlice_crop(:,:,2) < 250) &...
        (frameSlice_crop(:,:,3) > 140) & (frameSlice_crop(:,:,3) < 250)) &...
        ~((frameSlice_crop(:,:,1) < 120) & (frameSlice_crop(:,:,2) < 120) &...
        (frameSlice_crop(:,:,3) < 120)) &...
        ~((frameSlice_crop(:,:,1) > 10) & (frameSlice_crop(:,:,1) < 240) &...
        (frameSlice_crop(:,:,2) > 94) & (frameSlice_crop(:,:,2) < 160) &...
        (frameSlice_crop(:,:,3) < 160));
            
    % filters the black and white frame further, removing unwanted white 
    % patches
    img_arr_bw = imfilter(img_arr_bw, ones(10)/100);
   
   
    [cent(2,n),cent(1,n)] = Centroid(img_arr_bw); % Uses custom function 
    % centroid to return row and col of centroid
    % centroid returns the x and y coordinates of the centroid given a
    % binary image file
    
    t(n) = k/frameRate;
    
    centroid_plot_Xmin = 0;
    centroid_plot_Xmax = crop_x_size;
    centroid_plot_Ymin = 0;
    centroid_plot_Ymax = crop_y_size; 
       
    figure(1)
    img_sp = subplot(1,2,1); % binary image subplot
    imshow(img_arr_bw)
    title('Binary Image')
    
    cent_sp = subplot(1,2,2); % centroid subplot
    plot(cent(1,:),cent(2,:),'LineWidth',1.5)
    title('Centroid')
    hold on
    p = plot(cent(1,n),cent(2,n),'x','MarkerFaceColor','red','LineWidth',2);
    hold off
    axis([centroid_plot_Xmin centroid_plot_Xmax centroid_plot_Ymin centroid_plot_Ymax])
    pbaspect([crop_x_size crop_y_size 1])
    set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse')
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    drawnow
    
    
    n = n + 1;
end


%% Plot position, velocity, and acceleration
% Corrects centroid position, then calculates velocity and acceleration
for k = 1:size(cent,2)-1
    time_elapsed(k) = t(k+1) - t(1); % Time elapsed frameStart to frameStop
end

time_elapsed = [0 time_elapsed];

% Plot pos/vel/acc versus time
figure(2)
post_sp = subplot(3,1,1);
position_y = cent(2,:)*px2m; % Centroid y-coordinates position in meters
position_y = -position_y; % Flips orientation of y
plot(time_elapsed,position_y)
xlabel('Time (s)')
ylabel('Position (m)')
title('Ball Position')

vel_sp = subplot(3,1,2);
% Velocity - Derivative (dy/dt) of position and time
vel = [0 diff(position_y)./diff(time_elapsed)];
plot(time_elapsed,vel)
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Ball Velocity')

acc_sp = subplot(3,1,3);
% Acceleration - Second Derivative (d2y/dt2) of position and time
acc = [0 diff(position_y,2)]./(diff(time_elapsed)).^2; 
plot(time_elapsed,[0 acc])
xlabel('Time (s)')
ylabel('Acceleration (m/s^{2})')
title('Ball Acceleration')

%% Calculate coefficient of restitution
% Find heights of bounces, compare heights, calculate average excluding
% last bounces without much height

% Calculates height of ball with y-origin at the bottom (ground)
actual_heights = (crop_y_size * px2m) + position_y;

actual_pks = findpeaks(actual_heights); % Peak heights from ground to ball.
actual_pks(1) = actual_heights(1); % Make first peak the position the ball was dropped from
%actual_pks = actual_pks + 0.3;
for n = 1:length(actual_pks)-1
    % Calclates the coefficient of restitution
    e_cr(n) = sqrt(actual_pks(n+1)/actual_pks(n));
end

avg_e_cr = mean(e_cr);

%% Plot peak heights as a stem plot (use stem() )
figure(3)
stem(actual_pks, 'LineWidth',2)
xlim([1 length(actual_pks)])
xlabel('Peak Number')
ylabel('Bounce Height (m)')
title('Bounce Height')

%% Plot for potential, kinetic & total energy
PE_y = PE(actual_heights);
KE_y = KE(vel);
TE_y = PE_y + KE_y;

figure(4)
plot(time_elapsed,PE_y,'k','LineWidth',1.5)
hold on
plot(time_elapsed,KE_y,'r','LineWidth',1.5)
hold on
plot(time_elapsed,TE_y,'b','LineWidth',1.5)
hold off
xlabel('Time (s)')
ylabel('Energy (J)')
title('Energy Plots')
legend('Potential','Kinetic','Total')


%% Output
fprintf('The average coefficient of restitution is %.4f .\n', avg_e_cr);
