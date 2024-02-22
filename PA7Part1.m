% Amarachi Akunna Onyekachi
% akunna1@live.unc.edu
% 06/13/2020
% PA7Part1.m
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
ballMass = 0.03;    % ball mass (kg)
ballDIn = 4;        % ball diameter (in)
ballDm = ballDIn*0.0254; % ball diameter (m)
ballDPx = 110; % ball diameter (px)
px2m = ballDm/ballDPx; % pixel to meter conversion factor
bounceCount = 5;

% Video information
vidFile = 'bounce.avi'; % with extension

% Load video
vid = VideoReader(vidFile);    % reads in .avi file
frameRate = vid.FrameRate;
vid_duration = vid.Duration;
total_frames = frameRate * vid_duration;
vid_width = vid.Width;
vid_height = vid.Height;
vid_bits = vid.BitsPerPixel;
frameStart = floor(0.411*frameRate);
frameStop = floor(2.431*frameRate); % the ball rolls off after 2.431s

% Cropping rectangle dimension
crop_x_start = 800; % x-coordinate to begin cropping image
crop_y_start = 0; % y-coordinate to begin cropping image
crop_x_size = 300; % shifts crop_x_start by crop_x_size along x-axis
crop_y_size = 720; % shifts crop_y_start by crop_y_size along y-axis

% Threshold
rTh = 160; % Ball color (red) threshold;

% Equations
PE =  @(h) ballMass * G * -h; % Potential energy formula
KE = @(v) 0.5 * ballMass * v.^2; % Kinteic Energy formula


%% Threshold Video
n = 1; % place holder for centroid index
% Step through each frame 
for k = frameStart:frameStop
    frameSlice = read(vid,k); % loads current frame into frameSlice
    
    % Crop image
    frameSlice_crop = imcrop(frameSlice, [crop_x_start,crop_y_start,crop_x_size,crop_y_size]); 

    % Changes image to only have red color
    img_arr_zeros = zeros(size(frameSlice_crop)); % makes copy of cropped
    % image black
    
    img_arr(:,:,1) = frameSlice_crop(:,:,1); % keeps red color matrix from cropped frame
    img_arr(:,:,2) = img_arr_zeros(:,:,1); % changes all green colors to black
    img_arr(:,:,3) = img_arr_zeros(:,:,1); % changes all blue colors to black
    
    img_arr(:,:,1) = (img_arr(:,:,1) > rTh).*255; % changes all red colors 
    % greater than 160 to 255 and the rest to 0
   
    
    img_arr_bw = imbinarize(img_arr(:,:,1)); % changes frame to black/white  
   
    
    [cent(2,n),cent(1,n)] = Centroid(img_arr_bw);  
    % Uses custom function centroid to return row and col of centroid
    % centroid returns the x and y coordinates of the centroid given a
    % binary image file
    
    t(n) = k/frameRate;
    
    centroid_plot_Xmin = 0;
    centroid_plot_Xmax = crop_x_size;
    centroid_plot_Ymin = crop_y_start;
    centroid_plot_Ymax = crop_y_size;
    % Display the thresholded image and plot centroid movement dynamically
    % Make sure image and plot are same size, and no change in axes from
    % plot to plot
    
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
% Corrects centroid position, then calculate velocity and acceleration
for k = 1:size(cent,2)-1
    dist(k) = sqrt((cent(2,k+1) - cent(2,k))^2 + (cent(1,k+1) - cent(1,k))^2) * px2m;
    time_elapsed(k) = t(k+1) - t(1); % Time elapsed frameStart to frameStop
end

time_elapsed = [0 time_elapsed];

% Plot pos/vel/acc versus time
figure(2)
post_sp = subplot(3,1,1);
position_y = cent(2,:)*px2m; % Centroid y-coordinates position in meters
position_y = -position_y;
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
actual_heights = ((crop_y_size - crop_y_start) * px2m) + position_y;

actual_pks = findpeaks(actual_heights); % Peak heights from ground to ball.
actual_pks(1) = actual_heights(1);

for n = 1:length(actual_pks)-1
    e_cr(n) = sqrt(actual_pks(n+1)/actual_pks(n));
end

avg_e_cr = mean(e_cr);

%% Plot peak heights as a stem plot (use stem() )
figure(3)
stem(actual_pks, 'LineWidth',2)
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
