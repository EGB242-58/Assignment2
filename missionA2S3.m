%%%% EGB242 Assignment 2, Section 3 %%
% This file is a template for your MATLAB solution to Section 3.
%
% Before starting to write code, generate your data with the startHere
% script as described in the assignment task.

%% Initialise workspace
clear all; close all;
load DataA2 imagesReceived;

% Begin writing your MATLAB solution below this line.

%% Section 3

% Grabbing the first image from the pixel stream and converting it
% to a 480 by 640 matrix (2D)
im1D = imagesReceived(1,:);
im2D = reshape(im1D, 480, 640);

%% 3.1 Displaying the first raw image

% Displaying im2D, the first image
figure();
imshow(im2D);
imwrite(im2D, 'Figures/raw_image.png')

%% 3.2 Creating and viewing the time and frequency vectors

samples = length(im1D);
fs = 1000; % pixels per second
T = samples/fs; % total time taken for image signal to be received

t = linspace(0,T, samples+1); t(end) = []; % time vector
f = linspace(-fs/2, fs/2, samples+1); f(end) = []; % frequency vector

% Plotting im1D image in the time domain
margin = 10; % margin for time axis to make start and end values more visible

figure
plot(t, im1D)
%title('Image Signal in the Time Domain')
xlabel('Time (s)')
xlim([-margin, T+margin])
ylabel('Amplitude')
ylim([-35, 30])
set(gca, 'FontSize', 18)
saveas(gcf, 'Figures/noisy_time.png')

% Converting im1D to the frequency domain with the Fast Fourier Transform
Im1D = fftshift(fft(im1D))/fs;

% Displaying Im1D in the frequency domain
figure
plot(f,abs(Im1D))
xlabel('Frequency (Hz)')
ylabel('Amplitude')
set(gca, 'FontSize', 18)
saveas(gcf, 'Figures/noisy_ft.png')
xlim([-50, 500])
ylim([0, 5])
saveas(gcf, 'Figures/noisy_ft_noise.png')
xlim([-10, 40])
ylim([0, 6])
saveas(gcf, 'Figures/noisy_ft_spikes.png')

%% 3.3 Filter modelling

% Component values
R1 = 1.2e3;     % Ohms
R2 = 1e3;       % Ohms
C1 = 10e-6;     % Farads
C2 = 4.7e-6;    % Farads
R = 820;        % Ohms
C = 1e-6;       % Farads

s = tf('s'); % helper variable for rational expressions of TFs

% Filter transfer functions
H_p1 = ( s*(R1/C2) ) / ( (s^2)*R1*R2 + s*(R1/C2+R2/C1) + 1/(C1*C2) );
H_p2 = ( 1/(C1*C2) ) / ( (s^2)*R1*R2 + s*(R1/C2+R2/C1) + 1/(C1*C2) );
H_a1 = ( 1/((R*C)^2) ) / ( s^2 + s*2/(R*C) + 1/((R*C)^2) );
H_a2 = ( s^2 ) / ( s^2 + s*2/(R*C) + 1/((R*C)^2) );

% Bode magnitude plots of each filter
% figure
% bodemag(H_p1, H_p2, H_a1, H_a2)
% legend('Passive filter 1', 'Passive filter 2', 'Active filter 1', ...
%     'Active filter 2', 'Location', 'south')
% title('')
% xlabel('Frequency', 'FontSize', 18)
% ylabel('Magnitude', 'FontSize', 18)

% Pole zero plots of each filter
fig = figure();
% Eliminate internal labels
plotoptions = pzoptions;
plotoptions.Title.String = '';
plotoptions.Title.FontSize = 0.0001;
plotoptions.XLabel.String = '';
plotoptions.XLabel.FontSize = 0.0001;
plotoptions.YLabel.String = '';
plotoptions.YLabel.FontSize = 0.0001;
subplot(2, 2, 1)
pzplot(H_p1, plotoptions);
subplot(2, 2, 2)
pzplot(H_p2, plotoptions)
subplot(2, 2, 3)
pzplot(H_a1, plotoptions)
subplot(2, 2, 4)
pzplot(H_a2, plotoptions)
% Common axis labels
hdl=axes(fig,'visible','off'); 
hdl.XLabel.Visible='on';
hdl.YLabel.Visible='on';
xlabel(hdl,'Real Axis (seconds^{-1})');
ylabel(hdl,'Imaginary Axis (seconds^{-1})');
saveas(gcf, 'Figures/pole_zero.png')

% Analyse chosen filter
%ltiview(H_a1)

%% 3.4 Applying the Filter

im1Dfiltered = lsim(H_a1, im1D, t);
im2Dfiltered = reshape(im1Dfiltered, 480, 640);

Im1Dfiltered = fftshift(fft(im1Dfiltered))/fs;

figure
plot(f,abs(Im1Dfiltered))
xlabel('Frequency (Hz)')
ylabel('Amplitude')
set(gca, 'FontSize', 18)
saveas(gcf, 'Figures/filtered_magnitude_spectrum.png')

figure
imshow(im2Dfiltered)
imwrite(im2Dfiltered, 'Figures/filtered_image.png')

%% 3.5 Apply to all images

for i = 1:4
im1D = imagesReceived(i,:);
im1Dfiltered = lsim(H_a1, im1D, t);
im2Dfiltered = reshape(im1Dfiltered, 480, 640);
imwrite(im2Dfiltered, ['Filtered_Images/filtered_image_' num2str(i) '.png'])
end

%% Extra fun code to normalise filtered images and see edge artefacts

im1Dfiltered = lsim(H_a2, im1D, t); % Replace H_a2 with desired tf
im1Dnorm = ((im1Dfiltered - mean(im1Dfiltered)) ./ (6 * std(im1Dfiltered))) + 0.5;
im2Dnorm = reshape(im1Dnorm, 480, 640);
figure
imshow(im2Dnorm)