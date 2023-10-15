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
imwrite(im2D, 'raw_image.png')

%% 3.2 Creating and viewing the time and frequency vectors

samples = length(imagesReceived);
fs = 1000; % pixels per second
T = samples/fs; % total time taken for image signal to be received

t = linspace(0,T, samples+1); t(end) = []; % time vector
f = linspace(-fs/2, fs/2, samples+1); f(end) = []; % frequency vector

% Plotting im1D image in the time domain
margin = 10; % margin for time axis to make start and end values more visible

figure
plot(t, imagesReceived(1,:))
title('Image Signal in the Time Domain')
xlabel('Time (s)')
xlim([-margin, T+margin])
ylabel('Amplitude')

% Converting im1D to the frequency domain with the Fast Fourier Transform
Im1D = fftshift(fft(im1D))/fs;

% Displaying Im1D in the frequency domain
figure
tiledlayout(2,1)

nexttile
plot(f,abs(Im1D))
title('Magnitude Spectrum of Image Signal')
xlabel('Frequency (Hz)')
ylabel('Amplitude')

nexttile
plot(f,angle(Im1D))
title('Phase Spectrum of Image Signal')
xlabel('Frequency (Hz)')
ylabel('Phase (rads)')

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
H_p1 = ( s*R1 ) / ( (s^2)*C2*R1*(R1+R2) + s*((C2/C1)*(R1+R2)+R1) + 1/C1 );
H_p2 = ( 1/C2 ) / ( s*( (s^2)*C1*R1*R2 + s*(R2+R1*C1/C2) + 1/C2 + R2 ) );
H_a1 = ( 1/((R*C)^2) ) / ( s^2 + s*2/(R*C) + 1/((R*C)^2) );
H_a2 = ( s^2 ) / ( s^2 + s*2/(R*C) + 1/((R*C)^2) );

%ltiview(H_a1)

Im1Dfiltered = Im1D .* H_a1;

figure
tiledlayout(2,1)

nexttile
plot(f,abs(Im1D))
title('Magnitude Spectrum of Image Signal')
xlabel('Frequency (Hz)')
ylabel('Amplitude')

nexttile
plot(f,angle(Im1D))
title('Phase Spectrum of Image Signal')
xlabel('Frequency (Hz)')
ylabel('Phase (rads)')

im1Dfiltered = lsim(H_a1, im1D, t);
im2Dfiltered = reshape(im1Dfiltered, 480, 640);
figure
imshow(im2Dfiltered)