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
%Grabbing the first image from the pixel stream and converting it
%to a 480 by 640 matrix (2D)
im1d = reshape(imagesReceived(1,:), 480, 640);

%% Displaying the first raw image
%displaying im1d, the first image
figure();
imshow(im1d);

%% Creating and viewing the time and frequency vectors

samples = length(imagesReceived); % the total length of each image
fs = 1000; %number of pixels per second
T = samples/fs; %The time it takes to get all the pixels needed for an image

t = linspace(0,T, samples+1); t(end) = []; % Created a time Vector
f = linspace(-fs/2, fs/2, samples+1); f(end) = []; % Created a frequency Vector

%Graph displaying the noise found in the matrix imagesRecieved
figure();
plot(t, imagesReceived(1,:));
xlabel('Time [s]')
xlim([0, 300]);
ylabel('Amplitude')
title('The noise in the Time domain of the matrix imagesRecieved')

% changing imagesReceived from being viewed in the time domain to the
% frequency domain by using Fast Fourier Transfer function
Image1 = fftshift(fft(imagesReceived(1,:)))/fs;

%displaying Image1 using the frequency domain
figure()
plot(f,abs(Image1) )
xlabel('Frequency')
ylabel('Amplitude')
title('The noise in the Frequency domain of the new created matrix Image1')

%% Finding a Filter to remove the noise
%Passive values
R1 = 1200; %Ohms
R2 = 1000; %Ohms
C1 = 10*10^-6; %Faradays
C2 = 4.7*10^-6; %Faradays

%????

%Active values
C = 1*10^-6; %Faradays
R = 820; %Ohms

s = 1j.*2.*pi.*f;
vout1 = 1/ (R.*C).^2;
vin1 = s.^2 + (2.*s)./(R.*C) + 1/(R.*C).^2;
vout2 = s.^2;
vin2 = s.^2 + (2.*s)./(R.*C) + 1/(R.*C).^2;
ActFilter1 = vout1./vin1;
ActFilter2 = vout2./vin2;
%plotting Frequency domain of clean image
CleanImage1 = Image1 .* ActFilter1;
figure, plot(f,abs(CleanImage1));
xlabel('Frequency');
ylabel('Amplitude');
title('Frequency Domain of the first Filtered Image');
%reshaping from string to 480 by 640 and displaying image

CleanIm1d = ifft(ifftshift(CleanImage1)) * fs; % was moved into Freq Domain 
%forgraphing, needs to be ing time for image

cleanIm1d_new = reshape(CleanIm1d, [480,640]);
figure;
imshow(real(cleanIm1d_new));
title('De-Noised Image 1');
figure;
plot(t,abs(CleanIm1d));
xlim([0,300]);
xlabel('Time');
ylabel('Amplitude');
title('Time Domain of the first Filtered Image');
% the clipping at a amplitude of zero tells us that there are values 
%cleanImg that are in the imaginary region.

% now for the rest of the images
Image2 = fftshift(fft(imagesReceived(2,:)))/fs;
Image3 = fftshift(fft(imagesReceived(3,:)))/fs;
Image4 = fftshift(fft(imagesReceived(4,:)))/fs;
%image 2
Clean2 = Image2 .* ActFilter1;
CleanIm2d = ifft(ifftshift(Clean2)) * fs;
cleanIm2d_new = reshape(CleanIm2d, [480,640]);
figure;
imshow(real(cleanIm2d_new));
title('De-Noised Image 2');
%image 3
Clean3 = Image3 .* ActFilter1;
CleanIm3d = ifft(ifftshift(Clean3)) * fs;
cleanIm3d_new = reshape(CleanIm3d, [480,640]);
figure;
imshow(real(cleanIm3d_new));
title('De-Noised Image 3');
%image 4
Clean4 = Image4 .* ActFilter1;
CleanIm4d = ifft(ifftshift(Clean4)) * fs;
cleanIm4d_new = reshape(CleanIm4d, [480,640]);
figure;
imshow(real(cleanIm4d_new));
title('De-Noised Image 4');



