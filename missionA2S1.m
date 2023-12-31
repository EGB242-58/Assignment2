%% EGB242 Assignment 2, Section 1 %%
% This file is a template for your MATLAB solution to Section 1.
%
% Before starting to write code, generate your data with the startHere
% script as described in the assignment task.

%% Initialise workspace
clear all; close all;
load DataA2 audioMultiplexNoisy fs sid;

% Begin writing your MATLAB solution below this line.
%1.1 Time and frequency domains of audioMultiplexNoisy

t0 = 0; T = 20; 
%time and frequency vectors
t = linspace(t0, t0+T, (fs*T) + 1);
t(end) =[];

fA = linspace(-fs/2, fs/2, (fs*T) + 1);
fA(end) = [];

audioMultiplexNoisyfft = fftshift(fft(audioMultiplexNoisy))/ fs;


figure; hold on;

plot(t, audioMultiplexNoisy);
title ('Time Domain of Multiplexed audio signal');
xlabel ('Time [sec]');
ylabel ('Amplitude');

figure;
plot(fA/1000,abs(audioMultiplexNoisyfft));
title('Frequency Domain of Multiplexed audio signal');
xlabel ('Frequency [kHz]');
ylabel ('Magnitude');


%%1.2 demodulate
fdemod1 = 8030; %frequency tuning in Hz,Audio signals within signal = 8030Hz, 24220Hz, 40250, 56050Hz, 72170Hz
                %noisy frequencies within Multiplex audio frequencies =
                %5365Hz, 26643Hz, 37865Hz, 58491Hz, 74434Hz = 2665, -2423,
                %2385, -2264.
low = 5000; %lowpass filter bounds Hz

audioMultiplexdemod1 = (audioMultiplexNoisy .* cos(2*pi*-fdemod1*t)) + (audioMultiplexNoisy .* cos(2*pi*fdemod1*t)); %demodulate using modulation property and linerity property.

audioSignal1 = lowpass(audioMultiplexdemod1,low,fs);

fdemod2 = 24220; 

audioMultiplexdemod2 = (audioMultiplexNoisy .* cos(2*pi*-fdemod2*t)) + (audioMultiplexNoisy .* cos(2*pi*fdemod2*t)); %demodulate using modulation property and linerity property.

audioSignal2 = lowpass(audioMultiplexdemod2,low,fs);

fdemod3 = 40250;

audioMultiplexdemod3 = (audioMultiplexNoisy .* cos(2*pi*-fdemod3*t)) + (audioMultiplexNoisy .* cos(2*pi*fdemod3*t)); %demodulate using modulation property and linerity property.

audioSignal3 = lowpass(audioMultiplexdemod3,low,fs);

fdemod4 = 56050;

audioMultiplexdemod4 = (audioMultiplexNoisy .* cos(2*pi*-fdemod4*t)) + (audioMultiplexNoisy .* cos(2*pi*fdemod4*t)); %demodulate using modulation property and linerity property.

audioSignal4 = lowpass(audioMultiplexdemod4,low,fs);


fdemod5 = 72170; 

audioMultiplexdemod5 = (audioMultiplexNoisy .* cos(2*pi*-fdemod5*t)) + (audioMultiplexNoisy .* cos(2*pi*fdemod5*t)); %demodulate using modulation property and linerity property.

audioSignal5 = lowpass(audioMultiplexdemod5,low,fs);


%sound(audioSignal#,fs); 1-5 To play any of the Signals

audioSignal1fft = fftshift(fft(audioSignal1))/fs;
audioSignal2fft = fftshift(fft(audioSignal2))/fs;
audioSignal3fft = fftshift(fft(audioSignal3))/fs;
audioSignal4fft = fftshift(fft(audioSignal4))/fs;
audioSignal5fft = fftshift(fft(audioSignal5))/fs;

%Time Domain Plot
figure;
subplot(5,1,1)
plot(t,audioSignal1);
title(sprintf('Time domain plot of the demodulated signal at %.2f kHz', fdemod1/1000));
xlabel ('Time [sec]');
ylabel ('Amplitude');

subplot(5,1,2)
plot(t,audioSignal2);
title(sprintf('Time domain plot of the demodulated signal at %.2f kHz', fdemod2/1000));
xlabel ('Time [sec]');
ylabel ('Amplitude');

subplot(5,1,3)
plot(t,audioSignal3);
title(sprintf('Time domain plot of the demodulated signal at %.2f kHz', fdemod3/1000));
xlabel ('Time [sec]');
ylabel ('Amplitude');

subplot(5,1,4)
plot(t,audioSignal4);
title(sprintf('Time domain plot of the demodulated signal at %.2f kHz', fdemod4/1000));
xlabel ('Time [sec]');
ylabel ('Amplitude');

subplot(5,1,5)
plot(t,audioSignal5);
title(sprintf('Time domain plot of the demodulated signal at %.2f kHz', fdemod5/1000));
xlabel ('Time [sec]');
ylabel ('Amplitude');

%Frequency Domain Plot

figure;
subplot(5,1,1)
plot(fA/1000,abs(audioSignal1fft));
title(sprintf('Frequency domain plot of the demodulated signal at %.2f kHz', fdemod1/1000));
xlabel ('Frequency [kHz]');
ylabel ('Amplitude');
xlim([-low/1000,low/1000]);

subplot(5,1,2)
plot(fA/1000,abs(audioSignal2fft));
title(sprintf('Frequency domain plot of the demodulated signal at %.2f kHz', fdemod2/1000));
xlabel ('Frequency [kHz]');
ylabel ('Amplitude');
xlim([-low/1000,low/1000]);

subplot(5,1,3)
plot(fA/1000,abs(audioSignal3fft));
title(sprintf('Frequency domain plot of the demodulated signal at %.2f kHz', fdemod3/1000));
xlabel ('Frequency [kHz]');
ylabel ('Amplitude');
xlim([-low/1000,low/1000]);

subplot(5,1,4)
plot(fA/1000,abs(audioSignal4fft));
title(sprintf('Frequency domain plot of the demodulated signal at %.2f kHz', fdemod4/1000));
xlabel ('Frequency [kHz]');
ylabel ('Amplitude');
xlim([-low/1000,low/1000]);

subplot(5,1,5)
plot(fA/1000,abs(audioSignal5fft));
title(sprintf('Frequency domain plot of the demodulated signal at %.2f kHz', fdemod5/1000));
xlabel ('Frequency [kHz]');
ylabel ('Amplitude');
xlim([-low/1000,low/1000]);


%% 1.3

ts = 1/fs;
impulse = [1/ts, zeros(1, (fs*T) -1)]; 

impulseResp1 = channel(sid, impulse, fs);


impulserespfreq = fftshift(fft(impulseResp1))/fs;

% plot the inpulse response of the channel
figure();
plot(fA, abs(impulserespfreq));
hold on;
plot(fA, abs(audioMultiplexNoisyfft))
title('Frequency domain of the inpulse response')
xlabel('Frequency')
ylabel('Magnitude')

%% 1.4
x1= ones(1, (T*fs));%normalisation function
AudioDenoisedfreq = audioMultiplexNoisyfft ./ impulserespfreq;

% create the de-noised audio signal in the time domain by performing an 
% inverse Fourier Transform

AudioDenoised = (ifft(ifftshift(AudioDenoisedfreq)) * fs);

audioDenoisedfreqdemod1 = (AudioDenoised .* cos(2*pi*-fdemod1*t)) + (AudioDenoised .* cos(2*pi*fdemod1*t))-x1;
audioDenoisedfreqdemod2 = (AudioDenoised .* cos(2*pi*-fdemod2*t)) + (AudioDenoised .* cos(2*pi*fdemod2*t))-x1; %demodulate using modulation property and linerity property.
audioDenoisedfreqdemod3 = (AudioDenoised .* cos(2*pi*-fdemod3*t)) + (AudioDenoised .* cos(2*pi*fdemod3*t))-x1; %demodulate using modulation property and linerity property.
audioDenoisedfreqdemod4 = (AudioDenoised .* cos(2*pi*-fdemod4*t)) + (AudioDenoised .* cos(2*pi*fdemod4*t))-x1; %demodulate using modulation property and linerity property.
audioDenoisedfreqdemod5 = (AudioDenoised .* cos(2*pi*-fdemod5*t)) + (AudioDenoised .* cos(2*pi*fdemod5*t))-x1;

AudioSignal1 = lowpass(audioDenoisedfreqdemod1,low,fs);
AudioSignal2 = lowpass(audioDenoisedfreqdemod2,low,fs);
AudioSignal3 = lowpass(audioDenoisedfreqdemod3,low,fs);
AudioSignal4 = lowpass(audioDenoisedfreqdemod4,low,fs);
AudioSignal5 = lowpass(audioDenoisedfreqdemod5,low,fs);
% listen to the de-noised audio
%sound(AudioDenoised, fs)

figure;
subplot(5,1,1)
plot(t,AudioSignal1);
title(sprintf('Time domain plot of the demodulated signal at %.2f kHz', fdemod1/1000));
xlabel ('Time [sec]');
ylabel ('Amplitude');

subplot(5,1,2)
plot(t,AudioSignal2);
title(sprintf('Time domain plot of the demodulated signal at %.2f kHz', fdemod2/1000));
xlabel ('Time [sec]');
ylabel ('Amplitude');

subplot(5,1,3)
plot(t,AudioSignal3);
title(sprintf('Time domain plot of the demodulated signal at %.2f kHz', fdemod3/1000));
xlabel ('Time [sec]');
ylabel ('Amplitude');

subplot(5,1,4)
plot(t,AudioSignal4);
title(sprintf('Time domain plot of the demodulated signal at %.2f kHz', fdemod4/1000));
xlabel ('Time [sec]');
ylabel ('Amplitude');

subplot(5,1,5)
plot(t,AudioSignal5);
title(sprintf('Time domain plot of the demodulated signal at %.2f kHz', fdemod5/1000));
xlabel ('Time [sec]');
ylabel ('Amplitude');

% plot the frequency domain of the de-noised audio
figure();
plot(fA, abs(AudioDenoisedfreq));
title('Frequency domain of the de-noised audio')
xlabel('Frequency')
ylabel('Magnitude')

% plot the time domain of the de-noised audio
figure();
plot(t, AudioDenoised);
title('Time domain of the de-noised audio')
xlabel ('Time (s)'); 
ylabel('Channel(t)');

%% 1.5
%remove frequency components to fully denoise the signal

AudioDenoisedfreq(fA==5365) = 0;  
AudioDenoisedfreq(fA==-5365) = 0;
AudioDenoisedfreq(fA==26643) = 0;
AudioDenoisedfreq(fA==-26643) = 0;
AudioDenoisedfreq(fA==37858) = 0;
AudioDenoisedfreq(fA==-37858) = 0;
AudioDenoisedfreq(fA==58491) = 0;
AudioDenoisedfreq(fA==-58491) = 0;
AudioDenoisedfreq(fA==74434) = 0;
AudioDenoisedfreq(fA==-74434) = 0;

%5365Hz, 26643Hz, 37865Hz, 58491Hz, 74434Hz


AudioDenoised2 = ifft(ifftshift(AudioDenoisedfreq)) * fs;


% fully de-noised audio in time domain
figure(14)
subplot(2, 1, 1);
plot(t, AudioDenoised2);
title('Time domain of the fully de-noised audio')
xlabel ('Time (s)'); 
ylabel('Channel(t)');

% fully de-noised audio in frequency domain
subplot(2, 1, 2);
plot(fA, abs(AudioDenoisedfreq));
title('Frequency domain of the fully de-noised audio');
xlabel('Frequency');
ylabel('Magnitude');


AudioDenoisedfreqdemod1 = (AudioDenoised2 .* cos(2*pi*-fdemod1*t)) + (AudioDenoised2 .* cos(2*pi*fdemod1*t))-x1; %demodulate using modulation property and linerity property.

FAudioSignal1 = lowpass(AudioDenoisedfreqdemod1,low,fs);

FAudioSignal1fft = fftshift(fft(FAudioSignal1))/fs;

AudioDenoisedfreqdemod2 = (AudioDenoised2 .* cos(2*pi*-fdemod2*t)) + (AudioDenoised2 .* cos(2*pi*fdemod2*t))-x1; %demodulate using modulation property and linerity property.
AudioDenoisedfreqdemod3 = (AudioDenoised2 .* cos(2*pi*-fdemod3*t)) + (AudioDenoised2 .* cos(2*pi*fdemod3*t))-x1; %demodulate using modulation property and linerity property.
AudioDenoisedfreqdemod4 = (AudioDenoised2 .* cos(2*pi*-fdemod4*t)) + (AudioDenoised2 .* cos(2*pi*fdemod4*t))-x1; %demodulate using modulation property and linerity property.
AudioDenoisedfreqdemod5 = (AudioDenoised2 .* cos(2*pi*-fdemod5*t)) + (AudioDenoised2 .* cos(2*pi*fdemod5*t))-x1; %demodulate using modulation property and linerity property.

FAudioSignal2 = lowpass(AudioDenoisedfreqdemod2,low,fs);
FAudioSignal3 = lowpass(AudioDenoisedfreqdemod3,low,fs);
FAudioSignal4 = lowpass(AudioDenoisedfreqdemod4,low,fs);
FAudioSignal5 = lowpass(AudioDenoisedfreqdemod5,low,fs);

FAudioSignal2fft = fftshift(fft(FAudioSignal2))/fs;
FAudioSignal3fft = fftshift(fft(FAudioSignal3))/fs;
FAudioSignal4fft = fftshift(fft(FAudioSignal4))/fs;
FAudioSignal5fft = fftshift(fft(FAudioSignal5))/fs;

%Time Domain Plot
figure;
subplot(5,1,1)
plot(t,FAudioSignal1);
title(sprintf('Time domain plot of the demodulated clean signal at %.2f kHz', fdemod1/1000));
xlabel ('Time [sec]');
ylabel ('Amplitude');

subplot(5,1,2)
plot(t,FAudioSignal2);
title(sprintf('Time domain plot of the demodulated clean signal at %.2f kHz', fdemod2/1000));
xlabel ('Time [sec]');
ylabel ('Amplitude');

subplot(5,1,3)
plot(t,FAudioSignal3);
title(sprintf('Time domain plot of the demodulated clean signal at %.2f kHz', fdemod3/1000));
xlabel ('Time [sec]');
ylabel ('Amplitude');

subplot(5,1,4)
plot(t,FAudioSignal4);
title(sprintf('Time domain plot of the demodulated clean signal at %.2f kHz', fdemod4/1000));
xlabel ('Time [sec]');
ylabel ('Amplitude');

subplot(5,1,5)
plot(t,FAudioSignal5);
title(sprintf('Time domain plot of the demodulated clean signal at %.2f kHz', fdemod5/1000));
xlabel ('Time [sec]');
ylabel ('Amplitude');

%Frequency Domain Plot

figure;
subplot(5,1,1)
plot(fA/1000,abs(FAudioSignal1fft));
title(sprintf('Frequency domain plot of the demodulated clean signal at %.2f kHz', fdemod1/1000));
xlabel ('Frequency [kHz]');
ylabel ('Amplitude');
xlim([-low/1000,low/1000]);

subplot(5,1,2)
plot(fA/1000,abs(FAudioSignal2fft));
title(sprintf('Frequency domain plot of the demodulated clean signal at %.2f kHz', fdemod2/1000));
xlabel ('Frequency [kHz]');
ylabel ('Amplitude');
xlim([-low/1000,low/1000]);

subplot(5,1,3)
plot(fA/1000,abs(FAudioSignal3fft));
title(sprintf('Frequency domain plot of the demodulated clean signal at %.2f kHz', fdemod3/1000));
xlabel ('Frequency [kHz]');
ylabel ('Amplitude');
xlim([-low/1000,low/1000]);

subplot(5,1,4)
plot(fA/1000,abs(FAudioSignal4fft));
title(sprintf('Frequency domain plot of the demodulated clean signal at %.2f kHz', fdemod4/1000));
xlabel ('Frequency [kHz]');
ylabel ('Amplitude');
xlim([-low/1000,low/1000]);

subplot(5,1,5)
plot(fA/1000,abs(FAudioSignal5fft));
title(sprintf('Frequency domain plot of the demodulated clean signal at %.2f kHz', fdemod5/1000));
xlabel ('Frequency [kHz]');
ylabel ('Amplitude');
xlim([-low/1000,low/1000]);

