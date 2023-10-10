%% EGB242 Assignment 2, Section 2 %%
% This file is a template for your MATLAB solution to Section 2.
%
% Before starting to write code, generate your data with the startHere
% script as described in the assignment task.

%% Initialise workspace
clear all; close all;

% Begin writing your MATLAB solution below this line.
%% 2.1 step response of the LTI system
T = 20;
samples = 10^4;

t = linspace(0,T,samples+1); t(end)=[];

%create the step input and response functions
step_input = max(0, min(t+1,1));
step_response = 2*(1-exp(-0.5*t));

%plotting both the step input and response function against time vector
figure;
plot(t, step_input, 'b');
hold on;
plot(t, step_response, 'r');
title('Step Response of the Camera System (Motor Only)');
xlabel('time [s]');
ylabel('yaw angle [rads]');
legend('step input', 'step response');

%% 2.2 create feed back in LTI
sys = tf(1, [1,0.5,1]);

%plot the system with step input and the time vector
figure();
lsim(sys, step_input, t);
title('Step Response of Camera System (Motor with Feedback)');
xlabel('time [s]');
ylabel('yaw angle [rads]');

%% 2.4 Gain and its affects
Kfwd01 = 0.1;
Kfwd02 = 0.2;
Kfwd05 = 0.5;
Kfwd1 = 1;
Kfwd2 = 2;

Kfb01 = 0.1;
Kfb02 = 0.2;
Kfb05 = 0.5;
Kfb1 = 1;
Kfb2 = 2;

gain1 = Kfwd1 * Kfb01;
gain2 = Kfwd1 * Kfb02;
gain3 = Kfwd1 * Kfb05;
gain4 = Kfwd1 * Kfb1;
gain5 = Kfwd1 * Kfb2;

gain6 = Kfwd01 * Kfb1;
gain7 = Kfwd02 * Kfb1;
gain8 = Kfwd05 * Kfb1;
gain9 = Kfwd2 * Kfb1;

sys1 = tf(Kfwd1, [1, 0.5, gain1]);
sys2 = tf(Kfwd1, [1, 0.5, gain2]);
sys3 = tf(Kfwd1, [1, 0.5, gain3]);
sys4 = tf(Kfwd1, [1, 0.5, gain4]);
sys5 = tf(Kfwd1, [1, 0.5, gain5]);

sys6 = tf(Kfb1, [1, 0.5, gain6]);
sys7 = tf(Kfb1, [1, 0.5, gain7]);
sys8 = tf(Kfb1, [1, 0.5, gain8]);
sys9 = tf(Kfb1, [1, 0.5, gain9]);

% Plot each of the systems with step_input and time vector
figure();
lsim(sys1, step_input, t);
hold on;
lsim(sys2, step_input, t);
lsim(sys3, step_input, t);
lsim(sys4, step_input, t);
lsim(sys5, step_input, t);hold off
title('Step Response of Camera System testing with different Kfb values');
xlabel('time [s]');
ylabel('yaw angle [rads]');
legend('Kfb = 0.1', 'Kfb = 0.2', 'Kfb = 0.5', 'Kfb = 1', 'Kfb = 2');

figure();
lsim(sys6, step_input, t);
hold on
lsim(sys7, step_input, t);
lsim(sys8, step_input, t);
lsim(sys9, step_input, t);hold off
title('Step Response of Camera System testing with different Kfwd values');
xlabel('time [s]');
ylabel('yaw angle [rads]');
legend('Kfwd = 0.1', 'Kfwd = 0.2', 'Kfwd = 0.5', 'Kfwd = 2');

%% 2.5 Finding Kfwd and Kfb values that will allow for an amplitude of 2pi 
% and a time to peak time equal to 13 seconds

%By testing different values for Kfwd and Kfb the values of 0.7 and 0.16
%respectivly was found
usedKfwd = 0.7;
usedKfb = 0.16;
usedgain = usedKfwd * usedKfb;

cameraTF = tf(usedKfwd, [1, 0.5, usedgain]);

%plotting the used Kfwd and Kfb values
figure;
lsim(cameraTF, step_input, t);
title('Step Response of Camera System');
xlabel('time [s]');
ylabel('yaw angle [rads]')

%% 2.6 testing the rover camera

%setting variables to be able to from from 30 to 210 degrees
startVoltage = 30/360;
endVoltage = 210/360;

% Test camera function
[startIm, finalIm] = cameraPan(startVoltage, endVoltage, cameraTF);
