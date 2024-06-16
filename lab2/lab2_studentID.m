close all; clear; clc;
addpath ./tasks;

% Environment Configurations
tx_node_number = 1;          % Number of Tx users
rx_node_number = 2;          % Number of Rx users
analog_antenna_number = 16;  % Number of Tx antennas of analog
digital_antenna_number = 2;  % Number of Tx antennas of digital
rx_antenna_number = 1;       % Number of Rx antennas

% Generate receivers with beam direction and distances
% [x y] = rand_rx_location_list(tx_beam_direction_idx)
rand_rx_location_list = [];
for tx_beam = 0:10:180
    tmp = [];
    for d = 50:50:500
        offset = -5 + 10 * rand();       % -5~5 degrees
        x = d * cosd(tx_beam + offset);  % Add a small random offset
        y = d * sind(tx_beam + offset);  % Add a small random offset
        tmp = [tmp x y];
    end
    rand_rx_location_list = [rand_rx_location_list; tmp];
end

% Define coordination and power for transmitter
origin = [0, 0];
tx_location = origin;
P_tx_dBm = 10;          % Transmission power of Tx (dBm)
N0_dBm = -95;           % Assume noise power is -90 dBm

% Randomly select receiver's coordination at d=200
rand_idx = randperm(numel(0:10:180), 2);
rx1_location = rand_rx_location_list(rand_idx(1), (7:8));
rx2_location = rand_rx_location_list(rand_idx(2), (7:8));
while true
    if rx1_location(2) > 0 && rx2_location(2) > 0
        break;
    end
    rand_idx = randperm(numel(0:10:180), 2);
    rx1_location = rand_rx_location_list(rand_idx(1), (7:8));
    rx2_location = rand_rx_location_list(rand_idx(2), (7:8));
end

% TODO: Implement Analog Beamforming functions in /tasks
% Hint: you can adjust input/output for reports
[rx1_SNR_dbm, rx2_INR_dbm] = analog_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, analog_antenna_number);

fprintf('Task1: SNR Calculation\n');
fprintf('Receiver1\tSNR: %f dBm\tReceiver2 INR: %f dBm\n', rx1_SNR_dbm, rx2_INR_dbm);

% TODO: Implement Digital Beamforming functions in /tasks
% Hint: you can adjust input/output for reports
[rx1_SNR_dbm, rx2_SNR_dbm, ori_rx1_SNR_dbm, ori_rx2_SNR_dbm] = digital_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, tx_node_number*digital_antenna_number, rx_node_number*rx_antenna_number);

fprintf('Task2: SNR Calculation\n');
fprintf('Receiver1\tSNR: %f dBm\t without precoding SNR: %f dBm\n', rx1_SNR_dbm, ori_rx1_SNR_dbm);
fprintf('Receiver2\tSNR: %f dBm\t without precoding SNR: %f dBm\n', rx2_SNR_dbm, ori_rx2_SNR_dbm);
