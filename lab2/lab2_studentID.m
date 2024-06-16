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

avg_SNR = [];
avg_INR = [];
for d_idx = 1:2:19
    % fprintf("===================================================================\n");
    % fprintf("distance: %d\n", (d_idx + 1) * 25);
    SNR = [];
    INR = [];
    for i = 1:10 % repeat 10 times
        [rx1_location, rx2_location] = generate_rx_location(rand_rx_location_list, d_idx);
        % TODO: Implement Analog Beamforming functions in /tasks
        % Hint: you can adjust input/output for reports
        [rx1_SNR_dbm, rx2_INR_dbm, P_rx_dbm] = analog_beamforming(P_tx_dBm, N0_dBm, 10, tx_location, rx1_location, rx2_location, analog_antenna_number);
        % fprintf('Task1: SNR Calculation\n');
        % fprintf('i = %d, Receiver1\tSNR: %f dBm\tReceiver2 INR: %f dBm\n', i, rx1_SNR_dbm, rx2_INR_dbm);
        SNR = [SNR rx1_SNR_dbm];
        INR = [INR rx2_INR_dbm];
    end
    avg_SNR = [avg_SNR mean(SNR)];
    avg_INR = [avg_INR mean(INR)];
    % fprintf("===================================================================\n");
end

figure
plot(50:50:500, avg_SNR, 50:50:500, avg_INR, 'LineWidth', 2);
legend('SNR', 'INR');
title('average SNR and INR over Distance');

antenna_num = 16;
data = [];
SNR = [];
INR = [];
for codebook_size = [10 5 2.5]
    P_rx_dbm = [];
    for i = 1:10
        [rx1_location, rx2_location] = generate_rx_location(rand_rx_location_list, 7);
        [rx1_SNR_dbm, rx2_INR_dbm, P_rx_dbm_] = analog_beamforming(P_tx_dBm, N0_dBm, codebook_size, tx_location, rx1_location, rx2_location, antenna_num);
        P_rx_dbm = [P_rx_dbm P_rx_dbm_];
        if codebook_size == 10
            SNR = [SNR rx1_SNR_dbm];
            INR = [INR rx2_INR_dbm];
        end
    end
    data = [data; P_rx_dbm];
end

figure
plot(1:10, data(1,:), 1:10, data(2,:), 1:10, data(3,:), 'LineWidth', 2);
legend('codebook size: 19', 'codebook size: 37', 'codebook size: 73', 'Location', 'southeast');
title('power received by receiver 1 over 10 trials with different codebook size')

figure
plot(1:10, SNR, 1:10, INR, 'LineWidth', 2);
legend('SNR', 'INR', 'Location', 'southeast');
title('SNR of receiver 1 and INR of receiver 2 over 10 trials');

% TODO: Implement Digital Beamforming functions in /tasks
% Hint: you can adjust input/output for reports
% [rx1_SNR_dbm, rx2_SNR_dbm, ori_rx1_SNR_dbm, ori_rx2_SNR_dbm] = digital_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, tx_node_number*digital_antenna_number, rx_node_number*rx_antenna_number);

rx1_SNR_dbms = [];
rx2_SNR_dbms = [];
ori_rx1_SNR_dbms = [];
ori_rx2_SNR_dbms = [];
h_eqs = [];
errors = [];
for d_idx = 1:2:19
    rx1_SNR_dbm_sum = 0;
    rx2_SNR_dbm_sum = 0;
    ori_rx1_SNR_dbm_sum = 0;
    ori_rx2_SNR_dbm_sum = 0;
    for i = 1:10
        [rx1_location, rx2_location] = generate_rx_location(rand_rx_location_list, d_idx);
        [rx1_SNR_dbm, rx2_SNR_dbm, ori_rx1_SNR_dbm, ori_rx2_SNR_dbm, H_eq, error] = digital_beamforming(P_tx_dBm, N0_dBm, tx_location, rx1_location, rx2_location, tx_node_number*digital_antenna_number, rx_node_number*rx_antenna_number);
        rx1_SNR_dbm_sum = rx1_SNR_dbm_sum + rx1_SNR_dbm;
        rx2_SNR_dbm_sum = rx2_SNR_dbm_sum + rx2_SNR_dbm;
        ori_rx1_SNR_dbm_sum = ori_rx1_SNR_dbm_sum + ori_rx1_SNR_dbm;
        ori_rx2_SNR_dbm_sum = ori_rx2_SNR_dbm_sum + ori_rx2_SNR_dbm;
        if d_idx == 7
            h_eqs = [h_eqs H_eq];
            errors = [errors error];
        end
    end
    rx1_SNR_dbms = [rx1_SNR_dbms rx1_SNR_dbm_sum/10];
    rx2_SNR_dbms = [rx2_SNR_dbms rx2_SNR_dbm_sum/10];
    ori_rx1_SNR_dbms = [ori_rx1_SNR_dbms ori_rx1_SNR_dbm_sum/10];
    ori_rx2_SNR_dbms = [ori_rx2_SNR_dbms ori_rx2_SNR_dbm_sum/10];
end

figure
plot(50:50:500, rx1_SNR_dbms, 50:50:500, rx2_SNR_dbms, 50:50:500, ori_rx1_SNR_dbms, 50:50:500, ori_rx2_SNR_dbms, 'LineWidth', 2);
legend('receiver 1 with precoding', 'receiver 2 with precoding', 'receiver 1 without precoding', 'receiver 2 without precoding');
title('average SNR over Distance');

figure
plot(h_eqs, errors, 'o');
title('H_eq vs. error');

% fprintf('Task2: SNR Calculation\n');
% fprintf('Receiver1\tSNR: %f dBm\t without precoding SNR: %f dBm\n', rx1_SNR_dbm, ori_rx1_SNR_dbm);
% fprintf('Receiver2\tSNR: %f dBm\t without precoding SNR: %f dBm\n', rx2_SNR_dbm, ori_rx2_SNR_dbm);

function [rx1_location, rx2_location] = generate_rx_location(rand_rx_location_list, distance_index)
    rand_idx = randperm(numel(0:10:180), 2);
    rx1_location = rand_rx_location_list(rand_idx(1), (distance_index:distance_index+1));
    rx2_location = rand_rx_location_list(rand_idx(2), (distance_index:distance_index+1));
    while true
        if rx1_location(2) > 0 && rx2_location(2) > 0
            break;
        end
        rand_idx = randperm(numel(0:10:180), 2);
        rx1_location = rand_rx_location_list(rand_idx(1), (distance_index:distance_index+1));
        rx2_location = rand_rx_location_list(rand_idx(2), (distance_index:distance_index+1));
    end
end