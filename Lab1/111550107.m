throughput_x = zeros(1, 12);
throughput_BPSK = zeros(1, 12);
throughput_QPSK = zeros(1, 12);
throughput_16QAM = zeros(1, 12);
throughput_64QAM = zeros(1, 12);
for distance = 50:50:600
    throughput_x(distance / 50) = distance;
    [Prx_watt, Prx_dbm, ...
    SNR_BPSK, SNR_QPSK, SNR_16QAM, SNR_64QAM, SNR_BPSK_db, SNR_QPSK_db, SNR_16QAM_db, SNR_64QAM_db, ...
    BER_BPSK, BER_QPSK, BER_16QAM, BER_64QAM, ...
    throughput_BPSK_theory, throughput_QPSK_theory, throughput_16QAM_theory, throughput_64QAM_theory, ...
    throughput_BPSK_, throughput_QPSK_, throughput_16QAM_, throughput_64QAM_, optimal_scheme] = Lab1(distance);
    throughput_BPSK(distance / 50) = throughput_BPSK_;
    throughput_QPSK(distance / 50) = throughput_QPSK_;
    throughput_16QAM(distance / 50) = throughput_16QAM_;
    throughput_64QAM(distance / 50) = throughput_64QAM_;
end
figure; hold on;
line1 = plot(throughput_x, throughput_BPSK, 'color', '#D3BBB7', 'LineWidth', 2); label1 = "BPSK";
line2 = plot(throughput_x, throughput_QPSK, 'color', '#B6BBBE', 'LineWidth', 2); label2 = "QPSK";
line3 = plot(throughput_x, throughput_16QAM, 'color', '#9297AB', 'LineWidth', 2); label3 = "16QAM";
line4 = plot(throughput_x, throughput_64QAM, 'color', '#D6BBBE', 'LineWidth', 2); label4 = "64QAM";
legend([line1;line2;line3;line4], label1, label2, label3, label4);

disp("=========================================================")
disp("=========================================================")
theo_optimal_scheme = theoretical_rate_selection();

function [Prx_watt, Prx_dbm, ...
        SNR_BPSK, SNR_QPSK, SNR_16QAM, SNR_64QAM, SNR_BPSK_db, SNR_QPSK_db, SNR_16QAM_db, SNR_64QAM_db, ...
        BER_BPSK, BER_QPSK, BER_16QAM, BER_64QAM, ...
        throughput_BPSK_theory, throughput_QPSK_theory, throughput_16QAM_theory, throughput_64QAM_theory, ...
        throughput_BPSK, throughput_QPSK, throughput_16QAM, throughput_64QAM, optimal_scheme] = Lab1(dis)
    % task 1
    % dis = 5;
    disp("distance: " + dis);
    wavelength = 3e8 / 2.4e9;
    rx_gain = 1;
    tx_gain = 1;
    Ptx_dbm = 10;
    Ptx_watt = 10 ^ (Ptx_dbm / 10.0) / 1000;
    Prx_watt = Ptx_watt * rx_gain * tx_gain * (wavelength / 4 / pi / dis) ^ 2;
    Prx_dbm = 10 * log10(Prx_watt * 1000);

    disp(["Prx: " + Prx_watt + " Watt / " + Prx_dbm + " dBm"])
    % disp(["Ptx: " + Ptx_watt + " Watt / " + Ptx_dbm + " dBm"])

    % task 2
    rng(0);
    h_a = randn(1);
    h_b = randn(1);
    h = h_a + h_b * 1i;
    % Scale h to the receiving power Prx
    h_power = h_a^2 + h_b^2;
    h = h / sqrt(h_power) * sqrt(Prx_watt);

    % task 3
    msg = randi(2, [1, 300000]) - 1;
    msg_size = 300000;
    % BPSK modulation
    msg_BPSK = msg;
    for i = 1:msg_size
        if msg(i) == 0
            msg_BPSK(i) = 1;
        else
            msg_BPSK(i) = -1;
        end
    end
    % QPSK modulation
    msg_QPSK = zeros(1, msg_size / 2);
    for i = 1:2:msg_size
        if msg(i) == 0 && msg(i + 1) == 0
            msg_QPSK((i + 1) / 2) = (1 + 1i) / sqrt(2);
        elseif msg(i) == 0 && msg(i + 1) == 1
            msg_QPSK((i + 1) / 2) = (1 - 1i) / sqrt(2);
        elseif msg(i) == 1 && msg(i + 1) == 0
            msg_QPSK((i + 1) / 2) = (-1 + 1i) / sqrt(2);
        else
            msg_QPSK((i + 1) / 2) = (-1 - 1i) / sqrt(2);
        end
    end

    % 16QAM modulation
    msg_16QAM = zeros(1, msg_size / 4);
    % table_(i) => the modulation result of i-1
    table_16 = [3+3i, 3+1i, 3-1i, 3-3i, 1+3i, 1+1i, 1-1i, 1-3i, -1+3i, -1+1i, -1-1i, -1-3i, -3+3i, -3+1i, -3-1i, -3-3i];
    table_16 = table_16 / sqrt(10.0);
    for i = 1:4:msg_size
        index = 8 * msg(i) + 4 * msg(i + 1) + 2 * msg(i + 2) + msg(i + 3) + 1;
        msg_16QAM((i + 3) / 4) = table_16(index);
    end

    % 64QAM modulation
    msg_64QAM = zeros(1, msg_size / 6);
    table_64 = zeros(1, 64);
    real_ = 7;
    % build table
    for i = 0:7
        imaginary = 7;
        for j = 0:7
            table_64(8 * i + j + 1) = real_ + imaginary * 1i;
            imaginary = imaginary - 2;
        end
        real_ = real_ - 2;
    end
    table_64 = table_64 / sqrt(42.0);
    % modulation
    for i = 1:6:msg_size
        index = 32 * msg(i) + 16 * msg(i + 1) + 8 * msg(i + 2) + 4 * msg(i + 3) + 2 * msg(i + 4) + msg(i + 5) + 1;
        msg_64QAM((i + 5) / 6) = table_64(index);
    end

    % task 4
    N0_dbm = -90;
    N0_watt = 10 ^ (N0_dbm / 10.0) / 1000;
    N_real = sqrt(N0_watt / 2) * randn(1, msg_size);
    N_imaginary = sqrt(N0_watt / 2) * randn(1, msg_size);
    N = N_real + N_imaginary * 1i;

    sum_ = 0;
    for i = 1:msg_size
        sum_ = sum_ + real(N(i)) ^ 2 + imag(N(i)) ^ 2;
    end

    y_BPSK = h * msg_BPSK + N;
    y_QPSK = h * msg_QPSK + N(1:msg_size/2);
    y_16QAM = h * msg_16QAM + N(1:msg_size/4);
    y_64QAM = h * msg_64QAM + N(1:msg_size/6);

    % disp("channel: " + h);
    % disp("generated noise: ");
    % disp(N(1:5));
    % disp("transmitted message: ");
    % disp(msg_BPSK(1:5))
    % disp("received message: ");
    % disp(y_BPSK(1:5))
    % disp(y_QPSK(1:5));
    % disp(y_16QAM(1:5));
    % disp(y_64QAM(1:5));

    % task 5
    x_BPSK = y_BPSK / h;
    x_QPSK = y_QPSK / h;
    x_16QAM = y_16QAM / h;
    x_64QAM = y_64QAM / h;

    x_BPSK_hat = zeros(1, msg_size);
    x_QPSK_hat = zeros(1, msg_size/2);
    x_16QAM_hat = zeros(1, msg_size/4);
    x_64QAM_hat = zeros(1, msg_size/6);

    bit_BPSK = zeros(1, msg_size);
    bit_QPSK = zeros(1, msg_size);
    bit_16QAM = zeros(1, msg_size);
    bit_64QAM = zeros(1, msg_size);

    % BPSK demodulation
    for i = 1:msg_size
        if real(x_BPSK(i)) > 0
            x_BPSK_hat(i) = 1;
            bit_BPSK(i) = 0;
        else
            x_BPSK_hat(i) = -1;
            bit_BPSK(i) = 1;
        end
    end

    % QPSK demodulation
    for i = 1:msg_size/2
        %i = typecast(i_, 'uint16');
        if real(x_QPSK(i)) > 0
            if imag(x_QPSK(i)) > 0
                x_QPSK_hat(i) = (1+1i) / sqrt(2);
                bit_QPSK(2*i-1) = 0;
                bit_QPSK(2*i) = 0;
            else
                x_QPSK_hat(i) = (1-1i) / sqrt(2);
                bit_QPSK(2*i-1) = 0;
                bit_QPSK(2*i) = 1;
            end
        else
            if imag(x_QPSK(i)) > 0
                x_QPSK_hat(i) = (-1+1i) / sqrt(2);
                bit_QPSK(2*i-1) = 1;
                bit_QPSK(2*i) = 0;
            else
                x_QPSK_hat(i) = (-1-1i) / sqrt(2);
                bit_QPSK(2*i-1) = 1;
                bit_QPSK(2*i) = 1;
            end
        end
    end

    % 16QAM demodulation
    % table_16 contain all the constellation points
    % table_16(i+1) => the modulation result of i / j = 1 => constellation point 0
    for i = 1:msg_size/4
        min_dis = 10000;
        for j = 1:16
            dis = real(x_16QAM(i) - table_16(j)) ^ 2 + imag(x_16QAM(i) - table_16(j)) ^ 2;
            if dis < min_dis
                min_dis = dis;
                x_16QAM_hat(i) = table_16(j);
                bit_16QAM(4 * i - 3) = fix((j-1) / 8);
                bit_16QAM(4 * i - 2) = fix((j - 1 - bit_16QAM(4 * i - 3) * 8) / 4);
                bit_16QAM(4 * i - 1) = mod(fix((j-1)/2), 2);
                bit_16QAM(4 * i) = mod(j-1, 2);
            end
        end
    end

    % 64QAM demodulation
    for i = 1:msg_size/6
        min_dis = 10000;
        min_idx = -1;
        for j = 1:64
            dis = real(x_64QAM(i) - table_64(j)) ^ 2 + imag(x_64QAM(i) - table_64(j)) ^ 2;
            if dis < min_dis
                min_dis = dis;
                min_idx = j;
            end
        end
        x_64QAM_hat(i) = table_64(min_idx);
        k = min_idx-1;
        idx = 6 * i;
        for j = 1:6
            bit_64QAM(idx) = mod(k, 2);
            k = fix(k / 2);
            idx = idx - 1;
        end
    end

    % disp("msg_BPSK: ");
    % disp(msg_BPSK(1:5))
    % disp("send over air / y_BPSK = h * msg_BPSK + N")
    % disp(y_BPSK(1:5))
    % disp("receive / x_BPSK = y_BPSK / h")
    % disp(x_BPSK(1:5))
    % disp("demodulate: ")
    % disp(x_BPSK_hat(1:5))

    % disp("========================================")
    % disp(["x_BPSK: ", x_BPSK(1:5)])
    % disp(["x_BPSK_hat: ", x_BPSK_hat(1:5)])
    % disp(["bit_BPSK: ", bit_BPSK(1:5)])
    % disp("========================================")
    % disp(["x_QPSK: ", x_QPSK(1:5)])
    % disp(["x_QPSK_hat: ", x_QPSK_hat(1:5)])
    % disp(["bit_QPSK: ", bit_QPSK(1:10)])
    % disp("========================================")
    % disp(["x_16QAM: ", x_16QAM(1:5)])
    % disp(["x_16QAM_hat: ", x_16QAM_hat(1:5)])
    % disp(["bit_16QAM: ", bit_16QAM(1:20)])
    % disp("========================================")
    % disp(["x_64QAM: ", x_64QAM(1:5)])
    % disp(["x_64QAM_hat: ", x_64QAM_hat(1:5)])
    % disp(["bit_64QAM: ", bit_64QAM(1:30)])
    % disp("========================================")

    % task 6
    noise_BPSK = x_BPSK - msg_BPSK;
    noise_QPSK = x_QPSK - msg_QPSK;
    noise_16QAM = x_16QAM - msg_16QAM;
    noise_64QAM = x_64QAM - msg_64QAM;

    noise_power_BPSK_watt = 0;
    noise_power_QPSK_watt = 0;
    noise_power_16QAM_watt = 0;
    noise_power_64QAM_watt = 0;

    for i = 1:msg_size
        noise_power_BPSK_watt = noise_power_BPSK_watt + real(noise_BPSK(i)) ^ 2 + imag(noise_BPSK(i)) ^ 2;
        if i <= msg_size / 2
            noise_power_QPSK_watt = noise_power_QPSK_watt + real(noise_QPSK(i)) ^ 2 + imag(noise_QPSK(i)) ^ 2;
        end
        if i <= msg_size / 4
            noise_power_16QAM_watt = noise_power_16QAM_watt + real(noise_16QAM(i)) ^ 2 + imag(noise_16QAM(i)) ^ 2;
        end
        if i <= msg_size / 6
            noise_power_64QAM_watt = noise_power_64QAM_watt + real(noise_64QAM(i)) ^ 2 + imag(noise_64QAM(i)) ^ 2;
        end
    end
    noise_power_BPSK_watt = noise_power_BPSK_watt / msg_size;
    noise_power_QPSK_watt = noise_power_QPSK_watt / (msg_size / 2);
    noise_power_16QAM_watt = noise_power_16QAM_watt / (msg_size / 4);
    noise_power_64QAM_watt = noise_power_64QAM_watt / (msg_size / 6);

    noise_power_BPSK_dbm = 10 * log10(noise_power_BPSK_watt * 1000);
    noise_power_QPSK_dbm = 10 * log10(noise_power_QPSK_watt * 1000);
    noise_power_16QAM_dbm = 10 * log10(noise_power_16QAM_watt * 1000);
    noise_power_64QAM_dbm = 10 * log10(noise_power_64QAM_watt * 1000);

    disp("empirical noise power(BPSK): " + noise_power_BPSK_watt + " Watt / " + noise_power_BPSK_dbm + " dBm");
    disp("empirical noise power(QPSK): " + noise_power_QPSK_watt + " Watt / " + noise_power_QPSK_dbm + " dBm");
    disp("empirical noise power(16QAM): " + noise_power_16QAM_watt + " Watt / " + noise_power_16QAM_dbm + " dBm");
    disp("empirical noise power(64QAM): " + noise_power_64QAM_watt + " Watt / " + noise_power_64QAM_dbm + " dBm");

    % task 6
    signal_power_BPSK_watt = 0;
    signal_power_QPSK_watt = 0;
    signal_power_16QAM_watt = 0;
    signal_power_64QAM_watt = 0;
    for i = 1:msg_size
        signal_power_BPSK_watt = signal_power_BPSK_watt + real(x_BPSK(i)) ^ 2 + imag(x_BPSK(i)) ^ 2;
        if i <= msg_size / 2
            signal_power_QPSK_watt = signal_power_QPSK_watt + real(msg_QPSK(i)) ^ 2 + imag(msg_QPSK(i)) ^ 2;
        end
        if i <= msg_size / 4
            signal_power_16QAM_watt = signal_power_16QAM_watt + real(msg_16QAM(i)) ^ 2 + imag(msg_16QAM(i)) ^ 2;
        end
        if i <= msg_size / 6
            signal_power_64QAM_watt = signal_power_64QAM_watt + real(msg_64QAM(i)) ^ 2 + imag(msg_64QAM(i)) ^ 2;
        end
    end
    signal_power_BPSK_watt = signal_power_BPSK_watt / msg_size;
    signal_power_QPSK_watt = signal_power_QPSK_watt / (msg_size / 2);
    signal_power_16QAM_watt = signal_power_16QAM_watt / (msg_size / 4);
    signal_power_64QAM_watt = signal_power_64QAM_watt / (msg_size / 6);

    SNR_BPSK = signal_power_BPSK_watt / noise_power_BPSK_watt;
    SNR_QPSK = signal_power_QPSK_watt / noise_power_QPSK_watt;
    SNR_16QAM = signal_power_16QAM_watt / noise_power_16QAM_watt;
    SNR_64QAM = signal_power_64QAM_watt / noise_power_64QAM_watt;

    SNR_BPSK_db = 10 * log10(SNR_BPSK);
    SNR_QPSK_db = 10 * log10(SNR_QPSK);
    SNR_16QAM_db = 10 * log10(SNR_16QAM);
    SNR_64QAM_db = 10 * log10(SNR_64QAM);

    % disp("signal power(BPSK): " + signal_power_BPSK_watt + " Watt / " + 10 * log10(signal_power_BPSK_watt * 1000) + " dBm");
    % disp("signal power(QPSK): " + signal_power_QPSK_watt + " Watt / " + 10 * log10(signal_power_QPSK_watt * 1000) + " dBm");
    % disp("signal power(16QAM): " + signal_power_16QAM_watt + " Watt / " + 10 * log10(signal_power_16QAM_watt * 1000) + " dBm");
    % disp("signal power(64QAM): " + signal_power_64QAM_watt + " Watt / " + 10 * log10(signal_power_64QAM_watt * 1000) + " dBm");
    disp("SNR_BPSK: " + SNR_BPSK_db + " dB / " + SNR_BPSK + " Watt");
    disp("SNR_QPSK: " + SNR_QPSK_db + " dB / " + SNR_QPSK + " Watt");
    disp("SNR_16QAM: " + SNR_16QAM_db + " dB / " + SNR_16QAM + " Watt");
    disp("SNR_64QAM: " + SNR_64QAM_db + " dB / " + SNR_64QAM + " Watt");

    % task 7
    % calculate the bit error rate
    error_BPSK = zeros(1, msg_size);
    error_QPSK = zeros(1, msg_size);
    error_16QAM = zeros(1, msg_size);
    error_64QAM = zeros(1, msg_size);
    for i = 1:msg_size
        if bit_BPSK(i) ~= msg(i)
            error_BPSK(i) = 1;
        end
        if bit_QPSK(i) ~= msg(i)
            error_QPSK(i) = 1;
        end
        if bit_16QAM(i) ~= msg(i)
            error_16QAM(i) = 1;
        end
        if bit_64QAM(i) ~= msg(i)
            error_64QAM(i) = 1;
        end
    end
    BER_BPSK = sum(error_BPSK) / msg_size;
    BER_QPSK = sum(error_QPSK) / msg_size;
    BER_16QAM = sum(error_16QAM) / msg_size;
    BER_64QAM = sum(error_64QAM) / msg_size;

    disp("empirical BER (BPSK): " + BER_BPSK);
    disp("empirical BER (QPSK): " + BER_QPSK);
    disp("empirical BER (16QAM): " + BER_16QAM);
    disp("empirical BER (64QAM): " + BER_64QAM);

    % PDR
    sample_duration = 3.2 * (10 ^ (-6));
    % bit_per_pkt = 500 * 8;
    bit_per_pkt = 2000;
    PDR_BPSK_theory = (1 - BER_BPSK) ^ bit_per_pkt;
    PDR_QPSK_theory = (1 - BER_QPSK) ^ bit_per_pkt;
    PDR_16QAM_theory = (1 - BER_16QAM) ^ bit_per_pkt;
    PDR_64QAM_theory = (1 - BER_64QAM) ^ bit_per_pkt;
    throughput_BPSK_theory = PDR_BPSK_theory * 1 / sample_duration;
    throughput_QPSK_theory = PDR_QPSK_theory * 2 / sample_duration;
    throughput_16QAM_theory = PDR_16QAM_theory * 4 / sample_duration;
    throughput_64QAM_theory = PDR_64QAM_theory * 6 / sample_duration;

    pkt_cnt_BPSK = 0;
    pkt_cnt_QPSK = 0;
    pkt_cnt_16QAM = 0;
    pkt_cnt_64QAM = 0;
    total_pkt = msg_size / bit_per_pkt;
    for i = 1:msg_size / bit_per_pkt
        if sum(error_BPSK((i - 1) * bit_per_pkt + 1:i * bit_per_pkt)) == 0
            pkt_cnt_BPSK = pkt_cnt_BPSK + 1;
        end
        if sum(error_QPSK((i - 1) * bit_per_pkt + 1:i * bit_per_pkt)) == 0
            pkt_cnt_QPSK = pkt_cnt_QPSK + 1;
        end
        if sum(error_16QAM((i - 1) * bit_per_pkt + 1:i * bit_per_pkt)) == 0
            pkt_cnt_16QAM = pkt_cnt_16QAM + 1;
        end
        if sum(error_64QAM((i - 1) * bit_per_pkt + 1:i * bit_per_pkt)) == 0
            pkt_cnt_64QAM = pkt_cnt_64QAM + 1;
        end
    end
    PDR_BPSK = pkt_cnt_BPSK / total_pkt;
    PDR_QPSK = pkt_cnt_QPSK / total_pkt;
    PDR_16QAM = pkt_cnt_16QAM / total_pkt;
    PDR_64QAM = pkt_cnt_64QAM / total_pkt;
    throughput_BPSK = PDR_BPSK * 1 / sample_duration;
    throughput_QPSK = PDR_QPSK * 2 / sample_duration;
    throughput_16QAM = PDR_16QAM * 4 / sample_duration;
    throughput_64QAM = PDR_64QAM * 6 / sample_duration;

    % bugs
    disp("empirical throughput (BPSK): " + throughput_BPSK + " bps");
    disp("empirical throughput (QPSK): " + throughput_QPSK + " bps");
    disp("empirical throughput (16QAM): " + throughput_16QAM + " bps");
    disp("empirical throughput (64QAM): " + throughput_64QAM + " bps");
    disp("theoretical throughput (BPSK): " + throughput_BPSK_theory + " bps");
    disp("theoretical throughput (QPSK): " + throughput_QPSK_theory + " bps");
    disp("theoretical throughput (16QAM): " + throughput_16QAM_theory + " bps");
    disp("theoretical throughput (64QAM): " + throughput_64QAM_theory + " bps");

    optimal_scheme = "BPSK";
    if throughput_QPSK > throughput_BPSK
        optimal_scheme = "QPSK";
    elseif throughput_16QAM > throughput_QPSK
        optimal_scheme = "16QAM";
    elseif throughput_64QAM > throughput_16QAM
        optimal_scheme = "64QAM";
    end
    disp("optimal modulation scheme: " + optimal_scheme);
    disp("========================================================= ")

    % task 8
    % plot the constellation diagram
    BPSK_x = zeros(1, msg_size);
    BPSK_y = zeros(1, msg_size);
    QPSK_x = zeros(1, msg_size / 2);
    QPSK_y = zeros(1, msg_size / 2);
    QAM16_x = zeros(1, msg_size / 4);
    QAM16_y = zeros(1, msg_size / 4);
    QAM64_x = zeros(1, msg_size / 6);
    QAM64_y = zeros(1, msg_size / 6);
    for i = 1:msg_size
        BPSK_x(i) = real(x_BPSK(i));
        BPSK_y(i) = imag(x_BPSK(i));
        if i <= msg_size / 2
            QPSK_x(i) = real(x_QPSK(i));
            QPSK_y(i) = imag(x_QPSK(i));
        end
        if i <= msg_size / 4
            QAM16_x(i) = real(x_16QAM(i));
            QAM16_y(i) = imag(x_16QAM(i));
        end
        if i <= msg_size / 6
            QAM64_x(i) = real(x_64QAM(i));
            QAM64_y(i) = imag(x_64QAM(i));
        end
    end
    % f_BPSK = figure;
    % f_QPSK = figure;
    % f_16QAM = figure;
    % f_64QAM = figure;
    % figure(f_BPSK);
    % scatter(BPSK_x, BPSK_y, 1, 'b', 'filled');
    % figure(f_QPSK);
    % scatter(QPSK_x, QPSK_y, 1, 'b', 'filled');
    % figure(f_16QAM);
    % scatter(QAM16_x, QAM16_y, 1, 'b', 'filled');
    % figure(f_64QAM);
    % scatter(QAM64_x, QAM64_y, 1, 'b', 'filled');
end

function [optimal_scheme] = theoretical_rate_selection()
    SNR_BER_table = load("SNR_BER.mat");
    % SNR_BER_table.SNR_BER(1, :) => BPSK, the second index is the SNR
    % SNR_BER_table.SNR_BER(2, :) => QPSK, the second index is the SNR
    distance = [50:50:600];
    pkt_size = [100 2000 4000];
    Ptx_dbm = 10;
    Ptx = 10 ^ (Ptx_dbm / 10.0) / 1000;
    wavelength = 3e8 / 2.4e9;
    noise_power_dbm = -90;
    optimal_scheme = zeros(3, 12);
    for i = 1:3
        for j = 1:12
            Prx = Ptx * (wavelength / 4 / pi / distance(j)) ^ 2;
            Prx_dbm = 10 * log10(Prx * 1000);
            SNR = Prx_dbm - noise_power_dbm;
            
            BER = [SNR_BER_table.SNR_BER(1, fix(SNR)) SNR_BER_table.SNR_BER(2, fix(SNR))...
                SNR_BER_table.SNR_BER(3, fix(SNR)) SNR_BER_table.SNR_BER(4, fix(SNR))];
            PDR = (1 - BER) .^ pkt_size(i);
            throughput = PDR .* pkt_size(i);
            optimal_scheme(i, j) = 1;
            for k = 2:4
                if throughput(k) >= throughput(optimal_scheme(i, j))
                    optimal_scheme(i, j) = k;
                end
            end
            % disp("distance: " + distance(j) + " / pkt_size: " + pkt_size(i));
            % disp("SNR: " + SNR + " / Prx_dbm: " + Prx_dbm + " / noise_power_dbm: " + noise_power_dbm);
            % disp("BER: " + BER);
            % disp("PDR: " + PDR);
            % disp("throughput: " + throughput);
        end
    end
    schemes = ["BPSK" "QPSK" "16QAM" "64QAM"];
    disp("theoretical optimal scheme: ");
    disp(schemes(optimal_scheme));
end