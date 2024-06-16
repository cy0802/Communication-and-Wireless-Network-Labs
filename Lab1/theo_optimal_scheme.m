SNR_BER_table = load("SNR_BER.mat");
% SNR_BER_table.SNR_BER(1, :) => BPSK, the second index is the SNR
% SNR_BER_table.SNR_BER(2, :) => QPSK, the second index is the SNR
distance = [50:50:600];
pkt_size = [100 2000 4000];
Ptx_dbm = 10;
Ptx = 10 ^ (Ptx_dbm / 10.0) / 1000;
wavelength = 3e8 / 2.4e9;
noise_power_dbm = -90;
noise_power = 10 ^ (noise_power_dbm / 10.0) / 1000;
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
        disp("distance: " + distance(j) + " / pkt_size: " + pkt_size(i));
        disp("SNR: " + SNR + " / Prx_dbm: " + Prx_dbm + " / noise_power_dbm: " + noise_power_dbm);
        disp("BER: " + BER);
        disp("PDR: " + PDR);
        disp("throughput: " + throughput);
    end
end
schemes = ["BPSK" "QPSK" "16QAM" "64QAM"];
disp(schemes(optimal_scheme));
