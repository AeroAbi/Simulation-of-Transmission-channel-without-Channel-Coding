%1 Modulator/demodulator
% Parameters
clc;
close all;
clear all;

N=1e6;                % Number of bits to simulate
L=6;
SNRdB=0:1:4;             %With Noise
%SNRdB=100;              %Absence of Noise

% Generate random bits
bits = randi([0 1], 1, N);

%Modulator
M= 4;                %QPSK   
Ns=8 ;               %Oversampling factor-sampling_freq/symbol_rate

%mapping
bits_I=bits(1:2:length(bits));
bits_Q=bits(2:2:length(bits));
symbols =1-2*bits_I+1i*(1-2*bits_Q);

%shaping filter-roll off factor α = 0.35.
diracsM=kron(symbols,[1 zeros(1,Ns-1)]);
ht = rcosdesign(0.35, L, Ns);          
signalt=filter(ht,1,diracsM);           %Transmitted signal
signalpower=mean(abs(signalt).^2);      %Signal Power
figure
plot(signalt,'*');
plot(ht)
title('impulse response Plot');


% Generating the AWGN noise


% Create an array for storing the BER values
BER_values = zeros(size(SNRdB));

% Simulating for each SNR value
for eachSNRdB_index = 1:numel(SNRdB)
    % Calculating the SNR for the given SNR_dB value
    SNR_value = 10.^(0.1 * SNRdB(eachSNRdB_index));
  
    %noise
    noise_power=(signalpower*Ns)/(2*log2(M)*SNR_value);    %noise power
    %complex noise
        noisevar=randn(1, length(signalt));
        noisevar2=randn(1, length(signalt));
        noise=noisevar + 1i*noisevar2;
        awgn_noise = sqrt(noise_power) * noise;
        signalt2 = signalt + awgn_noise;                %Signal with noise
        mean(abs(awgn_noise));
       
    %Setting up the receiver
    hr=rcosdesign(0.35, L, Ns);                   %hr=ht
    signalr=filter(hr,1,signalt2);
  
    
    %eyediagram
    eye=reshape(signalr,Ns,length(signalr)/Ns);
    plot(real(eye));
    xlabel('symbol period Ts');
    ylabel('height');
    title('Eye Diagram @SNR:',SNRdB(eachSNRdB_index));
    
    % sampling
    down_sampled_signal = signalr(1:Ns:end);
    %scatterplot
    figure
    plot(down_sampled_signal,'*')
    title('Scatter Plot @SNR:',SNRdB(eachSNRdB_index));
    
    %Decision Device
    I=real(down_sampled_signal);
    Decidedsym_I=sign(I);                                        %Real part
    Q=imag(down_sampled_signal);
    Decidedsym_Q=sign(Q);                                        %Imaginary part
    estimated_symbols=Decidedsym_I+1i*Decidedsym_Q;      
    estimated_symbols(1:L)=[];                                   %to shorten the length of symbols       
    symbols(length(estimated_symbols)+1:end)=[];
    SER =length(find(estimated_symbols~=symbols))/length(symbols);   %Symbol Error Rate
    figure
    plot(SER,'*')
    title('SER Plot');
       
    %BER
    decidedbits=zeros(1,2*length(estimated_symbols));        
    decidedbits(1:2:end)=real(estimated_symbols);
    decidedbits(2:2:end)=imag(estimated_symbols);
    bits(length(decidedbits)+1:end)=[];                          %to shorten the length of bits-remove the end of bits vector
    BER =length(find(decidedbits~=1-2*bits))/length(decidedbits);
    BER_values(eachSNRdB_index)=BER;
    fprintf("BER @SNR: %f = %f\n", SNRdB,BER);
end

BER_theory = berawgn(SNRdB, 'psk', 4,'nondiff');
figure;
semilogy(SNRdB, BER_theory, 'LineWidth', 2, 'color', 'red');
hold on;
grid on;
semilogy(SNRdB, BER_values, 'LineWidth', 2, 'color', 'black', 'LineStyle', '--', 'Marker', '*');
legend('Theoretical BER', 'Simulated BER');
title('Bit Error Rate Without Coding');
xlabel('Eb/No')
ylabel('BER')
