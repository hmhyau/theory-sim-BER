clear, close all
clc

% Declare predefined variables for AWGN channel and 4-PAM
N = 10^6;
d = 1;
variance = [1.25 0.789 0.498 0.314 0.198 0.125 0.0789 0.0395];
L = 2;
numSymbols = 4;
Es = abs(2*(d^2 + (3*d)^2))/4;
EbN0 = (5/4)*(d^2./variance);
EbN0db = (10.*log10(EbN0));

% Allocate dynamic array for number of bit errors and simulated BER
numbitErr = zeros(1, length(variance));
simErrBER = zeros(1, length(variance));

for varIndex = 1:length(variance)
    % Generate equally probable bits and map them into different signal levels
    for i = 1:N
        RandomBits = randi([0,1], 1, 2);
        if RandomBits == [0 0]
            Tx(:,i) = -3*d;
        elseif RandomBits == [0 1] 
            Tx(:,i) = -d;
        elseif RandomBits == [1 1]
            Tx(:,i) = d;
        elseif RandomBits == [1 0]
            Tx(:,i) = 3*d;
        end
    end
    
    % Pass the transmitted symbols into AWGN channel
    awgnNoise = randn(1, N) * sqrt(variance(1,varIndex));
    NoisySignal = Tx  + awgnNoise;
    
    % Optimum detector implementation: Using mid-point of symbol as
    % detection criterion
    Rx(find(NoisySignal >= 2*d)) = 3*d;
    Rx(find(NoisySignal < 2*d & NoisySignal >= 0)) = d;
    Rx(find(NoisySignal < 0 & NoisySignal >= -2*d)) = -d;
    Rx(find(NoisySignal < -2*d)) = -3*d;
    
    % Calculate number of bit errors
    for i = 1:N
        if (abs((Rx(:,i)-Tx(:,i))) == 2)
            numbitErr(varIndex) = (numbitErr(varIndex) + 1);
        elseif (abs((Rx(:,i)-Tx(:,i))) == 4)
            numbitErr(varIndex) = (numbitErr(varIndex) + 2);
        end
    end
end

% Plotting constellation diagram of 4-PAM
% x = linspace(-3*d, 3*d, 4);
% scatterplot(x);
% grid on
% axis([-4 4 -1 1]);
% title('Constellation diagram of 4-PAM');

% Number of bits = 2 * Number of symbols
simErrBER = numbitErr/(L*N);
% theoryBER = 0.75*erfc(sqrt(2*EbN0/5))*(1/L); % Alternative expression for
% theoretical BER
% theoryBER = 0.75*erfc(sqrt(d^2./(2.*variance)))*(1/L)
theoryBER = 0.75*erfc(sqrt((2/5).*(EbN0)))*(1/L)

% Plot both theoretical BER and simulated BER onto the same graph
semilogy(EbN0db, theoryBER, 'rx-');
hold on
semilogy(EbN0db, simErrBER, 'bo');
grid on
legend('TheoryBER', 'simulatedBER')
title('Comparison between Theoretical and Simulated BER');