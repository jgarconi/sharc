% Leia os dados do arquivo
data = dlmread('[SYS] CDF of system interference power from IMT DL.txt');

% Separe as colunas de dados
x = data(:, 1);
y = data(:, 2) * 100;

%x = flipud(x);

% Crie o gráfico
figure;
semilogy(x, y); % Use semilogy for logarithmic y-axis
title('IMT Downlink');
xlabel('Interference [dBW]'); % Interference on the x-axis
ylabel('Probability of Interference < X [%]'); % Probability on the y-axis
grid on;

% Exiba o gráfico
print -dpng 'resultados_grafico.png'; % Salve o gráfico como uma imagem PNG
