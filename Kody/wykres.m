% Dane wejściowe
A2 = [0 5 -1; -1 5 -1; -1 5 0]; % Macierz A
b2 = [2 8 3]; % Wektor wyrazów wolnych
omega = 1.5;
max_iter = 100;

% Przedziały dokładności
tolerances = logspace(-1, -10, 10); % 10 wartości dokładności w skali logarytmicznej
iterations = zeros(size(tolerances)); % Miejsce na wyniki liczby iteracji

% Oblicz liczbę iteracji dla każdej dokładności
for i = 1:length(tolerances)
    tol = tolerances(i);
    [~, iterations(i), ~] = sor_iteration_matrix(A2, b2, omega, tol, max_iter);
end

% Tworzenie wykresu
figure;
semilogx(tolerances, iterations, '-o', 'LineWidth', 2);
grid on;
xlabel('Dokładność (tolerancja)', 'FontSize', 12);
ylabel('Liczba iteracji', 'FontSize', 12);
title('Zależność liczby iteracji od wymaganej dokładności', 'FontSize', 14);

disp("Macierz:")
disp([5 -1 0; -1 5 -1; 0 -1 5])