function x = sor_classical(A, b, omega, tol, max_iter)
    % Funkcja rozwiązuje układ równań liniowych za pomocą klasycznej metody BSOR (Backward Successive Over-Relaxation).
    % Parametry wejściowe:
    %   A      - macierz współczynników układu równań liniowych (rozmiar n x 3),
    %            gdzie A(:,1) to dolna przekątna, A(:,2) to główna przekątna,
    %            A(:,3) to górna przekątna (A jest macierzą trójdiagonalną).
    %   b      - wektor wyrazów wolnych (rozmiar n x 1).
    %   omega  - parametr relaksacji (0 < omega < 2), który wpływa na szybkość zbieżności metody.
    %   tol    - tolerancja błędu, określająca warunek zatrzymania iteracji (jeśli różnica między kolejnymi przybliżeniami jest mniejsza od tol, algorytm kończy obliczenia).
    %   max_iter - maksymalna liczba iteracji, po której algorytm zakończy się, nawet jeśli nie osiągnięto tolerancji. 
    % Parametry wyjściowe:
    %   x      - wektor rozwiązania układu równań (rozmiar n x 1), który zawiera przybliżone wartości zmiennych układu. 
    % Opis działania:
    % 1. Funkcja wydobywa dolną, główną i górną przekątną macierzy A.
    %    - lower_diag to dolna przekątna,
    %    - main_diag to główna przekątna,
    %    - upper_diag to górna przekątna.
    % 2. Tworzony jest wektor x, który początkowo zawiera same zera (rozpoczynamy od przybliżenia zerowego).
    % 3. W pętli for wykonywane są iteracje metody BSOR. Dla każdej zmiennej i (od n do 1) obliczany jest nowy
    %    przybliżony wynik x_i^(k+1) z uwzględnieniem poprzednich i następnych elementów układu równań.
    %    Formuła aktualizacji to:
    %    x(i) = (1 - omega) * x(i) + (omega / main_diag(i)) * (b(i) - sum_lower - sum_upper),
    %    gdzie sum_lower i sum_upper to składniki wynikające z sąsiednich elementów macierzy A.
    % 4. Po każdej iteracji sprawdzany jest warunek zatrzymania. Jeśli blad wzgledny jest mniejszy niz zadana tolerancja,
    %    algorytm zatrzymuje się
    % 5. Jeśli po wykonaniu maksymalnej liczby iteracji rozwiązanie nie spełnia warunku tolerancji, proces
    %    zostaje zakończony.
    % Zwracany jest wektor x, który zawiera przybliżone rozwiązanie układu równań.
    % Przykład użycia:
    % A = [a1 b1 c1; a2 b2 c2; ...];  % Przykładowa macierz trójdiagonalna
    % b = [b1; b2; ...];              % Wektor wyrazów wolnych
    % omega = 1.25;                   % Parametr relaksacji
    % tol = 1e-6;                     % Tolerancja błędu
    % max_iter = 1000;                % Maksymalna liczba iteracji
    % x = sor_classical(A, b, omega, tol, max_iter);
    
    % Wydobycie elementów trójdiagonalnych macierzy A
    n = length(b);
    lower_diag = A(:, 1);  % Dolna przekątna
    main_diag = A(:, 2);   % Główna przekątna
    upper_diag = A(:, 3);  % Górna przekątna


    % Inicjalizacja wektora x (rozpoczynamy od zera)
    x = zeros(n, 1);
    
    for iter = 1:max_iter
        x_old = x; % Zapisujemy poprzednią wersję x
        
        for i = 1:n
            % Obliczanie x_i^(k+1) za pomocą metody BSOR
            sum_lower = 0;
            if i > 1
                sum_lower = lower_diag(i) * x(i-1);
            end
            
            sum_upper = 0;
            if i < n
                sum_upper = upper_diag(i) * x(i+1);
            end
            
            x(i) = (1 - omega) * x(i) + (omega / main_diag(i)) * (b(i) - sum_lower - sum_upper);
        end
        
        % Sprawdzamy warunek zatrzymania
        if norm(x - x_old, inf) < tol
            break;
        end
    end
end
