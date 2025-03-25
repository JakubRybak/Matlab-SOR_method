function [x, liczba_iteracji, macierz_iteracyjna] = sor_iteration_matrix(A, b, omega, tol, max_iter)
    % Funkcja rozwiązuje układ równań liniowych za pomocą metody BSOR
    % (Backward Successive Over-Relaxation), (Używając macierzy iteracyjnej).
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
    %   liczba_iteracji 
    %   macierz_iteracyjna
    % Opis działania:
    % 1. Funkcja najpierw wydobywa dolną, główną i górną przekątną macierzy A.
    % 2. Następnie tworzy macierze D, L, U, które odpowiadają za rozkład macierzy A:
    %    - D to macierz diagonalna (zawiera elementy z głównej przekątnej),
    %    - L to macierz dolna (zawiera elementy z dolnej przekątnej),
    %    - U to macierz górna (zawiera elementy z górnej przekątnej).
    % 3. Wyznacza macierz iteracyjną B_SOR oraz wektor c_SOR, które są używane do aktualizacji przybliżenia rozwiązania.
    % 4. Funkcja iteracyjnie aktualizuje wektor rozwiązania x, wykorzystując wzór metody BSOR:
    %    x_new = B_SOR * x + c_SOR
    % 5. Warunek zatrzymania sprawdza, czy blad wzgledny jest mniejszy od zadanej tolerancji.
    %    Jeśli tak, iteracje są zatrzymywane.
    % 6. Jeśli po wykonaniu maksymalnej liczby iteracji rozwiązanie nie spełnia warunku tolerancji, proces zostaje zakończony.
    % Zwracany jest wektor x, który zawiera przybliżone rozwiązanie układu równań.
    % Przykład użycia:
    % A = [a1 b1 c1; a2 b2 c2; ...];  % Przykładowa macierz trójdiagonalna
    % b = [b1; b2; ...];              % Wektor wyrazów wolnych
    % omega = 1.25;                   % Parametr relaksacji
    % tol = 1e-6;                     % Tolerancja błędu
    % max_iter = 1000;                % Maksymalna liczba iteracji
    % x = sor_iteration_matrix(A, b, omega, tol, max_iter);
    
    n = length(b);
    liczba_iteracji = inf;
    % Wydobycie przekątnych
    lower_diag = A(:, 1);  % Dolna przekątna
    main_diag = A(:, 2);   % Główna przekątna
    upper_diag = A(:, 3);  % Górna przekątna
    
    % Tworzenie macierzy D, L, U
    D = diag(main_diag);   % Macierz diagonalna
    L = diag(lower_diag(2:end), -1);  % Macierz dolna
    U = diag(upper_diag(1:end-1), 1); % Macierz górna

    % Wyznaczamy macierz iteracji B_SOR oraz wektor c_SOR
    B_SOR = inv(D + omega * L) * ((1 - omega) * D - omega * U);
    macierz_iteracyjna = B_SOR;
    c_SOR = omega * inv(D + omega * L) * b';
    
    % Inicjalizacja wektora x
    x = zeros(n, 1);
    
    % Iteracyjne rozwiązanie
    for iter = 1:max_iter
        x_new = B_SOR * x + c_SOR;
        
        % Sprawdzamy warunek zatrzymania
        if norm(x_new - x) < tol * norm(x)
            liczba_iteracji = iter;
            x = x_new;
            break;
        end
        
        x = x_new;
    end
end
