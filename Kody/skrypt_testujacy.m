% Przygotowanie przykładowych danych
A1 = [0 4 -1; -1 4 -1; -1 4 0];  % Macierz przechowywana jak n x 3
Macierz1 = [4 -1 0; -1 4 -1; 0 -1 4];
b1 = [5 3 2];  % Wektor wyrazów wolnych

omega1 = 1.5;  % Parametr relaksacji
tol1 = 1e-5;    % Tolerancja błędu
max_iter1 = 100;  % Maksymalna liczba iteracji

% Wywołanie funkcji dla klasycznej metody SOR\
disp("-----------------------------------------")
disp("Pierwszy przypadek:")
[x1, liczba_iteracji1, macierz_iteracyjna1] = sor_iteration_matrix(A1, b1, omega1, tol1, max_iter1);
rozwiazanie1 = linsolve(Macierz1, b1');
wskaznik1 = norm(inv(Macierz1)) * norm(Macierz1);
blad_wzgledny1 = norm(rozwiazanie1 - x1)/norm(rozwiazanie1);
disp("Macierz: ")
disp(Macierz1)
disp("Omega, dokładność, maksymalna liczba iteracji: ");
fprintf('omega = %.2f, tolerancja = %.1e, max liczba iteracji = %d\n', omega1, tol1, max_iter1);
disp('Rozwiązanie za pomocą metody BSOR:');
disp(x1);
disp("Roziwązanie dokładne: ")
disp(rozwiazanie1)
disp("Błąd względny:")
disp(blad_wzgledny1)
disp("Liczba iterecji:")
disp(liczba_iteracji1)
disp("Promien spektralny macierzy iteracyjnej:")
disp(max(abs(eig(macierz_iteracyjna1))))
disp("Wskaznik uwarunkowania:")
disp(wskaznik1)
disp(" ")
disp(" ")
disp(" ")

% Przygotowanie przykładowych danych
A2 = [0 4 -2; -2 5 -2; 2 4 0];  % Macierz przechowywana jak n x 3
Macierz2 = [4 -2 0; -2 5 -2; 0 -2 4];
b2 = [2 8 3];  % Wektor wyrazów wolnych

omega2 = 1.25;  % Parametr relaksacji
tol2 = 1e-6;    % Tolerancja błędu
max_iter2 = 100;  % Maksymalna liczba iteracji

% Wywołanie funkcji dla klasycznej metody SOR
disp("Drugi przypadek:")
[x2, liczba_iteracji2, macierz_iteracyjna2] = sor_iteration_matrix(A2, b2, omega2, tol2, max_iter2);
rozwiazanie2 = linsolve(Macierz2, b2');
wskaznik2 = norm(inv(Macierz2)) * norm(Macierz2);
blad_wzgledny2 = norm(rozwiazanie2 - x2)/norm(rozwiazanie2);
disp("Macierz: ")
disp(Macierz2)
disp("Omega, dokładność, maksymalna liczba iteracji: ");
fprintf('omega = %.2f, tolerancja = %.1e, max liczba iteracji = %d\n', omega2, tol2, max_iter2);
disp('Rozwiązanie za pomocą metody BSOR:');
disp(x2);
disp("Roziwązanie dokładne: ")
disp(rozwiazanie2)
disp("Błąd względny:")
disp(blad_wzgledny2)
disp("Liczba iterecji:")
disp(liczba_iteracji2)
disp("Promien spektralny macierzy iteracyjnej:")
disp(max(abs(eig(macierz_iteracyjna2))))
disp("Wskaznik uwarunkowania:")
disp(wskaznik2)
disp(" ")
disp(" ")
disp(" ")


% Przygotowanie przykładowych danych
A3 = [0 1 -1; -1 2 -1; -1 1 0];  % Macierz przechowywana jak n x 3
Macierz3 = [1 -1 0; -1 2 -1; 0 -1 5];
b3 = [7 11 1];  % Wektor wyrazów wolnych

omega3 = 1.35;  % Parametr relaksacji
tol3 = 1e-5;    % Tolerancja błędu
max_iter3 = 100;  % Maksymalna liczba iteracji

% Wywołanie funkcji dla klasycznej metody SOR
disp("Trzeci przypadek:")
[x3, liczba_iteracji3, macierz_iteracyjna3] = sor_iteration_matrix(A3, b3, omega3, tol3, max_iter3);
rozwiazanie3 = linsolve(Macierz3, b3');
wskaznik3 = norm(inv(Macierz3)) * norm(Macierz3);
blad_wzgledny3 = norm(rozwiazanie3 - x3)/norm(rozwiazanie3);
disp("Macierz: ")
disp(Macierz3)
disp("Omega, dokładność, maksymalna liczba iteracji: ");
fprintf('omega = %.2f, tolerancja = %.1e, max liczba iteracji = %d\n', omega3, tol3, max_iter3);
disp('Rozwiązanie za pomocą metody BSOR:');
disp(x3);
disp("Roziwązanie dokładne: ")
disp(rozwiazanie3)
disp("Błąd względny:")
disp(blad_wzgledny3)
disp("Liczba iterecji:")
disp(liczba_iteracji3)
disp("Promien spektralny macierzy iteracyjnej:")
disp(max(abs(eig(macierz_iteracyjna3))))
disp("Wskaznik uwarunkowania:")
disp(wskaznik3)
disp(" ")
disp(" ")
disp(" ")

% Przygotowanie przykładowych danych
A4 = [0 1 0; 0.999 1 0.999; 0.999 1 0.999];  % Macierz przechowywana jak n x 3
Macierz4 = [1 0.999 0; 0.999 1 0.999; 0 0.999 1];
b4 = [1 2 3];  % Wektor wyrazów wolnych

omega4 = 1.35;  % Parametr relaksacji
tol4 = 1e-5;    % Tolerancja błędu
max_iter4 = 100;  % Maksymalna liczba iteracji

% Wywołanie funkcji dla klasycznej metody SOR
disp("Czwarty przypadek:")
[x4, liczba_iteracji4, macierz_iteracyjna4] = sor_iteration_matrix(A4, b4, omega4, tol4, max_iter4);
rozwiazanie4 = linsolve(Macierz4, b4');
wskaznik4 = norm(inv(Macierz4)) * norm(Macierz4);
blad_wzgledny4 = norm(rozwiazanie4 - x4)/norm(rozwiazanie4);
disp("Macierz: ")
disp(Macierz4)
disp("Omega, dokładność, maksymalna liczba iteracji: ");
fprintf('omega = %.2f, tolerancja = %.1e, max liczba iteracji = %d\n', omega4, tol4, max_iter4);
disp('Rozwiązanie za pomocą metody BSOR:');
disp(x4);
disp("Roziwązanie dokładne: ")
disp(rozwiazanie4)
disp("Błąd względny:")
disp(blad_wzgledny4)
disp("Liczba iterecji:")
disp(liczba_iteracji4)
disp("Promien spektralny macierzy iteracyjnej:")
disp(max(abs(eig(macierz_iteracyjna4))))
disp("Wskaznik uwarunkowania:")
disp(wskaznik4)
disp(" ")
disp(" ")
disp(" ")

% Przygotowanie przykładowych danych
A5 = [0 0.1 -1; -1 0.1 -1; -1 0.1 0];  % Macierz przechowywana jak n x 3
Macierz5 = [0.1 -1 0; -1 0.1 -1; 0 -1 0.1];
b5 = [1 1 1];  % Wektor wyrazów wolnych

omega5 = 1.35;  % Parametr relaksacji
tol5 = 1e-5;    % Tolerancja błędu
max_iter5 = 100;  % Maksymalna liczba iteracji

% Wywołanie funkcji dla klasycznej metody SOR
disp("Piąty przypadek:")
[x5, liczba_iteracji5, macierz_iteracyjna5] = sor_iteration_matrix(A5, b5, omega5, tol5, max_iter5);
rozwiazanie5 = linsolve(Macierz5, b5');
wskaznik5 = norm(inv(Macierz5)) * norm(Macierz5);
blad_wzgledny5 = norm(rozwiazanie5 - x5)/norm(rozwiazanie5);
disp("Macierz: ")
disp(Macierz5)
disp("Omega, dokładność, maksymalna liczba iteracji: ");
fprintf('omega = %.2f, tolerancja = %.1e, max liczba iteracji = %d\n', omega5, tol5, max_iter5);
disp('Rozwiązanie za pomocą metody BSOR:');
disp(x5);
disp("Roziwązanie dokładne: ")
disp(rozwiazanie5)
disp("Błąd względny:")
disp(blad_wzgledny5)
disp("Liczba iterecji:")
disp(liczba_iteracji5)
disp("Promien spektralny macierzy iteracyjnej:")
disp(max(abs(eig(macierz_iteracyjna5))))
disp("Wskaznik uwarunkowania:")
disp(wskaznik5)
disp(" ")
disp(" ")
disp(" ")

% Przygotowanie przykładowych danych
A6 = [0 4 -1; -1 4 -1; -1 4 0];  % Macierz przechowywana jak n x 3
Macierz6 = [4 -1 0; -1 4 -1; 0 -1 4];
b6 = [1 1 1];  % Wektor wyrazów wolnych

omega6 = 1.35;  % Parametr relaksacji
tol6 = 1e-5;    % Tolerancja błędu
max_iter6 = 100;  % Maksymalna liczba iteracji

% Wywołanie funkcji dla klasycznej metody SOR
disp("Szósty przypadek:")
[x6, liczba_iteracji6, macierz_iteracyjna6] = sor_iteration_matrix(A6, b6, omega6, tol6, max_iter6);
rozwiazanie6 = linsolve(Macierz6, b6');
wskaznik6 = norm(inv(Macierz6)) * norm(Macierz6);
blad_wzgledny6 = norm(rozwiazanie6 - x6)/norm(rozwiazanie6);
disp("Macierz: ")
disp(Macierz6)
disp("Omega, dokładność, maksymalna liczba iteracji: ");
fprintf('omega = %.2f, tolerancja = %.1e, max liczba iteracji = %d\n', omega6, tol6, max_iter6);
disp('Rozwiązanie za pomocą metody BSOR:');
disp(x6);
disp("Roziwązanie dokładne: ")
disp(rozwiazanie6)
disp("Błąd względny:")
disp(blad_wzgledny6)
disp("Liczba iterecji:")
disp(liczba_iteracji6)
disp("Promien spektralny macierzy iteracyjnej:")
disp(max(abs(eig(macierz_iteracyjna6))))
disp("Wskaznik uwarunkowania:")
disp(wskaznik6)
disp(" ")
disp(" ")
disp(" ") 



n = 20;

main_diag = 2 * ones(n, 1);     
upper_diag = 1 * ones(n-1, 1); 
lower_diag = -1 * ones(n-1, 1);  %

A = diag(main_diag) + diag(upper_diag, 1) + diag(lower_diag, -1);


% Wektory dla przekątnych
main_diag = 2 * ones(n, 1);           
upper_diag = [1 * ones(n-1, 1); 0];   
lower_diag = [0; -1 * ones(n-1, 1)];  

tridiagonal_form = [lower_diag, main_diag, upper_diag];


% Przygotowanie przykładowych danych
A7 = tridiagonal_form;  % Macierz przechowywana jak n x 3
Macierz7 = A;
b7 = [1 3 5 7 8 5 4 7 8 1 9 2 3 4 2 5 8 9 1 2];  % Wektor wyrazów wolnych

omega7 = 1.35;  % Parametr relaksacji
tol7 = 1e-5;    % Tolerancja błędu
max_iter7 = 100;  % Maksymalna liczba iteracji

% Wywołanie funkcji dla klasycznej metody SOR
disp("Siódmy przypadek:")
[x7, liczba_iteracji7, macierz_iteracyjna7] = sor_iteration_matrix(A7, b7, omega7, tol7, max_iter7);
rozwiazanie7 = linsolve(Macierz7, b7');
wskaznik7 = norm(inv(Macierz7)) * norm(Macierz7);
blad_wzgledny7 = norm(rozwiazanie7 - x7)/norm(rozwiazanie7);
disp("Omega, dokładność, maksymalna liczba iteracji: ");
fprintf('omega = %.2f, tolerancja = %.1e, max liczba iteracji = %d\n', omega7, tol7, max_iter7);
disp("Błąd względny:")
disp(blad_wzgledny7)
disp("Liczba iterecji:")
disp(liczba_iteracji7)
disp("Promien spektralny macierzy iteracyjnej:")
disp(max(abs(eig(macierz_iteracyjna7))))
disp("Wskaznik uwarunkowania:")
disp(wskaznik7)
disp(" ")
disp(" ")
disp(" ") 


n = 100;

main_diag = 3 * ones(n, 1);     
upper_diag = 2 * ones(n-1, 1); 
lower_diag = -3 * ones(n-1, 1);  %

A = diag(main_diag) + diag(upper_diag, 1) + diag(lower_diag, -1);


% Wektory dla przekątnych
main_diag = 3 * ones(n, 1);           
upper_diag = [2 * ones(n-1, 1); 0];   
lower_diag = [0; -3 * ones(n-1, 1)];  

tridiagonal_form = [lower_diag, main_diag, upper_diag];

% Przygotowanie przykładowych danych
A8 = tridiagonal_form;  % Macierz przechowywana jak n x 3
Macierz8 = A;
b8 = ones(1, 100);  % Wektor wyrazów wolnych

omega8 = 1.35;  % Parametr relaksacji
tol8 = 1e-5;    % Tolerancja błędu
max_iter8 = 100;  % Maksymalna liczba iteracji

% Wywołanie funkcji dla klasycznej metody SOR
disp("Ósmy przypadek:")
[x8, liczba_iteracji8, macierz_iteracyjna8] = sor_iteration_matrix(A8, b8, omega8, tol8, max_iter8);
rozwiazanie8 = linsolve(Macierz8, b8');
wskaznik8 = norm(inv(Macierz8)) * norm(Macierz8);
blad_wzgledny8 = norm(rozwiazanie8 - x8)/norm(rozwiazanie8);
disp("Omega, dokładność, maksymalna liczba iteracji: ");
fprintf('omega = %.2f, tolerancja = %.1e, max liczba iteracji = %d\n', omega8, tol8, max_iter8);
disp("Błąd względny:")
disp(blad_wzgledny8)
disp("Liczba iterecji:")
disp(liczba_iteracji8)
disp("Promien spektralny macierzy iteracyjnej:")
disp(max(abs(eig(macierz_iteracyjna8))))
disp("Wskaznik uwarunkowania:")
disp(wskaznik8)
disp(" ")
disp(" ")
disp(" ") 