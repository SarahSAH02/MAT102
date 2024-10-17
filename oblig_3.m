load ('Arbeidskrav3.mat');

[n, m] = size(X1);

meanX = mean(X1); % gjennomsnitt
stdX = std(X1); % standardavvik

% Deloppgave a: Preprosesser datamatrisa av å trekke gjennomsnittet og dele på standardavvik
for j = 1:m
    X1(:,j) = X1(:,j) - meanX(j);
    if stdX(j) ~= 0
        X1(:,j) = X1(:,j) / stdX(j);
    end
end
X_pre = X1; % Den preprosesserte matrisen

disp(X_pre);
disp('Preprosessering fullført.');

% Deloppgave b: Forklaring på hvorfor det er viktig å preprosessere matrise
% Standarisering er viktig for å sammenligne ulike data fordi variabler kan ha ulike skalaer.
% Dersom vi ikke standardiserer, kan variabler med større skala dominere PCA resultatene.

% Deloppgave c: Finn de to første prinsipal komponentene ved NIPALS-metoden
a = 2; % Antall prinsipal komponenter
T = zeros(n, a); % Matriser for scores (T) og loadings (P)
P = zeros(m, a); 

% PCA ved hjelp av NIPALS-metoden
for i = 1:a
    t_old = rand(n, 1);  
    p = X_pre' * t_old; % Bruk den preprosesserte matrisen X_pre
    p = p / norm(p); 
    t = X_pre * p;
    
    while norm(t - t_old) > 1e-5
        p = X_pre' * t;
        p = p / norm(p);
        t_old = t;
        t = X_pre * p;
    end
    
    T(:, i) = t;  % Scores for komponent i
    P(:, i) = p;  % Loadings for komponent i
    
    X_pre = X_pre - t * p';  % Oppdaterer X_pre ved å fjerne komponent i
end

% Plot de to første komponentene
figure
scatter(T(:,1), T(:,2))  % Plotter de to første komponentene
text(T(:,1) + 0.01, T(:,2), objNames1)

% Deloppgave d: Kommenter på smakspanelets konsistens
% De to innstillingene med 5a og 5b er ganske nære på plottet, dermed kan vi si at smakspanelet er konsistente.

% Deloppgave e: Beregn total varians og andel forklart varians for PC1 og PC2
% Beregn total varians og forklart varians
total_variasjon = trace(X1' * X1); % Total varians i dataene
forklart_variasjon = trace(T' * T); % Forklart varians basert på scores

% Beregn andelen forklart varians
andel_forklart = forklart_variasjon / total_variasjon * 100; % I prosent

% Skriv ut resultatet
fprintf('Andel av variasjon forklart er %.2f%%\n', andel_forklart);

% Andel av variasjon forklart er 98.36%.

% Deloppgave f: Plot loadingene for de to første komponentene
figure
scatter(P(:,1), P(:,2))
text(P(:,1) + 0.01, P(:,2), varNames1)

%Deloppgave g
% matrise X2
[n2, m2] = size(X2); 

% Beregne gjennomsnitt og standardavvik
meanX2 = mean(X2); 
stdX2 = std(X2); 

% Preprosesser datamatrise 
for j = 1:
    X2(:, j) = X2(:, j) - meanX2(j); % Trekke gjennomsnittet
    if stdX2(j) ~= 0
        X2(:, j) = X2(:, j) / stdX2(j); % Deler på standardavviket
    end
end
X_pre2 = X2; 

a = 2; % Antall prinsipal komponenter
T2 = zeros(n2, a); % Matriser for scores (T)

% NIPALS algoritmen for PCA
for i = 1:a
    t_old = rand(n2, 1);  
    p = X_pre2' * t_old; % Initialisering
    p = p / norm(p); 
    t = X_pre2 * p;
    
    % Oppdatering
    while norm(t - t_old) > 1e-5
        p = X_pre2' * t; % Oppdatering
        p = p / norm(p);
        t_old = t;
        t = X_pre2 * p;
    end
    
    T2(:, i) = t; % Lagre scores for komponent i
    X_pre2 = X_pre2 - t * p'; % Fjerne komponenten fra dataene
end

figure
scatter(T2(:, 1), T2(:, 2), 'filled'); 
hold on;
