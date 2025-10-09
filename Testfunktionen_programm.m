%%Einzelne Polstelle
syms z
f = 1/(z - (1 + 1j));  % Definiere die Funktion
poles = solve(1/ f == 0, z);  % Setze den Nenner gleich null und löse nach z auf
disp(poles)  % Ausgabe der Polstelle

%% mehrere, nicht identische Polstellen
syms z
f = 1 / ((z - (1 + 1j)) * (z - (-2 + 2j)) * (z - (3 - 1j)));  
denom = simplify(1 / f);  % Den Nenner extrahieren
poles = solve(denom == 0, z);  % Nenner = 0 setzen und nach z auflösen
disp(poles)  % Ausgabe der Polstellen

%% mehrere identische Polstellen
syms z
f = 1 / ((z - (1 + 1j))^3 * (z - (-2 + 2j))^2);  % Funktion mit mehrfachen Polen
poles = solve((z - (1 + 1j))^3 * (z - (-2 + 2j))^2 == 0, z);  % Nenner = 0 setzen
disp('Polstellen:');
disp(poles);


%% Integralbeispiel Residuum

r = 2; % gegebener Wert für r


f = @(phi) exp(r * exp(1j * phi)) .* exp(1j * phi) ./ (r^2 * exp(2j * phi) + 1);
I = 2j*integral(f, 0, 2*pi, 'ArrayValued', true);

disp(I);

%% Riemannsheets Wurzelfunktion

clc; clear; close all;

% Erstellen des Gitters in Polarkoordinaten
theta = linspace(-pi, pi, 200); % Winkel von -pi bis pi
r = linspace(0.1, 2, 100); % Radius (kein Nullwert wegen sqrt(0))
[Theta, R] = meshgrid(theta, r);

% Umwandlung in kartesische Koordinaten
X = R .* cos(Theta);
Y = R .* sin(Theta);

% Berechnung der Riemannschen Fläche für sqrt(z)
Z1 = sqrt(R) .* exp(1i * Theta / 2); % Erstes Blatt
Z2 = sqrt(R) .* exp(1i * (Theta / 2 + pi)); % Zweites Blatt

% Real- und Imaginärteil für die 3D-Darstellung
Z1_Re = real(Z1);
Z2_Re = real(Z2);

% Erstellen der 3D-Darstellung
figure;
hold on;
surf(X, Y, Z1_Re, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
surf(X, Y, Z2_Re, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
colormap cool;
xlabel('Re(z)');
ylabel('Im(z)');
zlabel('Re(sqrt(z))');
title('Riemannsche Fläche für sqrt(z)');

% Achsen und Ansicht anpassen
view(3);
grid on;
axis tight;
hold off;

saveas(gcf, 'C:\Users\KNECHTT\Desktop\Sta2\riemann_sqrt.png');