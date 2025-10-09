% Parameter des Rechtecks

a = 2;  % Breite des Rechtecks
b = 2;  % Höhe des Rechtecks
Genauigkeit = 1e-5; %min Größe des Endbereichs
Genauigkeit_simplex = 1e-10;
Startbereich = [-a-1e-7, -b-1e-7, 2*a, 2*b];
p = 0;

Aufteilung(Startbereich, Genauigkeit, Genauigkeit_simplex, p)  

