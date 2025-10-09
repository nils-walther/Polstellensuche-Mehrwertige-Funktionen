function result = Cauchysch(Startbereich, p)
if p <=3
Toleranz = 1e-16;  %Toleranz für Hüllintegral, kleine Werte ungleich Null können aufgrund ungenauigkeiten entstehen
else 
    Toleranz = 1e-3;
end
result = false;
num_points = 100000;

a = Startbereich(1);  % x-Koordinate der unteren linken Ecke
b = Startbereich(2);  % y-Koordinate der unteren linken Ecke
c = Startbereich(3);  % Breite des Rechtecks
d = Startbereich(4);  % Höhe des Rechtecks

% Unterkante (von links nach rechts)
x1 = linspace(a, a+c, num_points);
y1 = b * ones(1, num_points);

% Rechte Kante (von unten nach oben)
x2 = (a+c) * ones(1, num_points);
y2 = linspace(b, b+d, num_points);

% Oberkante (von rechts nach links)
x3 = linspace(a+c, a, num_points);
y3 = (b+d) * ones(1, num_points);

% Linke Kante (von oben nach unten)
x4 = a * ones(1, num_points);
y4 = linspace(b+d, b, num_points);

% Geschlossene Kurve z entlang der Bereichskantne
z = [x1 + 1i*y1, x2 + 1i*y2, x3 + 1i*y3, x4 + 1i*y4];


f = @(z) 1 ./ (z-1);
%f = @(z) z.^2 + 1;
%f = @(z) 1./(z - (1 + 1j));
%f = @(z) 1 ./ ((z - (1 + 1j)) .* (z - (-2 + 2j)) .* (z - (3 - 1j)));
%f = @(z) 1 ./ ((z - (1 + 1j)).^3 .* (z - (-2 + 2j)).^2);
%f = @(z) z./(z+(1+1j));
%f = @(z) (1-z)./((1-z).*(z-(2+2j)));
%f = @(z) sqrt(z);
%f = @(z) sqrt((z - 10).^2 - 121);
%f = @(z) z;

dz_real = gradient(real(z));
dz_imag = gradient(imag(z));
dz = dz_real + 1i * dz_imag;
integrand = f(z) .* dz;
integral = trapz(1:numel(z), integrand);

% Toleranzbedingung : Bei einfacher Umschließung von Pol- und
% Nullstellen gilt Hüllintegral einer != 0, genaue Berechnung über
% \frac{1}{2pi i}Hüllint\frac{f´(z)}{f(z)}dz = N(0) - N(inf)
if abs(integral) - Toleranz > 0
    result = true;
end

end
