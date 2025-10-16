%% ============================   Main    =================================
% ==        Polstellensuche mit Cauchysch und Newton-Verfahren           ==
% =========================================================================

% Beispiele
syms z;
f_symbol = 1 / ((z - (1+1i)) * (z - (-2+2i))); % symbolische Funktion f(z)
f = @(z) 1 ./ ((z - (1+1i)) .* (z - (-2+2i))); % f(z): Funktion
g = @(z) (z - (1+1i)) .* (z - (-2+2i)); % g(z): Nenner der Funktion
g_ableitung = @(z) 2*z + 1 - 3i; % Ableitung von g(z)

x_start = -5;
y_start = -5;
breite_start = 10;
hoehe_start = 10;
Startbereich = [x_start, y_start, breite_start, hoehe_start];


Eingrenzung(f, g, g_ableitung, Startbereich)



%% ------------------------------------------------------------------------
% Eingenzung

% in quadrate einteilen
% cauchysch aufrufen quadrate prüfen
% newton einsetzen für genaue ermittlung und ausgabe

function Eingrenzung(f_symbol, g, g_ableitung, Startbereich)

    x_start = Startbereich(1);
    y_start = Startbereich(2);
    breite_gesamt = Startbereich(3);
    hoehe_gesamt = Startbereich(4);

    anzahl_x = 2; % Für ein 2x2 Gitter -> 4 Quadrate
    anzahl_y = 2;
    
    schrittweite_x = breite_gesamt / anzahl_x;
    schrittweite_y = hoehe_gesamt / anzahl_y;



end


%% ------------------------------------------------------------------------
% Cauchysch
function pol_vorhanden = Cauchysch(f_symbol, Startbereich)

    syms t z
    a = Startbereich(1);
    b = Startbereich(2);
    c = Startbereich(3);
    d = Startbereich(4);


    gamma = [ (a+c*t)+1i*b, (a+c)+1i*(b+d*t), (a+c*(1-t))+1i*(b+d), a+1i*(b+d*(1-t)) ];

    integral = 0;
    for k = 1:4
        dz = diff(gamma(k), t);
        integral = integral + int(subs(f_symbol, z, gamma(k))*dz, t, 0, 1);
    end

    pol_vorhanden = double(integral) ~= 0;
end


%% ------------------------------------------------------------------------
% Newton und Ausgabe
function Newton_Verfahren(g, g_ableitung, z_start)
    z_n = z_start; % Punkt im quadrat welches aufgelöst worden ist
    
    Genauigkeit_Newton = 1e-9;
    while (abs(g(z_n)/g_ableitung(z_n))) > Genauigkeit_Newton
        g_wert = g(z_n);
        g_ableitung_wert = g_ableitung(z_n);
        z_n = z_n - (g_wert/g_ableitung_wert);
    end
    
    fprintf("Singularität bei %.1 und j%.1") real(z_n), imag(z_n);
end