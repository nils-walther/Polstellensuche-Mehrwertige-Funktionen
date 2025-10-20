%% ============================   Main    =================================
% ==        Polstellensuche mit rekursiver Eingrenzung und Newton        ==
% =========================================================================

% Beispiele:
% f = @(z) 1 ./ ((z - (1 + 1i)) .* (z - (-2 + 2i)));
% g = @(z) (z - (1 + 1i)) .* (z - (-2 + 2i));
% g_ableitung = @(z) 2*z + 1 - 3i;

f = @(z) 1 ./ (z.^3 - 1);
g = @(z) z.^3 - 1;
g_ableitung = @(z) 3*z.^2;


x_start = -5;
y_start = -5;
breite_start = 10;
hoehe_start = 10;
Genauigkeit_Eingrenzung = 1e-3; % Mindestgröße des Quadrats, bevor Newton startet

Startbereich = [x_start, y_start, breite_start, hoehe_start];
p = 0; % Starttiefe der Rekursion

fprintf('Starte Polstellensuche...\n');
Eingrenzung(f, g, g_ableitung, Startbereich, Genauigkeit_Eingrenzung, p);
fprintf('Suche beendet.\n');


%% ------------------------------------------------------------------------
% Eingrenzung (Rekursive Funktion)
function Eingrenzung(f, g, g_ableitung, Bereich, Genauigkeit, p)
    p = p + 1; % Rekursionstiefe erhöhen
    
    pol_im_Bereich = Cauchysch(f, Bereich, p);
    
    if pol_im_Bereich
        breite = Bereich(3);
        hoehe = Bereich(4);
        
        % Prüfung Bereich klein genug für Newton
        if breite <= Genauigkeit || hoehe <= Genauigkeit
            % Startpunkt für Newton (Mittelpunkt des kleinen Quadrats)
            x_mitte = Bereich(1) + breite/2;
            y_mitte = Bereich(2) + hoehe/2;
            z_start_newton = x_mitte + 1i*y_mitte;
            
            Newton_Verfahren(g, g_ableitung, z_start_newton);
            return;
        end
        
        halbe_breite = breite / 2;
        halbe_hoehe = hoehe / 2;
        x = Bereich(1);
        y = Bereich(2);

        sub_Bereiche = [
            x, y, halbe_breite, halbe_hoehe; % Unten links
            x + halbe_breite, y, halbe_breite, halbe_hoehe; % Unten rechts
            x, y + halbe_hoehe, halbe_breite, halbe_hoehe; % Oben links
            x + halbe_breite, y + halbe_hoehe, halbe_breite, halbe_hoehe % Oben rechts
        ];

        for i = 1:4
            Eingrenzung(f, g, g_ableitung, sub_Bereiche(i, :), Genauigkeit, p);
        end
    end
end


%% ------------------------------------------------------------------------
% Cauchysch (Numerische Version)
function pol_vorhanden = Cauchysch(f, Startbereich, p)
    if p <= 3; Toleranz = 1e-16; else; Toleranz = 1e-3; end
    pol_vorhanden = false;
    num_points = 1000; % Anzahl der Punkte pro Kante

    a = Startbereich(1); b = Startbereich(2);
    c = Startbereich(3); d = Startbereich(4);

    x_kontur = [linspace(a, a+c, num_points), (a+c)*ones(1,num_points), linspace(a+c, a, num_points), a*ones(1,num_points)];
    y_kontur = [b*ones(1,num_points), linspace(b, b+d, num_points), (b+d)*ones(1,num_points), linspace(b+d, b, num_points)];
    z_kontur = x_kontur + 1i * y_kontur;

    dz = gradient(z_kontur);
    integrand = f(z_kontur) .* dz;
    integral_wert = trapz(integrand);

    if abs(integral_wert) > Toleranz
        pol_vorhanden = true;
    end
end


%% ------------------------------------------------------------------------
% Newton und Ausgabe
function Newton_Verfahren(g, g_ableitung, z_start)
    z_n = z_start;
    Genauigkeit_Newton = 1e-9;
    max_iterationen = 50; % Sicherheitsabbruch
    iteration = 0;
    
    while (abs(g(z_n))) > Genauigkeit_Newton && iteration < max_iterationen
        g_wert = g(z_n);
        g_ableitung_wert = g_ableitung(z_n);
        
        if abs(g_ableitung_wert) < 1e-12
            fprintf('Ableitung nahe Null, Newton bricht ab.\n');
            break;
        end
        
        z_n = z_n - (g_wert / g_ableitung_wert);
        iteration = iteration + 1;
    end
    
    fprintf('Singularität gefunden bei: %.6f + j(%.6f)\n', real(z_n), imag(z_n));
end