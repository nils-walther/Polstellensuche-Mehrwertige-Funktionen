%% ============================   Main    =================================
% ==        Polstellensuche mit rekursiver Eingrenzung und fsolve        ==
% =========================================================================
% --- Definitionen der Funktionen in der w-Ebene
f_w = @(w) 1 ./ (w.^4 - 16);
g_w = @(w) w.^4 - 16; % Nennerfunktion
% --- Parameter für die Suche
x_start = -5;
y_start = -5;
breite_start = 10;
hoehe_start = 10;
Genauigkeit_Eingrenzung = 1e-3;
Startbereich = [x_start, y_start, breite_start, hoehe_start];

gefundene_z_pole = [];

fprintf('Starte Polstellensuche...\n');
gefundene_z_pole = Eingrenzung(f_w, g_w, Startbereich, Genauigkeit_Eingrenzung, 0, gefundene_z_pole); 
fprintf('Suche beendet.\n\n');

% --- Finale Ausgabe der Ergebnisse
if ~isempty(gefundene_z_pole)
    % fprintf('--- Eindeutige Pole in der z-Ebene ---\n');
    % disp(gefundene_z_pole);
else
    fprintf('Keine Pole im Suchbereich gefunden.\n');
end


%% =======================   Lokale Funktionen   =========================
% -------------------------------------------------------------------------
function gefundene_z = Eingrenzung(f_w, g_w, Bereich, Genauigkeit, p, gefundene_z)
    p = p + 1; % Rekursionsstufe
    pol_im_Bereich = Cauchysch(f_w, Bereich, p);
    if p == 1
        pol_im_Bereich = true; % Erste Teilung erzwingen. Verhindert, dass sich symmetrische Pole im Integral gegenseitig aufheben und so fälschlicherweise ein leeres Gebiet signalisieren.
    end
    if pol_im_Bereich
        breite = Bereich(3);
        hoehe = Bereich(4);
        
        if breite <= Genauigkeit || hoehe <= Genauigkeit
            x_mitte = Bereich(1) + breite/2;
            y_mitte = Bereich(2) + hoehe/2;
            w_start = x_mitte + 1i*y_mitte;
            
            w_pole = finde_singularitaet_fsolve(g_w, w_start);
            z_pole = resubstitution(w_pole); % Resubstitution
            
            toleranz_doppelt = 1e-4;
            if isempty(gefundene_z) || ~any(abs(gefundene_z - z_pole) < toleranz_doppelt)
                gefundene_z = [gefundene_z, z_pole];
                fprintf('z-Pol: %+.6f %+.6fi\n', real(z_pole), imag(z_pole));
            end
            return;
        end
        
        halbe_breite = breite / 2;
        halbe_hoehe = hoehe / 2;
        x = Bereich(1);
        y = Bereich(2);
        ungenauigkeit = 1e-9;
        sub_Bereiche = [
            x + ungenauigkeit, y + ungenauigkeit, halbe_breite, halbe_hoehe;  % Unten links
            x + halbe_breite, y + ungenauigkeit, halbe_breite, halbe_hoehe;  % Unten rechts
            x + ungenauigkeit, y + halbe_hoehe, halbe_breite, halbe_hoehe;  % Oben links
            x + halbe_breite, y + halbe_hoehe, halbe_breite, halbe_hoehe % Oben rechts
        ];
        for i = 1:4
            gefundene_z = Eingrenzung(f_w, g_w, sub_Bereiche(i, :), Genauigkeit, p, gefundene_z);
        end
    end
end


% -------------------------------------------------------------------------
function w_n = finde_singularitaet_fsolve(g_w_komplex, w_start_komplex)
    fun = @(v) [real(g_w_komplex(v(1) + 1i*v(2))); 
                imag(g_w_komplex(v(1) + 1i*v(2)))];
    w_start_xy = [real(w_start_komplex), imag(w_start_komplex)];
    options = optimoptions('fsolve', 'Display', 'none', 'FunctionTolerance', 1e-9); % Optionen (verhindert, dass fsolve seine eigene Ausgabe anzeigt)
    w_n_xy = fsolve(fun, w_start_xy, options);
    w_n = w_n_xy(1) + 1i * w_n_xy(2);
end


% -------------------------------------------------------------------------
% Cauchysch (Numerische Variante des Cauchy-Integrals)
function pol_vorhanden = Cauchysch(f, Startbereich, p)
    if p <= 3; Toleranz = 1e-19; else; Toleranz = 1e-3; end
    
    num_points = 800; % Anzahl der Punkte pro Kante
    a = Startbereich(1); b = Startbereich(2);
    c = Startbereich(3); d = Startbereich(4);
    
    % Konturweg definieren
    weg1 = linspace(a, a+c, num_points);
    weg2 = linspace(b, b+d, num_points);
    x_kontur = [weg1, (a+c)*ones(1,num_points), fliplr(weg1), a*ones(1,num_points)];
    y_kontur = [b*ones(1,num_points), weg2, (b+d)*ones(1,num_points), fliplr(weg2)];
    z_kontur = x_kontur + 1i * y_kontur;
    
    % Numerische Integration
    dz = diff([z_kontur z_kontur(1)]);
    integral_wert = sum(f(z_kontur) .* dz);
    pol_vorhanden = abs(integral_wert) > Toleranz;
end


% -------------------------------------------------------------------------
% Resubstitution
function z_pole = resubstitution(w_pole)
    z_pole = w_pole.^2;
end