function Polstelle = Polstellenbestimmung(f, Startpunkt, tolerance)

    % Hilfsfunktion zur Optimierung (Betrag der inversen Funktion)
    g = @(z) abs(1 ./ (f(z) + eps));  % eps zur Vermeidung von Division durch Null

    % Optimierungsoptionen mit stabilerem Verhalten
    options = optimset('Display', 'iter', 'TolX', tolerance, 'TolFun', tolerance);

    try
        % Numerische Suche nach der Polstelle direkt in der komplexen Ebene
        Polstelle = fminsearch(@(z) g(z(1) + 1i * z(2)), [real(Startpunkt), imag(Startpunkt)], options);
        
        % Umwandlung in komplexe Zahl
        Polstelle = Polstelle(1) + 1i * Polstelle(2);

        % Validierung des Ergebnisses
        if abs(f(Polstelle)) > 10 * tolerance
            warning('Erkannte Polstelle ist möglicherweise ungenau.');
        end
    catch
        warning('Optimierung fehlgeschlagen. Rückgabe von NaN.');
        Polstelle = NaN;
    end

end
