function Aufteilung(Bereich, Genauigkeit, Genauigkeit_simplex, p)
f = @(z) 1 ./ (z-1);
%f = @(z) 1./(z - (1 + 1j));
%f = @(z) 1 ./ ((z - (1 + 1j)) .* (z - (-2 + 2j)) .* (z - (3 - 1j)));
%f = @(z) 1 ./ ((z - (1 + 1j)).^3 .* (z - (-2 + 2j)).^2);
%f = @(z) z./(z+(1+1j));
%f = @(z) (1-z)./((1-z).*(z-(2+2j)));
%f = @(z) sqrt(z);
% f = @(z) sqrt((z - 10).^2 - 121);
%f = @(z) z;
x = Bereich(1);
y = Bereich(2);
width = Bereich(3);
height = Bereich(4);
p = p+1;
Bedingung = Cauchysch(Bereich, p);


Re = x + width/2;
Im = y + height/2;
% Prüfe, ob das Viereck klein genug ist, um nicht weiter zu unterteilen
if Bedingung == true
    if width <= Genauigkeit || height <= Genauigkeit
        
        Startpunkt = Re + 1i*Im;
        Polstelle = Polstellenbestimmung(f, Startpunkt, Genauigkeit_simplex);
        PolstelleRe = real(Polstelle);
        PolstelleIm = imag(Polstelle);
        fprintf('Singularitätsposition: %d + j%d \n', PolstelleRe, PolstelleIm)
        p
        return;
        
    end
end


if Bedingung == true
    % Viertelbereiche berechnen
    half_width = width / 2;
    half_height = height / 2;

    sub_Bereich = [
        x, y, half_width, half_height;
        x + half_width, y, half_width, half_height;
        x, y + half_height, half_width, half_height;
        x + half_width, y + half_height, half_width, half_height
        ];

    % Rekursiv jedes Viertel prüfen und ggf. weiter unterteilen
    for i = 1:4
        Aufteilung(sub_Bereich(i, :), Genauigkeit, Genauigkeit_simplex, p);
    end

end
end



