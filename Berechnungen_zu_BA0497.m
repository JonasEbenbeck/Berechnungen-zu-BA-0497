
%                                  1. Vorbereitende Berechnungen und Festlegungen

% ____________________________________________________________________________________________________________________

%                                               Parametereingabe
                          
% --------------------------------------------------------------------------------------------------------------------
% Geometrische Größen der Schneckenwelle (Welle 1):

% Welle und Hülsen:
d_I1 = 0.026;                   % Außendurchmesser eines Absatzes von Hülse 1 (siehe Bild 3.7) [m]
d_II1 = 0.020;                  % Außendurchmesser eines Absatzes von Hülse 1 (siehe Bild 3.7) [m]
d_III1 = 0.013;                 % Außendurchmesser der Hülse 2 (siehe Bild 3.7) [m]
d_IV1 = 0.010;                  % Durchmesser des Wellenabsatzes (siehe Bild 3.7) [m]
d_V1 = 0.038;                   % Durchmesser des Wellenabsatzes (siehe Bild 3.7) [m]
d_VI1 = 0.020;                  % Durchmesser des Wellenabsatzes (siehe Bild 3.7) [m]
d_VII1 = 0.015;                 % Durchmesser des Wellenabsatzes (siehe Bild 3.7) [m]
l_I1 = 0.008;                   % Länge eines Absatzes von Hülse 1 (siehe Bild 3.7) [m]
l_II1 = 0.015;                  % Länge eines Absatzes von Hülse 1 (siehe Bild 3.7) [m]
l_III1 = 0.040;                 % Gesamtlänge der Schnecke (siehe Bild 3.7) [m]
l_IV1 = 0.015;                  % Angenommene Länge des Eingriffsbereichs der Verzahnung (siehe Bild 3.7) [m]
l_V1 = 0.011;                   % Angenäherte Länge des Wellenabsatzes (siehe Bild 3.7) [m]
l_VI1 = 0.0335;                 % Hilfslänge (siehe Bild 3.7) [m]
l_VII1 = 0.0225;                % Hilfslänge (siehe Bild 3.7) [m]
l_VIII1 = 0.043;                % Halbe Länge des in der Umgebung rotierenden Zylinders (analog zu l_VI2 in Bild 3.10) [m]
A_51 = 0.00897;                 % Wärmeübertragende Fläche zur Berechnung des Wärmeübergangswiderstands zwischen k19 und Umgebung (k0) (analog zu A_62 in Bild 3.10) [m^2]

% Schnecke:
d_a1 = 0.024635;                % Kopfkreisdurchmesser der Schnecke [m]
d_f1 = 0.019010;                % Fußkreisdurchmesser der Schnecke [m] 
s_n1 = 0.002743;                % Zahndicke am Teilkreis (Schnecke) [m]
s_a1 = 0.00073;                 % Zahndicke am Kopfkreis (Schnecke) [m]
z_1 = 2;                        % Anzahl der Zähne der Schnecke [-]
m_n1 = 0.00125;                 % Normalmodul der Schnecke [m]

% Lager C:
d_C = 0.017;                    % Lagerinnendurchmesser für Lager C [m]
D_C = 0.035;                    % Lageraußendurchmesser für Lager C [m]
B_C = 0.010;                    % Lagerbreite für Lager C [m]

% Lager D1:
d_D1 = 0.020;                   % Lagerinnendurchmesser für Lager D1 [m]
D_D1 = 0.047;                   % Lageraußendurchmesser für Lager D1 [m]
B_D1 = 0.014;                   % Lagerbreite für Lager D1 [m]

% Lager D2:
d_D2 = 0.020;                   % Lagerinnendurchmesser für Lager D2 [m]
D_D2 = 0.047;                   % Lageraußendurchmesser für Lager D2 [m]
B_D2 = 0.014;                   % Lagerbreite für Lager D2 [m]

% --------------------------------------------------------------------------------------------------------------------
% Geometrische Größen der Schraubradwelle (Welle 2):

% Welle, Hülsen und Passstifte:
d_I2 = 0.017;                   % Durchmesser des Wellenabsatzes (siehe Bild 3.9) [m]
d_III2 = 0.036;                 % Durchmesser des Wellenabsatzes (siehe Bild 3.9) [m]
d_IV2 = 0.0252;                 % Außendurchmesser der Hülse (siehe Bild 3.9) [m]
l_I2 = 0.018;                   % Länge der Hülse (siehe Bild 3.9) [m]
l_II2 = 0.015;                  % Länge des Wellenabsatzes (siehe Bild 3.9) [m]
l_III2 = 0.0055;                % Länge des Wellenabsatzes (siehe Bild 3.9) [m]
l_IV2 = 0.029;                  % Hilfslänge (siehe Bild 3.9) [m]
l_V2 = 0.028;                   % Hilfslänge (siehe Bild 3.10) [m]
l_VI2 = 0.042;                  % Hilfslänge (siehe Bild 3.10) [m]
d_P2 = 0.004;                   % Durchmesser Passstifte [m]
l_P2 = 0.025;                   % Länge der Passstifte [m]
A_62 = 0.010746;                % Wärmeübertragende Fläche zur Berechnung des Wärmeübergangswiderstands zwischen k4 und Umgebung (k0) (siehe Bild 3.10) [m^2]

% Schraubrad:
d_II2 = 0.0285;                 % Lochkreisdurchmesser der Passstiftbohrungen (siehe Bild 3.8) [m]
d_V2 = 0.022;                   % Bohrungsdurchmesser des Schraubrads (siehe Bild 3.8) [m]       
d_VIII2 = 0.035;                % Hilfsdurchmesser (siehe Bild 3.8) [m]
d_IX2 = 0.042;                  % Hilfsdurchmesser (siehe Bild 3.8) [m]
d_Knoten = 0.0385;              % Durchmesser des Knotenpunkts Radkörper (siehe Bild 3.8) [m]
d_f2 = 0.048590;                % Fußkreisdurchmesser des Schraubrads (siehe Bild 3.8) [m]
d_a2 = 0.054215;                % Kopfkreisdurchmesser des Schraubrads (siehe Bild 3.8) [m]
b_R2 = 0.010;                   % Breite des Schraubrads (siehe Bild 3.9) [m]
b_S2 = 0.004;                   % Stegbreite des Schraubrads [m]
s_n2 = 0.002168;                % Zahndicke am Teilkreis (Schraubrad) [m]
A_Zahnstirn2 = 4.53*10^(-6);    % Zahnstirnfläche einer Seite eines Zahns am Schraubrad [m^2]
A_Zahnkopf2 = 4.53*10^(-6);     % Zahnkopffläche eines Zahns am Schraubrad [m^2]
z_2 = 40;                       % Anzahl der Zähne des Schraubrads [-]


% Lager A:
d_A = 0.017;                    % Lagerinnendurchmesser für Lager A [m]
D_A = 0.040;                    % Lageraußendurchmesser für Lager A [m]
B_A = 0.012;                    % Lagerbreite für Lager A [m]

% Lager B1:
d_B1 = 0.017;                   % Lagerinnendurchmesser für Lager B1 [m]
D_B1 = 0.040;                   % Lageraußendurchmesser für Lager B1 [m]
B_B1 = 0.012;                   % Lagerbreite für Lager B1 [m]

% Lager B2:
d_B2 = 0.017;                   % Lagerinnendurchmesser für Lager B2 [m]
D_B2 = 0.040;                   % Lageraußendurchmesser für Lager B2 [m]
B_B2 = 0.012;                   % Lagerbreite für Lager B2 [m]

% --------------------------------------------------------------------------------------------------------------------
% Geometrische Größen des Gehäusemodells und der Grundplatten:

b_Gehaeuse = 0.099;             % Breite des Gehäusemodells [m]
l_Gehaeuse = 0.138;             % Länge des Gehäusemodells [m]
h_Gehaeuse = 0.250;             % Höhe des Gehäusemodells [m]
d_PlatteI = 0.030;              % Dicke der oberen Platte (Winkelverstellung) [m]
d_PlatteII = 0.100;             % Dicke der unteren Platte (Grundplatte) [m]
l_II0 = 0.250;                  % Hilfslänge zur Bestimmung von A_40 (siehe Bild 3.13) [m]

% --------------------------------------------------------------------------------------------------------------------
% Geometrische Größen zur Bestimmung der dimensionslosen Kennzahlen:

L_K12 = 0.16666;                % Charakteristische Länge zur Bestimmung von Re_K12 [m]

% --------------------------------------------------------------------------------------------------------------------
% Weitere Kenngrößen:

epsilon_n = 1.5;                % Überdeckung [-]
a = 0.0369;                     % Achsabstand [m]
i = 20;                         % Übersetzung [-]
n_1 = 12000;                    % Drehzahl Schneckenwelle [U/min] (Nicht größer als 12000 U/min)
n_2 = n_1/i;                    % Drehzahl Schraubradwelle [U/min]
VPunkt_OelA = 0.25;             % Eingespritzter Ölvolumenstrom bei Lager A [l/min]
VPunkt_OelB1 = 0.125;           % Eingespritzter Ölvolumenstrom bei Lager B1 [l/min]
VPunkt_OelB2 = 0.125;           % Eingespritzter Ölvolumenstrom bei Lager B2 [l/min]
VPunkt_OelC = 0.25;             % Eingespritzter Ölvolumenstrom bei Lager C [l/min]
VPunkt_OelD1 = 0.125;           % Eingespritzter Ölvolumenstrom bei Lager D1 [l/min]
VPunkt_OelD2 = 0.125;           % Eingespritzter Ölvolumenstrom bei Lager D2 [l/min]
A_innen = 0.037990;             % Innere Gehäuseoberfläche ohne der Fläche der Scheibe [m^2]
VPunkt_zuVZ = 1.5/60000;        % Insgesamt zur Verzahnung eingespritzter Ölvolumenstrom [m^3/s]
VPunkt_zuRF = 0.01/60000;       % Eingespritzter Ölvolumenstrom, der mit der Flanke des Schraubrads wechselwirkt [m^3/s]
VPunkt_zuSF = 0.01/60000;       % Eingespritzter Ölvolumenstrom, der mit der Flanke der Schnecke wechselwirkt [m^3/s]
T_0 = 20;                       % Umgebungstemperatur [°C]
T_zu = 80;                      % Temperatur des zur Verzahnung eingespritzten Ölvolumenstroms (vor sämtlichen Wechselwirkungen) [°C]

% --------------------------------------------------------------------------------------------------------------------
% Näherungsweise temperaturunabhängige Stoffwerte:

lambda_Schraubrad = 0.95;       % Wärmeleitfähigkeit des Schraubradwerkstoffs [W/mK] [Vic24, S.2]
lambda_Schnecke = 0.95;         % Wärmeleitfähigkeit des Schneckenwerkstoffs [W/mK] [Vic24, S.2]
lambda_Huelse5 = 50;            % Wärmeleitfähigkeit der Hülse 5 (siehe Abbildung) [W/mK] [CNC21, S.1]
lambda_Huelse2 = 50;            % Wärmeleitfähigkeit der Hülse 2 (siehe Abbildung) [W/mK] [CNC21, S.1]
lambda_Passstift = 50;          % Wärmeleitfähigkeit des Passstiftwerkstoffs [W/mK] (Annahme)
lambda_Welle1 = 39.8;           % Wärmeleitfähigkeit des Schneckenwellenwerkstoffs [W/mK] [thy19b, S.1]
lambda_Welle2 = 39.8;           % Wärmeleitfähigkeit des Schraubradwellenwerkstoffs [W/mK] [thy19b, S.1]
lambda_A = 33;                  % Wärmeleitfähigkeit des Lagerwerkstoffs von Lager A [W/mK] [thy19a, S.1]
lambda_B1 = 33;                 % Wärmeleitfähigkeit des Lagerwerkstoffs von Lager B1 [W/mK] [thy19a, S.1]
lambda_B2 = 33;                 % Wärmeleitfähigkeit des Lagerwerkstoffs von Lager B2 [W/mK] [thy19a, S.1]
lambda_C = 33;                  % Wärmeleitfähigkeit des Lagerwerkstoffs von Lager C [W/mK] [thy19a, S.1]
lambda_D1 = 33;                 % Wärmeleitfähigkeit des Lagerwerkstoffs von Lager D1 [W/mK] [thy19a, S.1]
lambda_D2 = 33;                 % Wärmeleitfähigkeit des Lagerwerkstoffs von Lager D2 [W/mK] [thy19a, S.1]
lambda_Oel = 0.125;             % Wärmeleitfähigkeit des Öls (bei 100°C angenommen) [W/mK] [She06, S.1]
lambda_Luft = 0.02808;          % Wärmeleitfähigkeit von trockener Luft bei 1bar und 50°C [W/mK] [Ver13, S.197]
epsilon_Stahl = 0.52;           % Emissionsgrad von poliertem Stahlguss (Annahme) [Ver13, S.1087]
lambda_PlatteI = 50;            % Wärmeleitfähigkeit der oberen Platte (Winkelverstellung) [W/mK] (Annahme)
lambda_PlatteII = 50;           % Wärmeleitfähigkeit der unteren Platte (Grundplatte) [W/mK] (Annahme)
Pr_Luft = 0.7045;               % Prandtl-Zahl von trockener Luft bei 1bar und 50°C [-] [Ver13, S.197]

% --------------------------------------------------------------------------------------------------------------------
% Temperaturabhängige Stoffwerte:

% Prandtl-Zahl von Öl [-] als Funktion der Temperatur in °C [She06, S.1]
% Gültigkeitsbereich: -20°C bis 200°C
p_PrOel = polyfit([-20 0 20 40 100 125 150 175 200 225],[2225 562 215 111 32 21 17 14 12 11],9);
Pr_Oel = @(x) p_PrOel(1)*x.^9 + p_PrOel(2)*x.^8 + p_PrOel(3)*x.^7 + p_PrOel(4)*x.^6 + p_PrOel(5)*x.^5 + p_PrOel(6)*x.^4 + p_PrOel(7)*x.^3 + p_PrOel(8)*x.^2 + p_PrOel(9).*x + p_PrOel(10);

% Kinematische Viskosität von trockener Luft [m^2/s] als Funktion der Temperatur in °C [Ver13, S.197]
% Gültigkeitsbereich 0°C bis 100°C
p_nyLuft = polyfit([0 10 20 30 40 50 60 70 80 90 100],[135 144 153.2 162.6 172.3 182.2 192.2 202.5 213 223.7 234.6]*10^(-7),9);
ny_Luft = @(x) p_nyLuft(1)*x.^9 + p_nyLuft(2)*x.^8 + p_nyLuft(3)*x.^7 + p_nyLuft(4)*x.^6 + p_nyLuft(5)*x.^5 + p_nyLuft(6)*x.^4 + p_nyLuft(7)*x.^3 + p_nyLuft(8)*x.^2 + p_nyLuft(9).*x + p_nyLuft(10);

% Kinematische Viskosität von Öl [m^2/s] als Funktion der Temperatur in °C [She06, S.1]
% Gültigkeitsbereich 0°C bis 150°C
p_nyOel = polyfit([0 10 20 30 40 60 80 100 125 150 170 180 190 200 225 250],[47.4 25.8 17.28 11.2 8.36 4.83 3.11 2.19 1.51 1.14 0.91 0.81 0.77 0.73 0.59 0.52]*10^(-6),9);
ny_Oel = @(x) p_nyOel(1)*x.^9 + p_nyOel(2)*x.^8 + p_nyOel(3)*x.^7 + p_nyOel(4)*x.^6 + p_nyOel(5)*x.^5 + p_nyOel(6)*x.^4 + p_nyOel(7)*x.^3 + p_nyOel(8)*x.^2 + p_nyOel(9).*x + p_nyOel(10);

% Spezifische Wärmekapazität von Öl [J/kgK] als Funktion der Temperatur in °C [She06, S.1]
% Gültigkeitsbereich 0°C bis 200°C
p_cOel = polyfit([0 20 40 60 80 100 125 150 175 200],[1781 1853 1925 1995 2073 2140 2230 2319 2409 2499],3);
c_Oel = @(x) p_cOel(1)*x.^3 + p_cOel(2)*x.^2 + p_cOel(3)*x + p_cOel(4);

% Dichte von Öl [kg/m^3] als Funktion der Temperatur in °C [She06, S.1]
% Gültigkeitsbereich 0°C bis 200°C
p_rohOel = polyfit([0 20 40 60 80 100 125 150 175 200],[903 890 877 864 851 838 822 805 789 773],3);
roh_Oel = @(x) p_rohOel(1)*x.^3 + p_rohOel(2)*x.^2 + p_rohOel(3)*x +p_rohOel(4);

% Temperaturleitfähigkeit von Öl [m^2/s] als Funktion der Temperatur in °C
% Gütigkeitsbereich 0°C bis 200°C
a_Oel = @(x) lambda_Oel./(roh_Oel(x).*c_Oel(x));

% Thermischer Ausdehnungskoeffizient von Luft [1/K] (modelliert als ideales Gas) als Funktion der Temperatur in °C [Pol09, S.288]
% Gültigkeitsbereich 0°C bis 200°C
beta_Luft = @(x) 1/(x+273.15);

% ____________________________________________________________________________________________________________________

%                                 Temperaturunabhängige thermische Widerstände

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L1 [Suc12, S.60-62]:

alpha_E = 2*acos((2*a-d_a2)/(d_a1));
A_Ersatz = (d_a1^2)/8*(alpha_E-sin(alpha_E));
R_ZR = s_n2/(2*lambda_Schraubrad*A_Ersatz);
R_L1 = R_ZR/epsilon_n;
G_L1 = 1/R_L1;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L2 [Suc12, S.60-62]:

R_ZS = s_n1/(2*lambda_Schnecke*A_Ersatz);
R_L2 = R_ZS/epsilon_n;
G_L2 = 1/R_L2;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L3, R_L4 und R_L18 [Gei14, S.69-71; Suc12, S.60]:

R_12 = ((d_a2-d_f2)/4)/(lambda_Schraubrad*s_n2*b_R2);
R_12 = R_12/epsilon_n;
R_22 = log(d_f2/d_IX2)/(2*pi*lambda_Schraubrad*b_R2);
R_32 = log(d_IX2/d_Knoten)/(2*pi*lambda_Schraubrad*b_S2);
R_42 = log(d_Knoten/d_VIII2)/(2*pi*lambda_Schraubrad*b_S2);
R_52 = log(d_VIII2/d_II2)/(2*pi*lambda_Schraubrad*b_R2);
R_62 = log(d_II2/d_V2)/(2*pi*lambda_Schraubrad*b_R2);
R_72 = log(d_V2/d_I2)/(2*pi*lambda_Huelse5*b_R2);

R_L3 = R_12+R_22+R_32;
G_L3 = 1/R_L3;

R_L4 = R_62+R_72;
G_L4 = 1/R_L4;

R_L18 = R_42+R_52;
G_L18 = 1/R_L18;

% Hinweis zur Notation: Bei den Widerständen R_12 bis R_72 handelt es sich um Hilfswiderstände, die zur Berechnung 
% von R_L3, R_L4 und R_L18 verwendet wurden. Der Index der Hilfswiderstände setzt sich aus zwei Zahlen zusammen. 
% Die zweite Zahl gibt an, ob der Hilfswiderstand zur Berechnung eines Widerstands der Schraubradwelle samt 
% Schraubrad (2) oder der Schneckenwelle inklusive Schnecke (1) verwendet wird. Die erste Zahl des Index ist lediglich
% eine fortlaufende Nummerierung. So ist der Hilfswiderstand R_52 der fünfte Hilfswiderstand, der zur Berechnung
% eines Widerstands der Schraubradwelle oder des Schraubrads verwendet wird. Gleiche Indizierungslogik wird auch 
% bei den im Verlauf des Codes vorkommenden Hilfsflächen (z.B. A_12) verwendet.

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L5 [Pol09, S.67]:

A_Pges = 10*pi*(0.5*d_P2)^2;                            % Summe der Querschittsflächen aller Passstifte im Schraubrad
R_L5 = (0.5*(b_R2+l_II2))/(lambda_Passstift*A_Pges);
G_L5 = 1/R_L5;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L6 [Gei14, S.71]:

R_L6 = log(d_II2/d_I2)/(2*pi*lambda_Welle2*l_II2);
G_L6 = 1/R_L6;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L7 [Suc12, S.54]:

R_L7 = (0.5*(b_R2+l_II2))/(lambda_Welle2*pi*(0.5*d_I2)^2);
G_L7 = 1/R_L7;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L8 [Suc12, S.54]:

R_82 = (0.5*b_R2)/(lambda_Welle2*pi*(0.5*d_I2)^2);
R_92 = (0.5*l_I2)/(lambda_Welle2*pi*(0.5*d_IV2)^2);

R_L8 = R_82+R_92;
G_L8 = 1/R_L8;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L9 [Gei14, S.72 ;Suc12, S.54]:

R_102 = (0.5*l_I2)/(lambda_Welle2*pi*(0.5*d_IV2)^2);
R_112 = (0.5*B_A)/(lambda_Welle2*pi*(0.5*d_I2)^2);

R_L9 = R_102+R_112;
G_L9 = 1/R_L9;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L10 und R_L11 [Gei14 S.76-77]:

d_mA = (d_A+D_A)/2;

R_L10 = log(d_mA/d_A)/(2*pi*lambda_A*B_A);
G_L10 = 1/R_L10;

R_L11 = log(D_A/d_mA)/(2*pi*lambda_A*B_A);
R_L11 = 5*R_L11;                                   % Berücksichtigung der eigentlichen Inhomogenität der Gehäusetemperatur [Gei14, S.76]
G_L11 = 1/R_L11;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L12 und R_L13 [Suc12, S.54]:

R_L12 = (0.5*l_II2+l_III2+0.5*B_B1)/(lambda_Welle2*pi*(0.5*d_I2)^2);
G_L12 = 1/R_L12;

R_L13 = (0.5*B_B1+(l_IV2-B_B1-B_B2)+0.5*B_B2)/(lambda_Welle2*pi*(0.5*d_I2)^2);
G_L13 = 1/R_L13;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L14 und R_L16 [Gei14, S.76-77]:

d_mB1 = (d_B1+D_B1)/2;

R_L14 = log(d_mB1/d_B1)/(2*pi*lambda_B1*B_B1);
G_L14 = 1/R_L14;

R_L16 = log(D_B1/d_mB1)/(2*pi*lambda_B1*B_B1);  
R_L16 = 5*R_L16;                                   % Berücksichtigung der eigentlichen Inhomogenität der Gehäusetemperatur [Gei14, S.76]
G_L16 = 1/R_L16;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L15 und R_L17 [Gei14, S.76-77]:

d_mB2 = (d_B2+D_B2)/2;

R_L15 = log(d_mB2/d_B2)/(2*pi*lambda_B2*B_B2);
G_L15 = 1/R_L15;

R_L17 = log(D_B2/d_mB2)/(2*pi*lambda_B2*B_B2);
R_L17 = 5*R_L17;                                   % Berücksichtigung der eigentlichen Inhomogenität der Gehäusetemperatur [Gei14, S.76]
G_L17 = 1/R_L17;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L19 [Gei14, S.69-71; Suc12, S.60]:

R_101 = ((d_a1-d_f1)*0.25)/(lambda_Schnecke*s_n1*b_R2);
R_101 = R_101/epsilon_n;
R_111 = log(d_f1/d_III1)/(2*pi*lambda_Schnecke*l_IV1);
R_141 = log(d_III1/d_IV1)/(2*pi*lambda_Huelse2*l_IV1);

R_L19 = R_101+R_111+R_141;
G_L19 = 1/R_L19;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L20 [Suc12, S.54]:

R_11 = (l_III1-0.5*l_IV1)/(lambda_Welle1*pi*(0.5*d_IV1)^2);
R_21 = (0.5*l_V1)/(lambda_Welle1*pi*(0.5*d_V1)^2);

R_L20 = R_11+R_21;
G_L20 = 1/R_L20;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L21 [Suc12, S.54]:

R_31 = (0.5*l_V1)/(lambda_Welle1*pi*(0.5*d_V1)^2);
R_41 = (0.5*B_D1)/(lambda_Welle1*pi*(0.5*d_VI1)^2);

R_L21 = R_31+R_41;
G_L21 = 1/R_L21; 

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L22 [Suc12, S.54]:

R_L22 = (l_VI1-0.5*B_D1-0.5*B_D2)/(lambda_Welle1*pi*(0.5*d_VI1)^2);
G_L22 = 1/R_L22;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L24 und R_L26 [Gei14, S.76-77]:

d_mD1 = (d_D1+D_D1)/2;

R_L24 = log(d_mD1/d_D1)/(2*pi*lambda_D1*B_D1);
G_L24 = 1/R_L24;

R_L26 = log(D_D1/d_mD1)/(2*pi*lambda_D1*B_D1);
R_L26 = 5*R_L26;                                   % Berücksichtigung der eigentlichen Inhomogenität der Gehäusetemperatur [Gei14, S.76]
G_L26 = 1/R_L26;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L25 und R_L27 [Gei14, S.76-77]:

d_mD2 = (d_D2+D_D2)/2;

R_L25 = log(d_mD2/d_D2)/(2*pi*lambda_D2*B_D2);
G_L25 = 1/R_L25;

R_L27 = log(D_D2/d_mD2)/(2*pi*lambda_D2*B_D2);
R_L27 = 5*R_L27;                                   % Berücksichtigung der eigentlichen Inhomogenität der Gehäusetemperatur [Gei14, S.76]
G_L27 = 1/R_L27;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L28 [Suc12, S.54]:

R_51 = (0.5*l_IV1)/(lambda_Welle1*pi*(0.5*d_IV1)^2);
R_61 = (0.5*(l_I1+l_II1))/(lambda_Welle1*pi*(0.5*d_II1)^2);

R_L28 = R_51+R_61;
G_L28 = 1/R_L28; 

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L29 [Suc12, S.54]:

R_71 = (0.5*(l_I1+l_II1)-l_I1)/(lambda_Welle1*pi*(0.5*d_II1)^2);
R_81 = l_I1/(lambda_Welle1*pi*(0.5*d_I1)^2);
R_91 = (0.5*B_C)/(lambda_Welle1*pi*(0.5*d_C)^2);

R_L29 = R_71+R_81+R_91;
G_L29 = 1/R_L29;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_L30 und R_L31 [Gei14, S.76-77]:

d_mC = (d_C+D_C)/2;

R_L30 = log(d_mC/d_C)/(2*pi*lambda_C*B_C);
G_L30 = 1/R_L30;

R_L31 = log(D_C/d_mC)/(2*pi*lambda_C*B_C);
R_L31 = 5*R_L31;                                  % Berücksichtigung der eigentlichen Inhomogenität der Gehäusetemperatur [Gei14, S.76]
G_L31 = 1/R_L31;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K1 [Gei14, S.74]:

A_12 = 2*pi*(0.5*d_IV2)*l_I2;
alpha_K1 = 0.587*n_2^0.7*(d_IV2*1000)^0.4;

R_K1 = 1/(A_12*alpha_K1);
G_K1 = 1/R_K1;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K2 [Gei14, S.75-76; Sch20, S.86-87]:

R_K2 = 1/(28.6*VPunkt_OelA);
G_K2 = 1/R_K2;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K4 [Gei14, S.75-76; Sch20, S.86-87]:

R_K4 = 1/(28.6*VPunkt_OelB1);
G_K4 = 1/R_K4;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K5 [Gei14, S.75-76; Sch20, S.86-87]:

R_K5 = 1/(28.6*VPunkt_OelB2);
G_K5 = 1/R_K5;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K6 [Gei14, S.74]:

A_22 = 2*pi*(0.5*d_III2)*l_II2;
alpha_K6 = 0.587*n_2^0.7*(d_III2*1000)^0.4;

R_K6 = 1/(A_22*alpha_K6);
G_K6 = 1/R_K6;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K7 [Gei14, S.74]:

A_11 = 2*pi*(0.5*d_I1)*l_I1;
alpha_121 = 0.587*n_1^0.7*(d_I1*1000)^0.4;
R_121 = 1/(A_11*alpha_121);

A_21 = 2*pi*(0.5*d_II1)*l_II1;
alpha_131 = 0.587*n_1^0.7*(d_II1*1000)^0.4;
R_131 = 1/(A_21*alpha_131);

R_K7 = 1/(1/R_121+1/R_131);
G_K7 = 1/R_K7;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K8 [Gei14, S.74]:

A_31 = 2*pi*(0.5*d_V1)*l_V1;
alpha_K8 = 0.587*n_1^0.7*(d_V1*1000)^0.4;

R_K8 = 1/(A_31*alpha_K8);
G_K8 = 1/R_K8;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K9 [Gei14, S.75-76; Sch20, S.86-87]:

R_K9 = 1/(28.6*VPunkt_OelD1);
G_K9 = 1/R_K9;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K10 [Gei14, S.75-76; Sch20, S.86-87]:

R_K10 = 1/(28.6*VPunkt_OelD2);
G_K10 = 1/R_K10;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K11 [Gei14, S.75-76; Sch20, S.86-87]:

R_K11 = 1/(28.6*VPunkt_OelC);
G_K11 = 1/R_K11;

% ____________________________________________________________________________________________________________________

%                                 Weitere vorbereitende Definitionen und Rechnungen

% --------------------------------------------------------------------------------------------------------------------
% Startvektor der Potentialdifferenzen (Temperaturdifferenzen) mit von Null verschiednenen Einträgen ausschließlich
% an den Stellen, die zur Berechnung der temperaturabhängigen Widerstanden für die erste Iteration benötigt werden:

DeltaT = zeros(29,1);
DeltaT(1) = 40;                  % Startwert der Potentialdifferenz zwischen Gehäuse und Umgebung
DeltaT(2) = 50;                  % Stratwert der Potentialdiffernez zwischen Schmierstoff und Umgebung
DeltaT(14) = 50;                 % Stratwert der Potentialdiffernez zwischen Radkörper und Umgebung
DeltaT(15) = 90;                 % Stratwert der Potentialdiffernez zwischen der Zahnmitte am Schraubrad und Umgebung
DeltaT(16) = 110;                % Startwert der Potentialdifferenz zwischen dem eigespritzten Schmierstoff nach der Schraubradflanke und Umgebung
DeltaT(17) = 180;                % Startwert der Potentialdifferenz zwischen der Schraubradflanke und Umgebung
DeltaT(27) = 80;                 % Stratwert der Potentialdiffernez zwischen der Zahnmitte an der Schnecke und Umgebung
DeltaT(28) = 110;                % Startwert der Potentialdifferenz zwischen dem eigespritzten Schmierstoff nach der Schneckenflanke und Umgebung
DeltaT(29) = 180;                % Startwert der Potentialdifferenz zwischen der Schneckenflanke und Umgebung

% --------------------------------------------------------------------------------------------------------------------
% Temperaturdifferenz zwischen der Temperatur des zur Verzahnung eingespritzen Ölvolumenstroms und der
% Umgebungstemperatur:

DeltaT_zu0 = T_zu-T_0;  

% --------------------------------------------------------------------------------------------------------------------
% Bestimmung des zur Verzahnung eingespritzten Ölvolumenstroms, welcher nicht mit den Zahnflanken wechselwirkt
% und somit direkt in den Knoten Schmieröl fließt:

VPunkt_zuOel = VPunkt_zuVZ-(VPunkt_zuRF+VPunkt_zuSF);

% --------------------------------------------------------------------------------------------------------------------
% Für den Enthalpiestrom HPunkt_ab (siehe Abbildung) gilt:

HPunkt_ab = 0;

% ____________________________________________________________________________________________________________________         


%                                   2. Iterative Berechnung der Temperaturverteilung


delta_DeltaT = inf;
tol = 0.1;
niter = 0;

while delta_DeltaT > tol

% ____________________________________________________________________________________________________________________

%                                          Bestimmung der Verlustleistungen

% --------------------------------------------------------------------------------------------------------------------

P_VVSF = 60;                    % Verzahnungsverluste NACHPRFÜFEN
P_VVRF = 60;
P_VLA = 10;
P_VLB1 = 10;
P_VLB2 = 10;
P_VLC = 10;
P_VLD1 = 10;
P_VLD2 = 10;

% ____________________________________________________________________________________________________________________

%                                      Temperaturabhängige thermische Widerstände

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K3 [Gei14, S.78]:

L_K3 = pi*((d_f2+d_IV2)/2);
v_K3 = pi*((d_f2+d_IV2)/2)*(n_2/60);
Re_K3 = (v_K3*L_K3)/ny_Oel((DeltaT(2)+T_0+DeltaT(14)+T_0)/2);
Nu_K3 = 0.664*Re_K3^(1/2)*Pr_Oel((DeltaT(2)+T_0+DeltaT(14)+T_0)/2)^(1/3);
alpha_K3 = (Nu_K3*lambda_Oel)/L_K3;
A_42 = pi*((0.5*d_f2)^2-(0.5*d_IV2)^2);
A_52 = pi*((0.5*d_f2)^2-(0.5*d_III2)^2);

R_K3 = 1/(alpha_K3*(A_42+A_52));
G_K3 = 1/R_K3;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K12 [Gei14, S.77-78]:

v_K12 = pi*d_a2*(n_2/60);
Re_K12 = (v_K12*L_K12)/ny_Oel((DeltaT(2)+T_0+DeltaT(1)+T_0)/2);
Nu_K12 = 0.664*Re_K12^(1/2)*Pr_Oel((DeltaT(2)+T_0+DeltaT(1)+T_0)/2)^(1/3);
alpha_K12 = (Nu_K12*lambda_Oel)/L_K12;

R_K12 = 1/(alpha_K12*A_innen);
G_K12 = 1/R_K12;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K13 (Gei14, S.22-25; Suc12, S.64)

L_10 = h_Gehaeuse;
Gr_10 = (9.81*beta_Luft((DeltaT(1)+T_0+T_0)/2)*DeltaT(1)*L_10^3)/(ny_Luft((DeltaT(1)+T_0+T_0)/2)^2);
Ra_10 = Gr_10*Pr_Luft;
f_1 = (1+(0.492/Pr_Luft)^(9/16))^(-16/9);
Nu_10 = (0.825+0.387*(Ra_10*f_1)^(1/6))^2;
alpha_10 = (Nu_10*lambda_Luft)/L_10;
R_10 = 1/(alpha_10*(2*b_Gehaeuse*h_Gehaeuse+2*l_Gehaeuse*h_Gehaeuse));

L_20 = (b_Gehaeuse*l_Gehaeuse)/(2*(b_Gehaeuse+l_Gehaeuse));
Gr_20 = (9.81*beta_Luft((DeltaT(1)+T_0+T_0)/2)*DeltaT(1)*L_20^3)/(ny_Luft((DeltaT(1)+T_0+T_0)/2)^2);
Ra_20 = Gr_20*Pr_Luft;
f_2 = (1+(0.322/Pr_Luft)^(11/20))^(-20/11);
Nu_20 = 0.766*(Ra_20*f_2)^(1/5);
alpha_20 = (Nu_20*lambda_Luft)/L_20;
R_20 = 1/(alpha_20*b_Gehaeuse*l_Gehaeuse);

alpha_30 = epsilon_Stahl*5.67*10^(-8)*(((DeltaT(1)+T_0+273.15)^4-(T_0+273.15)^4)/DeltaT(1));
R_30 = 1/(alpha_30*(2*b_Gehaeuse*h_Gehaeuse+2*l_Gehaeuse*h_Gehaeuse+b_Gehaeuse*l_Gehaeuse));

R_40 = d_PlatteI/(lambda_PlatteI*b_Gehaeuse*l_Gehaeuse);

R_50 = d_PlatteII/(lambda_PlatteII*b_Gehaeuse*l_Gehaeuse);

l_I0 = sqrt(b_Gehaeuse*l_Gehaeuse);
S_L60 = (2*pi)/(0.93*log((l_II0+l_I0)/(2*l_I0))-0.0502);
R_60 = 1/(lambda_PlatteI*d_PlatteI*S_L60);

alpha_70 = 5;
A_40 = l_II0^2;
R_70 = 1/(alpha_70*(A_40-b_Gehaeuse*l_Gehaeuse));

alpha_80 = 4;
R_80 = 1/(alpha_80*b_Gehaeuse*l_Gehaeuse);

R_K13 = 1/(1/R_10+1/R_20+1/R_30+1/(R_40+R_50+R_80)+1/(R_60+R_70));
G_K13 = 1/R_K13;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K14 [Gei14, S.74-75]:

% tau_2 = 1/(n_2/60);
% omega_2 = 2*pi*(n_2/60);
% h_Zahn2 = (d_a2-d_f2)/2;
% d_w2 = (d_a2+d_f2)/2;
% psi_2 = ((d_w2*a_Oel*(tau_2*omega_2)^2)/(2*ny_Oel(80)*h_Zahn2))^0.25;
% 
% if psi_2 <= 0.68
%     R_K14 = (2*pi*sqrt(a_Oel))/(1.14*b_R2*2*z_2*h_Zahn2*lambda_Oel*omega_2*sqrt(tau_2))
% elseif (0.68 < psi_2) & (psi_2 <= 1.5)
%     R_K14 = (2*pi*sqrt(a_Oel))/((1.55-0.6*psi_2)*2*b_R2*z_2*h_Zahn2*lambda_Oel*omega_2*sqrt(tau_2))
% else
%     psi_2 = 1.5;
%     R_K14 = (2*pi*sqrt(a_Oel))/((1.55-0.6*psi_2)*2*b_R2*z_2*h_Zahn2*lambda_Oel*omega_2*sqrt(tau_2))
% end
% 
% G_K14 = 1/R_K14;

G_K14 = 0.0000000001;

% Da die zuvor durch G_K14 beschriebene thermische Verbindung zwischen der Zahnflanke des Schraubrads und dem 
% Schmierstoff im aktuellen Stand des Modells nicht mehr berücksichtigt wird, wurde die Verbindung im Matlabcode
% durch einen extrem hohen Widerstand bzw. sehr kleinen Leitwert gekappt. Der Leitwert wurde dabei nicht direkt 
% zu Null gesetzt, um eine ungünstige Konditionierung der Matrix zu verhindern. Gleichzeitig wurde auch darauf 
% verzeichtet ihn gänzlich aus dem Modell (und damit aus der weiter unten definierten Matrix) zu entfernen, um 
% diese Verbindung in künftigen Optimierungen des Modells bei Bedarf leicht wieder integrieren zu können.

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K15 [Gei14, S.74-75]:

% gamma_1 = acos((0.5*d_a1-1.5*m_n1)/(0.5*d_a1));
% tau_1 = (1/(n_1/60))*((pi-gamma_1)/(2*pi));
% omega_1 = 2*pi*(n_1/60);
% h_Zahn1 = (d_a1-d_f1)/2;
% d_w1 = (d_a1+d_f1)/2;
% psi_1 = ((d_w1*a_Oel*(tau_1*omega_1)^2)/(2*ny_Oel(80)*h_Zahn1))^0.25;
% 
% if psi_1 <= 0.68
%     R_K15 = (2*pi*sqrt(a_Oel))/(1.14*b_R2*2*z_1*h_Zahn1*lambda_Oel*omega_1*sqrt(tau_1))
% elseif (0.68 < psi_1) & (psi_1 <= 1.5)
%     R_K15 = (2*pi*sqrt(a_Oel))/((1.55-0.6*psi_1)*2*b_R2*z_1*h_Zahn1*lambda_Oel*omega_1*sqrt(tau_1))
% else
%     psi_1 = 1.5;
%     R_K15 = (2*pi*sqrt(a_Oel))/((1.55-0.6*psi_1)*2*b_R2*z_1*h_Zahn1*lambda_Oel*omega_1*sqrt(tau_1))
% end
% 
% G_K15 = 1/R_K15;

G_K15 = 0.0000000001;

% Da die zuvor durch G_K15 beschriebene thermische Verbindung zwischen der Zahnflanke der Schnecke und dem 
% Schmierstoff im aktuellen Stand des Modells nicht mehr berücksichtigt wird, wurde die Verbindung im Matlabcode
% durch einen extrem hohen Widerstand bzw. sehr kleinen Leitwert gekappt. Der Leitwert wurde dabei nicht direkt 
% zu Null gesetzt, um eine ungünstige Konditionierung der Matrix zu verhindern. Gleichzeitig wurde auch darauf 
% verzeichtet ihn gänzlich aus dem Modell (und damit aus der weiter unten definierten Matrix) zu entfernen, um 
% diese Verbindung in künftigen Optimierungen des Modells bei Bedarf leicht wieder integrieren zu können.

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K16 [Gei14 S.77]:

L_K16 = s_n2;
v_K16 = pi*d_a2*(n_2/60);
Re_K16 = (v_K16*L_K16)/ny_Oel((DeltaT(2)+T_0+DeltaT(15)+T_0)/2);
Nu_K16 = 0.664*Re_K16^(1/2)*Pr_Oel((DeltaT(2)+T_0+DeltaT(15)+T_0)/2)^(1/3);
alpha_K16 = (Nu_K16*lambda_Oel)/L_K16;
A_32 = 2*A_Zahnstirn2+A_Zahnkopf2;

R_K16 = 1/(alpha_K16*A_32);
G_K16 = 1/R_K16;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K17 [Gei14, S.77]:

L_K17 = b_R2;
v_K17 = pi*d_a1*(n_1/60);
Re_K17 = (v_K17*L_K17)/ny_Oel((DeltaT(2)+T_0+DeltaT(27)+T_0)/2);
Nu_K17 = 0.664*Re_K17^(1/2)*Pr_Oel((DeltaT(2)+T_0+DeltaT(27)+T_0)/2)^(1/3);
alpha_K17 = (Nu_K17*lambda_Oel)/L_K17;
A_41 = s_a1*b_R2;

R_K17 = 1/(alpha_K17*A_41);
G_K17 = 1/R_K17;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K18 [Suc12, S.55-57; Kre58, S.324]

R_122 = (0.5*B_B2+l_V2)/(lambda_Welle2*pi*(0.5*d_I2)^2);

d_VI2 = A_62/(2*pi*l_VI2);
R_132 = l_VI2/(lambda_Welle2*pi*(0.5*d_VI2)^2);

Re_142 = (2*pi*(n_2/60)*pi*d_VI2^2)/ny_Luft(0.25*DeltaT(4)+T_0);
Gr_142 = (9.81*beta_Luft(0.25*DeltaT(4)+T_0)*0.5*DeltaT(4)*d_VI2^3)/ny_Luft(0.25*DeltaT(4)+T_0)^2;
Nu_142 = 0.11*((0.5*Re_142^2+Gr_142)*Pr_Luft)^0.35;
alpha_142 = (Nu_142*lambda_Luft)/d_VI2;
R_142 = 1/(alpha_142*A_62);

R_K18 = R_122+R_132+R_142;
G_K18 = 1/R_K18;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_K19 [Suc12, S.55-57; Bec63, S.1054]

R_151 = (0.5*B_D2+l_VII1)/(lambda_Welle1*pi*(0.5*d_VII1)^2);

d_VIII1 = A_51/(2*pi*l_VIII1);
R_161 = l_VIII1/(lambda_Welle1*pi*(0.5*d_VIII1)^2);

Re_171 = (2*pi*(n_1/60)*pi*d_VIII1^2)/ny_Luft(0.25*DeltaT(19)+T_0);
Nu_171 = 0.076*Re_171^0.7;
alpha_171 = (Nu_171*lambda_Luft)/d_VIII1;
R_171 = 1/(alpha_171*A_51);

R_K19 = R_151+R_161+R_171;
G_K19 = 1/R_K19;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_HPunkt1 [Gei14, S.54-55]:

VPunkt_1716 = VPunkt_zuRF;
R_HPunkt1 = 1/(VPunkt_1716*roh_Oel((DeltaT(17)+T_0+DeltaT(16)+T_0)/2)*c_Oel((DeltaT(17)+T_0+DeltaT(16)+T_0)/2));
G_HPunkt1 = 1/R_HPunkt1;

% -------------------------------------------------------------------------------------------------------------------- 
% Berechnung von R_HPunkt2 [Gei14, S.54-55]:

VPunkt_162 = VPunkt_zuRF;
R_HPunkt2 = 1/(VPunkt_162*roh_Oel((DeltaT(16)+T_0+DeltaT(2)+T_0)/2)*c_Oel((DeltaT(16)+T_0+DeltaT(2)+T_0)/2));
G_HPunkt2 = 1/R_HPunkt2;

% --------------------------------------------------------------------------------------------------------------------
% Berechnung von R_HPunkt3 [Gei14, S.54-55]:

VPunkt_2928 = VPunkt_zuSF;
R_HPunkt3 = 1/(VPunkt_2928*roh_Oel((DeltaT(29)+T_0+DeltaT(28)+T_0)/2)*c_Oel((DeltaT(29)+T_0+DeltaT(28)+T_0)/2));
G_HPunkt3 = 1/R_HPunkt3;

% -------------------------------------------------------------------------------------------------------------------- 
% Berechnung von R_HPunkt4 [Gei14, S.54-55]:

VPunkt_282 = VPunkt_zuSF;
R_HPunkt4 = 1/(VPunkt_282*roh_Oel((DeltaT(28)+T_0+DeltaT(2)+T_0)/2)*c_Oel((DeltaT(28)+T_0+DeltaT(2)+T_0)/2));
G_HPunkt4 = 1/R_HPunkt4;

% ____________________________________________________________________________________________________________________

%                                     Eingabe und Lösung des Gleichungssystems A*U=b 

%                         1                                                                                                                                          2                                                                                                                              3                           4                            5                         6                          7                       8                        9                       10                        11                       12                          13                          14                       15                         16                                                                           17                                                                                  18                        19                       20                      21                      22                       23                      24                     25                       26                      27                        28                                                                                 29      
A = [-G_L11-G_L31-G_K12-G_K13-G_L27-G_L26-G_L16-G_L17                                                                                                              G_K12                                                                                                                          G_L17                         0                          G_L16                       0                          0                     G_L11                      0                        0                         0                        0                           0                           0                        0                          0                                                                            0                                                                                 G_L27                       0                      G_L26                     0                       0                      G_L31                     0                      0                        0                       0                         0                                                                                  0
                        G_K12                                       -G_K1-G_K2-G_K3-G_K16-G_K14-G_K15-G_K17-G_K11-G_K7-G_K12-G_K8-G_K10-G_K9-G_HPunkt4-G_HPunkt2-G_K6-G_K4-G_K5-VPunkt_zuOel*roh_Oel((DeltaT_zu0+T_0+DeltaT(2)+T_0)/2)*c_Oel((DeltaT_zu0+T_0+DeltaT(2)+T_0)/2)                     G_K5                         0                          G_K4                        0                          0                      G_K2                      0                       G_K1                       0                      G_K6                          0                         G_K3                     G_K16                   G_HPunkt2                                                                       G_K14                                                                               G_K10                       0                      G_K9                      0                     G_K8                     G_K11                     0                    G_K7                       0                     G_K17                  G_HPunkt4                                                                             G_K15
                        G_L17                                                                                                                                       G_K5                                                                                                                    -G_L15-G_K5-G_L17                 G_L15                          0                         0                          0                       0                        0                        0                         0                        0                           0                           0                        0                          0                                                                            0                                                                                   0                         0                        0                       0                       0                        0                       0                      0                        0                       0                         0                                                                                  0   
                          0                                                                                                                                          0                                                                                                                            G_L15                -G_L13-G_L15-G_K18                    0                       G_L13                        0                       0                        0                        0                         0                        0                           0                           0                        0                          0                                                                            0                                                                                   0                         0                        0                       0                       0                        0                       0                      0                        0                       0                         0                                                                                  0
                        G_L16                                                                                                                                       G_K4                                                                                                                            0                           0                   -G_L14-G_K4-G_L16                G_L14                        0                       0                        0                        0                         0                        0                           0                           0                        0                          0                                                                            0                                                                                   0                         0                        0                       0                       0                        0                       0                      0                        0                       0                         0                                                                                  0 
                          0                                                                                                                                          0                                                                                                                              0                         G_L13                       G_L14               -G_L12-G_L13-G_L14                G_L12                     0                        0                        0                         0                        0                           0                           0                        0                          0                                                                            0                                                                                   0                         0                        0                       0                       0                        0                       0                      0                        0                       0                         0                                                                                  0
                          0                                                                                                                                          0                                                                                                                              0                           0                            0                       G_L12                -G_L7-G_L12-G_L6                0                        0                        0                        G_L7                    G_L6                          0                           0                        0                          0                                                                            0                                                                                   0                         0                        0                       0                       0                        0                       0                      0                        0                       0                         0                                                                                  0
                        G_L11                                                                                                                                       G_K2                                                                                                                            0                           0                            0                         0                          0               -G_L10-G_K2-G_L11              G_L10                      0                         0                        0                           0                           0                        0                          0                                                                            0                                                                                   0                         0                        0                       0                       0                        0                       0                      0                        0                       0                         0                                                                                  0
                          0                                                                                                                                          0                                                                                                                              0                           0                            0                         0                          0                     G_L10                 -G_L9-G_L10                  G_L9                       0                        0                           0                           0                        0                          0                                                                            0                                                                                   0                         0                        0                       0                       0                        0                       0                      0                        0                       0                         0                                                                                  0
                          0                                                                                                                                         G_K1                                                                                                                            0                           0                            0                         0                          0                       0                       G_L9               -G_L8-G_K1-G_L9                 G_L8                      0                           0                           0                        0                          0                                                                            0                                                                                   0                         0                        0                       0                       0                        0                       0                      0                        0                       0                         0                                                                                  0
                          0                                                                                                                                          0                                                                                                                              0                           0                            0                         0                         G_L7                     0                        0                       G_L8                -G_L4-G_L8-G_L7                 0                          G_L4                         0                        0                          0                                                                            0                                                                                   0                         0                        0                       0                       0                        0                       0                      0                        0                       0                         0                                                                                  0
                          0                                                                                                                                         G_K6                                                                                                                            0                           0                            0                         0                         G_L6                     0                        0                        0                         0                 -G_L6-G_L5-G_K6                   G_L5                         0                        0                          0                                                                            0                                                                                   0                         0                        0                       0                       0                        0                       0                      0                        0                       0                         0                                                                                  0
                          0                                                                                                                                          0                                                                                                                              0                           0                            0                         0                          0                       0                        0                        0                        G_L4                     G_L5                  -G_L18-G_L4-G_L5                 G_L18                      0                          0                                                                            0                                                                                   0                         0                        0                       0                       0                        0                       0                      0                        0                       0                         0                                                                                  0        
                          0                                                                                                                                         G_K3                                                                                                                            0                           0                            0                         0                          0                       0                        0                        0                         0                        0                          G_L18                -G_L3-G_L18-G_K3                G_L3                        0                                                                            0                                                                                   0                         0                        0                       0                       0                        0                       0                      0                        0                       0                         0                                                                                  0
                          0                                                                                                                                        G_K16                                                                                                                            0                           0                            0                         0                          0                       0                        0                        0                         0                        0                           0                         G_L3                -G_L1-G_L3-G_K16                  0                                                                           G_L1                                                                                 0                         0                        0                       0                       0                        0                       0                      0                        0                       0                         0                                                                                  0
                          0                                                                                                                                      G_HPunkt2                                                                                                                          0                           0                            0                         0                          0                       0                        0                        0                         0                        0                           0                           0                        0                -G_HPunkt1-G_HPunkt2                                                               G_HPunkt1                                                                               0                         0                        0                       0                       0                        0                       0                      0                        0                       0                         0                                                                                  0
                          0                                                                                                                                        G_K14                                                                                                                            0                           0                            0                         0                          0                       0                        0                        0                         0                        0                           0                           0                       G_L1                     G_HPunkt1               -VPunkt_zuRF*roh_Oel((DeltaT_zu0+T_0+DeltaT(17)+T_0)/2)*c_Oel((DeltaT_zu0+T_0+DeltaT(17)+T_0)/2)-G_L1-G_K14-G_HPunkt1                       0                         0                        0                       0                       0                        0                       0                      0                        0                       0                         0                                                                                  0
                        G_L27                                                                                                                                      G_K10                                                                                                                            0                           0                            0                         0                          0                       0                        0                        0                         0                        0                           0                           0                        0                          0                                                                            0                                                                          -G_L27-G_K10-G_L25               G_L25                      0                       0                       0                        0                       0                      0                        0                       0                         0                                                                                  0
                          0                                                                                                                                          0                                                                                                                              0                           0                            0                         0                          0                       0                        0                        0                         0                        0                           0                           0                        0                          0                                                                            0                                                                                 G_L25              -G_L22-G_K19-G_L25                0                     G_L22                     0                        0                       0                      0                        0                       0                         0                                                                                  0
                        G_L26                                                                                                                                      G_K9                                                                                                                             0                           0                            0                         0                          0                       0                        0                        0                         0                        0                           0                           0                        0                          0                                                                            0                                                                                   0                         0                -G_L24-G_L26-G_K9             G_L24                     0                        0                       0                      0                        0                       0                         0                                                                                  0
                          0                                                                                                                                          0                                                                                                                              0                           0                            0                         0                          0                       0                        0                        0                         0                        0                           0                           0                        0                          0                                                                            0                                                                                   0                       G_L22                     G_L24           -G_L21-G_L24-G_L22             G_L21                      0                       0                      0                        0                       0                         0                                                                                  0
                          0                                                                                                                                        G_K8                                                                                                                             0                           0                            0                         0                          0                       0                        0                        0                         0                        0                           0                           0                        0                          0                                                                            0                                                                                   0                         0                        0                     G_L21            -G_L20-G_L21-G_K8                 0                       0                      0                      G_L20                     0                         0                                                                                  0
                        G_L31                                                                                                                                      G_K11                                                                                                                            0                           0                            0                         0                          0                       0                        0                        0                         0                        0                           0                           0                        0                          0                                                                            0                                                                                   0                         0                        0                       0                       0               -G_L30-G_K11-G_L31             G_L30                    0                        0                       0                         0                                                                                  0
                          0                                                                                                                                          0                                                                                                                              0                           0                            0                         0                          0                       0                        0                        0                         0                        0                           0                           0                        0                          0                                                                            0                                                                                   0                         0                        0                       0                       0                      G_L30               -G_L29-G_L30               G_L29                      0                       0                         0                                                                                  0
                          0                                                                                                                                        G_K7                                                                                                                             0                           0                            0                         0                          0                       0                        0                        0                         0                        0                           0                           0                        0                          0                                                                            0                                                                                   0                         0                        0                       0                       0                        0                     G_L29            -G_L28-G_K7-G_L29              G_L28                     0                         0                                                                                  0
                          0                                                                                                                                          0                                                                                                                              0                           0                            0                         0                          0                       0                        0                        0                         0                        0                           0                           0                        0                          0                                                                            0                                                                                   0                         0                        0                       0                     G_L20                      0                       0                    G_L28             -G_L19-G_L28-G_L20             G_L19                       0                                                                                  0
                          0                                                                                                                                        G_K17                                                                                                                            0                           0                            0                         0                          0                       0                        0                        0                         0                        0                           0                           0                        0                          0                                                                            0                                                                                   0                         0                        0                       0                       0                        0                       0                      0                      G_L19             -G_L2-G_L19-G_K17                 0                                                                                 G_L2
                          0                                                                                                                                      G_HPunkt4                                                                                                                          0                           0                            0                         0                          0                       0                        0                        0                         0                        0                           0                           0                        0                          0                                                                            0                                                                                   0                         0                        0                       0                       0                        0                       0                      0                        0                       0               -G_HPunkt3-G_HPunkt4                                                                      G_HPunkt3
                          0                                                                                                                                        G_K15                                                                                                                            0                           0                            0                         0                          0                       0                        0                        0                         0                        0                           0                           0                        0                          0                                                                            0                                                                                   0                         0                        0                       0                       0                        0                       0                      0                        0                      G_L2                   G_HPunkt3                     -VPunkt_zuSF*roh_Oel((DeltaT_zu0+T_0+DeltaT(29)+T_0)/2)*c_Oel((DeltaT_zu0+T_0+DeltaT(29)+T_0)/2)-G_L2-G_K15-G_HPunkt3];     



b = [                                                          0
       HPunkt_ab-VPunkt_zuOel*roh_Oel((DeltaT_zu0+T_0+DeltaT(2)+T_0)/2)*c_Oel((DeltaT_zu0+T_0+DeltaT(2)+T_0)/2)*DeltaT_zu0
                                                            -P_VLB2
                                                               0
                                                            -P_VLB1
                                                               0
                                                               0
                                                            -P_VLA
                                                               0
                                                               0
                                                               0
                                                               0
                                                               0
                                                               0
                                                               0
                                                               0
        -P_VVRF-VPunkt_zuRF*roh_Oel((DeltaT_zu0+T_0+DeltaT(17)+T_0)/2)*c_Oel((DeltaT_zu0+T_0+DeltaT(17)+T_0)/2)*DeltaT_zu0
                                                           -P_VLD2
                                                               0
                                                           -P_VLD1
                                                               0
                                                               0
                                                            -P_VLC
                                                               0
                                                               0
                                                               0
                                                               0
                                                               0
       -P_VVSF-VPunkt_zuSF*roh_Oel((DeltaT_zu0+T_0+DeltaT(29)+T_0)/2)*c_Oel((DeltaT_zu0+T_0+DeltaT(29)+T_0)/2)*DeltaT_zu0];
    
   
DeltaT_neu = A\b;
delta_DeltaT = norm(DeltaT_neu-DeltaT);
DeltaT = DeltaT_neu;
niter = niter+1;

end

% ____________________________________________________________________________________________________________________

%                                            Ausgabe der gesuchten Werte

T = DeltaT+T_0                           % Lösungsvektor mit den Temperaturen aller Knotenpunkte in °C

T_Gehaeuse = T(1);                       % Gehäusetemperatur in °C (Knoten 1)
T_Schmierstoff = T(2);                   % Schmierstofftemperatur in °C (Knoten 2)
T_LagerA = T(8);                         % Temperatur von Lager A in °C (Knoten 8)
T_LagerB1 = T(5);                        % Temperatur von Lager B1 in °C (Knoten 5)
T_LagerB2 = T(3);                        % Temperatur von Lager B2 in °C (Knoten 3)
T_LagerC = T(23);                        % Temperatur von Lager C in °C (Knoten 23)
T_LagerD1 = T(20);                       % Temperatur von Lager D1 in °C (Knoten 20)
T_LagerD2 = T(18);                       % Temperatur von Lager D2 in °C (Knoten 18)
T_RF = T(17);                            % Temperatur der Schraubradflanke in °C (Knoten 17)
T_SF = T(29);                            % Temperatur der Schneckenflanke in °C (Knoten 29)
T_ZahnmitteRad = T(15);                  % Temperatur der Zahnmitte des Schraubrads in °C (Knoten 15)
T_ZahnmitteSchnecke = T(27);             % Temperatur der Zahnmitte der Schnecke in °C (Knoten 27)
T_Radkoerper = T(14);                    % Temperatur des Schraubradkörpers in °C (Knoten 14)

