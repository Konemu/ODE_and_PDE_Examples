clear;
% Numerische Lösung von y* = k*y mit Eulerverfahren
% Setup
k = 1;          % Wachstums-/Zerfallskonstante
x_0 = 0;        % Stelle des Anfangswerts = linke Intervallgrenze
y_0 = 1;        % Anfangswert
x_max = 2;      % Rechte Intervallgrenze

f = @(y) k.*y ;         % Rechte Seite der DGL y' = k*y
sol = @(x) exp(k.*x);   % Exakte Lösung
n = [10 50 100];   % Zahl der Raumpunkte nach Diskretisierung

% Iteration und Plot für gegebene n
for j = 1:length(n)
    [a,b] = fwdEuler(f, n(j), x_max/n(j), y_0, x_0, x_max);
    plot(a, b)
    hold on
end

% Plot der exakten Lösung
plot(a, sol(a))
hold off;
title({'Numerische L\"osung von $d_x y = k y$ mit dem Euler-Verfahren','f\"ur verschiedene Schrittweiten und $y(0) = 1$'},"Interpreter","latex","FontSize",14)
xlabel('$x$',"Interpreter","latex","FontSize",16)
ylabel('$y(x)$',"Interpreter","latex","FontSize",16)
legend('$n = 10$', '$n = 50$', '$n = 100$', 'Exakt',"Interpreter","latex","Location","Northwest")
%%
clear;
% Numerische Lösung gedämpfter Oszillator mit Euler sowie Runge-Kutta-Verfahren
% x'' = - 2*d*x' - w^2*x 

% Setup
t_0 = 0;                        % Zeitintervallgrenzen
t_max = 25;
x_0 = 1;                        % Anfangswerte
v_0 = 0;
n = 50;                         % Schrittzahl
h = (t_max - t_0) / n;          % Schrittweite
d = 0.1;                        % Dämpfung
w2 = 1;                         % Eigenfrequenz

t = linspace(t_0, t_max, n);    % Zeitdiskretisierung
x = zeros(1,n);                 % Init x-Array
v = zeros(1,n);                 % Init v-Array
x(1) = x_0;                     % Init Anfangswerte
v(1) = v_0;

% Exakte Lösung (Schwingfall):
w_d = sqrt(w2 - d^2);
x_exakt = @(t) exp(-d.*t) .* ( x_0 .* cos(w_d.*t) + (x_0.*d + v_0)./w_d .*sin(w_d.*t) );

% Rechte Seiten des DGL-Systems 1. Ordnung
% v' = -2d**v - w^2*x
% x' = v
f_v = @(x, v) -w2 * x - 2*d*v;
f_x = @(x, v) v;

% Euler zum Vergleich:
x_Euler = x;
v_Euler = v;
for i = 2:n
    v_Euler(i) = v_Euler(i-1) + h*f_v(x_Euler(i-1), v_Euler(i-1));
    x_Euler(i) = x_Euler(i-1) + h*v_Euler(i);
end

% Plot Euler
%plot(t, x_Euler)
%hold on;
t_exakt = linspace(t_0, t_max, 1000);
%plot(t_exakt, x_exakt(t_exakt));
%plot(t, v_Euler)
%hold off;
%legend('Euler', 'Exakt')
%title({'Numerische L\"osung von $\ddot{x} = -2d\dot{x} - \omega_0^2 x$ mit dem','Euler-Verfahren und $x(0) = 1$, $\dot{x}(0) = 0$, $dt = 0.5$'},"Interpreter","latex","FontSize",14)
%xlabel('$t$',"Interpreter","latex","FontSize",16)
%ylabel('$x(t)$',"Interpreter","latex","FontSize",16)
%legend('Euler', 'Exakt',"Interpreter","latex","Location","Northeast")


% Runge-Kutta:
for i = 2:n
    k1v = f_v(x(i-1), v(i-1));                          % Koeffizienten
    k2v = f_v(x(i-1), v(i-1) + h/2 * k1v);
    k3v = f_v(x(i-1), v(i-1) + h/2 * k2v);
    k4v = f_v(x(i-1), v(i-1) + h   * k3v);
    
    v(i) = v(i-1) + h/6 * (k1v + 2*k2v + 2*k3v + k4v);  % Schritt
    
    % Da die Rechte Seite der Gleichung für x x-unabhängig ist,
    % entspricht der Runge-Kutta-Schritt für x dem Euler-Verfahren:
    x(i) = x(i-1) + h*v(i);                             
end

% Plot Runge-Kutta
plot(t, x)
hold on;
plot(t, x_Euler)
plot(t_exakt, x_exakt(t_exakt));
%plot(t, v)
hold off;
title({'Numerische L\"osung von $\ddot{x} = -2d\dot{x} - \omega_0^2 x$ mit dem','Runge-Kutta- und Euler-Verfahren und $x(0) = 1$, $\dot{x}(0) = 0$, $dt =0.5$'},"Interpreter","latex","FontSize",12)
xlabel('$t$',"Interpreter","latex","FontSize",16)
ylabel('$x(t)$',"Interpreter","latex","FontSize",16)
legend('Runge-Kutta', 'Euler', 'Exakt',"Interpreter","latex","Location","Northeast")
%%

clear;
% Numerische Betrachtung der Strahlungsbilanz einer dünnen
% Atmosphärenschicht
% init: Konstanten
h   = 8300;              % Dicke der Troposphäre (m)
rho = 1.2;              % Dichte der Troposphäre (kg m^-3)
C   = 1000;             % Massenspez. Wärmekapazität (J kg^-1 K^-1)
a   = 0.3;              % Albedo
I   = 1367;             % mittlere Sonnenintensität bei der Erde (W m^-2)
e   = 0.6;              % Emissivität der Atmosphäre
s   = 5.67e-8;          % Stefan-Boltzmann-Konstante (W m^-2 K^-4)

% init: Diskretisierung
t_0 = 0;
t_max = 365;                    % 365 Tage (d)
dt = 35 ;                        % dt = 1 Tag (d)
n = round( (t_max - t_0)/dt );
T_0 = 300;                      % Anfangsbedingung T(t_0) = T_0 = 300 K

% Gleichgewichtstemperatur:
T_equilibrium = ((1-a)*I/(4*e*s))^(1/4)

% Rechte Seite, ||Einheit K/d||
f = @(t, T) 1/(4*h*rho*C) * ( (1-a)*I - 4 * e * s * T^4 ) * 86400;


t = linspace(t_0, t_max, n);
T = zeros(1, n);
T(1) = T_0;

for i = 2:n
    T_curr = T(i-1);
    k1 = f(i*dt, T_curr);
    k2 = f(i*dt + dt/2, T_curr + dt/2*k1);
    k3 = f(i*dt + dt/2, T_curr + dt/2*k2);
    k4 = f(i*dt + dt, T_curr + dt*k3);
    
    T(i) = T_curr + dt/6 * ( k1 + 2*k2 + 2*k3 + k4 );    
end

% Analytische Näherung
tau = (h*rho*C)/(4*e*s*T_equilibrium^3) / 86400;
T_analytical = T_equilibrium + (T_0 - T_equilibrium) * exp(-linspace(t_0, t_max, 365)/tau);

clf;
plot(t, T);
hold on;
plot(linspace(t_0, t_max, 365), T_analytical);
hold off;
legend('RK, $dt = 35$ d', 'analytische N\"aherung', 'interpreter', 'latex')
xlabel('Zeit $t$ (d)', "Interpreter","latex", "FontSize",14)
ylabel('Temperatur $T$ (K)', "Interpreter","latex", "FontSize",14)
title('Numerische L\"osung der Strahlungsbilanzgleichung', "Interpreter","latex")

%%
function [xVec, solVec] = fwdEuler(f, n, h, y_0, x_0, x_max) 
% Forwärts-Eulerverfahren
% f: Linke Seite der DGL, n: Zahl der Raumpunkte nach Diskretisierung
% h: Schrittweite, y_0: Anfangswert
% x_0, x_max: Intervallgrenzen
    xVec = x_0:h:x_max;
    
    solVec = zeros(1,n+1);
    solVec(1) = y_0;
    for i = 2:n+1
        solVec(i) = solVec(i-1) + h * f(solVec(i-1));
    end
end

