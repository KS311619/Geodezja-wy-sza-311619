%Kornel Samociuk 311619
clear;

%Gwiazdozbiór Raka, gwiazda Delta Cancri
alfa = 8.744749722; %rektascencja (h)
delta = 18.15431; %deklinacja (stopnie)

%Współrzędne na półkuli N
phiN = 52.220072;
lambdaN = 21.012115673748383;

%Współrzędne na półkuli S
phiS = -33.8548157;
lambdaS = 151.2164539;

%Wspórzędne w okolicach równika
phiR = 0.390002;
lambdaR = 9.454001;

%Data i godzina
p = 96; %sprawdzenie pozycji co 15 min

time = zeros(p,1);
for i = 2:p
    time(i,1) = time(i-1,1) + 1/4; 
end
y = 2021;
m = 12;
d = 29;


%kąty godzinne
tN = katgodz(y, m, d, time, lambdaN, alfa);
tS = katgodz(y, m, d, time, lambdaS, alfa);
tR = katgodz(y, m, d, time, lambdaR, alfa);

%azymuty
lN = zeros(p,1);
mN = zeros(p,1);
lS = zeros(p,1);
mS = zeros(p,1);
lR = zeros(p,1);
mR = zeros(p,1);

for i = 1:p %liczniki i mianowniki azymutów
    lN(i,1) = -cosd(delta) * sind(tN(i,1));
    mN(i,1) = cosd(phiN)*sind(delta) - sind(phiN)*cosd(delta)*cosd(tN(i,1));
    lS(i,1) = -cosd(delta) * sind(tS(i,1));
    mS(i,1) = cosd(phiS)*sind(delta) - sind(phiS)*cosd(delta)*cosd(tS(i,1));
    lR(i,1) = -cosd(delta) * sind(tR(i,1));
    mR(i,1) = cosd(phiR)*sind(delta) - sind(phiR)*cosd(delta)*cosd(tR(i,1));
end

    azymN = atand(lN./mN);
    azymS = atand(lS./mS);
    azymR = atand(lR./mR);

for i = 1:p %azym360
    azymN(i,1) = azym360(azymN(i,1), lN(i,1), mN(i,1));
    azymS(i,1) = azym360(azymS(i,1), lS(i,1), mS(i,1));
    azymR(i,1) = azym360(azymR(i,1), lR(i,1), mR(i,1));
end

%odległość zenitalna i wysokość
zN = zeros(p,1);
zS = zeros(p,1);
zR = zeros(p,1);

hN = zeros(p,1);
hS = zeros(p,1);
hR = zeros(p,1);

for i = 1:p
    zN(i,1) = acosd(sind(phiN)*sind(delta) + cosd(phiN)*cosd(delta)*cosd(tN(i,1)));
    zS(i,1) = acosd(sind(phiS)*sind(delta) + cosd(phiS)*cosd(delta)*cosd(tS(i,1)));
    zR(i,1) = acosd(sind(phiR)*sind(delta) + cosd(phiR)*cosd(delta)*cosd(tR(i,1)));

    hN(i,1) = asind(sind(phiN)*sind(delta) + cosd(phiN)*cosd(delta)*cosd(tN(i,1)));
    hS(i,1) = asind(sind(phiS)*sind(delta) + cosd(phiS)*cosd(delta)*cosd(tS(i,1)));
    hR(i,1) = asind(sind(phiR)*sind(delta) + cosd(phiR)*cosd(delta)*cosd(tR(i,1)));
end

%transformacja współrzędnych 
wxN = zeros(p,1);
wyN = zeros(p,1);
wzN = zeros(p,1);

wxS = zeros(p,1);
wyS = zeros(p,1);
wzS = zeros(p,1);

wxR = zeros(p,1);
wyR = zeros(p,1);
wzR = zeros(p,1);

for i = 1:p
    wxN(i,1) = sind(zN(i,1)) * cosd(azymN(i,1));
    wyN(i,1) = sind(zN(i,1)) * sind(azymN(i,1));
    wzN(i,1) = cosd(zN(i,1));

    wxS(i,1) = sind(zS(i,1)) * cosd(azymS(i,1));
    wyS(i,1) = sind(zS(i,1)) * sind(azymS(i,1));
    wzS(i,1) = cosd(zS(i,1));

    wxR(i,1) = sind(zR(i,1)) * cosd(azymR(i,1));
    wyR(i,1) = sind(zR(i,1)) * sind(azymR(i,1));
    wzR(i,1) = cosd(zR(i,1));
end

%rysowanie tras
figure;
[X, Y, Z] = sphere(50);
Z(Z < 0) = 0;
S = surf(X,Y,Z);
axis equal
S.FaceAlpha = 0.05;
hold on
scatter3(wxN,wyN,wzN, 'blue', 'filled');
grid on;
title('Trasa Delta Cancri na nieboskłonie półkuli północnej');
xlabel('x');
ylabel('y');
zlabel('z');

figure;
[X, Y, Z] = sphere(50);
Z(Z < 0) = 0;
S = surf(X,Y,Z);
axis equal
S.FaceAlpha = 0.05;
hold on
scatter3(wxS,wyS,wzS, 'blue', 'filled');
title('Trasa Delta Cancri na nieboskłonie półkuli południowej');
xlabel('x');
ylabel('y');
zlabel('z');

figure;
[X, Y, Z] = sphere(50);
Z(Z < 0) = 0;
S = surf(X,Y,Z);
axis equal
S.FaceAlpha = 0.05;
hold on
scatter3(wxR,wyR,wzR, 'blue', 'filled');
title('Trasa Delta Cancri na nieboskłonie w okolicach równika');
xlabel('x');
ylabel('y');
zlabel('z');

%wykresy h od time
figure;
plot(time,hN)
grid on;
title('Wykres wysokości od czasu na nieboskłonie półkuli północnej');
xlabel('czas');
ylabel('wysokość');

figure;
plot(time,hS)
grid on;
title('Wykres wysokości od czasu na nieboskłonie półkuli południowej');
xlabel('czas');
ylabel('wysokość');

figure;
plot(time,hR)
grid on;
title('Wykres wysokości od czasu na nieboskłonie w okolicach równika');
xlabel('czas');
ylabel('wysokość');

%wykresy azym od time
figure;
plot(time,azymN)
grid on;
title('Wykres azymutu od czasu na nieboskłonie półkuli północnej');
xlabel('czas');
ylabel('azymut');

figure;
plot(time,azymS)
grid on;
title('Wykres azymutu od czasu na nieboskłonie półkuli południowej');
xlabel('czas');
ylabel('azymut');

figure;
plot(time,azymR)
grid on;
title('Wykres azymutu od czasu na nieboskłonie w okolicach równika');
xlabel('czas');
ylabel('azymut');

%momenty wschodu i zachodu gwiazdy
for i = 2:p
    if wzN(i,1) * wzN(i-1,1) < 0
        if wzN(i-1,1) < 0
            dawnN = i * 0.25;
        else 
            duskN = i * 0.25;
        end
    end
end

for i = 2:p
    if wzS(i,1) * wzS(i-1,1) < 0
        if wzS(i-1,1) < 0
            dawnS = i * 0.25;
        else 
            duskS = i * 0.25;
        end
    end
end

for i = 2:p
    if wzR(i,1) * wzR(i-1,1) < 0
        if wzR(i-1,1) < 0
            dawnR = i * 0.25;
        else 
            duskR = i * 0.25;
        end
    end
end

%funkcje
function [t] = katgodz(y,d,m,h,lambda,alfa)
    jd = juliandate(datetime(y,m,d)); %dni
    g = GMST(jd); %stopnie
    UT1 = h*1.002737909350795; %godziny
    %obliczenie czasu gwiazdowego(w stopniach)
    S = UT1*15 + lambda + g;
    %obliczenie kąta godzinnego(w stopniach)
    t = S - alfa*15;
end 

function g = GMST(JD)
    T = (JD - 2451545)/36525;
    g = 280.46061837 + 360.98564736629 * (JD - 2451545.0) + 0.000387933*T.^2 - T.^3/38710000;
    g = mod(g,360);
end

function az = azym360(a, l, m)
p = 0;
    if l < 0
        a = a + 360;
        p = 1;
    end
    if m < 0
        if p == 1
            a = a - 360;
        end
        a = a + 180;
    end
    az = a;
end