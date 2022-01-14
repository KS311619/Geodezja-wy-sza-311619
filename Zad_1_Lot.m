%Kornel Samociuk 311619
clear;
a = 6378137;
e2 = 0.00669437999013;

%wsp samolotu geo
macierzDane = load('dane.txt');
phi = macierzDane(:,1);
lambda = macierzDane(:,2);
h = macierzDane(:,3);

%wsp lotniska geo
philot = 42.0029;
lambdalot = -87.9168;
hlot = 114;

%wsp lotniska xyz
[Nlot, xlot, ylot, zlot] = geo2xyz(a, e2, philot, lambdalot, hlot);

%tworzenie macierzy dla wsp samolotu xyz
Nairp = zeros(252,1);
xairp = zeros(252,1);
yairp = zeros(252,1);
zairp = zeros(252,1);

%wsp samolotu xyz
for i = 1:252
[Nairp(i,1), xairp(i,1), yairp(i,1), zairp(i,1)] = geo2xyz(a, e2, phi(i,1), lambda(i,1), h(i,1));
end

%macierz delt
delta = [xairp - xlot, yairp - ylot, zairp - zlot]';

%macierz do transpozycji
mdt = [-sind(philot)*cosd(lambdalot), -sind(lambdalot), cosd(philot)*cosd(lambdalot);
       -sind(philot)*sind(lambdalot), cosd(lambdalot), cosd(philot)*sind(lambdalot);
       cosd(philot), 0, sind(philot)];

%macierz neu
neu = mdt' * delta;

%wsp samolotu względem lotniska (neu)
n = neu(1,:)';
e = neu(2,:)';
u = neu(3,:)';

%azymut, odległość skośna i odległość zenitalna
azym = atand(e./n);
s = sqrt(n.^2 + e.^2 + u.^2);
zenit = acosd(u./s);

%wyznaczanie widoczności samolotu na podstawie u
for i = 1:252
    if u(i,1) > 0
        disp('Samolot pojawi się nad horyzontem dla wsp.:')
        fprintf('x: %d', xairp(i,1)); fprintf('\n');
        fprintf('y: %d', yairp(i,1)); fprintf('\n');
        fprintf('z: %d', zairp(i,1)); fprintf('\n');
        break
    end
end

%wyznaczanie widoczności samolotu na podstawie z 
for i = 1:252
    if zenit(i,1) < 90
        disp('Samolot pojawi się nad horyzontem dla wsp.:')
        fprintf('x: %d', xairp(i,1)); fprintf('\n');
        fprintf('y: %d', yairp(i,1)); fprintf('\n');
        fprintf('z: %d', zairp(i,1)); fprintf('\n');
        break
    end
end

%rysowanie trasy lotu (odkomentować przy użyciu)
geoscatter(phi,lambda,5,'.r');

%wizualizacja lotu nue (w razie użycia odkomentować)
figure;
plot3(n,e,u)
grid on
title('Lot w układzie topocentrycznym (n,e,u)');
xlabel('n');
ylabel('e');
zlabel('u');

%wizualizacja lotu xyz (w razie użycia odkomentować)
figure;
plot3(xairp,yairp,zairp); 
grid on
title('Lot w układzie geodezyjnym (x,y,z)');
xlabel('x');
ylabel('y');
zlabel('z');

function[N,x,y,z] = geo2xyz(a,e2,p,l,h)
N = a/sqrt(1 - e2 * sind(p)^2);
x = (N + h) * cosd(p) * cosd(l);
y = (N + h) * cosd(p) * sind(l);
z = (N * (1-e2) + h) * sind(p);
end