%Kornel Samociuk 311619
clear;

%Dane do GRS80
ag = 6378137; %metry
e2g = 0.00669437999013; %liczba

%Dane do Krasowskiego
ak = 6378245; %metry
e2k = 0.0066934215520398155; %liczba
X0 = -33.4297; %metry
Y0 = 146.5746; %metry
Z0 = 76.2865; %metry
Ex = -0.35867 / 3600; %stopnie
Ey = -0.05283 / 3600; %stopnie
Ez = 0.84354 / 3600; %stopnie
kappa = 0.8407728 * 10^(-6); %liczba

%Dane wierzchołków GRS80 (stopnie)
fiA = 50.25; %B
lamA = 20.75; %L
fiB = 50.0;
lamB = 20.75;
fiC = 50.25;
lamC = 21.25;
fiD = 50.0 ;
lamD = 21.25;

%Punkt średniej szerokości (stopnie)
fiS = (fiA + fiD) / 2;
lamS = (lamA + lamD) / 2;

%Punkt środkowy (stopnie)
[fiSr, lamSr] = Kivioj(fiA, lamA, fiD, lamD, ag, e2g);

%wsp BHL GRS80 na XYZ GRS80
[xA, yA, zA] =  BLH2XYZ(fiA, lamA, 0, ag, e2g);
[xB, yB, zB] =  BLH2XYZ(fiB, lamB, 0, ag, e2g);
[xC, yC, zC] =  BLH2XYZ(fiC, lamC, 0, ag, e2g);
[xD, yD, zD] =  BLH2XYZ(fiD, lamD, 0, ag, e2g);
[xS, yS, zS] =  BLH2XYZ(fiS, lamS, 0, ag, e2g);
[xSr, ySr, zSr] =  BLH2XYZ(fiSr, lamSr, 0, ag, e2g);

%wsp XYZ GRS80 na XYZ Krasowskiego
[KxA, KyA, KzA] = Transform(xA, yA, zA, X0, Y0, Z0, Ex, Ey, Ez, kappa);
[KxB, KyB, KzB] = Transform(xB, yB, zB, X0, Y0, Z0, Ex, Ey, Ez, kappa);
[KxC, KyC, KzC] = Transform(xC, yC, zC, X0, Y0, Z0, Ex, Ey, Ez, kappa);
[KxD, KyD, KzD] = Transform(xD, yD, zD, X0, Y0, Z0, Ex, Ey, Ez, kappa);
[KxS, KyS, KzS] = Transform(xS, yS, zS, X0, Y0, Z0, Ex, Ey, Ez, kappa);
[KxSr, KySr, KzSr] = Transform(xSr, ySr, zSr, X0, Y0, Z0, Ex, Ey, Ez, kappa);

%wsp XYZ Krasowskiego na BHL Krasowskiego
[KbA, KlA, KhA] = Hirvonen(KxA, KyA, KzA, ak, e2k);
[KbB, KlB, KhB] = Hirvonen(KxB, KyB, KzB, ak, e2k);
[KbC, KlC, KhC] = Hirvonen(KxC, KyC, KzC, ak, e2k);
[KbD, KlD, KhD] = Hirvonen(KxD, KyD, KzD, ak, e2k);
[KbS, KlS, KhS] = Hirvonen(KxS, KyS, KzS, ak, e2k);
[KbSr, KlSr, KhSr] = Hirvonen(KxSr, KySr, KzSr, ak, e2k);

%Transformacja na elipsoidę Krasowskiego (metry -> metry)
function [X, Y, Z] = Transform(Xp, Yp, Zp, X0, Y0, Z0, Ex, Ey, Ez, kappa)
    mXYZ0 = [X0; Y0; Z0];
    mXYZp = [Xp; Yp; Zp];
    mPar = [kappa, deg2rad(Ez), deg2rad(-Ey);
            deg2rad(-Ez), kappa, deg2rad(Ex);
            deg2rad(Ey), deg2rad(-Ex), kappa];
    mXYZ = (mXYZp + mPar*mXYZp + mXYZ0)';

    X = mXYZ(1,1);
    Y = mXYZ(1,2);
    Z = mXYZ(1,3);
end

%Zamiana wsp BLH na XYZ (stopnie -> metry)
function [X, Y, Z] = BLH2XYZ(B, L, H, a, e2)
    B = deg2rad(B);
    L = deg2rad(L);
    N = a / (sqrt(1 - e2 * sin(B)^2));
    X = (N + H) * cos(B) * cos(L);
    Y = (N + H) * cos(B) * sin(L);
    Z = (N * (1 - e2) + H) * sin(B);
end

%Algorytm Hirvonena (metry -> stopnie)
function [B, L, H] = Hirvonen(X, Y, Z, a, e2)
    r = sqrt(X^2 + Y^2);
    Bo = atan((Z/r) * (1-e2)^(-1));

    while 1
         N = a / (sqrt(1 - e2 * sin(Bo)^2));
         H = (r/cos(Bo)) - N;
         Bn = atan((Z/r) * (1-e2*(N/(N + H)))^(-1));

         if abs(Bn - Bo) < deg2rad(0.00005 / 3600)
            break
        else 
            Bo = Bn;
         end
    end

    N = a / (sqrt(1 - e2 * sin(Bn)^2));
    B = rad2deg(Bn);
    L = rad2deg(atan(Y/X));
    H = (r/cos(Bn)) - N;

end

%Algorytm Vincentego (stopnie -> stopnie)
function [sAD, Aad, Ada] = Vincent(fiA, lamA, fiD, lamD, a, e2)
    b = a * sqrt(1-e2); %metry
    f = 1 - (b/a); %liczba
    dellab = lamD - lamA; %stopnie

    Ua = atand((1-f) * tand(fiA)); %stopnie
    Ud = atand((1-f) * tand(fiD)); %stopnie

    L = zeros(1,1);
    L(1) = dellab; %stopnie

    i = 1;

    while 1
        sinsig = sqrt(((cosd(Ud)*sind(L(1,i))) * (cosd(Ud)*sind(L(1,i)))) + ...
            ((cosd(Ua)*sind(Ud) - sind(Ua)*cosd(Ud)*cosd(L(1,i))) * (cosd(Ua)*sind(Ud) - sind(Ua)*cosd(Ud)*cosd(L(1,i))))); %liczba

        cossig = sind(Ua)*sind(Ud) + cosd(Ua)*cosd(Ud)*cosd(L(1,i)); %liczba

        sigma = atand(sinsig/cossig); %stopnie

        sinalfa = cosd(Ua)*cosd(Ud)*sind(L(1,i)) / sind(sigma); %liczba

        cos2alfa = 1 - (sinalfa * sinalfa); %liczba

        cos2sigm = cossig - (2*sind(Ua)*sind(Ud)) / cos2alfa; %liczba

        C = (f/16)*cos2alfa * (4 + f*(4 - 3*cos2alfa)); %liczba

        L(1,i+1) = dellab + (1-C)*f*sinalfa * (sigma + C*sinsig * (cos2sigm + C*cossig * (-1 + 2*cos2sigm))); %stopnie
    
        if abs(L(1,i+1) - L(1,i)) < (0.000001/3600)
            Lend = L(1,i+1); %stopnie
            break
        else 
            i = i+1;
        end
    end

    u2 = (((a*a) - (b*b)) / (b*b)) * cos2alfa; %metry^2
    A = 1 + (u2/16384) * (4096 + u2*(-768 + u2*(320 - 175*u2))); %metry^2
    B = (u2/1024) * (256 + u2*(-128 + u2*(74 - 47*u2))); %metry^2

    delsig = B * sinsig * (cos2sigm + (B/4)*(cossig*(-1 + 2*(cos2sigm*cos2sigm)) ...
        - (B/6)*cos2sigm*(-3 + 4*(sinsig*sinsig)) * (-3 + 4*(cos2sigm*cos2sigm)))); %radiany

    sAD = b * A * (deg2rad(sigma) - delsig); %metry
    
    Aad = atand( (cosd(Ud)*sind(Lend)) / (cosd(Ua)*sind(Ud) - sind(Ua)*cosd(Ud)*cosd(Lend)));
    Ada = atand( (cosd(Ua)*sind(Lend)) / (cosd(Ua)*sind(Ud)*cosd(Lend) - sind(Ua)*cosd(Ud)));

    Aad = azym360(Aad, (cosd(Ud)*sind(Lend)), (cosd(Ua)*sind(Ud) - sind(Ua)*cosd(Ud)*cosd(Lend)));
    Ada = azym360(Ada, (cosd(Ua)*sind(Lend)), (cosd(Ua)*sind(Ud)*cosd(Lend) - sind(Ua)*cosd(Ud))) + 180;
end

%Algorytm Kivioja (stopnie -> stopnie)
function [FiS, LambdaS, AzymS, AzymSo] = Kivioj(fiA, lamA, fiD, lamD, a, e2)
    
    [sAD, Aad] = Vincent(fiA, lamA, fiD, lamD, a, e2);
    
    n = 22;
    ds = sAD/(2*n);

    Fi = deg2rad(fiA);
    Lambda = deg2rad(lamA);
    Azym = deg2rad(Aad);

    for i = 1:n
        M = (a*(1-e2))/(sqrt((1-e2*(sin(Fi)^2))^3)); %metry
        N = a / sqrt(1 - e2*(sin(Fi)^2)); %metry
        DelFi = (ds * cos(Azym)) / M; %radiany 
        DelAzym = ((sin(Azym)*tan(Fi)) / N)*ds; %radiany
        Fim = Fi + (DelFi/2); %radiany
        Azymm = Azym + (DelAzym/2); %radiany
        Mm = (a*(1-e2)) / sqrt((1 - e2*(sin(Fim))^2)^3); %metry
        Nm = a / sqrt(1 - e2*(sin(Fim))^2); %metry
        DelFi = (ds*cos(Azymm)) / Mm; %radiany
        DelLambda = (ds*sin(Azymm)) / (Nm*cos(Fim)); %radiany
        DelAzym = (sin(Azymm)*tan(Fim)*ds) / Nm; %radiany
        Fi = Fi + DelFi; %radiany
        Lambda = Lambda + DelLambda; %radiany
        Azym = Azym + DelAzym; %radiany
    end

    FiS = rad2deg(Fi);
    LambdaS = rad2deg(Lambda);
    AzymS = rad2deg(Azym);
    AzymSo = AzymS + 180;

end

%przekształcanie azymutów na ćwiartki (stopnie -> stopnie)
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