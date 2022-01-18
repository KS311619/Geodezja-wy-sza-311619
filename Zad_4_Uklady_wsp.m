%Kornel Samociuk 311619
clear;
a = 6378137; %metry
e2 = 0.00669437999013; %liczba

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
[fiSr, lamSr] = Kivioj(fiA, lamA, fiD, lamD, a, e2);

%wsp GRS80 na Gaussa-Krugera (metry)
[GKxA, GKyA, mgkA] = GRS80_2_GK(fiA, lamA, a, e2, 1992);
[GKxB, GKyB, mgkB] = GRS80_2_GK(fiB, lamB, a, e2, 1992);
[GKxC, GKyC, mgkC] = GRS80_2_GK(fiC, lamC, a, e2, 1992);
[GKxD, GKyD, mgkD] = GRS80_2_GK(fiD, lamD, a, e2, 1992);
[GKxS, GKyS, mgkS] = GRS80_2_GK(fiS, lamS, a, e2, 1992);
[GKxSr, GKySr, mgkSr] = GRS80_2_GK(fiSr, lamSr, a, e2, 1992);

%wsp Gausas-Krugera na PL-1992 (metry)
[PL92xA, PL92yA, m92A] = GK_2_1992(GKxA, GKyA, mgkA);
[PL92xB, PL92yB, m92B] = GK_2_1992(GKxB, GKyB, mgkB);
[PL92xC, PL92yC, m92C] = GK_2_1992(GKxC, GKyC, mgkC);
[PL92xD, PL92yD, m92D] = GK_2_1992(GKxD, GKyD, mgkD);
[PL92xS, PL92yS, m92S] = GK_2_1992(GKxS, GKyS, mgkS);
[PL92xSr, PL92ySr, m92Sr] = GK_2_1992(GKxSr, GKySr, mgkSr);

%wsp GRS80 na PL-2000 (metry)
[PL20xA, PL20yA, m20A] = GRS80_2_2000(fiA, lamA, a, e2, 2000);
[PL20xB, PL20yB, m20B] = GRS80_2_2000(fiB, lamB, a, e2, 2000);
[PL20xC, PL20yC, m20C] = GRS80_2_2000(fiC, lamC, a, e2, 2000);
[PL20xD, PL20yD, m20D] = GRS80_2_2000(fiD, lamD, a, e2, 2000);
[PL20xS, PL20yS, m20S] = GRS80_2_2000(fiS, lamS, a, e2, 2000);
[PL20xSr, PL20ySr, m20Sr] = GRS80_2_2000(fiSr, lamSr, a, e2, 2000);

%wsp Gaussa-Krugera na GRS80 (stopnie)
[BgkA, LgkA] = GK_2_GRS80(GKxA, GKyA, 19, a, e2);
[BgkB, LgkB] = GK_2_GRS80(GKxB, GKyB, 19, a, e2);
[BgkC, LgkC] = GK_2_GRS80(GKxC, GKyC, 19, a, e2);
[BgkD, LgkD] = GK_2_GRS80(GKxD, GKyD, 19, a, e2);
[BgkS, LgkS] = GK_2_GRS80(GKxS, GKyS, 19, a, e2);
[BgkSr, LgkSr] = GK_2_GRS80(GKxSr, GKySr, 19, a, e2);

%wsp PL-1992 na GRS80 (stopnie)
[B92A, L92A] = PL1992_2_GRS80(PL92xA, PL92yA, a, e2);
[B92B, L92B] = PL1992_2_GRS80(PL92xB, PL92yB, a, e2);
[B92C, L92C] = PL1992_2_GRS80(PL92xC, PL92yC, a, e2);
[B92D, L92D] = PL1992_2_GRS80(PL92xD, PL92yD, a, e2);
[B92S, L92S] = PL1992_2_GRS80(PL92xS, PL92yS, a, e2);
[B92Sr, L92Sr] = PL1992_2_GRS80(PL92xSr, PL92ySr, a, e2);

%wsp PL-2000 na GRS80 (stopnie)
[B20A, L20A] = PL2000_2_GRS80(PL20xA, PL20yA, a, e2);
[B20B, L20B] = PL2000_2_GRS80(PL20xB, PL20yB, a, e2);
[B20C, L20C] = PL2000_2_GRS80(PL20xC, PL20yC, a, e2);
[B20D, L20D] = PL2000_2_GRS80(PL20xD, PL20yD, a, e2);
[B20S, L20S] = PL2000_2_GRS80(PL20xS, PL20yS, a, e2);
[B20Sr, L20Sr] = PL2000_2_GRS80(PL20xSr, PL20ySr, a, e2);

%Obliczanie pól dla kolejnych układów (metry^2)
PGRS80 = PoleBL(fiA, lamA, fiB, lamC, a, e2);
Pgk = Pole(GKxA, GKyA, GKxB, GKyB, GKxC, GKyC, GKxD, GKyD);
P92 = Pole(PL92xA, PL92yA, PL92xB, PL92yB, PL92xC, PL92yC, PL92xD, PL92yD);
P20 = Pole(PL20xA, PL20yA, PL20xB, PL20yB, PL20xC, PL20yC, PL20xD, PL20yD);

%Obliczanie zniekształceń dla kolejnych układów 
[kgkA, khagkA, m2gkA] = Zniekszt(mgkA);
[kgkB, khagkB, m2gkB] = Zniekszt(mgkB);
[kgkC, khagkC, m2gkC] = Zniekszt(mgkC);
[kgkD, khagkD, m2gkD] = Zniekszt(mgkD);
[kgkS, khagkS, m2gkS] = Zniekszt(mgkS);
[kgkSr, khagkSr, m2gkSr] = Zniekszt(mgkSr);

[k92A, kha92A, m2_92A] = Zniekszt(m92A);
[k92B, kha92B, m2_92B] = Zniekszt(m92B);
[k92C, kha92C, m2_92C] = Zniekszt(m92C);
[k92D, kha92D, m2_92D] = Zniekszt(m92D);
[k92S, kha92S, m2_92S] = Zniekszt(m92S);
[k92Sr, kha92Sr, m2_92Sr] = Zniekszt(m92Sr);

[k20A, kha20A, m2_20A] = Zniekszt(m20A);
[k20B, kha20B, m2_20B] = Zniekszt(m20B);
[k20C, kha20C, m2_20C] = Zniekszt(m20C);
[k20D, kha20D, m2_20D] = Zniekszt(m20D);
[k20S, kha20S, m2_20S] = Zniekszt(m20S);
[k20Sr, kha20Sr, m2_20Sr] = Zniekszt(m20Sr);

%Zniekształcenie
function[k, kha, m2] = Zniekszt(m_u)
    k = 1 - m_u;
    k = round(-k*1000, 3);
    kha = 1 - m_u^2;
    kha = round(-kha*10000, 6);
    m2 = round(m_u^2, 6);
end

%Obliczanie pola (stopnie -> metry^2)
function [P] = PoleBL(fiA, lamA, fiB, lamB, a, e2)
    b = a * sqrt(1-e2); %metry
    fiA = deg2rad(fiA); 
    lamA = deg2rad(lamA);
    fiB = deg2rad(fiB);
    lamB = deg2rad(lamB);
    phiA = (sin(fiA) / (1-e2 *sin(fiA)^2)) + (1/(2*sqrt(e2))) * log((1+(sqrt(e2)*sin(fiA))) / (1-(sqrt(e2)*sin(fiA))));
    phiB = (sin(fiB) / (1-e2 *sin(fiB)^2)) + (1/(2*sqrt(e2))) * log((1+(sqrt(e2)*sin(fiB))) / (1-(sqrt(e2)*sin(fiB))));
    P = (((b^2)*(lamB-lamA))/2) * (phiA - phiB);
end

%Obliczanie pola (metry -> metry^2)
function [P] = Pole(xA, yA, xB, yB, xC, yC, xD, yD)
    Pol = [xA yA; xB yB; xD yD; xC yC];
    Polin = polyshape(Pol);
    P = area(Polin);
end

%PL-2000 na GRS80 (metry -> stopnie)
function [B, L] = PL2000_2_GRS80(X20, Y20, a, e2)
    L0 = 0; 
    p = mod(Y20,1000000);
    nr = (Y20 - p) / 1000000;
   
    if nr == 5
        L0 = 15;
    elseif nr == 6
        L0 = 18;
    elseif nr == 7
        L0 = 21;
    elseif nr == 8
        L0 = 24;
    end
    
    Xgk = X20 / 0.999923;
    Ygk = (Y20 - 500000 - nr*1000000) / 0.999923;

    [B, L] = GK_2_GRS80(Xgk, Ygk, L0, a, e2);
end

%PL-1992 na GRS80 (metry -> stopnie)
function [B, L] = PL1992_2_GRS80(X92, Y92, a, e2)
    Xgk = (X92 + 5300000) / 0.9993;
    Ygk = (Y92 - 500000) / 0.9993;
    [B, L] = GK_2_GRS80(Xgk, Ygk, 19, a, e2);
end

%Gauss-Kruger na GRS80 (metry -> stopnie)
function [B, L] = GK_2_GRS80(Xgk, Ygk, L0, a, e2)
    A0 = 1 - (e2/4) - (3*(e2^2)/64) - (5*(e2^3)/256);
    A2 = (3/8) * (e2 + ((e2^2)/4) + (15*(e2^3))/128);
    A4 = (15/256) * (e2^2 + (3*(e2^3))/4);
    A6 = 35*(e2^3)/3072;

    B1 = zeros(1,1);
    B1(1,1) = Xgk/(a*A0);
    delta = a * (A0*B1(1,1) - A2*sin(2*B1(1,1)) + A4*sin(4*B1(1,1)) - A6*sin(6*B1(1,1)));

    i = 2;

    while 1
        B1(1,i) = B1(1,i-1) + (Xgk - delta)/(a*A0);
        delta = a * (A0*B1(1,i) - A2*sin(2*B1(1,i)) + A4*sin(4*B1(1,i)) - A6*sin(6*B1(1,i)));
        if abs(B1(1,i) - B1(1,i-1)) < (0.000001/3600)
            Bend = B1(1,i);
            break;
        else 
            i = i+1;
        end
    end
    
    ep2 = 0.00673949677548;
    N = a / sqrt(1 - e2*sin(Bend)^2);
    M = (a*(1-e2))/(sqrt((1-e2*(sin(Bend)^2))^3));
    t = tan(Bend);
    eta2 = ep2 * cos(Bend)^2;
    
    B = Bend - ((Ygk^2)*t)/(2*M*N) * (1 - ((Ygk^2)/(12*N^2)) * (5 + 3*(t^2) + eta2 - (9*eta2*(t^2)) - 4*(eta2^2)) +...
        (Ygk^4)/(360*N^4) * (61 + 90*t^2 + 45*t^4));
    L = deg2rad(L0) + (Ygk/(N*cos(Bend))) * (1 - (Ygk^2)/(6*N^2) * (1 + 2*t^2 + eta2) + (Ygk^4)/(120*N^4) * ...
        (5 + 28*t^2 + 24*t^4 + 6*eta2 + 8*eta2*t^2));

    B = rad2deg(B);
    L = rad2deg(L);
end

%GRS80 na PL-2000 (stopnie -> metry)
function [X20, Y20, m20] = GRS80_2_2000(B, L, a, e2, uklad)
    [Xgk, Ygk, mgk, nr] = GRS80_2_GK(B, L, a, e2, uklad);
    X20 = 0.999923 * Xgk;
    Y20 = 0.999923 * Ygk + nr*1000000 + 500000;
    m20 = 0.999923 * mgk;
end

%Gauss-Kruger na PL-1992 (metry -> metry)
function [X92, Y92, m92] = GK_2_1992(Xgk, Ygk, mgk)
    m92 = mgk * 0.9993;
    X92 = Xgk * 0.9993 - 5300000;
    Y92 = Ygk * 0.9993 + 500000;
end

%GRS80 na Gaussa-Krugera (stopnie -> metry)
function [Xgk, Ygk, mgk, nr] = GRS80_2_GK(B, L, a, e2, uklad)
    L0 = 0;
    nr = 0;

    if uklad == 1992
        L0 = deg2rad(19.0); %radiany
    elseif uklad == 2000
        if (L >= 13.5) && (L < 16.5)
            L0 = deg2rad(15.0);
            nr = 5;
        elseif (L >= 16.5) && (L < 19.5)
            L0 = deg2rad(18.0);
            nr = 6;
        elseif (L >= 19.5) && (L < 22.5)
            L0 = deg2rad(21.0);
            nr = 7;
        elseif (L >= 22.5) && (L < 25.5)
            L0 = deg2rad(24.0);
            nr = 8;
        end
    end

    B = deg2rad(B);
    L = deg2rad(L);
    ep2 = 0.00673949677548;
    N = a / sqrt(1 - e2*sin(B)^2);
    t = tan(B);
    eta2 = ep2 * cos(B)^2;
    l = L - L0; %radiany
    
    %współczynniki
    A0 = 1 - (e2/4) - (3*(e2^2)/64) - (5*(e2^3)/256); 
    A2 = (3/8) * (e2 + ((e2^2)/4) + (15*(e2^3))/128);
    A4 = (15/256) * (e2^2 + (3*(e2^3))/4);
    A6 = 35*(e2^3)/3072;

    delta = a * (A0*B - A2*sin(2*B) + A4*sin(4*B) - A6*sin(6*B));

    Xgk = delta + ((l^2)/2)*N*sin(B)*cos(B) * (1+((l^2)/12)*(cos(B)^2)*(5-(t^2)+(9*eta2)+(4*(eta2^2))) +...
        ((l^4)/360)*(cos(B)^4)*(61-58*(t^2)+(t^4)+(270*eta2)-(330*eta2*(t^2))));
    Ygk = l*N*cos(B) * (1+((l^2)/6)*(cos(B)^2)*(1-(t^2)+eta2)+((l^4)/120)*(cos(B)^4)*(5-(18*(t^2))+(t^4)+ ...
        (14*eta2)-(58*eta2*(t^2))));
    mgk = 1 + ((l^2)/2) * (cos(B)^2) * (1+eta2) + ((l^4)/24) * (cos(B)^4) * (5-4*(t^2));
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