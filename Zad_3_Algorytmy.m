%Kornel Samociuk 311619
clear;
a = 6378137; %metry
e2 = 0.00669437999013; %liczba

%Dane wierzchołków (stopnie)
fiA = 50.25;
lamA = 20.75;
fiB = 50.0;
lamB = 20.75;
fiC = 50.25;
lamC = 21.25;
fiD = 50.0 ;
lamD = 21.25;

%Punkt średniej szerokości (stopnie)
fiS = (fiA + fiD) / 2;
lambdaS = (lamA + lamD) / 2;

%Algorytm Vincentego dla przekątnej AD (wynik w metrach i stopniach)
[sAD, Aad, Ada] = Vincent(fiA, lamA, fiD, lamD, a, e2);

%Algorytm Kivioja dla punktu A (wynik w stopniach)
[FiS, LambdaS, AzymS, AzymSo] = Kivioj(fiA, lamA, sAD, Aad, a, e2);

%Obliczanie różnicy odległości między punktem średniej szerokości i punkem środkowym
[sSredSzer_Srod, Aseso, Asose] = Vincent(fiS, lambdaS, FiS, LambdaS, a, e2);

%Obliczanie pola czworokąta
[P] = Pole(fiA, lamA, fiB, lamC, a, e2);

%Zaokrąglanie wyników
sAD = round(sAD, 3);
sSredSzer_Srod = round(sSredSzer_Srod, 3);
P = round(P, 6);

%Pole czworokąta
function [P] = Pole(fiA, lamA, fiB, lamB, a, e2)
    b = a * sqrt(1-e2); %metry
    fiA = deg2rad(fiA); 
    lamA = deg2rad(lamA);
    fiB = deg2rad(fiB);
    lamB = deg2rad(lamB);
    phiA = (sin(fiA) / (1-e2 *sin(fiA)^2)) + (1/(2*sqrt(e2))) * log((1+(sqrt(e2)*sin(fiA))) / (1-(sqrt(e2)*sin(fiA))));
    phiB = (sin(fiB) / (1-e2 *sin(fiB)^2)) + (1/(2*sqrt(e2))) * log((1+(sqrt(e2)*sin(fiB))) / (1-(sqrt(e2)*sin(fiB))));
    P = (((b^2)*(lamB-lamA))/2) * (phiA - phiB);
end

%Algorytm Vincentego 
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

%Algorytm Kivioja 
function [FiS, LambdaS, AzymS, AzymSo] = Kivioj(fiA, lamA, sAD, Aad, a, e2)
    n = 19;
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

%przekształcanie azymutów na ćwiartki
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