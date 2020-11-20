clear;
clear all;
f = 3.6 * 10^9;               %cz�stotliwo��
c = 3 * 10^8;                 %pr�dko�� �wiat�a
d_max = 2.8;                  %d�ugo�� pokoju
s_max = 2.6;                  %szeroko�� pokoju
h_max = 2.6;                  %wysoko�� pokoju
step1 = 0.005;                  %krok oblicze� dla zadania 1
nd = fix(d_max / step1);           %rozmiar tablicy
ns = fix(s_max / step1);
step2 = 0.01;                 %krok oblicze� dla zadania 2
d = [1:step1:d_max];          %zakres dla zadania 1
sx = [0:step2:(2 * s_max)];   %zakres dla zadania 2
s3 = [0:step1:s_max];         %zakres dla zadania 3
d3 = [0:step1:d_max];         %zakres dla zadania 3
ns3 = [0:step1:10 * s_max];         %zakres dla zadania 3
nd3 = [0:step1:10 * d_max];         %zakres dla zadania 3
a = -0.5;                       %wsp�czynnik odbicia
h_anten = 1;                  %wysoko�� anteny
TXs0 = s_max / 2;             %wsp�rz�dne nadajnika 
TXd0 = 0;                     %wsp�rz�dne nadajnika
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  LOS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fTXlosRX = 2 * pi * f * (d ./ c);

prpo = (1 ./ d) .* exp(-1i .* fTXlosRX);
prpolos = (abs(prpo)) .^ 2;
logprpolos = 10 * log10(prpolos);

% Tau = 1.5 / c;
% AvgTau = (sum(Tau .* prpolos)) / (sum(prpolos));
% SquareAvgTau = (sum((Tau .^2) .* (prpolos))) / (sum(prpolos));
% TauRMS = sqrt((SquareAvgTau) - (AvgTau .^ 2));
% Bc50 = 1 / (5 .* TauRMS);

figure
plot(d, prpolos)
plot(d, logprpolos)
grid
grid minor
title('Wzgl�dny spadek mocy dla �cie�ki LOS')
xlabel('Odleg�o�� [m]');
ylabel('Moc [dB]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  Pojedyncze odbicie + LOS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dsc = hypot(d, s_max);
dsu = hypot(d, (2 * (h_max - 1)));
dpo = hypot (d, (2 * h_anten));
dscp = (2 * d_max) - d;

fTXscRX = 2 * pi * f * (dsc ./ c); %�ciana
fTXsuRX = 2 * pi * f * (dsu ./ c); %sufit
fTXpoRX = 2 * pi * f * (dpo ./ c); %pod�oga
fTXscpRX = 2 * pi * f * (dscp ./ c); %�ciana w linii LOS

prposc = ((a ./ dsc) .* exp(-1i .* fTXscRX));
prposu = ((a ./ dsu) .* exp(-1i .* fTXsuRX));
prpopo = ((a ./ dpo) .* exp(-1i .* fTXpoRX));
prposcp = ((a ./ dscp) .* exp(-1i .* fTXscpRX));

sumaprpo1 = (abs(prpo + ...
    prposc + prposc + prposu + prpopo + prposcp)).^2;

logsumaprpo1 = (10 * log10(sumaprpo1));

figure
plot(d, logsumaprpo1)
grid
grid minor
title('Wzgl�dny spadek mocy dla sumy �cie�ki LOS i �cie�ek powsta�ych na skutek pojedynczego odbicia')
xlabel('Odleg�o�� [m]');
ylabel('Moc [dB]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  Pojedyncze, pod�wjne odbicie + LOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dscsc = hypot(d, (2 * s_max)); %�ciana �ciana
dsupo = hypot(d, (2 * h_max)); %sufit pod�oga / pod�oga sufit
dscscp = (3 * d_max) - d; % �ciana �ciana w lini LOS

fTXscscRX = 2 * pi * f * (dscsc ./ c); %�ciana
fTXsupoRX = 2 * pi * f * (dsupo ./ c); %sufit
fTXscscpRX = 2 * pi * f * (dscscp ./ c); %�ciana za plecami

prposcsc = (((a * a) ./ dscsc) .* exp(-1i .* fTXscscRX));
prposupo = (((a * a) ./ dsupo) .* exp(-1i .* fTXsupoRX));
prposcscp = (((a * a) ./ dscscp) .* exp(-1i .* fTXscscpRX));

sumaprpo2 = (abs(prpo + ...
    prposc + prposc + prposu + prpopo + prposcp + ...
    prposcsc + prposcsc + prposupo + prposupo + prposcscp)).^2;

logsumaprpo2 = (10 * log10(sumaprpo2));

figure
plot(d, logsumaprpo2)
grid
grid minor
title('Wzgl�dny spadek mocy dla sumy �cie�ki LOS i �cie�ek powsta�ych na skutek pojedynczego oraz podw�jnego odbicia')
xlabel('Odleg�o�� [m]');
ylabel('Moc [dB]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Pojedyncze, pod�wjne i potr�jne odbicie + LOS  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dscscsc = hypot(d, (3 * s_max)); %�ciana �ciana �ciana
dsuposu = hypot(d, (3 * h_max)); %sufit pod�oga sufit/ pod�oga sufit pod�oga
dscscscp = (4 * d_max) - d; % �ciana �ciana �ciana w lini LOS

fTXscscscRX = 2 * pi * f * (dscscsc ./ c); %�ciana
fTXsuposcRX = 2 * pi * f * (dsuposu ./ c); %sufit
fTXscscscpRX = 2 * pi * f * (dscscscp ./ c); %�ciana za plecami

prposcscsc = (((a * a * a) ./ dscscsc) .* exp(-1i .* fTXscscscRX));
prposuposu = (((a * a * a) ./ dsuposu) .* exp(-1i .* fTXsuposcRX));
prposcscscp = (((a * a * a) ./ dscscscp) .* exp(-1i .* fTXscscscpRX));

sumaprpo3 = (abs(prpo + ...
    prposc + prposc + prposu + prpopo + prposcp + ...
    prposcsc + prposcsc + prposupo + prposupo + prposcscp + ...
    prposcscsc + prposcscsc + prposuposu + prposuposu + prposcscscp)).^2;

logsumaprpo3 = (10 * log10(sumaprpo3));

figure
plot(d, logsumaprpo3)
grid
grid minor
title('Wzgl�dny spadek mocy dla sumy �cie�ki LOS i �cie�ek powsta�ych na skutek pojedynczego, podw�jnego oraz potr�jnego odbicia')
xlabel('Odleg�o�� [m]');
ylabel('Moc [dB]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%  Wsp�lny wykres  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
figure
hold on
plot(d, logprpolos)
plot(d, logsumaprpo1)
plot(d, logsumaprpo2)
plot(d, logsumaprpo3)
grid
grid minor
hold off
title('Por�wnanie wzgl�dnych spadek�w mocy dla kolejnych sum �cie�ek')
xlabel('Odleg�o�� [m]');
ylabel('Moc [dB]');
legend('LOS','LOS + 1','LOS + 1 + 2','LOS + 1 + 2 + 3')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   funkcja og�lna dla x odbi�   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 0:3
    if n == 0
        flos = -2 * pi * f * (d ./ c);
        prpolos = 1 ./ d;
        sumaprpo = prpolos;
    elseif n == 1
        dscn = hypot(d, (n * s_max));
        dsun = hypot(d, (2 * (h_max - 1)));
        dpon = hypot (d, (2 * h_anten));
        dscpn = abs(((n+1) * d_max) - d);

        fTXscRXn = 2 * pi * f * (dscn ./ c); %�ciana
        fTXsuRXn = 2 * pi * f * (dsun ./ c); %sufit
        fTXpoRXn = 2 * pi * f * (dpon ./ c); %pod�oga
        fTXscpRXn = 2 * pi * f * (dscpn ./ c); %�ciana w linii LOS

        prposcn = (((a^n) ./ dscn) .* exp(-1i .* fTXscRXn));
        prposun = (((a^n) ./ dsun) .* exp(-1i .* fTXsuRXn));
        prposun = (((a^n) ./ dpon) .* exp(-1i .* fTXpoRXn));
        prposcpn = (((a^n) ./ dscpn) .* exp(-1i .* fTXscpRXn));

        sumaprpo = sumaprpo + ((2 * prposcn) + prposun + prposun + prposcpn);
    else
        dscn = hypot(d, (n * s_max));
        dsun = hypot(d, (n * h_max));
        dscpn = abs(((n+1) * d_max) - d);

        fTXscRXn = 2 * pi * f * (dscn ./ c); %�ciana
        fTXsuRXn = 2 * pi * f * (dsun ./ c); %sufit
        fTXscpRXn = 2 * pi * f * (dscpn ./ c); %�ciana w linii LOS

        prposcn = (((a^n) ./ dscn) .* exp(-1i .* fTXscRXn));
        prposun = (((a^n) ./ dsun) .* exp(-1i .* fTXsuRXn));
        prposcpn = (((a^n) ./ dscpn) .* exp(-1i .* fTXscpRXn));

        sumaprpo = sumaprpo + ((2 * prposcn) + (2 * prposun) + prposcpn);
    end
end

abssumaprpo = (abs(sumaprpo)).^2;

logabssumaprpo = (10 * log10(abssumaprpo));

figure
plot(d, (10 * log10(abssumaprpo)))
grid
grid minor
title(['Wzgl�dny spadek mocy dla sumy �cie�ki LOS i �cie�ek powsta�ych na skutek do ', num2str(n) ,' odbi�'])
xlabel('Odleg�o�� [m]');
ylabel('Moc [dB]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   cz�� druga   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prpo2 = [];
logprpo2 = [];
absprpo2 = [];
gamma1 = hypot((0.8 * d_max), (0.5 * s_max));
for s = 0:step2:(2 * s_max)
    z = fix(((1000 * s) * (1 / (1000 * step2))) + 1);
    if z == 0 
        z = (s * (1 / step2)) + 1;
    end
    if s < (1.0625 * s_max)
        if s == (0.5 * s_max)
            fTXlosRXn = 2 * pi * f * ((0.9 * d_max) / c);
            prpolos = (1 / (0.9 * d_max)) * exp(-1i * fTXlosRXn);
            x = hypot((0.9 * d_max), (2 * (h_max - h_anten)));
            fTXpsRXn = 2 * pi * f * (x / c);
            prpos = ((a / x) * exp(-1i * fTXpsRXn));
            x = hypot((0.9 * d_max), (2 * h_anten));
            fTXpsRXn = 2 * pi * f * (x / c);
            prpop = ((a / x) * exp(-1i * fTXpsRXn));
            prpo2(z) = prpolos +  prpos + prpop;
            absprpo2(z) = (abs(prpo2(z)))^2;
            logprpo2(z) = (10 * log10(absprpo2(z)));
        elseif s < (0.5 * s_max)
            los = hypot((0.9 * d_max), ((0.5 * s_max) - s));
            fTXlosRXn = 2 * pi * f * (los / c);
            prpolos = (1 / (los)) * exp(-1i * fTXlosRXn);
            x = hypot((0.9 * d_max), ((0.5 * s_max) - s));
            y = hypot(x, (2 * (h_max - h_anten)));
            fTXpsRXn = 2 * pi * f * (y / c);
            prpos = ((a / y) * exp(-1i * fTXpsRXn));
            y = hypot(x, (2 * h_anten));
            fTXpsRXn = 2 * pi * f * (y / c);
            prpop = ((a / y) * exp(-1i * fTXpsRXn));
            prpo2(z) = prpolos + prpos + prpop;
            absprpo2(z) = (abs(prpo2(z)))^2;
            logprpo2(z) = (10 * log10(absprpo2(z)));
        elseif ((s > (0.5 * s_max)) && (s < (1.0625 * s_max)))
            los = hypot((0.9 * d_max), (s - (0.5 * s_max)));
            fTXlosRXn = 2 * pi * f * (los / c);
            prpolos = (1 / (los)) * exp(-1i * fTXlosRXn);
            x = hypot((0.9 * d_max), (s - (0.5 * s_max)));
            y = hypot(x, (2 * (h_max - h_anten)));
            fTXpsRXn = 2 * pi * f * (y / c);
            prpos = ((a / y) * exp(-1i * fTXpsRXn));
            y = hypot(x, (2 * h_anten));
            fTXpsRXn = 2 * pi * f * (y / c);
            prpop = ((a / y) * exp(-1i * fTXpsRXn));
            prpo2(z) = prpolos + prpos + prpop;
            absprpo2(z) = (abs(prpo2(z)))^2;
            logprpo2(z) = (10 * log10(absprpo2(z)));
        end
    elseif s >= (1.0625 * s_max)
        %gamma1 = hypot((0.8 * d_max), (0.5 * s_max));
        k(z) = (0.1 * d_max);
        j(z) = (s - s_max);
        gamma2(z) = hypot((0.1 * d_max), (s - s_max));
        if gamma2(z) == 0
            k = 0;
        end
        s1s2 = hypot((0.9 * d_max), (s - (0.5 * s_max)));
        OBtroj(z) = gamma1 + gamma2(z) + s1s2;
        p(z) = OBtroj(z) / 2;
        Ptroj(z) = sqrt(p(z) * (p(z) - gamma1) * (p(z) - gamma2(z)) * (p(z) - s1s2));
        h(z) = (2 * Ptroj(z)) / s1s2;
        lambda = c / f;
        v(z) = h(z) * sqrt((2 / lambda) * (s1s2 / (gamma1 * gamma2(z))));
        C(z) = 6.9 + (20 * log10(sqrt((v(z) - 0.1)^2 + 1) + v(z) - 0.1));
        logprpo2(z) = -C(z);
    end
    %absprpo2(z) = (logprpo2(z));
end

figure
plot(sx, logprpo2)
grid
grid minor
title('Wzgl�dny spadek mocy dla sumy �cie�ki LOS i �cie�ki powsta�ej na skutek jednego odbicia oraz dla cienia radarowego tylko dla dyfrakcji')
xlabel('Szeroko�� [m]');
ylabel('Moc [dB]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%   cz�� trzecia   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:ns
%     for j = 1:nd
%         prpo3(i:ns,j:nd)=zeros; 
%         logprpo3(i:ns,j:nd)=zeros; 
%         absprpo3(i:ns,j:nd)=zeros; 
%     end
% end
% prpo3;
% logprpo3;
% absprpo3;
for n = 0 : 3
    for s = 0:step1:(s_max)
        for d = 0:step1:(d_max)
            z = fix(((1000 * s) * (1 / (1000 * step1))) + 1);
            w = fix(((1000 * d) * (1 / (1000 * step1))) + 1);
            Z(z) = z;
            W(w) = w;
            if ((s == (0.5 * s_max)) && (d > 0))
                x = d;
                fTXlosRXn = 2 * pi * f * (x / c);
                prpolos = (1 / x) * exp(-1i * fTXlosRXn);
                
                x = hypot(d, s_max);
                fTXscRXn = 2 * pi * f * (x / c);
                prposc = ((a / x) * exp(-1i * fTXscRXn));
                
                y = hypot(d, (2 * h_anten));
                fTXpRXn = 2 * pi * f * (y / c);
                prpop = ((a / y) * exp(-1i * fTXpRXn));
                
                y = hypot(d, (2 * (h_max - h_anten)));
                fTXsRXn = 2 * pi * f * (y / c);
                prpos = ((a / y) * exp(-1i * fTXsRXn));
                
                x = hypot(d, (2 * s_max));
                fTXscscRXn = 2 * pi * f * (x / c);
                prposcsc = (((a * a) / x) * exp(-1i * fTXscscRXn));
                
                y = hypot(d, (2 * h_max));
                fTXpsRXn = 2 * pi * f * (y / c);
                prpops = ((a / y) * exp(-1i * fTXpsRXn));
                
                prpo3(z,w) = prpolos + (2 * prposc) + (2 * prposcsc) + prpop + prpos + (2 * prpops);
                absprpo3(z,w) = (abs(prpo3(z,w)))^2;
                logprpo3(z,w) = (10 * log10(absprpo3(z,w)));
            elseif ((s > (0.5 * s_max)) && (d == 0))
                x = s - (0.5 * s_max);
                fTXlosRXn = 2 * pi * f * (x / c);
                prpolos = (1 / x) * exp(-1i * fTXlosRXn);
                x = hypot(d, s_max);
                fTXpsRXn = 2 * pi * f * (x / c);
                prpops = ((a / x) * exp(-1i * fTXpsRXn));
                x = hypot(d, (2 * s_max));
                fTXpspsRXn = 2 * pi * f * (x / c);
                prpopsps = (((a * a) / x) * exp(-1i * fTXpspsRXn));
                prpo3(z,w) = prpolos + (2 * prpops) + (2 * prpopsps);
                absprpo3(z,w) = (abs(prpo3(z,w)))^2;
                logprpo3(z,w) = (10 * log10(absprpo3(z,w)));
            elseif ((s < (0.5 * s_max)) && (d == 0))
                x = (0.5 * s_max) - s;
                fTXlosRXn = 2 * pi * f * (x / c);
                prpolos = (1 / x) * exp(-1i * fTXlosRXn);
                
                x = hypot(d, s_max);
                fTXpsRXn = 2 * pi * f * (x / c);
                prpops = ((a / x) * exp(-1i * fTXpsRXn));
                x = hypot(d, (2 * s_max));
                fTXpspsRXn = 2 * pi * f * (x / c);
                prpopsps = (((a * a) / x) * exp(-1i * fTXpspsRXn));
                prpo3(z,w) = prpolos + (2 * prpops) + (2 * prpopsps);
                absprpo3(z,w) = (abs(prpo3(z,w)))^2;
                logprpo3(z,w) = (10 * log10(absprpo3(z,w)));
            elseif (s < (0.5 * s_max) && (d > 0))
                los = hypot(d, ((0.5 * s_max) - s));
                fTXlosRXn = 2 * pi * f * (los / c);
                prpolos = (1 / (los)) * exp(-1i * fTXlosRXn);
                
                x = hypot(d, ((0.5 * s_max) - s));
                
                y = hypot(x, (2 * (h_max - h_anten)));
                fTXpsRXn = 2 * pi * f * (y / c);
                prpos = ((a / y) * exp(-1i * fTXpsRXn));
                y = hypot(x, (2 * h_anten));
                fTXpsRXn = 2 * pi * f * (y / c);
                prpop = ((a / y) * exp(-1i * fTXpsRXn));
                
                y = hypot(d , (s_max - ((0.5 * s_max) - s)));
                fTXpsRXn = 2 * pi * f * (y / c);
                prposc = ((a / y) * exp(-1i * fTXpsRXn));
                
                y = hypot(d , ((2 * s_max) - ((0.5 * s_max) - s)));
                fTXpsRXn = 2 * pi * f * (y / c);
                prposcsc = ((a / y) * exp(-1i * fTXpsRXn));
                
                y = hypot(x, (2 * h_max));
                fTXpsRXn = 2 * pi * f * (y / c);
                prposp = (((a * a) / y) * exp(-1i * fTXpsRXn));
                fTXpsRXn = 2 * pi * f * (y / c);
                prpops = (((a * a) / y) * exp(-1i * fTXpsRXn));

                prpo3(z,w) = prpolos + prpop + prpos + prposp + prpops + (2 * prposc) + (2 * prposcsc);
                absprpo3(z,w) = (abs(prpo3(z,w)))^2;
                logprpo3(z,w) = (10 * log10(absprpo3(z,w)));
            elseif ((s > (0.5 * s_max)) && (d > 0))
                los = hypot(d, (s - (0.5 * s_max)));
                fTXlosRXn = 2 * pi * f * (los / c);
                prpolos = (1 / (los)) * exp(-1i * fTXlosRXn);
                
                x = hypot(d, (s - (0.5 * s_max)));
                
                y = hypot(x, (2 * (h_max - h_anten)));
                fTXpsRXn = 2 * pi * f * (y / c);
                prpos = ((a / y) * exp(-1i * fTXpsRXn));
                y = hypot(x, (2 * h_anten));
                fTXpsRXn = 2 * pi * f * (y / c);
                prpop = ((a / y) * exp(-1i * fTXpsRXn));
                
                y = hypot(d , (s_max - (s - (0.5 * s_max))));
                fTXpsRXn = 2 * pi * f * (y / c);
                prposc = ((a / y) * exp(-1i * fTXpsRXn));
                
                y = hypot(d , ((2 * s_max) - (s - (0.5 * s_max))));
                fTXpsRXn = 2 * pi * f * (y / c);
                prposcsc = ((a / y) * exp(-1i * fTXpsRXn));
                
                y = hypot(x, (2 * h_max));
                fTXpsRXn = 2 * pi * f * (y / c);
                prposp = (((a * a) / y) * exp(-1i * fTXpsRXn));
                fTXpsRXn = 2 * pi * f * (y / c);
                prpops = (((a * a) / y) * exp(-1i * fTXpsRXn));
                
                prpo3(z,w) = prpolos + prpop + prpos + prposp + prpops + (2 * prposc) + (2 * prposc);
                absprpo3(z,w) = (abs(prpo3(z,w)))^2;
                logprpo3(z,w) = (10 * log10(absprpo3(z,w)));
            end
        end
    end
    end
%los = hypot(d, s);
% % los = sqrt((d.^2) + (s.^2));
% % 
% % fTXlosRXn = 2 * pi * f * (los ./ c);
% % prpolos = (1 ./ (los)) .* exp(-1i .* fTXlosRXn);
% % %x = hypot((0.9 * d_max), (2 * ((0.5 * s_max) - s)));
% % %fTXpsRXn = 2 * pi * f * (x / c);
% % %prpops = ((a / x) * exp(-1i * fTXpsRXn));
% % prpo3 = prpolos;% + (2 * prpops);
% % absprpo3 = (abs(prpo3)).^2;
% % logprpo3 = (10 * log10(absprpo3));

%%%figure
%plot([s3 d3], logprpo3)
%plot3([Z W], logprpo3(Z,W))
%plot([Z, W],logprpo3(Z,W))
%%%plot(logprpo3)

figure
[X, Y] = meshgrid(d3,s3);

mesh(X, Y,logprpo3)
xlim([0 d_max])
ylim([0 s_max])
colorbar

title('Wzgl�dny spadek mocy dla ca�ego pokoju')
xlabel('D�ugo�� [m]');
ylabel('Szeroko�� [m]');
zlabel('Moc [dB]');