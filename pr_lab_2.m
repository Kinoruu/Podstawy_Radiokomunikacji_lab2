f = 3.6 * 10^9;               %czêstotliwoœæ
c = 3 * 10^8;                 %prêdkoœæ œwiat³a
d_max = 2.8;                  %d³ugoœæ pokoju
s_max = 2.6;                  %szerokoœæ pokoju
h_max = 2.6;                  %wysokoœæ pokoju
step1 = 0.02;                 %krok obliczeñ dla zadania 1
step2 = 0.02;                 %krok obliczeñ dla zadania 2
d = [1:step1:d_max];          %zakres dla zadania 1
sx = [0:step2:(2 * s_max)];   %zakres dla zadania 2
a = -0.5;                     %wspó³czynnik odbicia
%%%%%%%%%%%%%%%%%  LOS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flos = -2 * pi * f * (d ./ c);

prpo = 1 ./ d;

prpolos = (abs(prpo)) .^ 2;

semilogx(d, prpolos)
grid
grid minor
title('Wzglêdny spadek mocy dla œcie¿ki LOS')
xlabel('Odleg³oœæ [m]');
ylabel('Moc [dB]');
%%%%%%%%%  Pojedyncze odbicie + LOS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dsc = hypot(d, s_max);
dsu = hypot(d, h_max);
dscp = (2 * d_max) - d;

fTXscRX = -2 * pi * f * (dsc ./ c); %œciana
fTXsuRX = -2 * pi * f * (dsu ./ c); %sufit
fTXscpRX = -2 * pi * f * (dscp ./ c); %œciana w linii LOS

prposc = ((a ./ dsc) .* exp(-j .* fTXscRX));
prposu = ((a ./ dsu) .* exp(-j .* fTXsuRX));
prposcp = ((a ./ dscp) .* exp(-j .* fTXscpRX));

sumaprpo1 = (abs(prpo + ...
    prposc + prposc + prposu + prposu + prposcp)).^2;

figure
semilogx(d, sumaprpo1)
grid
grid minor
title('Wzglêdny spadek mocy dla sumy œcie¿ki LOS i œcie¿ek powsta³ych na skutek pojedynczego odbicia')
xlabel('Odleg³oœæ [m]');
ylabel('Moc [dB]');
%%%%%%%%%  Pojedyncze, podówjne odbicie + LOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dscsc = hypot(d, (2 * s_max)); %œciana œciana
dsupo = hypot(d, (2 * h_max)); %sufit pod³oga / pod³oga sufit
dscscp = (3 * d_max) - d; % œciana œciana w lini LOS

fTXscscRX = -2 * pi * f * (dscsc ./ c); %œciana
fTXsupoRX = -2 * pi * f * (dsupo ./ c); %sufit
fTXscscpRX = -2 * pi * f * (dscscp ./ c); %œciana za plecami

prposcsc = (((a * a) ./ dscsc) .* exp(-j .* fTXscscRX));
prposupo = (((a * a) ./ dsupo) .* exp(-j .* fTXsupoRX));
prposcscp = (((a * a) ./ dscscp) .* exp(-j .* fTXscscpRX));

sumaprpo2 = (abs(prpo + ...
    prposc + prposc + prposu + prposu + prposcp + ...
    prposcsc + prposcsc + prposupo + prposupo + prposcscp)).^2;

figure
semilogx(d, sumaprpo2)
grid
grid minor
title('Wzglêdny spadek mocy dla sumy œcie¿ki LOS i œcie¿ek powsta³ych na skutek pojedynczego oraz podwójnego odbicia')
xlabel('Odleg³oœæ [m]');
ylabel('Moc [dB]');
%%%%  Pojedyncze, podówjne i potrójne odbicie + LOS  %%%%%%%%%%%%%%%%%%%%%%
dscscsc = hypot(d, (3 * s_max)); %œciana œciana œciana
dsuposu = hypot(d, (3 * h_max)); %sufit pod³oga sufit/ pod³oga sufit pod³oga
dscscscp = (4 * d_max) - d; % œciana œciana œciana w lini LOS

fTXscscscRX = -2 * pi * f * (dscscsc ./ c); %œciana
fTXsuposcRX = -2 * pi * f * (dsuposu ./ c); %sufit
fTXscscscpRX = -2 * pi * f * (dscscscp ./ c); %œciana za plecami

prposcscsc = (((a * a * a) ./ dscscsc) .* exp(-j .* fTXscscscRX));
prposuposu = (((a * a * a) ./ dsuposu) .* exp(-j .* fTXsuposcRX));
prposcscscp = (((a * a * a) ./ dscscscp) .* exp(-j .* fTXscscscpRX));

sumaprpo3 = (abs(prpo + ...
    prposc + prposc + prposu + prposu + prposcp + ...
    prposcsc + prposcsc + prposupo + prposupo + prposcscp + ...
    prposcscsc + prposcscsc + prposuposu + prposuposu + prposcscscp)).^2;

figure
semilogx(d, sumaprpo3)
grid
grid minor
title('Wzglêdny spadek mocy dla sumy œcie¿ki LOS i œcie¿ek powsta³ych na skutek pojedynczego, podwójnego oraz potrójnego odbicia')
xlabel('Odleg³oœæ [m]');
ylabel('Moc [dB]');
%%%%%%%%%%%  Wspólny wykres  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
semilogx(d, prpolos)
semilogx(d, sumaprpo1)
semilogx(d, sumaprpo2)
semilogx(d, sumaprpo3)
grid
grid minor
hold off
title('Porównanie wzglêdnych spadeków mocy dla kolejnych sum œcie¿ek')
xlabel('Odleg³oœæ [m]');
ylabel('Moc [dB]');
legend('LOS','LOS + 1','LOS + 1 + 2','LOS + 1 + 2 + 3')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 0:3
    if n == 0
        flos = -2 * pi * f * (d ./ c);
        prpolos = 1 ./ d;
        sumaprpo = prpolos;
    else
        dscn = hypot(d, (n * s_max));
        dsun = hypot(d, (n * h_max));
        dscpn = abs(((n+1) * d_max) - d);

        fTXscRXn = -2 * pi * f * (dscn ./ c); %œciana
        fTXsuRXn = -2 * pi * f * (dsun ./ c); %sufit
        fTXscpRXn = -2 * pi * f * (dscpn ./ c); %œciana w linii LOS

        prposcn = (((a^n) ./ dscn) .* exp(-j .* fTXscRXn));
        prposun = (((a^n) ./ dsun) .* exp(-j .* fTXsuRXn));
        prposcpn = (((a^n) ./ dscpn) .* exp(-j .* fTXscpRXn));

        sumaprpo = sumaprpo + ((2 * prposcn) + (2 * prposun) + prposcpn);
    end
end

abssumaprpo = (abs(sumaprpo)).^2;

figure
semilogx(d, abssumaprpo)
grid
grid minor
title(['Wzglêdny spadek mocy dla sumy œcie¿ki LOS i œcie¿ek powsta³ych na skutek do ', num2str(n) ,' odbiæ'])
xlabel('Odleg³oœæ [m]');
ylabel('Moc [dB]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      czêœæ druga
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prpo2 = [];
absprpo2 = [];
for s = 0:step2:(2 * s_max)
    z = fix((s * (1 / step2)) + 1);
    if s < (1.0625 * s_max)
        if s == (0.5 * s_max)
            prpolos = 1 / (0.9 * d_max);
            x = hypot((0.9 * d_max), s);
            fTXpsRXn = -2 * pi * f * (x / c);
            prpops = ((a / x) * exp(-j * fTXpsRXn));
            prpo2(z) = prpolos + (2 * prpops);
        elseif s < (0.5 * s_max)
            los = hypot((0.9 * d_max), ((0.5 * s_max) - s));
            prpolos = 1 / (los);
            x = hypot((0.9 * d_max), (2 * ((0.5 * s_max) - s)));
            fTXpsRXn = -2 * pi * f * (x / c);
            prpops = ((a / x) * exp(-j * fTXpsRXn));
            prpo2(z) = prpolos + (2 * prpops);
        elseif (s > (0.5 * s_max)) & (s < (s_max))
            los = hypot((0.9 * d_max), (s - (0.5 * s_max)));
            prpolos = 1 / (los);
            x = hypot((0.9 * d_max), (2 * (s - (0.5 * s_max))));
            fTXpsRXn = -2 * pi * f * (x / c);
            prpops = ((a / x) * exp(-j * fTXpsRXn));
            prpo2(z) = prpolos + (2 * prpops);
        elseif (s < (1.0625 * s_max)) & (s >= (s_max))
            los = hypot((0.9 * d_max), (s - (0.5 * s_max)));
            prpolos = 1 / (los);
            x = hypot((0.9 * d_max), (2 * (s - (0.5 * s_max))));
            fTXpsRXn = -2 * pi * f * (x / c);
            prpops = ((a / x) * exp(-j * fTXpsRXn));
            prpo2(z) = prpolos + (2 * prpops);
        end
    elseif s > (1.0625 * s_max)
        %{
        gamma1 = hypot((0.8 * d_max), (0.5 * s_max));
        gamma2 = hypot((0.1 * d_max), (s - s_max));
        s1s2 = hypot((0.9 * d_max), (s - (0.5 * s_max)));
        OBtroj = gamma1 + gamma2 + s1s2;
        p = OBtroj / 2;
        Ptroj = sqrt(p * (p - gamma1) * (p - gamma2) * (p - s1s2));
        h = (2 * Ptroj) / s1s2;
        lambda = c / f;
        v = h * sqrt((2 / lambda) * (s1s2 / (gamma1 * gamma2)));
        C = 6.9 + (20 * log(sqrt((v - 0.1)^2 + 1) + v - 0.1));
        prpo2(z) = C;
%}
        prpo2(z) = 0;
    end
    absprpo2(z) = (abs(prpo2(z)))^2;
end

figure
semilogx(sx, absprpo2)
grid
grid minor
title('Wzglêdny spadek mocy dla sumy œcie¿ki LOS i œcie¿ki powsta³ej na skutek jednego odbicia oraz dla cienia radarowego tylko dla dyfrakcji')
xlabel('Szerokoœæ [m]');
ylabel('Moc [dB]');