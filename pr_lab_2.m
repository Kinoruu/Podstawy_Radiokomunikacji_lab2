f = 3.6 * 10^9;               %czêstotliwoœæ
c = 3 * 10^8;                 %prêdkoœæ œwiat³a
d_max = 2.8;                  %d³ugoœæ pokoju
s_max = 2.6;                  %szerokoœæ pokoju
h_max = 2.6;                  %wysokoœæ pokoju
step1 = 0.0001;                 %krok obliczeñ dla zadania 1
step2 = 0.02;                 %krok obliczeñ dla zadania 2
d = [1:step1:d_max];          %zakres dla zadania 1
sx = [0:step2:(2 * s_max)];   %zakres dla zadania 2
a = -1;                     %wspó³czynnik odbicia
h_anten = 1;                  %wysokoœæ anteny
%%%%%%%%%%%%%%%%%  LOS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fTXlosRX = 2 * pi * f * (d ./ c);

prpo = (1 ./ d) .* exp(-j .* fTXlosRX);

prpolos = (abs(prpo)) .^ 2;

plot(d, prpolos)
grid
grid minor
title('Wzglêdny spadek mocy dla œcie¿ki LOS')
xlabel('Odleg³oœæ [m]');
ylabel('Moc [dB]');
%%%%%%%%%  Pojedyncze odbicie + LOS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dsc = hypot(d, s_max);
dsu = hypot(d, (2 * (h_max - 1)));
dpo = hypot (d, (2 * h_anten));
dscp = (2 * d_max) - d;

fTXscRX = 2 * pi * f * (dsc ./ c); %œciana
fTXsuRX = 2 * pi * f * (dsu ./ c); %sufit
fTXpoRX = 2 * pi * f * (dpo ./ c); %sufit
fTXscpRX = 2 * pi * f * (dscp ./ c); %œciana w linii LOS

prposc = ((a ./ dsc) .* exp(-j .* fTXscRX));
prposu = ((a ./ dsu) .* exp(-j .* fTXsuRX));
prpopo = ((a ./ dsu) .* exp(-j .* fTXpoRX));
prposcp = ((a ./ dscp) .* exp(-j .* fTXscpRX));

sumaprpo1 = (abs(prpo + ...
    prposc + prposc + prposu + prpopo + prposcp)).^2;

logsumaprpo1 = (10 * log10(sumaprpo1));

figure
plot(d, logsumaprpo1)
grid
grid minor
title('Wzglêdny spadek mocy dla sumy œcie¿ki LOS i œcie¿ek powsta³ych na skutek pojedynczego odbicia')
xlabel('Odleg³oœæ [m]');
ylabel('Moc [dB]');
%%%%%%%%%  Pojedyncze, podówjne odbicie + LOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dscsc = hypot(d, (2 * s_max)); %œciana œciana
dsupo = hypot(d, (2 * h_max)); %sufit pod³oga / pod³oga sufit
dscscp = (3 * d_max) - d; % œciana œciana w lini LOS

fTXscscRX = 2 * pi * f * (dscsc ./ c); %œciana
fTXsupoRX = 2 * pi * f * (dsupo ./ c); %sufit
fTXscscpRX = 2 * pi * f * (dscscp ./ c); %œciana za plecami

prposcsc = (((a * a) ./ dscsc) .* exp(-j .* fTXscscRX));
prposupo = (((a * a) ./ dsupo) .* exp(-j .* fTXsupoRX));
prposcscp = (((a * a) ./ dscscp) .* exp(-j .* fTXscscpRX));

sumaprpo2 = (abs(prpo + ...
    prposc + prposc + prposu + prpopo + prposcp + ...
    prposcsc + prposcsc + prposupo + prposupo + prposcscp)).^2;

logsumaprpo2 = (10 * log10(sumaprpo2));

figure
plot(d, logsumaprpo2)
grid
grid minor
title('Wzglêdny spadek mocy dla sumy œcie¿ki LOS i œcie¿ek powsta³ych na skutek pojedynczego oraz podwójnego odbicia')
xlabel('Odleg³oœæ [m]');
ylabel('Moc [dB]');
%%%%  Pojedyncze, podówjne i potrójne odbicie + LOS  %%%%%%%%%%%%%%%%%%%%%%
dscscsc = hypot(d, (3 * s_max)); %œciana œciana œciana
dsuposu = hypot(d, (3 * h_max)); %sufit pod³oga sufit/ pod³oga sufit pod³oga
dscscscp = (4 * d_max) - d; % œciana œciana œciana w lini LOS

fTXscscscRX = 2 * pi * f * (dscscsc ./ c); %œciana
fTXsuposcRX = 2 * pi * f * (dsuposu ./ c); %sufit
fTXscscscpRX = 2 * pi * f * (dscscscp ./ c); %œciana za plecami

prposcscsc = (((a * a * a) ./ dscscsc) .* exp(-j .* fTXscscscRX));
prposuposu = (((a * a * a) ./ dsuposu) .* exp(-j .* fTXsuposcRX));
prposcscscp = (((a * a * a) ./ dscscscp) .* exp(-j .* fTXscscscpRX));

sumaprpo3 = (abs(prpo + ...
    prposc + prposc + prposu + prpopo + prposcp + ...
    prposcsc + prposcsc + prposupo + prposupo + prposcscp + ...
    prposcscsc + prposcscsc + prposuposu + prposuposu + prposcscscp)).^2;

logsumaprpo3 = (10 * log10(sumaprpo3));

figure
plot(d, logsumaprpo3)
grid
grid minor
title('Wzglêdny spadek mocy dla sumy œcie¿ki LOS i œcie¿ek powsta³ych na skutek pojedynczego, podwójnego oraz potrójnego odbicia')
xlabel('Odleg³oœæ [m]');
ylabel('Moc [dB]');
%%%%%%%%%%%  Wspólny wykres  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
plot(d, prpolos)
plot(d, logsumaprpo1)
plot(d, logsumaprpo2)
plot(d, logsumaprpo3)
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
    elseif n == 1
        dscn = hypot(d, (n * s_max));
        dsun = hypot(d, (2 * (h_max - 1)));
        dpon = hypot (d, (2 * h_anten))
        dscpn = abs(((n+1) * d_max) - d);

        fTXscRXn = 2 * pi * f * (dscn ./ c); %œciana
        fTXsuRXn = 2 * pi * f * (dsun ./ c); %sufit
        fTXpoRXn = 2 * pi * f * (dpon ./ c);
        fTXscpRXn = 2 * pi * f * (dscpn ./ c); %œciana w linii LOS

        prposcn = (((a^n) ./ dscn) .* exp(-j .* fTXscRXn));
        prposun = (((a^n) ./ dsun) .* exp(-j .* fTXsuRXn));
        prposun = (((a^n) ./ dpon) .* exp(-j .* fTXpoRXn));
        prposcpn = (((a^n) ./ dscpn) .* exp(-j .* fTXscpRXn));

        sumaprpo = sumaprpo + ((2 * prposcn) + prposun + prposun + prposcpn);
    else
        dscn = hypot(d, (n * s_max));
        dsun = hypot(d, (n * h_max));
        dscpn = abs(((n+1) * d_max) - d);

        fTXscRXn = 2 * pi * f * (dscn ./ c); %œciana
        fTXsuRXn = 2 * pi * f * (dsun ./ c); %sufit
        fTXscpRXn = 2 * pi * f * (dscpn ./ c); %œciana w linii LOS

        prposcn = (((a^n) ./ dscn) .* exp(-j .* fTXscRXn));
        prposun = (((a^n) ./ dsun) .* exp(-j .* fTXsuRXn));
        prposcpn = (((a^n) ./ dscpn) .* exp(-j .* fTXscpRXn));

        sumaprpo = sumaprpo + ((2 * prposcn) + (2 * prposun) + prposcpn);
    end
end

abssumaprpo = (abs(sumaprpo)).^2;

logabssumaprpo = (10 * log10(abssumaprpo));

figure
plot(d, (10 * log10(abssumaprpo)))
grid
grid minor
title(['Wzglêdny spadek mocy dla sumy œcie¿ki LOS i œcie¿ek powsta³ych na skutek do ', num2str(n) ,' odbiæ'])
xlabel('Odleg³oœæ [m]');
ylabel('Moc [dB]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      czêœæ druga
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prpo2 = [];
logprpo2 = [];
absprpo2 = [];
for s = 0:step2:(2 * s_max)
    z = fix((s * (1 / step2)) + 1);
    if s < (1.0625 * s_max)
        if s == (0.5 * s_max)
            prpolos = 1 / (0.9 * d_max);
            x = hypot((0.9 * d_max), s);
            fTXpsRXn = 2 * pi * f * (x / c);
            prpops = ((a / x) * exp(-j * fTXpsRXn));
            prpo2(z) = prpolos + (2 * prpops);
            absprpo2(z) = (abs(prpo2(z)))^2
            logprpo2(z) = (10 * log10(absprpo2(z)));
        elseif s < (0.5 * s_max)
            los = hypot((0.9 * d_max), ((0.5 * s_max) - s));
            prpolos = 1 / (los);
            x = hypot((0.9 * d_max), (2 * ((0.5 * s_max) - s)));
            fTXpsRXn = 2 * pi * f * (x / c);
            prpops = ((a / x) * exp(-j * fTXpsRXn));
            prpo2(z) = prpolos + (2 * prpops);
            absprpo2(z) = (abs(prpo2(z)))^2
            logprpo2(z) = (10 * log10(absprpo2(z)));
        elseif (s > (0.5 * s_max)) & (s < (s_max))
            los = hypot((0.9 * d_max), (s - (0.5 * s_max)));
            prpolos = 1 / (los);
            x = hypot((0.9 * d_max), (2 * (s - (0.5 * s_max))));
            fTXpsRXn = 2 * pi * f * (x / c);
            prpops = ((a / x) * exp(-j * fTXpsRXn));
            prpo2(z) = prpolos + (2 * prpops);
            absprpo2(z) = (abs(prpo2(z)))^2
            logprpo2(z) = (10 * log10(absprpo2(z)));
        elseif (s < (1.0625 * s_max)) & (s >= (s_max))
            los = hypot((0.9 * d_max), (s - (0.5 * s_max)));
            prpolos = 1 / (los);
            x = hypot((0.9 * d_max), (2 * (s - (0.5 * s_max))));
            fTXpsRXn = 2 * pi * f * (x / c);
            prpops = ((a / x) * exp(-j * fTXpsRXn));
            prpo2(z) = prpolos + (2 * prpops);
            absprpo2(z) = (abs(prpo2(z)))^2
            logprpo2(z) = (10 * log10(absprpo2(z)));
        end
    elseif s >= (1.0625 * s_max)
        gamma1 = hypot((0.8 * d_max), (0.5 * s_max));
        gamma2 = hypot((0.1 * d_max), (s - s_max));
        s1s2 = hypot((0.9 * d_max), (s - (0.5 * s_max)));
        OBtroj = gamma1 + gamma2 + s1s2;
        p = OBtroj / 2;
        Ptroj = sqrt(p * (p - gamma1) * (p - gamma2) * (p - s1s2));
        h = (2 * Ptroj) / s1s2;
        lambda = c / f;
        v = h * sqrt((2 / lambda) * (s1s2 / (gamma1 * gamma2)));
        C = 6.9 + (20 * log10(sqrt((v - 0.1)^2 + 1) + v - 0.1));
        logprpo2(z) = -C;
    end
    %absprpo2(z) = (logprpo2(z));
end

figure
plot(sx, logprpo2)
grid
grid minor
title('Wzglêdny spadek mocy dla sumy œcie¿ki LOS i œcie¿ki powsta³ej na skutek jednego odbicia oraz dla cienia radarowego tylko dla dyfrakcji')
xlabel('Szerokoœæ [m]');
ylabel('Moc [dB]');