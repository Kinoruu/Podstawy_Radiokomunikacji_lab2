f = 3.6 * 10^9;      %cz�stotliwo��
c = 3 * 10^8;        %pr�dko�� �wiat�a
d_max = 2.8;         %d�ugo�� pokoju
s_max = 2.6;         %szeroko�� pokoju
h_max = 2.6;         %wysoko�� pokoju
step1 = 0.02;         %krok oblicze� dla zadania 1
step2 = 0.02;         %krok oblicze� dla zadania 2
d = [1:step1:d_max];    %zakres dla zadania 1
sx = [0:step2:(2 * s_max)];    %zakres dla zadania 2
a = -0.5;            %wsp�czynnik odbicia
%%%%%%%%%%%%%%%%%  LOS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flos = -2 * pi * f * (d ./ c);

prpo = 1 ./ d;

prpolos = (abs(prpo)) .^ 2;

semilogx(d, prpolos)
grid
grid minor
title('Wzgl�dny spadek mocy dla �cie�ki LOS')
xlabel('Odleg�o�� [m]');
ylabel('Moc [dB]');
%%%%%%%%%  Pojedyncze odbicie + LOS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dsc = hypot(d, s_max);
dsu = hypot(d, h_max);
dscp = (2 * d_max) - d;

fTXscRX = -2 * pi * f * (dsc ./ c); %�ciana
fTXsuRX = -2 * pi * f * (dsu ./ c); %sufit
fTXscpRX = -2 * pi * f * (dscp ./ c); %�ciana w linii LOS

prposc = ((a ./ dsc) .* exp(-j .* fTXscRX));
prposu = ((a ./ dsu) .* exp(-j .* fTXsuRX));
prposcp = ((a ./ dscp) .* exp(-j .* fTXscpRX));

sumaprpo1 = (abs(prpo + ...
    prposc + prposc + prposu + prposu + prposcp)).^2;

figure
semilogx(d, sumaprpo1)
grid
grid minor
title('Wzgl�dny spadek mocy dla sumy �cie�ki LOS i �cie�ek powsta�ych na skutek pojedynczego odbicia')
xlabel('Odleg�o�� [m]');
ylabel('Moc [dB]');
%%%%%%%%%  Pojedyncze, pod�wjne odbicie + LOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dscsc = hypot(d, (2 * s_max)); %�ciana �ciana
dsupo = hypot(d, (2 * h_max)); %sufit pod�oga / pod�oga sufit
dscscp = (3 * d_max) - d; % �ciana �ciana w lini LOS

fTXscscRX = -2 * pi * f * (dscsc ./ c); %�ciana
fTXsupoRX = -2 * pi * f * (dsupo ./ c); %sufit
fTXscscpRX = -2 * pi * f * (dscscp ./ c); %�ciana za plecami

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
title('Wzgl�dny spadek mocy dla sumy �cie�ki LOS i �cie�ek powsta�ych na skutek pojedynczego oraz podw�jnego odbicia')
xlabel('Odleg�o�� [m]');
ylabel('Moc [dB]');
%%%%  Pojedyncze, pod�wjne i potr�jne odbicie + LOS  %%%%%%%%%%%%%%%%%%%%%%
dscscsc = hypot(d, (3 * s_max)); %�ciana �ciana �ciana
dsuposu = hypot(d, (3 * h_max)); %sufit pod�oga sufit/ pod�oga sufit pod�oga
dscscscp = (4 * d_max) - d; % �ciana �ciana �ciana w lini LOS

fTXscscscRX = -2 * pi * f * (dscscsc ./ c); %�ciana
fTXsuposcRX = -2 * pi * f * (dsuposu ./ c); %sufit
fTXscscscpRX = -2 * pi * f * (dscscscp ./ c); %�ciana za plecami

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
title('Wzgl�dny spadek mocy dla sumy �cie�ki LOS i �cie�ek powsta�ych na skutek pojedynczego, podw�jnego oraz potr�jnego odbicia')
xlabel('Odleg�o�� [m]');
ylabel('Moc [dB]');
%%%%%%%%%%%  Wsp�lny wykres  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
semilogx(d, prpolos)
semilogx(d, sumaprpo1)
semilogx(d, sumaprpo2)
semilogx(d, sumaprpo3)
grid
grid minor
hold off
title('Por�wnanie wzgl�dnych spadek�w mocy dla kolejnych sum �cie�ek')
xlabel('Odleg�o�� [m]');
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

        fTXscRXn = -2 * pi * f * (dscn ./ c); %�ciana
        fTXsuRXn = -2 * pi * f * (dsun ./ c); %sufit
        fTXscpRXn = -2 * pi * f * (dscpn ./ c); %�ciana w linii LOS

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
title(['Wzgl�dny spadek mocy dla sumy �cie�ki LOS i �cie�ek powsta�ych na skutek do ', num2str(n) ,' odbi�'])
xlabel('Odleg�o�� [m]');
ylabel('Moc [dB]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      cz�� druga
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
            prpo2(z) = prpolos + prpops;
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
title('Wzgl�dny spadek mocy dla sumy �cie�ki LOS i �cie�ki powsta�ej na skutek jednego odbicia oraz dla cienia radarowego tylko dla dyfrakcji')
xlabel('Szeroko�� [m]');
ylabel('Moc [dB]');