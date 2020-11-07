d_max = 2.8;         %d�ugo�� pokoju
s_max = 2.6;         %szeroko�� pokoju
h_max = 2.6;         %wysoko�� pokoju
d = [1:0.02:2.6];    %zakres
a = -0.5;            %wsp�czynnik odbicia
%%%%%%%%%%%%%%%%%  LOS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flos = -2 * pi * 3.6 * 10^9 *(d./3*10^8);

prpo = 1./d;

prpolos = (abs(prpo)).^2;

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

fTXscRX = -2 * pi * 3.6 * 10^9 *(dsc./3*10^8); %�ciana
fTXsuRX = -2 * pi * 3.6 * 10^9 *(dsu./3*10^8); %sufit
fTXscpRX = -2 * pi * 3.6 * 10^9 *(dscp./3*10^8); %�ciana w linii LOS

prposc = ((a./dsc).*exp(-j.*fTXscRX));
prposu = ((a./dsu).*exp(-j.*fTXsuRX));
prposcp = ((a./dscp).*exp(-j.*fTXscpRX));

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
dscsc = hypot(d, 2 * s_max); %�ciana �ciana
dsupo = hypot(d, 2 * h_max); %sufit pod�oga / pod�oga sufit
dscscp = (3 * d_max) - d; % �ciana �ciana w lini LOS

fTXscscRX = -2 * pi * 3.6 * 10^9 *(dscsc./3*10^8); %�ciana
fTXsupoRX = -2 * pi * 3.6 * 10^9 *(dsupo./3*10^8); %sufit
fTXscscpRX = -2 * pi * 3.6 * 10^9 *(dscscp./3*10^8); %�ciana za plecami

prposcsc = (((a*a)./dscsc).*exp(-j.*fTXscscRX));
prposupo = (((a*a)./dsupo).*exp(-j.*fTXsupoRX));
prposcscp = (((a*a)./dscscp).*exp(-j.*fTXscscpRX));

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
dscscsc = hypot(d, 3 * s_max); %�ciana �ciana �ciana
dsuposu = hypot(d, 3 * h_max); %sufit pod�oga sufit/ pod�oga sufit pod�oga
dscscscp = (4 * d_max) - d; % �ciana �ciana �ciana w lini LOS

fTXscscscRX = -2 * pi * 3.6 * 10^9 *(dscscsc./3*10^8); %�ciana
fTXsuposcRX = -2 * pi * 3.6 * 10^9 *(dsuposu./3*10^8); %sufit
fTXscscscpRX = -2 * pi * 3.6 * 10^9 *(dscscscp./3*10^8); %�ciana za plecami

prposcscsc = (((a*a*a)./dscscsc).*exp(-j.*fTXscscscRX));
prposuposu = (((a*a*a)./dsuposu).*exp(-j.*fTXsuposcRX));
prposcscscp = (((a*a*a)./dscscscp).*exp(-j.*fTXscscscpRX));

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
        flos = -2 * pi * 3.6 * 10^9 *(d./3*10^8);
        prpolos = 1./d;
        sumaprpo = prpolos;
    else
        dscn = hypot(d, n * s_max);
        dsun = hypot(d, n * h_max);
        dscpn = abs(((n+1) * d_max) - d);

        fTXscRXn = -2 * pi * 3.6 * 10^9 * (dscn ./ 3 * 10^8); %�ciana
        fTXsuRXn = -2 * pi * 3.6 * 10^9 * (dsun ./ 3 * 10^8); %sufit
        fTXscpRXn = -2 * pi * 3.6 * 10^9 * (dscpn ./3 * 10^8); %�ciana w linii LOS

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


