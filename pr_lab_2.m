d_max = 2.8;         %d�ugo�� pokoju
s_max = 2.6;         %szeroko�� pokoju
h_max = 2.6;         %wysoko�� pokoju
d = [1:0.05:2.6];    %zakres
a = -0.5;            %wsp�czynnik odbicia
%%%%%%%%%%%%%%%%%  LOS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flos = -2 * pi * 3.6 * 10^9 *(d./3*10^8);

prpo = 1./d;

prpolos = (abs(prpo)).^2;

semilogx(d, prpolos)
grid
grid minor
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 0:3
    if n == 0
        flos = -2 * pi * 3.6 * 10^9 *(d./3*10^8);
        prpolos = 1./d;
        prpo = prpolos;
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

        prpo = ((2 * prposcn) + (2 * prposun) + prposcpn);
    end
    sumaprpo = sumaprpo + prpo;
end

abssumaprpo = (abs(sumaprpo)).^2;

figure
semilogx(d, abssumaprpo)
grid
grid minor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


