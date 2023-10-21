clear % zbrišem shranjene spremenljivke v workspace-u
clc % pobrišem komandno okno

% generiram imputno okno, kamor vpišem št. naključnih točk
st_tock = input("Število naključnih točk: ");

% generiram kote fi od 0 do 2pi(100 točk)
fi = linspace(0,2*pi,100);
%definiram anonimno fukcijo, ki izračuna točke na loku krožnice
lok_krog = @(fi) [sin(fi'),cos(fi')];
%koordinate krožnice
kroznica = lok_krog(fi);

% kličem funkcijo mcc_pi-shranim koordinate točk znotraj kroga in kvadrata
[ktkrog,ktkvad] = mcc_pi(st_tock);
%vrednosti pi-ja in napake shranim
[pi_value,pi_napaka] = area_pi(ktkrog,ktkvad);

% izpišem vrednost pi in napake
fprintf('Vrednost pi: %.4f\nNapaka pi: %.4f\n',pi_value,pi_napaka);

% vizualizacija
figure;
hold on;

% točke v kvadratu
plot(ktkvad(:,1),ktkvad(:,2),'b*');
% točke v krogu
plot(ktkrog(:,1),ktkrog(:,2),'ro');
% krožnica
plot(kroznica(:,1),kroznica(:,2),'k-','LineWidth',2);
% enakomerni razmik x in y osi
axis equal
% označba grafa
title("Metoda Monte Carlo")
xlabel("x-koordinate točk")
ylabel("y-koordiate točk")
legend('točke znotraj kvadrata', 'točke znotraj kroga','krožnica');


hold off;

% definiram programsko funkcijo area_pi
function [pi_value,pi_napaka] = area_pi(ktkrog,ktkvad)
    pi_value = 4 * size(ktkrog)/size(ktkvad);
    pi_napaka = pi_value - pi;
end