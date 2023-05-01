clear all
close all
% JacobiCzas = csvread("JacobiCzas.csv")./1000;
% GaussCzas = csvread("GaussCzas.csv")./1000;
% LUCzas = csvread("LUCzas.csv")./1000;
% Rozmiary = csvread("Rozmiary.csv");
resBJacobi = csvread("resBJacobi.csv");
resBGauss = csvread ("resBGauss.csv");
resCJacobi = csvread("resCJacobi.csv");
resCGauss = csvread("resCGauss.csv");


% figure("Name","Metoda LU",'NumberTitle','off');
% plot(Rozmiary, LUCzas);
% ylabel("Czas rozwiązania [s]");
% xlabel("Rozmiar macierzy");
% title("Czas rozwiazania rownania macierzowego metoda LU");
% saveas(gcf, 'LU_Czas.png');
% 
% figure("Name","Czas rozwiazania rownania macierzowego metodami iteracyjnymi",'NumberTitle','off');
% plot(Rozmiary, GaussCzas);
% ylabel("Czas rozwiązania [s]");
% xlabel("Rozmiar macierzy");
% title("Czas rozwiazania rownania macierzowego metodami iteracyjnymi");
% hold on
% plot(Rozmiary, JacobiCzas);
% legend("Gauss-Seidl", "Jacobi");
% saveas(gcf, 'Jacobi_Gauss_Czas.png');

%%%%%%%%%% NORMY BLEDOW %%%%%%%%%%%
figure("Name","Norma bledu rezydualnego w zad B",'NumberTitle','off');
semilogy(resBJacobi);
ylabel("Norma bledu reydualnego");
xlabel("nr iteracji");
title("Norma bledu rezydualnego w zad B");
hold on
semilogy(resBGauss);
legend("Jacobi", "Gauss");
saveas(gcf, 'res_B.png');

%%%%% NORMY BLEDOW NIEZBIEZNYCH %%%%%%

figure("Name","Norma bledu rezydualnego w zad C",'NumberTitle','off');
semilogy(resCJacobi);
ylabel("Norma bledu reydualnego");
xlabel("nr iteracji");
title("Norma bledu resydualnego w zad C");
hold on
semilogy(resCGauss);
legend("Jacobi", "Gauss");
saveas(gcf, 'res_C.png');

