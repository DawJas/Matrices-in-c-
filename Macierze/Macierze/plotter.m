clear all
close all
JacobiCzas = csvread("JacobiCzas.csv");
GaussCzas = csvread("GaussCzas.csv");
LUCzas = csvread("LUCzas.csv");
Rozmiary = csvread("Rozmiary.csv");


figure("Name","Metoda LU",'NumberTitle','off');
semilogy(Rozmiary, LUCzas)
ylabel("czas rozwiązania [ms]");
xlabel("rozmiar macierzy");
title("czas rozwiazania rownania macierzowego metoda LU");

figure("Name","Metoda Gaussa-Seidla",'NumberTitle','off');
semilogy(Rozmiary, GaussCzas)
ylabel("czas rozwiązania [ms]");
xlabel("rozmiar macierzy");
title("czas rozwiazania rownania macierzowego metoda Gaussa-Seidla");

figure("Name","Metoda Jacobiego",'NumberTitle','off');
semilogy(Rozmiary, JacobiCzas)
ylabel("czas rozwiązania [ms]");
xlabel("rozmiar macierzy");
title("czas rozwiazania rownania macierzowego metoda Jacobiego");