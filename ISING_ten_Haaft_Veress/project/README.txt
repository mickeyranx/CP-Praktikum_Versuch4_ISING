
-------------------------------------------------
1.                  Code
-------------------------------------------------
-die Methoden für die Simulationen sind aufgeteilt in main.cpp (mit main Methode) und Simulation.cpp (mit header file Simulation.h)
-in der Regel haben wir in Simulation.cpp die Algorithmen umgesetzt (z.B Draw, Sweep, Berechnung der Energieänderung, etc.)
  und in main.cpp die äußere Logik (z.B Thermalisieren -> Messen -> Observablen berechnen)

-die Methodenaufrufe für Aufgaben 1 und 2 sind auskommentiert in der Main methode
-für Aufgaben 4b und 4c haben wir seperatre Methoden geschrieben
-ansonsten lässt sich eine Simulation starten über single_measurements() wo dann die beispielsweise die Startparameter gewählt werden können

-------------------------------------------------
2.     Kompilieren und Starten in Windows
-------------------------------------------------
-Da wir in Visual Studio (nicht VS Code) gearbeitet ist dort bereits Compiler, Build tools und weiteres dabei
-Es wird aber Visual Studio (und C++ verion mindestens 20) benötigt um die Projektmappe (".sln"-Datei) zu öffnen

-------------------------------------------------
2.     Kompilieren und Starten in Linux 
-------------------------------------------------
-alternativ kann man in Linux oder WSL, sofern man make und g++ (Version die C++20 kompielieren kann) hat 
   das Projekt über die makefile kompilieren und builden

bash ->
make
bash ->
./run_ising
wenn "hello world" ausgegeben ist wurde alles richtig kompiliert


