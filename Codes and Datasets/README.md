# Optimization for Data Science Project

Contenuto cartella:

- RadialKer.m -> funzione kernel radiale
- Polyker.m -> funzione kernel polinomiale
- main_energy.m  -> caricamento dataset e chiamata alle funzioni main.m e mainL1.m
- main_housing.m -> caricamento dataset e chiamata alle funzioni main.m e mainL1.m
- main.m -> definizione matrice Q problema e chiamata funzione algoritmo risolutivo PGM.m e funzione Quadprog e generazione grafici errore
- PGM.m -> chiamata projection_routine e stampa statistiche per ogni iterazione
- projection_routine.m -> funzione che definisce la direzione con eventuale chiamata alla funzione AP.m ed aggiorna l'iterata corrente tramite chiamata a funzione exact_line_search.m   
- AP.m -> funzione che definisce la proiezione dell'antigradiente per la risoluzione del problema CQKnP
- compute_direction.m -> funzione che definisce la direzione basandosi sui vincoli del punto corrente e il breakpoint
- exact_line_search.m -> funzione che definisce il passo adattivo alfa
- mainL1.m -> definizione matrice Q problema e chiamata funzione algoritmo risolutivo con stepsize=1/L1 PGML1.m e funzione Quadprog e generazione grafici errore
- PGML1.m -> chiamata projection_routine_L1 e stampa statistiche per ogni iterazione
- projection_routine_L1.m -> funzione che definisce la direzione con eventuale chiamata alla funzione AP.m ed aggiorna l'iterata corrente con passo costante alfa=1/L1    