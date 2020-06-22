# exercises_LSN

01.1 
	compila: make
	esegui: ./main.exe
	cancella file output: ./clean.sh
	cancella file oggetto: make clean

01.2
        compila: make
	esegui: ./main.exe
	cancella file output: ./clean.sh
	cancella file oggetto: make clean

01.3
        compila: make
	esegui: ./main.exe
	cancella file output: ./clean.sh
	cancella file oggetto: make clean
        si possono cambiare gli input in input.dat

02.1
        compila: make
	esegui: ./main.exe
	cancella file output: ./clean.sh
	cancella file oggetto: make clean
        si possono cambiare gli input in input.dat

02.2
        compila: make
	esegui: ./main.exe
	cancella file output: ./clean.sh
	cancella file oggetto: make clean
        si possono cambiare gli input in input.dat

03.1
        compila: make
	esegui: ./main.exe
	cancella file output: ./clean.sh
	cancella file oggetto: make clean
        si possono cambiare gli input in parameters.dat

04
Il programma è nella cartella MolecularDynamics_NVE
compila: g++ -O3 -o MolDyn_NVE.exe MolDyn_NVE.cpp
esegui: ./MolDyn_NVE.exe
cancella file output: ./clean.sh
cancella file oggetto: make clean
per equilibrare la simulazione: ./equil_solid.sh
				./equil_liquid.sh
				./equil_gas.sh			
In ogni cartella relativa agli esercizi vi sono gli script per runnare il programma e 		copiare risultati
per runnare una simulazione già equilibrata:	./solid.sh
					   	./liquid.sh
						./gas.sh
I risultati sono nella cartella results results	

04.1
	./noequil_solid.sh
	./noequil_liquid.sh
	./noequil_gas.sh

05.1
	./1s_uniform.sh
	./2p_uniform.sh
	./1s_gauss.sh
	./2p_gauss.sh
	(o tutti insieme con il comando ./all.sh)

	compila: make
	esegui: ./main.exe
	cancella: make clean
	si possono cambiare gli input in input.dat
	si può cambiare la densità di probabilità da campionare in distribution.cpp
	Attenzione: se si cambia la densità di probabilità nel file distribution.cpp 
	è nescessario ricompilare il programma

06.1/ISING_1D
	./metropolis.sh
	./gibbs.sh

	compila: make
	esegui: ./Monte_Carlo_ISING_1D.exe
	cancella file output: ./clean.sh
	cancella file oggetto: make clean
	si possono cambiare gli input in input.dat
	Gli script sopra riportati eseguono il programma tante volte a diverse temperature.

07
Il programma è nella cartella MCNVT/Monte_Carlo_NVT
compila: make
esegui: ./Monte_Carlo_NVT.exe
cancella file output: ./clean.sh
cancella file oggetto: make clean
si possono cambiare gli input in input.dat
In ogni cartella relativa agli esercizi vi sono gli script per runnare il programma e 		copiare risultati

07.1
	ATTENZIONE: programmi abbastanza lunghi ~ qualche minuto
	./script_solid.sh 	
	./script_liquid.sh
	./script_gas.sh

07.4
	ATTENZIONE: programmi abbastanza lunghi ~ qualche minuto
	./script_solid.sh 	
	./script_liquid.sh 
	./script_gas.sh

08.2
	compila: make
	esegui: ./Variational_MC.exe
	cancella file output: ./clean.sh
	cancella file oggetto: make clean
	si possono cambiare gli input in input.dat
	ATTENZIONE: è necessario eliminare i file output prima di rifare una simulazione perché 	non vengono sovrascritti. 	

08.2_simulated_annealing
	ATTENZIONE: il programma non è pulitissimo, è solo servito per trovare il valore dei 		parametri mu e sigma che minimizzano l'energia.
	compila: make
	esegui: ./SA.exe
	cancella file output: ./clean.sh
	cancella file oggetto: make clean
	si possono cambiare gli input in input.dat
	si può cambiare l'annealing scheudule in schedule.dat

08.3/QMC_1D
	compila: make
	esegui: ./qmc1d
	cancella file output: ./clean.sh
	cancella file oggetto: make clean
	ATTENZIONE: qui è stata fatta una esecuzione manuale seguendo gli step:
	- in qmc1d.cpp scegliere nelle funzioni variationalWaveFunction e
	variationalWaveFunction_second i valori delle funzioni da cui parte l'algoritmo
	- scegliere se usare PIMC o PIGS e copiare il corrispondente file input in input.dat
	(nel caso di PIMC scegliere la temperatura inserendola in input.dat)
	- RICOMPILARE se è cambiata VariationalWaveFunction
	- eseguire il programma e salvare i file di output

09.1
	./circ.sh
	./square.sh	
	
	compila: make
	esegui: ./TSP.exe
	cancella file output: ./clean.sh
	cancella file oggetto: make clean
	si possono cambiare i dati in input su input.dat
	si possono cambiare le coordinate delle città su coordinates.dat
	si possono generare nuove coordinate random delle città con il comando:
	python cities.py
	nel file cities.py è possibile cambiare il seed del generatore random

	Per far partire una simulazione con un seed diverso:	./circ_diffseed.sh
								./square_diffseed.sh
								   

10.1
	./circ.sh
	./square.sh	
	
	compila: make
	esegui: ./TSP.exe
	cancella file output: ./clean.sh
	cancella file oggetto: make clean
	si possono cambiare i dati in input su input.dat
	si possono cambiare le coordinate delle città su coordinates.dat
	si può cambiare l'annealing schedule in schedule.dat
	il programma generate_schedule.py è stato utilizzato per generare una schedule.
	si possono generare nuove coordinate random delle città con il comando:
	python cities.py
	nel file cities.py è possibile cambiare il seed del generatore random 

10.2
	./circ.sh
	./square.sh	
	
	compila: make
	esegui: mpiexec -np 4 TSP.exe
	cancella file output: ./clean.sh
	cancella file oggetto: make clean
	si possono cambiare i dati in input su input.dat
	si possono cambiare le coordinate delle città su coordinates.dat
	si possono generare nuove coordinate random delle città con il comando:
	python cities.py
	nel file cities.py è possibile cambiare il seed del generatore random 

11 e 12 esercitazioni svolte direttamente su jupyter
