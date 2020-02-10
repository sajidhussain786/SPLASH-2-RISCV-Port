# SPLASH-2-RISCV-Port
The RISCV Port of SPLASH-2 benchmark suit



# 1--> Compile
	   in Makefile.config
	   set CC := to your riscv-gcc compiler
	   set BASEDIR :=  to the directory where this makefile is located. 



# 2--> Run

	Detail run instructions are given in the Readme file in the directory of each benchmark. Different sets of intputs are available to be used. For testing i only used the recommended inputs for running Splash-3 in a  Qemu-riscv64 simulator. I have used 4 processor core in each run and the exetuion time is provided in each README file just for example. The inputs used are based on the original Splash-2 characterization paper by Woo et al. [1].

      For the quick run use the following commands. 
      \# represents the number of processor cores/threads and RISCV-SIM represent your REISCV Simulator. 
	In my case #=4 and RISCV-SIM = qemu-riscv64

	APPS:
		RISCV-SIM ./BARNES < inputs/n16384-p#
		RISCV-SIM ./FMM < inputs/input.#.16384
		RISCV-SIM ./OCEAN -p# -n258
		RISCV-SIM ./RADIOSITY -p # -ae 5000 -bf 0.1 -en 0.05 -room -batch
		RISCV-SIM ./RAYTRACE -p# -m64 inputs/car.env
		RISCV-SIM ./VOLREND # inputs/head 8
		RISCV-SIM ./WATER-NSQUARED < inputs/n512-p#
		RISCV-SIM ./WATER-SPATIAL < inputs/n512-p#
	KERNELS:
		RISCV-SIM ./CHOLESKY -p# < inputs/tk15.O
		RISCV-SIM ./FFT -p# -m16
		RISCV-SIM ./LU -p# -n512
		RISCV-SIM ./RADIX -p# -n1048576




[1] Steven Cameron Woo, Moriyoshi Ohara, Evan Torrie, Jaswinder Pal Singh, and
Anoop Gupta. 1995. The SPLASH-2 programs: characterization and methodological
considerations. In Proceedings of the 22nd annual international symposium on
Computer architecture (ISCA '95). ACM, New York, NY, USA, 24-36.
DOI=http://dx.doi.org/10.1145/223982.223990
