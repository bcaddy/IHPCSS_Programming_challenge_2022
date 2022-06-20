# IHPCSS'22 Challenge #

You are taking part to the [International High-Performance Computing Summer School](https://ss22.ihpcss.org) programming challenge? That's where it starts!

<img src="https://github.com/capellil/IHPCSS_Programming_challenge_2022/blob/main/images/in_a_nutshell.png">

## Table of contents ##

* [What is the challenge?](#what-is-the-challenge)
* [Is the programming challenge for me?](#is-the-programming-challenge-for-me)
* [What is this repository for?](#what-is-this-repository-for)
* [How do I get set up?](#how-do-i-get-set-up)
  * [Download](#download)
  * [Compilation](#compilation)
  * [Submit](#submit)
  * [Verify](#verify)
* [What kind of optimisations are not allowed?](#what-kind-of-optimisations-are-not-allowed)
* [Send your solution to the competition](#send-your-solution-to-the-competition)
* [Whom do I talk to?](#whom-do-i-talk-to)
* [Acknowledgments](#acknowledgments)

## What is the challenge? ##

This challenge introduces a simple problem: take a metal plate, put a flame below it and simulate the temperature propagation across that metal plate. For simplicity, we assume that the area touched by the flame will always be at the same temperature (arbitrarily put at 50 degrees; yes it is a tiny flame). The rest of the metal plate however will start heating up as the heat propagates outwards the lighter zone.

We represent the metal plate as a 2D grid, and every iteration we calculate the temperature variation in each grid cell by averaging the temperatures from the previous iteration in the neighbouring cells. The challenge is that you have 30 seconds to process as many iterations as possible. You therefore have to optimise the code with the techniques you will learn in the IHPCSS because the faster your iterations, the more you can execute under 30 seconds.

[Go back to table of contents](#table-of-contents)
## Is the programming challenge for me? ##
Short answer? Yes. From experience there are two audiences, equally welcome:
1) **Participants who want to do it exclusively for fun, looking primarily for no pressure or judgement**. If you consider yourself in this category, know that:
- there is no registration needed to participate
- there is no list of who participates or not
- there is no obligation to even participate at all
- you work on it when you want, if at some point you no longer have time or no longer want to, you simply stop
- there is no deadline by which teams must be officially declared (until you actually submit, if you submit)
- sending your code at the end is optional
- if you do send it, because you still want to see if you happen to have the fastest code, there is no ranking or name disclosed beyond those of the winning team. In other words, either your name is disclosed and it means your team won, or nobody will know you even submitted
- you are allowed to send your code multiple times, for instance if found additional optimisations or bug fixes after submitting already. However, only the last submission received is considered
2) **Participants who are in a bit more of a competition mindset**. If you consider yourself in this category, know that in addition to the above:
- indeed the code has certain areas easy to optimise, but there are also inefficiencies in it that are sneakier to spot. Maybe you will be able to find the one nobody else found? ... or missed the one other teams found?
- C, FORTRAN-90, MPI, OpenMP, OpenACC etc... That's a lot of technologies and standards right there, are you sure your team is the best in your category?
- if you want extra challenge, you can always participate to both tracks: trying to develop the fastest CPU code, and the fastest GPU code.
- should there be a tie between two or more teams, the submission time will be the tiebreaker: the team in tie having submitted first will win. So, developing the fastest code is no garantee of winning, being the quickest at developing the fastest code is.
- and of course there are trophies to bring home, if your team is the fastest that is ;)

In 2021, the programming challenge got an overall feedback score of 88%. It therefore seems that both audiences had great fun doing it, why not giving it a try? :D

[Go back to table of contents](#table-of-contents)
## What is this repository for? ##

* You will find here everything you need; source codes, makefiles, documentation, scripts, tests... From compilation to submission, through run and verifications, you are covered.
* You will find all those for both tracks:
  * Classic parallel programming: MPI + OpenMP
  * Accelerator: MPI + OpenACC
* Each of these is available in C and FORTRAN 90
* It provides you with a pre-setup experimental protocol; it makes sure contestants compete in the same conditions and allows to compare experiments fairly.

[Go back to table of contents](#table-of-contents)
## How do I get set up? ##
### Download ###
All you have to do is clone this repository: ```git clone https://github.com/capellil/IHPCSS_Programming_challenge_2022.git```.

Note that you are strongly encouraged to work on the source files provided instead of making copies. To keep it short, you will discover in the sections below that multiple scripts have been written to make your life (much) easier (makefile, submitting to compute nodes etc..). However, these scripts are designed to work with the files provided, not arbitrary copies you could make.

[Go back to table of contents](#table-of-contents)
### Compilation ###
Due to some surprises from the `nvhpc` module, there are two scripts available; one that will happily compile all CPU codes, and one for GPU codes.
- ```./compile_cpu_versions``` for the CPU ones
- ```./compile_gpu_versions``` for the GPU ones
 
If the right module is not loaded, it will complain and will give you the command to issue before trying again.

What happens behind the scene?

As you will quickly see, there is one folder for C source codes, one for FORTRAN source codes. Inside, each version has a specific file:

| Model | C version | FORTRAN version |
|-------|-----------|-----------------|
| MPI + OpenMP | cpu.c | cpu.F90 |
| MPI + OpenACC | gpu.c | gpu.F90 |

And of course, you modify the file corresponding to the combination you want to work on. No need to make a copy, work on the original file, everything is version controlled remember.

[Go back to table of contents](#table-of-contents)
### Submit ###
(***Note**: Jobs submitted with this script will use the corresponding reservation queue for big jobs.*)

A script has been written for you to easily submit your work to Bridges2 via SLURM:

```./submit.sh LANGUAGE IMPLEMENTATION SIZE OUTPUT_FILE```

The parameters are always the same:
* LANGUAGE = ```c``` | ```f```
* IMPLEMENTATION = ```cpu``` | ```gpu```
* SIZE = ```small``` | ```big```

How does it work? As you have probably seen, there is a ```slurm_scripts``` folder. It contains two SLURM submission scripts for each version (OpenMP + MPI, OpenACC + MPI etc...): one for the small grid, one for the big grid. That allows each SLURM script to be tailored (number of nodes, type of nodes, walltime...) for the implementation and size demanded.

Examples:
* to submit the C version of the CPU code on the small dataset: ```./submit.sh c cpu small my_output.txt```
* to submit the FORTRAN version of the GPU code on the big dataset: ```./submit.sh f gpu big my_output.txt```

[Go back to table of contents](#table-of-contents)
### Verify ###
The correctness of your code will be evaluated using the temperature change observed throughout iterations. Once you have a file containing the output of your program, you can check the correctness by using the ```verify.sh``` as follows:

```./verify.sh LANGUAGE IMPLEMENTATION SIZE FILE_TO_VERIFY```

The parameters are always the same:
* LANGUAGE = ```c``` | ```f```
* IMPLEMENTATION = ```cpu``` | ```gpu```
* SIZE = ```small``` | ```big```

This will automatically fetch the corresponding reference file in the ```reference``` folder and compare all temperature changes with the ones contained in the file you specified. If your executable has run more iterations than the original executable, only the iterations in common are compared.

Examples:
* to verify the C version of the CPU code on the small dataset: ```./verify.sh c cpu small my_output.txt```
* to verify the FORTRAN version of the GPU code on the big dataset: ```./verify.sh f gpu big my_output.txt```

[Go back to table of contents](#table-of-contents)
## What kind of optimisations are not allowed? ##

The idea of this programming challenge is for you to practice on the different techniques you will learn. To keep you focused on those, certain optimisations that are not within the scope of this challenge are not allowed:

* Changing the compilation process (that is: using different compilers, compiler flags, external libraries etc...). The point in this challenge is not for you to read hundreds of pages of documentation to find an extra flag people may have missed.
* Changing the running process; provided scripts already use a sensible configuration. Again, you only have a few days; the objective of the hybrid challenge is for you to play with the code, not spend hours defining the best MPI / OpenMP ratio for instance.
* Changing the submission process; such as using more nodes for instance.
* Changing the algorithm; yes it is a naive one but it exposes good algorithmic characteristics for you to practice what you have learned in OpenMP, MPI and OpenACC.
* Reducing the amount of work to be done such as ignoring the cells whose value will be zero during the entire simulation or avoid the reduction every non-printing iteration.
* Removing the printing phase from the iteration loop or changing the frequency at which it prints.
* Bypassing the buffer copy using a pointer swap.
* Decreasing the accuracy of the calculations by switching from doubles to floats.


If you are not sure about whether a certain optimisation is allowed or not, simply ask :)

[Go back to table of contents](#table-of-contents)
## Send your solution to the competition ##
If you want to see how far you got, and maybe if you even happened to have developed the fastest code, you can send your solution for it to be assessed and evaluated as part of the competition.

Simply send an email **by Thursday 23rd of June 2022 11:59PM** to CAPELLI Ludovic (email address in the slack channel) containing:
* The full name of each team member (no more than 3 per team), to know who we need to congratulate if your team wins :)
* The source file of the version you optimised. Typically, it will mean:
  * ```cpu.c``` (or ```cpu.F90```) if you focused on CPU using the MPI + OpenMP version.
  * ```gpu.c``` (or ```gpu.F90```) if you focused on GPU using the MPI + OpenACC version.

**Note**: there is no need to send the ```makefile``` / ```submit.sh``` scripts and so on, send just the source file of the version you optimised. Your code will be compiled and run using the original ```makefile``` / ```submit.sh``` scripts provided, on the ```big``` grid. This way, every participant has their code compiled & run in strictly identical conditions.

[Go back to table of contents](#table-of-contents)
## Whom do I talk to? ##

Me or any of the IHPCSS staff. You can find me and my email address in the slack channel called "Programming challenge". (Click on the link posted in the general slack channel about the programming challenge, it will take you to the programming challenge channel so you can join.)

[Go back to table of contents](#table-of-contents)
## Acknowledgments ##
* [John Urbanic](https://www.psc.edu/staff/urbanic)
* [David Henty](https://www.epcc.ed.ac.uk/about/staff/dr-david-henty)

[Go back to table of contents](#table-of-contents)
