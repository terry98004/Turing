gcc -Wall -Wextra -pedantic-errors -std=gnu17 -pthread -c turing.c 
gcc -Wall -Wextra -pedantic-errors -std=gnu17 -pthread -c CompTuring.c 
gcc --static -pthread -o turing turing.o CompTuring.o  -L. -l:libhgt.a -l:libmpfr.a -l:libgmp.a
