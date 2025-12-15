gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c turing.c 
gcc -Wall -Wextra -pedantic-errors -std=c99 -posix -c CompTuring.c 
gcc --static -o turing turing.o CompTuring.o  -L. -l:libhgt.a -l:libquadmath.a -l:libmpfr.a -l:libgmp.a
