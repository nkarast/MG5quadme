
## MAIN LIB

for f in higgs*.cc common.cc; do
    g++ $(root-config --cflags) -fPIC -c -o ${f//.cc/.o} ${f}; 
done

## RAMBO

rootcint -f rambo_dict.cc -c rambo.h rambo_linkdef.h

g++ $(root-config --cflags) -fPIC -c -o rambo_dict.o rambo_dict.cc 
g++ $(root-config --cflags) -fPIC -c -o rambo.o rambo.cc
g++ $(root-config --libs) -fPIC -shared -o librambo.so rambo_dict.o rambo.o

## PDF libraries

rootcint -f pdf_dict.cc -c pdf.h pdf_linkdef.h
g++ $(root-config --cflags) -fPIC -c -o pdf_dict.o pdf_dict.cc 

g++ $(root-config --cflags) -fPIC -c -o pdf.o pdf.cc

gfortran -O3 -m64 -c -fPIC -ffixed-line-length-132 -ffast-math -fstrength-reduce -fexpensive-optimizations -fno-second-underscore -ext-names -c cteq6.f -o cteq6.o
gfortran -O3 -m64 -c -fPIC -ffixed-line-length-132 -ffast-math -fstrength-reduce -fexpensive-optimizations -fno-second-underscore -ext-names -c ct10pdf.f -o ct10pdf.o

g++ $(root-config --libs) -fPIC -lgfortran -shared -o libpdf.so pdf_dict.o pdf.o common.o cteq6.o ct10pdf.o

## NUMERICAL CHECKS

# g++ $(root-config --cflags) -fPIC -I$PWD -c -o tests/check_me.o tests/check_me.cc common.o
# g++ $(root-config --libs) -o tests/check_me.exe higgs*.o rambo.o tests/check_me.o common.o

# g++ $(root-config --cflags) -fPIC -I$PWD -c -o tests/check_vars.o tests/check_vars.cc common.o
# g++ $(root-config --libs) -o tests/check_vars.exe tests/check_vars.o rambo.o common.o

# g++ $(root-config --cflags) -fPIC -I$PWD -c -o tests/check_pdf.o tests/check_pdf.cc
# g++ $(root-config --libs) -lgfortran -L$PWD -lpdf -o tests/check_pdf.exe tests/check_pdf.o

## TEST ANSI C COMPAT

# for f in higgs*.c transformations.c; do
#     gcc -std=c99 -c -o cobj/${f//.cc/.o} ${f}; 
# done
