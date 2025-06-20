#!/bin/bash\n
TIME=time.dat
NOW=$(date +%m%d%H%M%S)
touch $TIME
echo 300000 100 $NOW > $TIME
INIT=init$NOW\n
ALPHA=0.001
SIND=-0.5
SIZE=1
DDZE=0.1
AMOD=0.01
MS=1
RMIN=0.1
RMAX=100.0
STAR=1
RDZEI=1.5
RDZEO=24.
GRID=1500
MD=0.01
gcc init_dust.ver.6.c -o $INIT -lm
DRIFT=drift$NOW$ALPHA$SIND
gcc -O3 dustdrift1Dver16_rdf.c -o $DRIFT -lm\n
DIR3=_o$SIZE$ALPHA$SIND\n
RUN=run$NOW
DIR=$RUN$DIR3
mkdir $DIR
mv $DRIFT $TIME $INIT $DIR
cd $DIR\n
./$INIT -ri $RMIN -ro $RMAX -index $SIND -rdzei $RDZEI -drdzei $DDZE -rdzeo $RDZEO -drdzeo $DDZE -alpha $ALPHA -amod $AMOD -n $GRID -md $MD -onesize $SIZE -m0 $STAR
./$DRIFT -n 5000. -twopop 0 -growth 0 &
