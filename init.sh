#!/bin/bash


echo "Ez a program a drift_16.ver programhoz general kezdeti parametereket!\n
Kerlek ird be, hogy hany evig fusson a szimulacio, majd uss egy [ENTER]-t!\n"

read RUNNING

#echo "A szimulacio $RUNNING evig fut\n"
echo "Ird be, hogy hany evenkent irja a program az adatokat, majd uss egy [ENTER]-t!\n"

read WRITEOUT
#echo "Az adatokat $WRITEOUT evenkent menti ki a program"

echo "Ird be, hogy mekkora legyen a kozponti csillag tomege, majd uss egy [ENTER]-t!\n"
read MS

echo "Ird be, hogy mekkora legyen a viszkozitas redukcio, majd uss egy [ENTER]-t!\n"
read AMOD

echo "Ird be, hogy mekkora legyen H atmenet AU-ban, majd uss egy [ENTER]-t!\n"
read DDZE

echo "Ird be, hogy mekkora legyen a reszecske kezdeti merete, majd uss egy [ENTER]-t!\n"
read SIZE

echo "Ird be, hogy milyen legyen SIGMA profil kitevoje, majd uss egy [ENTER]-t!\n"
read SIND

echo "Ird be, hogy mekkora legyen a kinetikus viszkozitas erteke, majd nyomj egy [ENTER]-t!\n"
read ALPHA



INITFILE="indit_am$AMOD.dze$DDZE.sind$SIND.size$SIZE.sh"
touch $INITFILE

echo "#!/bin/bash\n" >> $INITFILE
echo "TIME="time.dat"" >> $INITFILE
echo "NOW=\$(date +"%m%d%H%M%S")" >> $INITFILE
echo "touch \$TIME" >> $INITFILE
echo "echo "$RUNNING $WRITEOUT \$NOW" > \$TIME" >> $INITFILE

echo "INIT=init\$NOW\n" >> $INITFILE

echo "ALPHA=$ALPHA">> $INITFILE
echo "SIND=$SIND" >> $INITFILE
echo "SIZE=$SIZE" >> $INITFILE
echo "DDZE=$DDZE" >> $INITFILE
echo "AMOD=$AMOD" >> $INITFILE
echo "MS=$MS" >> $INITFILE
echo "RMIN=0.1" >> $INITFILE
echo "RMAX=100.0" >> $INITFILE
echo "STAR=$MS" >> $INITFILE
echo "RDZEI=1.5" >> $INITFILE
echo "RDZEO=24." >> $INITFILE
echo "GRID=1500" >> $INITFILE
echo "MD=0.01" >> $INITFILE

echo "gcc init_dust.ver.6.c -o \$INIT -lm" >> $INITFILE
echo "DRIFT=drift\$NOW\$ALPHA\$SIND" >> $INITFILE
echo "gcc -O3 dustdrift1Dver16_rdf.c -o \$DRIFT -lm\n" >> $INITFILE

echo "DIR3=_o\$SIZE\$ALPHA\$SIND\n" >> $INITFILE

echo "RUN=run\$NOW" >> $INITFILE
echo "DIR=\$RUN\$DIR3" >> $INITFILE
echo "mkdir \$DIR" >> $INITFILE
echo "mv \$DRIFT \$TIME \$INIT \$DIR" >> $INITFILE
echo "cd \$DIR\n" >> $INITFILE

echo "./\$INIT -ri \$RMIN -ro \$RMAX -index \$SIND -rdzei \$RDZEI -drdzei \$DDZE -rdzeo \$RDZEO -drdzeo \$DDZE -alpha \$ALPHA -amod \$AMOD -n \$GRID -md \$MD -onesize \$SIZE -m0 \$STAR" >> $INITFILE
echo "./\$DRIFT -n 5000. -twopop 0 -growth 0 &" >> $INITFILE


CURRENT=$(pwd)
echo "$CURRENT"

SINDDIR=r$SIND

if [ -d "$SINDDIR" ]; then 
	echo "van ilyen mappa: $SINDDIR"
else
	echo "nincs ilyen mappa, de letre hozok"
	mkdir $SINDDIR
fi

echo "Belepes a $SINDDIR mappaba!"
cd $SINDDIR


SIZEDIR=size$SIZE

if [ -d "$SIZEDIR" ]; then 
	echo "van ilyen mappa: $SIZEDIR"
else
	echo "nincs ilyen mappa, de letre hozok"
	mkdir $SIZEDIR
fi

echo "Belepes a $SIZEDIR mappaba!"
cd $SIZEDIR

SUNDIR=msun$MS

if [ -d "$SUNDIR" ]; then 
	echo "van ilyen mappa: $SUNDIR"
else
	echo "nincs ilyen mappa, de letre hozok"
	mkdir $SUNDIR
fi

cd $SUNDIR

CREATDIR=$(pwd) 


cd $CURRENT
echo $(pwd)
mv $INITFILE $CREATDIR
cp *.c $CREATDIR

cd $CREATDIR

nohup sh $INITFILE &