# --- Gnuplot Szkript: Porrészecskék sugárevolúciójának ábrázolása ---

# 1. Terminál és kimeneti fájl beállítása
set terminal pdfcairo size 10in,7.5in font "12"
set output 'dust_radius_evolution_fixed_y_range.pdf' # Új fájlnév a fix tengelytartomány jelzésére

# 2. Ábra címek és tengelyfeliratok
set xlabel "Idő [év]"
set ylabel "Sugár [AU]"
set title "Porrészecskék sugarának evolúciója (Fixált Y-tengely: T=0 tartomány)"

# 3. Színpaletta beállítása
set palette defined ( \
    0 "#440154", \
    0.25 "#41416E", \
    0.5 "#2C728E", \
    0.75 "#22A884", \
    1.0 "#FDE725" \
)

# Színskála (colorbox) beállítása: mutatja, melyik index melyik szín
set cblabel "Részecske Index"
# A cbrange automatikus lesz, ahogy korábban megbeszéltük, ami jó.
# set cbrange [0:23] # Ezt továbbra is elhagyjuk, vagy fixen beállíthatod pl. [0:1000]-re.
set colorbox user origin 0.9,0.1 size 0.03,0.8 # Elhelyezés és méret

# 4. Y-tengely tartományának meghatározása a T=0 adatok alapján
# Szűrjük az adatokat T=0-ra (1. oszlop értéke 0), és a 3. oszlop (Sugár) min/max értékét keressük.
# Az "every ::0::N_particles-1" rész biztosítja, hogy csak az első blokkot (T=0) olvassa be.
# Feltételezi, hogy a T=0 adatok az első blokkban vannak, és a fájl sorrendje idő szerint növekvő.
# Ehhez szükségünk van a T=0 időponthoz tartozó sorok számára.
# Alternatívaként használhatjuk a 'grep' parancsot is az időszűrőhöz.

# Használjuk a 'grep' módszert, mert az robusztusabb, ha az adatfájl nincs szigorúan blokkokra bontva dupla üres sorokkal:
stats '< grep "^0 " output_0002/LOGS/dust_particle_evolution.dat' using 3 nooutput

# A stats parancs eredményeit használjuk az yrange beállításához:
# STATS_min és STATS_max tartalmazza a T=0 időponthoz tartozó sugár min/max értékét.
# Kicsit kibővítjük a tartományt, hogy legyen "levegő" a pontok körül.
set yrange [STATS_min*0.95 : STATS_max*1.05] # 5% puffert adunk a min/max értékhez

# 5. Adatok ábrázolása
plot 'output_0002/LOGS/dust_particle_evolution.dat' using 1:3:2 with points pt 7 ps 0.8 lc palette notitle

# exit