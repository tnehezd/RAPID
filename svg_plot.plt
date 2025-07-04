# --- Gnuplot Szkript: Porrészecskék trajektóriáinak ábrázolása ---
# FONTOS: Ez a szkript feltételezi, hogy az adatokat előzőleg egy Python szkript
# (pl. prepare_data_for_gnuplot.py) Gnuplot-barát blokkos formátumba alakította át,
# és elmentette a 'output_0034/LOGS/dust_particle_evolution_prepared.dat' fájlba.
# Továbbá, a Python szkript létrehozza a 'max_particle_id.tmp' fájlt is.

# 1. Terminál és kimeneti fájl beállítása (PDF helyett SVG!)
set terminal svg size 1000,750
set output 'dust_particle_trajectories.svg' # A kimeneti SVG fájl

# 2. Ábra címek és tengelyfeliratok
set xlabel "Idő [év]"
set ylabel "Sugár [AU]"
set title "Porrészecskék trajektóriái (Részecske ID szerint színezve)"

# 3. Színpaletta beállítása
set palette defined ( \
    0 "#440154", \
    0.25 "#41416E", \
    0.5 "#2C728E", \
    0.75 "#22A884", \
    1.0 "#FDE725" \
)

# 4. A maximális részecske ID beolvasása a Python szkript által generált fájlból
max_particle_id = system("cat max_particle_id.tmp") + 0

# 5. Színskála (colorbox) beállítása: mutatja, melyik ID-hez melyik szín tartozik
set cblabel "Részecske Index"
set cbrange [0:max_particle_id]
set colorbox user origin 0.9,0.1 size 0.03,0.8

# 6. Y-tengely tartományának meghatározása a T=0 adatok alapján
stats '< grep "^0 " output_0002/LOGS/dust_particle_evolution.dat' using 3 nooutput
set yrange [STATS_min*0.95 : STATS_max*1.05]

# 7. Adatok ábrázolása
plot for [ID_index=0:int(max_particle_id)] \
     'output_0002/LOGS/dust_particle_evolution_prepared.dat' index ID_index \
     using 1:3 with lines \
     linecolor palette frac ID_index/real(max_particle_id) \
     notitle
