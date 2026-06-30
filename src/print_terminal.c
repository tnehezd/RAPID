#include <stdio.h>
#include <math.h>
#include <string.h>
#include "print_terminal.h"

#define ANSI_RESET   "\x1b[0m"
#define ANSI_BLUE    "\x1b[36m"
#define ANSI_GREEN   "\x1b[32m"
#define ANSI_YELLOW  "\x1b[33m"
#define ANSI_MAGENTA "\x1b[35m"
#define ANSI_GRAY    "\x1b[90m"

#define BOX_LINES 7

typedef enum {
    STATUS_IDLE = 0,
    STATUS_BOX  = 1
} StatusState;

static StatusState status_state = STATUS_IDLE;
static int live_box_printed_once = 0;

static void clearCurrentLine(void) {
    fprintf(stderr, "\r\033[2K");
}

static void moveCursorUp(int n) {
    if (n > 0) fprintf(stderr, "\033[%dA", n);
}

static void moveCursorDown(int n) {
    if (n > 0) fprintf(stderr, "\033[%dB", n);
}

static void renderBox(int step,
                      double deltat,
                      double current_time_years,
                      double output_time,
                      const char *mode,
                      double current_mass,
                      double target_mass,
                      double initial_mass,
                      double last_snapshot_time,
                      double interval)
{
    double mass_percentage = (initial_mass > 0.0) ? (current_mass / initial_mass) * 100.0 : 0.0;
    if (mass_percentage > 100.0) mass_percentage = 100.0;
    if (mass_percentage < 0.0) mass_percentage = 0.0;

    static double last_output_time = 0.0;

    double time_in_interval = current_time_years - last_snapshot_time;
    double time_progress = (interval > 0.0) ? (time_in_interval / interval) * 100.0 : 0.0;
    
    if (time_progress > 100.0) time_progress = 100.0;
    if (time_progress < 0.0) time_progress = 0.0;

    int bar_width = 20;
    int pos = (int)(time_progress * bar_width / 100.0);
    char bar_str[64];
    int p = 0;
    for (int i = 0; i < bar_width; i++) {
        p += sprintf(bar_str + p, "%c", (i < pos) ? '#' : '.');
    }

    // KIPUCOLVA: Kiszedve az összes \033[K a sorok végéről, ami törölte a múltat
    fprintf(stderr, ANSI_GRAY "+------------------------------------------------------------+" ANSI_RESET "\n");
    fprintf(stderr, " | " ANSI_BLUE "STEP:" ANSI_RESET "   %-8d | " ANSI_MAGENTA "Mode:" ANSI_RESET " %s\n", step, mode);
    fprintf(stderr, " | " ANSI_BLUE "TIME:" ANSI_RESET "   %.4e / %.4e yr  " ANSI_GRAY "[" ANSI_GREEN "%s" ANSI_GRAY "] " ANSI_GREEN "%5.1f%%" ANSI_RESET "\n",
            current_time_years, output_time, bar_str, time_progress);
    fprintf(stderr, " | " ANSI_BLUE "DT:" ANSI_RESET "     %.4e yr\n", deltat);
    fprintf(stderr, " | " ANSI_BLUE "MASS:" ANSI_RESET "   %.4e / %.4e M_Sun " ANSI_YELLOW "(Remaining: %5.2f%%)" ANSI_RESET "\n",
            current_mass, target_mass, mass_percentage);
    fprintf(stderr, ANSI_GRAY "+------------------------------------------------------------+" ANSI_RESET "\n");
}

void printStatus(int step,
                 double deltat,
                 double current_time_years,
                 double internal_time,
                 double output_time,
                 const char *mode,
                 int was_snapshot,
                 double current_mass,
                 double target_mass,
                 double initial_mass, 
                 double last_snapshot_time, 
                 double interval)
{
    (void)internal_time;

    if (was_snapshot) {
        if (status_state == STATUS_BOX) {
            moveCursorDown(BOX_LINES);
            clearCurrentLine();
        }

        fprintf(stderr, "\n[SAVE: %6d] Time: %.2e yr | dt: %.2e yr | HDF5 SNAPSHOT SAVED\n",
                step, current_time_years, deltat);
        fflush(stderr);

        status_state = STATUS_IDLE;
        live_box_printed_once = 0;
        return;
    }

    if (status_state == STATUS_BOX && live_box_printed_once) {
        moveCursorUp(BOX_LINES);
        clearCurrentLine();
    } else {
        if (live_box_printed_once) {
            fprintf(stderr, "\n");
        }
        status_state = STATUS_BOX;
        live_box_printed_once = 1;
    }

    renderBox(step, deltat, current_time_years, output_time, mode, current_mass, target_mass, initial_mass, last_snapshot_time, interval);
    fflush(stderr);
}