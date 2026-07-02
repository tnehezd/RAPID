#ifndef LOGGER_H
#define LOGGER_H

#include <stdio.h>

extern int g_log_level; 

// Main message logging macro (always visible)
#define MSG(fmt, ...) fprintf(stderr, fmt "\n", ##__VA_ARGS__)

// Warnings and errors (always visible, but on different channels)
#define LOG_WARN(fmt, ...)  fprintf(stderr, "WARNING [%s]: " fmt "\n",  __func__, ##__VA_ARGS__)
#define LOG_ERROR(fmt, ...) fprintf(stderr, "ERROR [%s]: " fmt "\n", __func__, ##__VA_ARGS__)

// Detailed logging (only visible at certain levels)
#define LOG_INFO(fmt, ...)  do { if (g_log_level >= 1) fprintf(stderr, "INFO [%s]: " fmt "\n",__func__, ##__VA_ARGS__); } while(0)
#define LOG_DEBUG(fmt, ...) do { if (g_log_level >= 2) fprintf(stderr, "DEBUG [%s]: " fmt "\n", __func__, ##__VA_ARGS__); } while(0)

#endif