#ifndef PHOTOEVAPORATION_WRAPPER_H
#define PHOTOEVAPORATION_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

// void* -ként adjuk át, így a C és a C++ is csont nélkül elfogadja anélkül, hogy ismernék a típust
void computePhotoevaporationSink(void *disk_opaque);

#ifdef __cplusplus
}
#endif

#endif // PHOTOEVAPORATION_WRAPPER_H