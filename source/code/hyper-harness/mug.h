#pragma once

// start inlining gkyl_alloc.h 

// start inlining gkyl_util.h 

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// random number generator
// start inlining pcg_basic.h 
/*
 * PCG Random Number Generation for C.
 *
 * Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * For additional information about the PCG random number generation scheme,
 * including its license and other licensing options, visit
 *
 *     http://www.pcg-random.org
 */

/*
 * This code is derived from the full C implementation, which is in turn
 * derived from the canonical C++ PCG implementation. The C++ version
 * has many additional features and is preferable if you can use C++ in
 * your project.
 */

#ifndef PCG_BASIC_H_INCLUDED
#define PCG_BASIC_H_INCLUDED 1

#include <inttypes.h>

#if __cplusplus
extern "C" {
#endif

struct pcg_state_setseq_64 {    // Internals are *Private*.
    uint64_t state;             // RNG state.  All values are possible.
    uint64_t inc;               // Controls which RNG sequence (stream) is
                                // selected. Must *always* be odd.
};
typedef struct pcg_state_setseq_64 pcg32_random_t;

// If you *must* statically initialize it, here's one.

#define PCG32_INITIALIZER   { 0x853c49e6748fea9bULL, 0xda3e39cb94b95bdbULL }

// pcg32_srandom(initstate, initseq)
// pcg32_srandom_r(rng, initstate, initseq):
//     Seed the rng.  Specified in two parts, state initializer and a
//     sequence selection constant (a.k.a. stream id)

void pcg32_srandom(uint64_t initstate, uint64_t initseq);
void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate,
                     uint64_t initseq);

// pcg32_random()
// pcg32_random_r(rng)
//     Generate a uniformly distributed 32-bit random number

uint32_t pcg32_random(void);
uint32_t pcg32_random_r(pcg32_random_t* rng);

// pcg32_boundedrand(bound):
// pcg32_boundedrand_r(rng, bound):
//     Generate a uniformly distributed number, r, where 0 <= r < bound

uint32_t pcg32_boundedrand(uint32_t bound);
uint32_t pcg32_boundedrand_r(pcg32_random_t* rng, uint32_t bound);

#if __cplusplus
}
#endif

#endif // PCG_BASIC_H_INCLUDED
// ended inlining pcg_basic.h 

#ifdef __cplusplus

// extern "C" guards needed when using code from C++
# define EXTERN_C_BEG extern "C" {
# define EXTERN_C_END }

#else

# define EXTERN_C_BEG 
# define EXTERN_C_END

#endif

// restrict keyword in C and C++ are different
#ifdef __cplusplus
# define GKYL_RESTRICT __restrict__
#else
# define GKYL_RESTRICT restrict
#endif

// Maximum configuration-space dimensions supported
#ifndef GKYL_MAX_CDIM
# define GKYL_MAX_CDIM 3
#endif

// Maximum velocity-space dimensions supported
#ifndef GKYL_MAX_VDIM
# define GKYL_MAX_VDIM 3
#endif

// Maximum dimensions supported
#ifndef GKYL_MAX_DIM
# define GKYL_MAX_DIM 7
#endif

// Maximum number of supported species
#ifndef GKYL_MAX_SPECIES
# define GKYL_MAX_SPECIES 16
#endif

// Maximum number of supported species
#ifndef GKYL_MAX_REACT
# define GKYL_MAX_REACT 3*GKYL_MAX_SPECIES
#endif

// Maximum number of supported sources
#ifndef GKYL_MAX_SOURCES
# define GKYL_MAX_SOURCES 4
#endif

// Maximum number of supported fdot multiplier types
#ifndef GKYL_MAX_FDOT_MUL
# define GKYL_MAX_FDOT_MUL 4
#endif

// Maximum number of supported charge states
#ifndef GKYL_MAX_CHARGE_STATE
# define GKYL_MAX_CHARGE_STATE 18
#endif

// Maximum number of supported densities for radiation
#ifndef GKYL_MAX_RAD_DENSITIES
# define GKYL_MAX_RAD_DENSITIES 26
#endif

// Maximum number of supported projection objects
#ifndef GKYL_MAX_PROJ
# define GKYL_MAX_PROJ 4
#endif

// Maximum number of ghost cells in each direction
#ifndef GKYL_MAX_NGHOST
# define GKYL_MAX_NGHOST 8
#endif

// Default alignment boundary
#ifndef GKYL_DEF_ALIGN
# define GKYL_DEF_ALIGN 64
#endif

// Maximum number of blocks
#ifndef GKYL_MAX_BLOCKS
# define GKYL_MAX_BLOCKS 12
#endif

// CUDA specific defines etc
#ifdef __NVCC__

#include <cuda_runtime.h>

#define GKYL_HAVE_CUDA

#define GKYL_CU_DH __device__ __host__
#define GKYL_CU_D __device__ 

// for directional copies
enum gkyl_cu_memcpy_kind {
  GKYL_CU_MEMCPY_H2H = cudaMemcpyHostToHost,
  GKYL_CU_MEMCPY_H2D = cudaMemcpyHostToDevice,
  GKYL_CU_MEMCPY_D2H = cudaMemcpyDeviceToHost,
  GKYL_CU_MEMCPY_D2D = cudaMemcpyDeviceToDevice
};

#define GKYL_DEFAULT_NUM_THREADS 256

// CUDA helper function to find CUDA errors
#define checkCuda(val)           __checkCudaErrors__ ( (val), #val, __FILE__, __LINE__ )
inline cudaError_t __checkCudaErrors__(cudaError_t code, const char *func, const char *file, int line)
{
  if (code) {
    fprintf(stderr, "CUDA error: %s (code=%u)  \"%s\" at %s:%d \n",
      cudaGetErrorString(code), (unsigned int)code, func, file, line);
    cudaDeviceReset();
    exit(EXIT_FAILURE);
  }
  return code;
}

#else

#undef GKYL_HAVE_CUDA
#define GKYL_CU_DH
#define GKYL_CU_D
#define checkCuda(val) 
// for directional copies
enum gkyl_cu_memcpy_kind {
  GKYL_CU_MEMCPY_H2H,
  GKYL_CU_MEMCPY_H2D,
  GKYL_CU_MEMCPY_D2H,
  GKYL_CU_MEMCPY_D2D,
};

#define GKYL_DEFAULT_NUM_THREADS 1

#endif // CUDA specific defines etc

// This funny looking macro allows getting a pointer to the 'type'
// struct that contains an object 'member' given the 'ptr' to the
// 'member' inside 'type'. (Did I just write this gobbledygook?!)
//
// See https://en.wikipedia.org/wiki/Offsetof
#define container_of(ptr, type, member)                                 \
    ((type *)((char *)(1 ? (ptr) : &((type *)0)->member) - offsetof(type, member)))

// Select type-specific compare function
#define gkyl_compare(a, b, eps)                 \
    _Generic((a),                               \
      float: gkyl_compare_float,                \
      double: gkyl_compare_double)              \
    (a, b, eps)

// a quick-and-dirty macro for testing (mostly) CUDA kernel code
#define GKYL_CU_CHECK(expr, cntr) do {                                  \
      if (!(expr)) {                                                    \
        *cntr += 1;                                                     \
        printf("%s failed! (%s:%d)\n", #expr, __FILE__, __LINE__);      \
      }                                                                 \
    } while (0)

// Computes length of string needed given a format specifier and data. Example:
//
// size_t len = gkyl_calc_strlen("%s-%d", "gkyl", 25);
// 
#define gkyl_calc_strlen(fmt, ...) snprintf(0, 0, fmt, __VA_ARGS__)

// Open file 'fname' with 'mode; into handle 'fp'. Handle is closed
// when block attached to with_file exits
#define with_file(fp, fname, mode)                                             \
  for (bool _break = (fp = fopen(fname, mode), (fp != NULL)); _break;          \
       _break = false, fclose(fp))

// Code

#define GKYL_MIN2(x,y) ((x)<(y) ? (x) : (y))
#define GKYL_MAX2(x,y) ((x)>(y) ? (x) : (y))
#define GKYL_SGN(b) (((b)>=0.) ? 1.0 : -1.0)

// GKYL_ALIGN_UP finds the integer that is closest to "a" that is
// multiple of "b"
#define GKYL_UDIV_UP(a, b) (((a) + (b) - 1) / (b))
#define GKYL_ALIGN_UP(a, b) (GKYL_UDIV_UP(a, b) * (b))

EXTERN_C_BEG

/**
 * Kernel floating-point op-counts
 */
struct gkyl_kern_op_count {
  size_t num_sum; // number of + and - operations
  size_t num_prod; // number of * and / operations
};

// String, integer pairs
struct gkyl_str_int_pair {
  const char *str;
  int val;
};

/**
 * Search @a pairs list for @a str and return the corresponding int
 * value. Return @a def if not found. The pair table must be NULL
 * terminated.
 *
 * @param pairs Pair list to search from. Final entry must be { 0, 0 }
 * @param str String to search for
 * @param def Default value to return
 * @return value corresponding to @a str, or @a def.
 */
int gkyl_search_str_int_pair_by_str(const struct gkyl_str_int_pair pairs[], const char *str, int def);

/**
 * Search @a pairs list for @a val and return the corresponding string
 * value. Return @a def if not found. The pair table must be NULL
 * terminated.
 *
 * @param pairs Pair list to search from. Final entry must be { 0, 0 }
 * @param val Value to search for
 * @param def Default value to return
 * @return value corresponding to @a val, or @a def.
 */
const char* gkyl_search_str_int_pair_by_int(const struct gkyl_str_int_pair pairs[], int val, const char *def);

/**
 * Time-trigger. Typical initialization is:
 * 
 * struct gkyl_tm_trigger tmt = { .dt = tend/nframe };
 */
struct gkyl_tm_trigger {
  int curr; // Current counter.
  double dt, tcurr; // Time-interval, current time.
};

/**
 * Check if the tcurr should trigger and bump internal counters if it
 * does. This only works if sequential calls to this method have the
 * tcurr monotonically increasing.
 *
 * @param tmt Time trigger object
 * @param tcurr Current time.
 * @return 1 if triggered, 0 otherwise
 */
int gkyl_tm_trigger_check_and_bump(struct gkyl_tm_trigger *tmt, double tcurr);

/**
 * Print error message to stderr and exit.
 *
 * @param msg Error message.
 */
void gkyl_exit(const char* msg);

/**
 * Compares two float numbers 'a' and 'b' to check if they are
 * sufficiently close by, where 'eps' is the relative tolerance.
 */
int gkyl_compare_float(float a, float b, float eps);

/**
 * Compares two double numbers 'a' and 'b' to check if they are
 * sufficiently close by, where 'eps' is the relative tolerance.
 */
int gkyl_compare_double(double a, double b, double eps);

/**
 * Copy (small) int arrays.
 *
 * @param n Number of elements to copy
 * @param inp Input array
 * @param out Output array
 */
GKYL_CU_DH
static inline void
gkyl_copy_int_arr(int n, const int* GKYL_RESTRICT inp, int* GKYL_RESTRICT out)
{
  for (int i=0; i<n; ++i) out[i] = inp[i];
}

/**
 * Copy (small) long arrays.
 *
 * @param n Number of elements to copy
 * @param inp Input array
 * @param out Output array
 */
GKYL_CU_DH
static inline void
gkyl_copy_long_arr(int n, const long* GKYL_RESTRICT inp, long* GKYL_RESTRICT out)
{
  for (int i=0; i<n; ++i) out[i] = inp[i];
}

/**
 * Copy (small) double arrays.
 *
 * @param n Number of elements to copy
 * @param inp Input array
 * @param out Output array
 */
GKYL_CU_DH
static inline void
gkyl_copy_double_arr(int n, const double* GKYL_RESTRICT inp, double* GKYL_RESTRICT out)
{
  for (int i=0; i<n; ++i) out[i] = inp[i];
}

/**
 *   Round a/b to nearest higher integer value
 */
GKYL_CU_DH
static inline int
gkyl_int_div_up(int a, int b)
{
  return (a%b != 0) ? (a/b+1) : (a/b);
}

/**
 *   Minmod limiter for choosing the minimal modification between 3 values (usually slopes)
 */
GKYL_CU_DH
static inline double 
gkyl_minmod(double a, double b, double c)
{
  double sa = GKYL_SGN(a);
  double sb = GKYL_SGN(b);
  double sc = GKYL_SGN(c);
  if( (sa==sb) && (sb==sc) ) {
    if (sa<0)
      return GKYL_MAX2(GKYL_MAX2(a,b),c);
    else
      return GKYL_MIN2(GKYL_MIN2(a,b),c);
  }
  else {
     return 0;
  }
}

/**
 * Gets wall-clock time in secs/nanoseconds.
 * 
 * @return Time object.
 */
struct timespec gkyl_wall_clock(void);

/**
 * Difference between two timespec objects.
 *
 * @param tstart Start time
 * @param tend End time 
 * @return Time object representing difference
 */
struct timespec gkyl_time_diff(struct timespec tstart, struct timespec tend);

/**
 * Difference between timespec object and "now", returned in seconds.
 * 
 * @param tm Timespec
 * @return Time in seconds
 */
double gkyl_time_diff_now_sec(struct timespec tm);

/**
 * Compute in secs time stored in timespec object.
 *
 * @param tm Timespec object
 * @return Time in seconds
 */
double gkyl_time_sec(struct timespec tm);

/**
 * Compute time in seconds since epoch.
 *
 * @return Time in seconds
 */
double gkyl_time_now(void);

/**
 * Initialize 32-bit RNG 
 *
 * @param nd_seed true if to use non-deterministic seed, or false for
 *   a determistic seed.
 */
pcg32_random_t gkyl_pcg32_init(bool nd_seed);

/**
 * Returns a unsigned 32-bit random number
 *
 * @param rng Pointer to RNG
 * @return Uniformly distributed 32-bit integer
 */
uint32_t gkyl_pcg32_rand_uint32(pcg32_random_t* rng);

/**
 * Returns a random number in [0,1), rounded to the nearest
 * 1/2^32. This may not be enough resolution for some applications.
 *
 * @param rng Pointer to RNG
 * @return Uniformly distributed double in [0,1)
 */
double gkyl_pcg32_rand_double(pcg32_random_t* rng);

/** 64-bit RNG: concatination of two 32-bit RNGs */
typedef struct { pcg32_random_t gen[2]; } pcg64_random_t;

/**
 * Initialize 64-bit RNG 
 *
 * @param nd_seed true if to use non-deterministic seed, or false for
 *   a determistic seed.
 */
pcg64_random_t gkyl_pcg64_init(bool nd_seed);

/**
 * Returns a unsigned 64-bit random number
 *
 * @param rng Pointer to RNG
 * @return Uniformly distributed 64-bit integer
 */
uint64_t gkyl_pcg64_rand_uint64(pcg64_random_t* rng);

/**
 * Returns a random number in [0,1), rounded to the nearest
 * 1/2^64.
 *
 * @param rng Pointer to RNG
 * @return Uniformly distributed double in [0,1)
 */
double gkyl_pcg64_rand_double(pcg64_random_t *rng);

/**
 * Check if file exists.
 *
 * @param fname Name of file
 * @return true if file exists, false otherwise
 */
bool gkyl_check_file_exists(const char *fname);

/**
 * Get file size in bytes.
 *
 * @param fname Name of file
 * @return file size in bytes
 */
int64_t gkyl_file_size(const char *fname);

/**
 * Read contents of file into a character buffer. The returned buffer
 * must be freed using gkyl_free.
 *
 * @param fname Name of file
 * @param sz On return this has the size of data read
 * @return Data in file as a char array.
 */
char* gkyl_load_file(const char *fname, int64_t *sz);

/**
 * Structure to store msgpack data
 */
struct gkyl_msgpack_data {
  size_t meta_sz; // size of data in bytes
  char *meta; // data pointer (pointer to bytes. NOT a NULL terminated string!)
};

// Msgpack element type.
enum gkyl_msgpack_elem_type { 
  GKYL_MP_BOOL, GKYL_MP_INT, GKYL_MP_UNSIGNED_INT,
  GKYL_MP_FLOAT, GKYL_MP_DOUBLE, GKYL_MP_STRING,
};

// Entry type into msgpack map.
struct gkyl_msgpack_map_elem {
  char *key; // name of element
  enum gkyl_msgpack_elem_type elem_type; // type of element
  union {
    // depending on elem_type one of following should be set    
    bool bval;
    unsigned int uval;
    int ival;
    float fval;
    double dval;
    char *cval; // null terminated string managed by caller
  };
};

// The following functions and the macro reduce the errors in
// constructing gkyl_msgpack_map_elem objects
static inline struct gkyl_msgpack_map_elem
gmpe_bval(char *key, bool val)
{
  return (struct gkyl_msgpack_map_elem) {
    .key = key,
    .elem_type = GKYL_MP_BOOL,
    .bval = val
  };
}
static inline struct gkyl_msgpack_map_elem
gmpe_uval(char *key, unsigned int val)
{
  return (struct gkyl_msgpack_map_elem) {
    .key = key,
    .elem_type = GKYL_MP_UNSIGNED_INT,
    .uval = val
  };
}
static inline struct gkyl_msgpack_map_elem
gmpe_ival(char *key, int val)
{
  return (struct gkyl_msgpack_map_elem) {
    .key = key,
    .elem_type = GKYL_MP_INT,
    .ival = val
  };
}
static inline struct gkyl_msgpack_map_elem
gmpe_fval(char *key, float val)
{
  return (struct gkyl_msgpack_map_elem) {
    .key = key,
    .elem_type = GKYL_MP_FLOAT,
    .fval = val
  };
}
static inline struct gkyl_msgpack_map_elem
gmpe_dval(char *key, double val)
{
  return (struct gkyl_msgpack_map_elem) {
    .key = key,
    .elem_type = GKYL_MP_DOUBLE,
    .dval = val
  };
}
static inline struct gkyl_msgpack_map_elem
gmpe_cval(char *key, char *val)
{
  return (struct gkyl_msgpack_map_elem) {
    .key = key,
    .elem_type = GKYL_MP_STRING,
    .cval = val
  };
}
#define GKYL_MSGPACK_MAP_ELEM(key, val) _Generic((val), \
      bool: gmpe_bval,                                  \
      int: gmpe_ival,                                   \
      unsigned int: gmpe_uval,                          \
      float: gmpe_fval,                                 \
      double: gmpe_dval,                                \
      char *: gmpe_cval)(key, val)

#undef GMPE_TAG

/**
 * Check if element list contains the element with name "key".
 *
 * @param nvals number of elements in elist_in.
 * @param elist_in List of elements.
 * @param key Name of the element to look for.
 * @return True is list has this element key, false otherwise.
 */
bool
gkyl_msgpack_map_elem_has_key(int nvals, const struct gkyl_msgpack_map_elem *elist,
  const char *key);

/**
 * Allocate a new list of MessagePack map elements by cloning
 * another one.
 *
 * @param nvals number of elements in elist_in.
 * @param elist_in List of elements to put in MessagePack.
 * @return New msgpack_map_elem object. Free with gkyl_msgpack_map_elem_release.
 */
struct gkyl_msgpack_map_elem* gkyl_msgpack_map_elem_clone(int nvals,
  const struct gkyl_msgpack_map_elem *elist_in);

/**
 * Allocate a new list of MessagePack map elements out of the union of one or
 * more map elem lists.
 *
 * @param numlist_union Number of lists.
 * @param nvals_union Number of values in the elist array, one for each list.
 * @param elist_union List of elements to insert into map, for each list.
 * @param elist_out_len Length of the map elem list produced.
 * @return New msgpack_map_elem object. Free with gkyl_msgpack_map_elem_release.
 */
struct gkyl_msgpack_map_elem* gkyl_msgpack_map_elem_union(int numlist_union, int *nvals_union,
  const struct gkyl_msgpack_map_elem **elist_union, int *elist_out_len);

/**
 * Update the type double value of an element in an element list.
 *
 * @param nvals number of elements in elist_in.
 * @param elist_in List of elements.
 * @param key Name of the element to update.
 * @param value Value to update element with.
 */
void gkyl_msgpack_map_elem_set_double(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key, double value);

/**
 * Update the type unsigned int value of an element in an element list.
 *
 * @param nvals number of elements in elist_in.
 * @param elist_in List of elements.
 * @param key Name of the element to update.
 * @param value Value to update element with.
 */
void gkyl_msgpack_map_elem_set_uint(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key, unsigned int value);

/**
 * Fetch the type double value of an element in an element list.
 *
 * @param nvals number of elements in elist_in.
 * @param elist_in List of elements.
 * @param key Name of the element to update.
 * @return Value of the specified element.
 */
double gkyl_msgpack_map_elem_get_double(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key);

/**
 * Fetch the type unsigned_int value of an element in an element list.
 *
 * @param nvals number of elements in elist_in.
 * @param elist_in List of elements.
 * @param key Name of the element to update.
 * @return Value of the specified element.
 */
unsigned int gkyl_msgpack_map_elem_get_uint(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key);

/**
 * Fetch the pointer to the string value of an element in an element list.
 *
 * @param nvals number of elements in elist_in.
 * @param elist_in List of elements.
 * @param key Name of the element to update.
 * @return Pointer to string value of the specified element.
 */
char* gkyl_msgpack_map_elem_get_string(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key);

/**
 * Free the memory allocated to store a string in an element of the given list.
 *
 * @param nvals number of elements in elist_in.
 * @param elist_in List of elements.
 * @param key Name of the element whose string value to release.
 */
void gkyl_msgpack_map_elem_release_string(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key);

/**
 * Release list of MessagePack map elements.
 *
 * @param nvals number of elements in elist_in.
 * @param elist_in List of elements to put in MessagePack.
 */
void gkyl_msgpack_map_elem_release(int nvals, struct gkyl_msgpack_map_elem *elist_in);

/**
 * Create a new msgpack data from list of map elements.
 *
 * @param nvals Number of values in the elist array
 * @param elist List of elements to insert into map
 * @return New msgpack_data object. Free using the release method
 */
struct gkyl_msgpack_data* gkyl_msgpack_create(int nvals, const struct gkyl_msgpack_map_elem *elist);

/**
 * Create a new msgpack data from union of a numlist lists of map elements.
 *
 * @param numlist_union Number of lists.
 * @param nvals_union Number of values in the elist array, one for each list.
 * @param elist_union List of elements to insert into map, for each list.
 * @return New msgpack_data object. Free using the release method
 */
struct gkyl_msgpack_data* gkyl_msgpack_create_union(int numlist_union,
  int *nvals_union, const struct gkyl_msgpack_map_elem **elist_union);

/**
 * Clone a msgpack.
 *
 * @param mdata_in MessagePack to clone.
 * @return New msgpack_data object. Free using the release method
 */
struct gkyl_msgpack_data* gkyl_msgpack_clone(struct gkyl_msgpack_data *mdata_in);

/**
 * Read a MessagePack for the elements in an element list.
 * NOTE: if one of the elements read is a string, its memory must be
 * freed when done using it (e.g. with gkyl_msgpack_map_elem_release_string);
 * 
 * @param mpack_in MessagePack to read.
 * @param nvals Number of values in element list elist.
 * @param elist Element list to populate.
 */
void gkyl_msgpack_to_map_elem_list(struct gkyl_msgpack_data* mpack_in, int nvals,
  struct gkyl_msgpack_map_elem *elist);

/**
 * Release data created by the gkyl_msgpack_create method.
 *
 * @param mdata Data object to free
 */
void gkyl_msgpack_data_release(struct gkyl_msgpack_data *mdata);

EXTERN_C_END
// ended inlining gkyl_util.h 

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

/**
 * Set the global flag to turn on memory allocation/deallocation
 * tracing.
 *
 * @param flag Flag to set
 */
void gkyl_mem_debug_set(bool flag);

/**
 * Set the global flag to turn on cuda memory allocation/deallocation
 * tracing.
 *
 * @param flag Flag to set
 */
void gkyl_cu_dev_mem_debug_set(bool flag);

// The following allocators are implemented as macros to allow
// debugging memory leaks. In general, it is best to use
// valgrind. However, valgrind can be very slow and also it gets
// confused with CUDA symbols in the executable (even for host
// code). Further, it seems that cuda-memcheck is not really tracking
// memory leaks.

#define gkyl_malloc(size)                                       \
    gkyl_malloc_(__FILE__, __LINE__, __FUNCTION__, size)
#define gkyl_calloc(num, size)                                  \
    gkyl_calloc_(__FILE__, __LINE__, __FUNCTION__, num, size)
#define gkyl_realloc(ptr, new_size)                                     \
    gkyl_realloc_(__FILE__, __LINE__, __FUNCTION__, ptr, new_size)
#define gkyl_free(ptr) \
    gkyl_free_(__FILE__, __LINE__, __FUNCTION__, ptr)

// Allocate memory that is aligned to given byte boundary. You must
// only use gkyl_aligned_realloc() and gkyl_aligned_free() methods to
// reallocate or free memory returned by this methods.

#define gkyl_aligned_alloc(align, size)                                 \
    gkyl_aligned_alloc_(__FILE__, __LINE__, __FUNCTION__, align, size)
#define gkyl_aligned_realloc(ptr, align, old_sz, new_sz)                       \
    gkyl_aligned_realloc_(__FILE__, __LINE__, __FUNCTION__, ptr, align, old_sz, new_sz)
#define gkyl_aligned_free(ptr) \
    gkyl_aligned_free_(__FILE__, __LINE__, __FUNCTION__, ptr)

// The following allocators have the same calling/return behavior as
// standard C allocators. However, an error is signaled if allocation
// fails.

void* gkyl_malloc_(const char *file, int line, const char *func, size_t size);
void* gkyl_calloc_(const char *file, int line, const char *func, size_t num, size_t size);
void* gkyl_realloc_(const char *file, int line, const char *func, void *ptr, size_t new_size);
void gkyl_free_(const char *file, int line, const char *func, void *ptr);

/**
 * Allocate memory that is aligned to given byte boundary. You must
 * only use gkyl_aligned_realloc() and gkyl_aligned_free() methods to
 * reallocate or free memory returned by this method.
 *
 * @param align Alignment boundary. Must be power of 2.
 * @param size Number of bytes to allocate.
 */
void* gkyl_aligned_alloc_(const char *file, int line, const char *func, size_t align, size_t size);

/**
 * Reallocate memory that is aligned to given byte boundary. The
 * pointer passed to this must be allocated by
 * gkyl_aligned_alloc(). Returned pointer must be freed using
 * gkyl_aligned_free() method.
 *
 * @param ptr Pointer to reallocate
 * @param align Alignment boundary. Must be power of 2.
 * @param old_sz Old size of memory.
 * @param new_sz New size of memory.
 */
void* gkyl_aligned_realloc_(const char *file, int line, const char *func,
  void *ptr, size_t align, size_t old_sz, size_t new_sz);

/**
 * Free memory allocated by gkyl_aligned_alloc().
 *
 * @param ptr Memory to free.
 */
void gkyl_aligned_free_(const char *file, int line, const char *func, void *ptr);

// Represents a sized chunk of memory
typedef struct gkyl_mem_buff_tag* gkyl_mem_buff;

/** Allocate new memory buffer with count bytes */
gkyl_mem_buff gkyl_mem_buff_new(size_t count);

/** Allocate new memory buffer on GPU with count bytes */
gkyl_mem_buff gkyl_mem_buff_cu_new(size_t count);

/** Resize the mem_buff */
gkyl_mem_buff gkyl_mem_buff_resize(gkyl_mem_buff mem, size_t count);

/** Get size of mem_buff */
size_t gkyl_mem_buff_size(gkyl_mem_buff mem);

/** Get pointer to data buffer */
char* gkyl_mem_buff_data(gkyl_mem_buff mem);

/** Free buffer */
void gkyl_mem_buff_release(gkyl_mem_buff mem);

// CUDA specific code (NV: Nvidia)

#define gkyl_cu_malloc(size)                                    \
    gkyl_cu_malloc_(__FILE__, __LINE__, __FUNCTION__, size)
#define gkyl_cu_free(ptr)                                       \
    gkyl_cu_free_(__FILE__, __LINE__, __FUNCTION__, ptr)

#define gkyl_cu_malloc_host(size)                               \
    gkyl_cu_malloc_host_(__FILE__, __LINE__, __FUNCTION__, size)
#define gkyl_cu_free_host(ptr)                                  \
    gkyl_cu_free_host_(__FILE__, __LINE__, __FUNCTION__, ptr)

/** Allocate memory on NV-GPU */
void* gkyl_cu_malloc_(const char *file, int line, const char *func, size_t size);

/** Free memory on device */
void gkyl_cu_free_(const char *file, int line, const char *func, void *ptr);

/** Allocate pinned host memory on NV-GPU */
void* gkyl_cu_malloc_host_(const char *file, int line, const char *func, size_t size);

/** Free pinned host memory on device */
void gkyl_cu_free_host_(const char *file, int line, const char *func, void *ptr);

/** Copy data between host/device */
void gkyl_cu_memcpy(void *dst, const void *src, size_t count, enum gkyl_cu_memcpy_kind kind);

/** Copy data between host/device */
#ifdef GKYL_HAVE_CUDA
void gkyl_cu_memcpy_async(void *dst, const void *src, size_t count, enum gkyl_cu_memcpy_kind kind, cudaStream_t stream);
#else
void gkyl_cu_memcpy_async(void *dst, const void *src, size_t count, enum gkyl_cu_memcpy_kind kind, int stream);
#endif

/** Set memory on device */
void gkyl_cu_memset(void *data, int val, size_t count);
// ended inlining gkyl_alloc.h 
// start inlining gkyl_app.h 

#include <stdbool.h>
// start inlining gkyl_array_rio.h 

#include <stdio.h>

// start inlining gkyl_array.h 

// start inlining gkyl_ref_count.h 

// Implementation essentially as presented by Chris Wellons:
// https://nullprogram.com/blog/2015/02/17/

// skipping file: gkyl_util.h 

/**
 * Object holding use count and pointer to destructor function.
 */
struct gkyl_ref_count {
  void (*free)(const struct gkyl_ref_count* );
  int count;
};

/**
 * Create a new ref object with specified function pointer that
 * deletes the container object.
 *
 * @param free Function pointer to the delete function
 * @return Ref object
 */
static inline struct gkyl_ref_count
gkyl_ref_count_init(void (*free)(const struct gkyl_ref_count* ))
{
  return (struct gkyl_ref_count) {
    .free = free,
    .count = 1,
  };
}

/**
 * Increment use count.
 *
 * @param ref Object to increment.
 */
static inline void
gkyl_ref_count_inc(const struct gkyl_ref_count *ref)
{
  ((struct gkyl_ref_count *)ref)->count++;
}

/**
 * Decrement use count. If use count is zero the free (destructor)
 * function is called.
 *
 * @param ref Object to decrement.
 */
static inline void
gkyl_ref_count_dec(const struct gkyl_ref_count *ref)
{
  if (--((struct gkyl_ref_count *)ref)->count == 0)
    ref->free(ref);
}
// ended inlining gkyl_ref_count.h 
// skipping file: gkyl_util.h 
// start inlining gkyl_elem_type.h 

#include <stdio.h>
#include <stdint.h>

// Type of element stored in array
enum gkyl_elem_type { GKYL_INT, GKYL_INT_64, GKYL_FLOAT, GKYL_DOUBLE, GKYL_USER };

// Array reduce operators
enum gkyl_array_op { GKYL_MIN, GKYL_MAX, GKYL_SUM };

// file types for raw IO of Gkeyll data
enum gkyl_file_type {
  GKYL_FIELD_DATA_FILE,
  GKYL_DYNVEC_DATA_FILE,
  GKYL_MULTI_RANGE_DATA_FILE,
  GKYL_BLOCK_TOPO_DATA_FILE,
  GKYL_MULTI_BLOCK_DATA_FILE
};
// ended inlining gkyl_elem_type.h 

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <math.h>

/**
 * Array object. This is an untype, undimensioned, reference counted
 * array object. All additional structure is provided else where,
 * mainly by the range object.
 */
struct gkyl_array {
  enum gkyl_elem_type type; // type of data stored in array
  size_t elemsz, ncomp; // size of elements, number of 'components'
  size_t size; // number of indices

  size_t esznc; // elemsz*ncomp
  void *data; // pointer to data
  uint32_t flags;  
  struct gkyl_ref_count ref_count;

  int nthreads, nblocks; // threads per block, number of blocks
  struct gkyl_array *on_dev; // pointer to itself or device data
#ifdef GKYL_HAVE_CUDA
  cudaStream_t iostream;
#else
  int iostream;
#endif
};

/**
 * Create new array. Delete using gkyl_array_release method.
 * 
 * @param type Type of data in array
 * @param ncomp Number of components at each index
 * @param size Number of indices 
 * @return Pointer to newly allocated array.
 */
struct gkyl_array* gkyl_array_new(enum gkyl_elem_type type, size_t ncomp, size_t size);

/**
 * Create new array from per-allocated memory. You must ensure that
 * the memory is sufficiently large for use in the array. Delete using
 * gkyl_array_release method.
 *
 * @param type Type of data in array
 * @param ncomp Number of components at each index
 * @param size Number of indices
 * @param buff Buffer to use for array data
 * @return Pointer to newly allocated array.
 */
struct gkyl_array *gkyl_array_new_from_buff(
  enum gkyl_elem_type type, size_t ncomp, size_t size, void *buff);

/**
 * Create new array with data on NV-GPU. Delete using
 * gkyl_array_release method.
 *
 * NOTE: the data member lives on GPU, but the struct lives on the
 * host.  However, the on_dev member for this cal is set to a device
 * clone of the host struct, and is what should be used to pass to
 * CUDA kernels which require the entire array struct on device.
 * 
 * @param type Type of data in array
 * @param ncomp Number of components at each index
 * @param size Number of indices 
 * @return Pointer to newly allocated array.
 */
struct gkyl_array* gkyl_array_cu_dev_new(enum gkyl_elem_type type, size_t ncomp, size_t size);

/**
 * Create new array with host-pinned data for use with NV-GPU. Delete using
 * gkyl_array_release method.
 *
 * @param type Type of data in array
 * @param ncomp Number of components at each index
 * @param size Number of indices 
 * @return Pointer to newly allocated array.
 */
struct gkyl_array* gkyl_array_cu_host_new(enum gkyl_elem_type type, size_t ncomp, size_t size);

/**
 * Returns true if array lives on NV-GPU.
 *
 * @param arr Array to check
 * @return true of array lives on NV-GPU, false otherwise
 */
bool gkyl_array_is_cu_dev(const struct gkyl_array *arr);

/**
 * Returns true if array uses external buffer for storage.
 *
 * @param arr Array to check
 * @return true if array uses external buffer for storage
 */
bool gkyl_array_is_using_buffer(const struct gkyl_array *arr);

/**
 * Copy into array: pointer to dest array is returned. 'dest' and
 * 'src' must not point to same data.
 *
 * @param dest Destination for copy.
 * @param src Source to copy from.
 * @return dest is returned
 */
struct gkyl_array* gkyl_array_copy(struct gkyl_array* dest,
  const struct gkyl_array* src);

/**
 * Copy into array using async methods for cuda arrays: pointer to dest array is returned. 'dest' and
 * 'src' must not point to same data.
 *
 * @param dest Destination for copy.
 * @param src Source to copy from.
 * @return dest is returned
 */
struct gkyl_array* gkyl_array_copy_async(struct gkyl_array* dest,
  const struct gkyl_array* src);

/**
 * Clone array: pointer to newly created array is returned.
 * 
 * @param arr Array to clone
 * @return Pointer to clone
 */
struct gkyl_array* gkyl_array_clone(const struct gkyl_array* arr);


/**
 * Fetches a pointer to the element stored at the index 'loc'.
 *
 * @param arr Array to fetch from
 * @param loc Element to fetch
 * @return Element at location 'loc'
 */

GKYL_CU_DH
static inline void*
gkyl_array_fetch(struct gkyl_array* arr, long loc)
{
  return ((char*) arr->data) + loc*arr->esznc;
}

/** Same as above, except fetches a constant pointer */
GKYL_CU_DH
static inline const void*
gkyl_array_cfetch(const struct gkyl_array* arr, long loc)
{
  return ((const char*) arr->data) + loc*arr->esznc;
}

/**
 * Acquire pointer to array. The pointer must be released using
 * gkyl_array_release method.
 *
 * @param arr Array to which a pointer is needed
 * @return Pointer to acquired array
 */
struct gkyl_array* gkyl_array_acquire(const struct gkyl_array* arr);

/**
 * Release pointer to array
 *
 * @param arr Array to release.
 */
void gkyl_array_release(const struct gkyl_array* arr);
// ended inlining gkyl_array.h 
// start inlining gkyl_range.h 

// skipping file: gkyl_util.h 
// start inlining gkyl_vargm.h 

/**
   These horrible looking set of macros allows choosing a macro based on
   number of arguments passed to it. It is ugly as hell and so I am
   putting it in its own file. Please do not use it unless you know what
   you are doing. See:

   https://stackoverflow.com/questions/11761703/overloading-macro-on-number-of-arguments/11763277#11763277
*/

// get number of arguments with __NARG__
#define __NARG__(...)  __NARG_I_(__VA_ARGS__,__RSEQ_N())
#define __NARG_I_(...) __ARG_N(__VA_ARGS__)
#define __ARG_N(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N,...) N
#define __RSEQ_N() 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0

#define _VFUNC_(name, n) name##n
#define _VFUNC(name, n) _VFUNC_(name, n)
// general definition for any function name
#define VFUNC(func, ...) _VFUNC(func, __NARG__(__VA_ARGS__)) (__VA_ARGS__)
// general definition for any function name (extra args)
#define VFUNC1(func, a, ...) _VFUNC(func, __NARG__(__VA_ARGS__)) (a, __VA_ARGS__)
// ended inlining gkyl_vargm.h 

#include <stdint.h>
#include <stdio.h>

/**
 * Series of indexing "functions" to compute linear index into range
 */
#define gkyl_ridx1(r, i1)                       \
    ((r).ac[0]+(i1)*(r).ac[1])
#define gkyl_ridx2(r, i1, i2)                   \
    ((r).ac[0]+((i1)*(r).ac[1]+(i2)*(r).ac[2]))
#define gkyl_ridx3(r, i1, i2, i3)                                       \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3]))
#define gkyl_ridx4(r, i1, i2, i3, i4)                                   \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3]+(i4)*(r).ac[4]))
#define gkyl_ridx5(r, i1, i2, i3, i4, i5)                               \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3])+((i4)*(r).ac[4]+(i5)*(r).ac[5]))
#define gkyl_ridx6(r, i1, i2, i3, i4, i5, i6)                           \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3])+((i4)*(r).ac[4]+(i5)*(r).ac[5]+(i6)*(r).ac[6]))
#define gkyl_ridx7(r, i1, i2, i3, i4, i5, i6, i7)                       \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3])+((i4)*(r).ac[4]+(i5)*(r).ac[5]+(i6)*(r).ac[6]) + (i7)*(r).ac[7])

/** Generic indexing: works for 1D-7D (VFUNC1 is defined-ed in
 * gkyl_vargm.h) */
#define gkyl_ridx(r, ...) VFUNC1(gkyl_ridx, r, __VA_ARGS__)

/** Indexing macro taking index defined as array of int */
#define gkyl_ridxn(r, idx) gkyl_range_idx(&(r), idx)

// Constants to represent lower/upper edges
enum gkyl_edge_loc { GKYL_LOWER_EDGE = 0, GKYL_UPPER_EDGE, GKYL_NO_EDGE };

// Direction and location of range
struct gkyl_range_dir_edge {
  int dir;
  enum gkyl_edge_loc eloc;
};  

/**
 * Range object, representing an N-dimensional integer index
 * set. Lower and upper limits are inclusive.
 */
struct gkyl_range {
  int ndim; // number of dimension
  int lower[GKYL_MAX_DIM]; // lower bound
  int upper[GKYL_MAX_DIM]; // upper bound (inclusive)
  long volume; // total volume of range
    
  // do not access directly
  uint32_t flags; // Flags for internal use
  int ilo[GKYL_MAX_DIM]; // for use in inverse indexer
  long ac[GKYL_MAX_DIM+1]; // coefficients for indexing
  long iac[GKYL_MAX_DIM+1]; // for use in sub-range inverse indexer
  long linIdxZero; // linear index of {0,0,...}
  int nsplit, tid; // number of splits, split ID

  // FOR CUDA ONLY
  int nthreads, nblocks; // CUDA kernel launch specifiers for range-based ops
};

/**
 * Iterator object into the range. You can read the 'idx' pointer but
 * must not modify it or any other members of this struct.
 */
struct gkyl_range_iter {
  int idx[GKYL_MAX_DIM]; // current index (do not modify)

  // do not access
  int is_first, ndim;
  long bumps_left;
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
};

/**
 * Skip-list based iterator object.
 */
struct gkyl_range_skip_iter {
  long delta; // number of contiguous elements
  struct gkyl_range range; // outer range for iteration
};

/**
 * Initialize new range object.
 *
 * @param rng Range object to initialize
 * @param ndim Dimension of range to create.
 * @param lower Lower indices of range
 * @param upper Upper indices of range
 */
void gkyl_range_init(struct gkyl_range *rng, int ndim,
  const int *lower, const int *upper);

/**
 * Create new range object from specified shape. This sets the lower
 * indices to [0,...] and upper indices to [shape[0]-1, ...].
 *
 * @param rng Range object to initialize
 * @param ndim Dimensiom of range to create.
 * @param shape Shape of region
 */

void gkyl_range_init_from_shape(struct gkyl_range *rng, int ndim,
  const int *shape);

/**
 * Create new range object from specified shape. This sets the lower
 * indices to [1,...] and upper indices to [shape[0], ...].
 *
 * @param rng Range object to initialize
 * @param ndim Dimensiom of range to create.
 * @param shape Shape of region
 */

void gkyl_range_init_from_shape1(struct gkyl_range *rng, int ndim,
  const int *shape);

/**
 * Create a new range which is a tensor product of @a a and @a b input
 * ranges.
 *
 * @param rng On output, rng = a X b
 * @param a First operand of tensor-product
 * @param b Second operand of tensor-product
 */
void gkyl_range_ten_prod(struct gkyl_range *rng, const struct gkyl_range *a,
  const struct gkyl_range *b);

/**
 * Create a new range that is the same shape as inp range, but the
 * indices are shifted in each direction by delta[dir]
 *
 * @param rng On output new shifted range
 * @param inp Input range to shift
 * @param delta Range indices are shifted by delta[dir] in each direction
 */
void gkyl_range_shift(struct gkyl_range *rng, const struct gkyl_range *inp,
  const int *delta);

/**
 * Create a new range that is the same shape as inp range, but lower
 * indices are reset to the specified ones.
 *
 * @param rng On output new reset range
 * @param inp Input range to reset
 * @param lower New lower indices
 */
void gkyl_range_reset_lower(struct gkyl_range *rng, const struct gkyl_range *inp,
  const int *lower);

/**
 * Shape in direction dir
 *
 * @param rng Range object
 * @param dir Direction to compute shape
 * @return Shape in direction dit
 */
GKYL_CU_DH
static inline int gkyl_range_shape(const struct gkyl_range *rng, int dir)
{
  return rng->upper[dir]-rng->lower[dir]+1;  
}

/**
 * Return 1 if range is a sub-range.
 *
 * @param rng Range object
 * @return 1 if true, 0 otherwise
 */
int gkyl_range_is_sub_range(const struct gkyl_range *rng);

/**
 * Return 1 if idx is inside the range.
 *
 * @param rng Range obkect
 * @return 1 if true, 0 otherwise
 */
int gkyl_range_contains_idx(const struct gkyl_range *rng, const int *idx);

/**
 * Create a sub-range from a given range. The sub-range must be fully
 * contained in the parent range or else it will be truncated. The
 * sub-range and the parent range will returns the same linear index
 * for a given index.
 *
 * @param rng New range object to initialize
 * @param bigrng Parent range object 
 * @param sublower Lower indices of sub-range
 * @param subupper Upper indices of sub-range
 */
void gkyl_sub_range_init(struct gkyl_range *rng,
  const struct gkyl_range *bigrng, const int *sublower, const int *subupper);

/**
 * Creates a new range that is a split of the given range. The only
 * place split matters is for iterators. Iterators for split-ranges
 * only walk over the set of indices owned by that split.
 *
 * @param rng Range object to split
 * @param nsplits Number of splits
 * @param tid Split ID [0, nsplits)
 * @return Split range
 */
struct gkyl_range gkyl_range_split(struct gkyl_range *rng, int nsplits, int tid);

/**
 * Return the number of elements looped over by iterator for this
 * range.
 *
 * @param rng Range object
 * @return number of elements looped over by iterator
 */
long gkyl_range_split_len(const struct gkyl_range *rng);

/**
 * Return range which has some directions removed by setting the index
 * in those directions to fixed values. The "deflated" range has lower
 * dimension than the parent 'rng' object. The indexing into the
 * returned lower dimensional range gives the same index as the
 * corresponding location in the parent range (with the missing
 * indices set to 'locDir[dir]'). 
 *
 * @param srng Deflated range.
 * @param rng Range object to deflate
 * @param remDir 'ndim' Array of flags: 0 to keep direction, 1 to remove
 * @param loc Index to set removed direction.
 */
void gkyl_range_deflate(struct gkyl_range* srng,
  const struct gkyl_range* rng, const int *remDir, const int *locDir);

/**
 * Return range which has 'dir' direction shortened to length
 * 'len', reducing the upper limit of 'range' in that direction.
 * The shortened range has the same dimensions and the same
 * start index in 'dir'.
 *
 * @param rng Shortened range.
 * @param range Range object to shorten
 * @param dir Direction to shorten
 * @param len Length of shortened direction
 */
void gkyl_range_shorten_from_above(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int len);

/**
 * Return range which has 'dir' direction shortened to length
 * 'len', increasing the lower limit of 'range' in that direction.
 * The shortened range has the same dimensions and the same
 * start index in 'dir'.
 *
 * @param rng Shortened range.
 * @param range Range object to shorten
 * @param dir Direction to shorten
 * @param len Length of shortened direction
 */
void gkyl_range_shorten_from_below(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int len);

/**
 * Return a new range that is an extension of the input range. The
 * lower index in dir is reduced by elo[dir] and upper index increased
 * by eup[dir].
 *
 * @param erng Extended range
 * @param rng Range to extend
 * @param elo Lower in dir is reduced by elo[dir]
 * @param eup Upper in dir is increased by eup[dir]
 */
void gkyl_range_extend(struct gkyl_range *erng, const struct gkyl_range *rng,
  const int *elo, const int *eup);

/**
 * Return a new range that is an extension of the input range. The
 * lower index in dir is reduced by elo[dir] and upper index increased
 * by eup[dir]. This method only extends the range in the directions
 * other than the input @a dir.
 *
 * @param erng Extended range
 * @param dir Direction to skip extension
 * @param rng Range to extend
 * @param elo Lower in dir is reduced by elo[dir]
 * @param eup Upper in dir is increased by eup[dir]
 */
void gkyl_range_perp_extend(struct gkyl_range *erng, int dir,
  const struct gkyl_range* rng, const int *elo, const int *eup);

/**
 * Return range in direction 'dir' which corresponds to the "lower
 * skin" cells.  Lower skin cells refer to the second inner-most layer
 * of cells on the lower end of the range.
 *
 * @param srng Skin range
 * @param range Range object to find lower skin cells
 * @param dir Direction to find lower skin cells in
 * @param nskin Number of skin cells
 */
void gkyl_range_lower_skin(struct gkyl_range* srng,
  const struct gkyl_range* range, int dir, int nskin);

/**
 * Return range in direction 'dir' which corresponds to the "upper
 * skin" cells.  Upper skin cells refer to the second inner-most layer
 * of cells on the upper end of the range.
 *
 * @param srng Skin range
 * @param range Range object to find upper skin cells
 * @param dir Direction to find upper skin cells in
 * @param nskin Number of skin cells
 */
void gkyl_range_upper_skin(struct gkyl_range* srng,
  const struct gkyl_range* range, int dir, int nskin);

/**
 * Create ghost and skin sub-ranges given parent *extended* range. The
 * skin and ghost ranges are sub-ranges of the parent range and DO NOT
 * include corners. For 2D, dir=1 and nghost = { 1, 1} skin and ghost
 * are the cells marked below ("S"kin, "G"ghost)
 *
 * Lower-edge:
 * +--+--+--+
 * |  |  |  |
 * +--+--+--+
 * |G |S |  |
 * +--+--+--+
 * |G |S |  |
 * +--+--+--+
 * |G |S |  |
 * +--+--+--+
 * |  |  |  |
 * +--+--+--+
 *
 * Upper-edge:
 * +--+--+--+
 * |  |  |  |
 * +--+--+--+
 * |  |S |G |
 * +--+--+--+
 * |  |S |G |
 * +--+--+--+
 * |  |S |G |
 * +--+--+--+
 * |  |  |  |
 * +--+--+--+
 *
 * @param skin On output, skin range
 * @param ghost On outout, ghost range
 * @param dir Direction in which skin/ghost are computed
 * @param edge Edge on which skin/ghost are computed
 * @param parent Range for which skin/ghost are computed
 * @param nghost Number of ghost cells in 'dir' are nghost[dir]
 */
void gkyl_skin_ghost_ranges(struct gkyl_range *skin, struct gkyl_range *ghost,
  int dir, enum gkyl_edge_loc edge, const struct gkyl_range *parent, const int *nghost);

/**
 * Create ghost and skin sub-ranges given parent *extended* range. The
 * skin and ghost ranges are sub-ranges of the parent range. The
 * ranges include the corners.  For 2D, dir=1 and nghost = { 1, 1}
 * skin and ghost are the cells marked below ("S"kin, "G"ghost)
 *
 * Lower-edge:
 * +--+--+--+
 * |G |S |  |
 * +--+--+--+
 * |G |S |  |
 * +--+--+--+
 * |G |S |  |
 * +--+--+--+
 * |G |S |  |
 * +--+--+--+
 * |G |S |  |
 * +--+--+--+
 *
 * Upper-edge:
 * +--+--+--+
 * |  |S |G |
 * +--+--+--+
 * |  |S |G |
 * +--+--+--+
 * |  |S |G |
 * +--+--+--+
 * |  |S |G |
 * +--+--+--+
 * |  |S |G |
 * +--+--+--+
 *
 * @param skin On output, skin range
 * @param ghost On outout, ghost range
 * @param dir Direction in which skin/ghost are computed
 * @param edge Edge on which skin/ghost are computed
 * @param parent Range for which skin/ghost are computed
 * @param nghost Number of ghost cells in 'dir' are nghost[dir]
 */
void gkyl_skin_ghost_with_corners_ranges(struct gkyl_range *skin, struct gkyl_range *ghost,
  int dir, enum gkyl_edge_loc edge, const struct gkyl_range *parent, const int *nghost);

/**
 * Compute intersection of two ranges. No sub-range information is
 * propagated to the new range object.
 * 
 * @param irng Intersection of r1 and r2
 * @param r1 Range to intersect
 * @param r2 Range to intersect
 * @return 1 if intersection is not-empty, 0 otherwise
 */
int gkyl_range_intersect(struct gkyl_range *irng, const struct gkyl_range *r1,
  const struct gkyl_range *r2);

/**
 * Compute intersection of two ranges. The intersection is a sub-range
 * of @a r1.
 * 
 * @param irng Intersection of r1 and r2. 
 * @param r1 Range to intersect. irng is sub-range of r1
 * @param r2 Range to intersect
 * @return 1 if intersection is not-empty, 0 otherwise
 */
int gkyl_sub_range_intersect(struct gkyl_range* irng,
  const struct gkyl_range *r1, const struct gkyl_range *r2);

/**
 * Check if range touches the lower edge of parent range in direction
 * dir.
 *
 * @param dir Direction to check
 * @param range Inner range
 * @param parent Parent range
 * @return true if range is on lower edge, false otherwise
 */
bool gkyl_range_is_on_lower_edge(int dir, const struct gkyl_range *range,
  const struct gkyl_range *parent);

/**
 * Check if range touches the upper edge of parent range in direction
 * dir.
 *
 * @param dir Direction to check
 * @param range Inner range
 * @param parent Parent range
 * @return true if range is on upper edge, false otherwise
 */
bool gkyl_range_is_on_upper_edge(int dir, const struct gkyl_range *range,
  const struct gkyl_range *parent);

/**
 * Check if @a targ range shares an edge with the @a base range. The
 * edges do not be fully shared but any edge overlap will be
 * checked.
 *
 * @param base Base range wrt which edge overlap is checked
 * @param targ Target range to check
 * @return direction and edge. Returned struct eloc is set
 *   to GKYL_NO_EDGE if ranges dont match.
 */
struct gkyl_range_dir_edge gkyl_range_edge_match(const struct gkyl_range *base,
  const struct gkyl_range *targ);
                                                     
/**
 * General indexing function. Returns linear index into the index
 * range mapped by 'range'.
 *
 * @param range Range object to index
 * @param idx Index for which to compute linear index
 */
GKYL_CU_DH
static inline long
gkyl_range_idx(const struct gkyl_range* range, const int *idx)
{
#define RI(...) gkyl_ridx(*range, __VA_ARGS__)
  switch (range->ndim) {
    case 0:
      return range->ac[0];
      break;    
    case 1:
      return RI(idx[0]); 
      break;
    case 2:
      return RI(idx[0], idx[1]);
      break;
    case 3:
      return RI(idx[0], idx[1], idx[2]);
      break;
    case 4:
      return RI(idx[0], idx[1], idx[2], idx[3]);
      break;
    case 5:
      return RI(idx[0], idx[1], idx[2], idx[3], idx[4]);
      break;
    case 6:
      return RI(idx[0], idx[1], idx[2], idx[3], idx[4], idx[5]);
      break;
    case 7:
      return RI(idx[0], idx[1], idx[2], idx[3], idx[4], idx[5], idx[6]);
      break;
  }
  return 0;
#undef RI
}

/**
 * Compute offset given relative index. So for a 2D range, idx[2] = {1,
 * 0} will compute the relative offset from {i, j} to {i+1, j} in the
 * mapping from indices to a linear integer space.
 *
 * @param range Range to find offset in.
 * @param idx Relative index for offset calculation
 * @return Relatice offset to idx.
 */
GKYL_CU_DH
static inline long
gkyl_range_offset(const struct gkyl_range* range, const int *idx)
{
  return gkyl_range_idx(range, idx) - range->linIdxZero;
}

/**
 * Inverse indexer, mapping a linear index to an N-dimension index
 * into 'range' object.
 *
 * @param range Range object to map into
 * @param loc Linear index in [0, range->volume)
 * @param idx On output, the N-dimensional index into 'range'
 */
GKYL_CU_DH
static inline void
gkyl_range_inv_idx(const struct gkyl_range *range, long loc, int *idx)
{
  long n = loc;
  for (int i=1; i<=range->ndim; ++i) {
    long quot = n/range->ac[i];
    long rem = n % range->ac[i];
    idx[i-1] = quot + range->ilo[i-1];
    n = rem;
  }
}

/**
 * Inverse indexer for use with a sub_range, mapping a linear index to
 * an N-dimension index into 'range' object.  Behavior is such that
 * loc = 0 gives idx = {0, 0, ...}.
 *
 * @param range Range object to map into
 * @param loc Linear index in [0, range->volume)
 * @param idx On output, the N-dimensional index into 'range'
 */
GKYL_CU_DH
static inline void
gkyl_sub_range_inv_idx(const struct gkyl_range *range, long loc, int *idx)
{
  long n = loc;
  for (int i=1; i<=range->ndim; ++i) {
    long quot = n/range->iac[i];
    long rem = n % range->iac[i];
    idx[i-1] = quot + range->lower[i-1];
    n = rem;
  }
}

/**
 * Create iterator. The returned iterator can be used in a 'while'
 * loop using gkyl_range_iter_next method to loop over the index set
 * spanned by the region.
 *
 * @param range Range object.
 * @return New iterator object for 'range'
 */
void gkyl_range_iter_init(struct gkyl_range_iter *iter,
  const struct gkyl_range* range);

/**
 * Create iterator, ignoring split information in range.
 *
 * @param range Range object.
 * @return New iterator object for 'range'
 */
void gkyl_range_iter_no_split_init(struct gkyl_range_iter *iter,
  const struct gkyl_range* range);

/**
 * Get next index into range. The iter->idx array holds the next
 * index. This should not be modified by the user!
 *
 * @param iter Iterator object. On exit, iter->idx has the next index
 * @return 1 if there are more indices remaining, 0 if done.
 */
int gkyl_range_iter_next(struct gkyl_range_iter *iter);

/**
 * Create skip-iterator. The returned iterator can be used in a nested
 * loop structure: a for loop inside a 'while' loop.
 *
 * @param range Range object.
 * @return New iterator object for 'range'
 */
void gkyl_range_skip_iter_init(struct gkyl_range_skip_iter *iter,
  const struct gkyl_range* range);

/**
 * Print range information to file object.
 *
 * @param range Range object to print
 * @param nm Name of range
 * @param fp File object to print range information
 */
void gkyl_print_range(const struct gkyl_range* range, const char *nm, FILE *fp);

/**
 * Compares two ranges: ranges are the same if they have the same
 * dimensions and lower and upper indices.
 *
 * @param r1 Range 1 to compare
 * @param r2 Range 2 to compare
 * @return true if ranges are same, false otherwise
 */
bool gkyl_range_compare(const struct gkyl_range* r1, const struct gkyl_range* r2);
// ended inlining gkyl_range.h 
// start inlining gkyl_rect_grid.h 

#include <stdio.h>
#include <stdbool.h>
#include <math.h>

// skipping file: gkyl_util.h 

/**
 * Rectangular grid object.
 */
struct gkyl_rect_grid {
  int ndim; // number of dimensions
  double lower[GKYL_MAX_DIM]; // lower-left corner
  double upper[GKYL_MAX_DIM]; // upper-right corner
  int cells[GKYL_MAX_DIM]; // number of cells    
  double dx[GKYL_MAX_DIM]; // cell spacing
  double cellVolume; // cell volume
};

/**
 * Create new grid object.
 *
 * @param grid Grid object to initialize.
 * @param ndim Dimension of grid
 * @param lower Coordinates of lower-left corner of grid
 * @param upper Coordinates of upper-right corner of grid
 * @param cells Number of cells in each direction
 */
void gkyl_rect_grid_init(struct gkyl_rect_grid *grid, int ndim,
  const double *lower, const double *upper, const int *cells);

/**
 * Find cell indices of point
 *
 * @param grid Grid object.
 * @param point The point to find the cell indices at.
 * @param pick_lower If point on cell boundary, pick lower cell if true, and upper if false.
 * @param known_index Any known indices of where the point is (<0 if not known).
 * @param cell_index Pointer to cell indices.
 * Asserts: point lies within cell(s) specified by knownIdx (if specified). 
 */
GKYL_CU_DH
void gkyl_rect_grid_find_cell(const struct gkyl_rect_grid *grid, const double *point,
  bool pick_lower, const int *known_index, int *cell_index);

/**
 * Get cell-center coordinates. Note that idx is a 1-based cell index,
 * i.e. the lower-left corner is (1,1,...).
 *
 * @param grid Grid object
 * @param idx Index of cell (lower-left corner has all index (1,1,...) )
 * @param xc On output, cell-center coordinates of cell 'idx'
 */
GKYL_CU_DH
static inline void
gkyl_rect_grid_cell_center(const struct gkyl_rect_grid *grid,
  const int *idx, double *xc)
{
  for (int i=0; i<grid->ndim; ++i)
    xc[i] = grid->lower[i]+(idx[i]-0.5)*grid->dx[i];
}

/**
 * Get coordinate of lower-left node. Note that idx is a 1-based cell
 * index, i.e. the lower-left corner is (1,1,...).
 *
 * @param grid Grid object
 * @param idx Index of cell (lower-left corner has all index (1,1,...) )
 * @param xn On output, coordinates of lower-left node
 */
GKYL_CU_DH
static inline void
gkyl_rect_grid_ll_node(const struct gkyl_rect_grid *grid,
  const int *idx, double *xc)
{
  for (int i=0; i<grid->ndim; ++i)
    xc[i] = grid->lower[i]+(idx[i]-1)*grid->dx[i];
}

/**
 * Get index extents in direction @a dir. The extents are inclusive.
 *
 * @param grid Grid object
 * @param dir Direction in which to get extents
 * @param ext On output, inclusive extents in direction @a dir.
 */
GKYL_CU_DH
static inline void
gkyl_rect_grid_extents(const struct gkyl_rect_grid *grid, int dir, int ext[2])
{
  ext[0] = 1; ext[1] = grid->cells[dir];
}

/**
 * Get index of point with coordinate @a xn
 *
 * @param grid Grid object
 * @param xn Coordinate of point in grid
 * @param idx On output, index of point in grid
 */
GKYL_CU_DH
static inline void
gkyl_rect_grid_coord_idx(const struct gkyl_rect_grid *grid,
  const double *xn, int *idx)
{
  for (int d=0; d<grid->ndim; ++d) {
    int ext[2]; gkyl_rect_grid_extents(grid, d, ext);
    double xlower = grid->lower[d], dx = grid->dx[d];
    idx[d] = ext[0] + (int) floor((xn[d]-xlower)/dx);
  }
}

/**
 * Compare grids
 *
 * @param grid1 Grid object to compare
 * @param grid2 Grid object to compare
 * @return true if the grids are the same, false otherwise
 */
bool gkyl_rect_grid_cmp(const struct gkyl_rect_grid *grid1, struct gkyl_rect_grid *grid2);

/**
 * Write grid data to file. File must be opened by caller of this
 * function. Data is written in binary format.
 *
 * @param grid Grid object to write
 * @param fp File handle to write to.
 */
void gkyl_rect_grid_write(const struct gkyl_rect_grid *grid, FILE *fp);

/**
 * Read data from file and initialize grid. File must be opened by
 * caller of this function.
 *
 * @param grid Grid object to initialize
 * @param fp File handle to read fromx.
 * @return True if read succeeded, false otherwise
 */
bool gkyl_rect_grid_read(struct gkyl_rect_grid *grid, FILE *fp);
// ended inlining gkyl_rect_grid.h 
// skipping file: gkyl_util.h 

// Read status flags
enum gkyl_array_rio_status {
  GKYL_ARRAY_RIO_SUCCESS = 0,
  GKYL_ARRAY_RIO_BAD_VERSION,
  GKYL_ARRAY_RIO_FOPEN_FAILED,
  GKYL_ARRAY_RIO_FREAD_FAILED,
  GKYL_ARRAY_RIO_DATA_MISMATCH,
  GKYL_ARRAY_RIO_META_FAILED
};

/**
 * Return character string corresponding to the status enum flag.
 *
 * @param status Status flag
 * @return string corresponding to flag
 */
const char* gkyl_array_rio_status_msg(enum gkyl_array_rio_status status);

// Array header data to write: this is for low-level control and is
// typically not something most users would ever encounter
struct gkyl_array_header_info {
  uint64_t file_type; // file type
  enum gkyl_elem_type etype; // element type
  uint64_t esznc; // elem sz * number of components
  uint64_t tot_cells; // total number of cells in grid
  uint64_t meta_size; // size in bytes of meta-data embedded in header
  char *meta; // meta-data as byte array
  uint64_t nrange; // number of ranges
};

/**
 * Read grid and array data header data from file. Note that only
 * HEADER is read and NOT the array data itself. If the header has
 * meta-data (meta_size > 0) then the meta char array must be freed
 * using gkyl_free.
 *
 * @param grid Grid object to read
 * @param hrd On output, Header data.
 * @param fname Name of output file (include .gkyl extension)
 * @return Status flag
 */
enum gkyl_array_rio_status gkyl_grid_sub_array_header_read(struct gkyl_rect_grid *grid,
  struct gkyl_array_header_info *hdr, const char *fname);

/**
 * Free header info if needed (only of meta_size > 0) does this call
 * actually free anything.
 *
 * @param info Header info to free
 */
void gkyl_array_header_info_release(struct gkyl_array_header_info *info);

/**
 * Write out grid and array data to file in .gkyl format so postgkyl
 * can understand it.
 *
 * @param grid Grid object to write
 * @param range Range describing portion of the array to output.
 * @param meta Meta-data to write. Set to NULL or 0 if no metadata
 * @param arr Array object to write
 * @param fname Name of output file (include .gkyl extension)
 * @return Status flag
 */
enum gkyl_array_rio_status gkyl_grid_sub_array_write(const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range, const struct gkyl_msgpack_data *meta,
  const struct gkyl_array *arr, const char *fname);

/**
 * Read grid and array data from file. The input array must be
 * pre-allocated and must be big enough to hold the read data.
 * 
 * @param grid Grid object to read
 * @param range Range describing portion of the array.
 * @param arr Array object to read
 * @param fname Name of input file
 * @return Status flag
 */
enum gkyl_array_rio_status gkyl_grid_sub_array_read(struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  struct gkyl_array *arr, const char* fname);

/**
 * Read grid and array data from file, creating a new array.
 * 
 * @param grid On outout, grid on which array is defined.
 * @param fname Name of input file
 * @return Newly created array object. NULL if failed
 */
struct gkyl_array *gkyl_grid_array_new_from_file(struct gkyl_rect_grid *grid,
  const char* fname);
// ended inlining gkyl_array_rio.h 

// Update status
struct gkyl_update_status {
  bool success; // status of update
  double dt_actual; // actual time-step taken
  double dt_suggested; // suggested stable time-step
};

// Status of restart
struct gkyl_app_restart_status {
  enum gkyl_array_rio_status io_status; // status of the file read
  int frame; // frame number of file read
  double stime; // simulation time at which data was read
};

// inputs from user for specifying range and communicator to use
struct gkyl_app_comm_low_inp {
  struct gkyl_range local_range; // local range over which App operates
  struct gkyl_comm *comm; // communicator to use
};

// BC for blocks
struct gkyl_block_physical_bcs {
  int bidx; // block index
  int dir;  // direction in which BC is specified
  enum gkyl_edge_loc edge; // which edge this BC is for
  int bc_type; // BC code
};

// Parallelization-related inputs.
struct gkyl_app_parallelism_inp {
  bool use_gpu; // Run on the GPU(s).
  int cuts[3]; // Number of subdomain in each dimension.
  struct gkyl_comm *comm; // Communicator to use.
};

// Boundary conditions on particles
enum gkyl_species_bc_type {
  GKYL_SPECIES_COPY = 0, // copy BCs
  GKYL_SPECIES_SKIP, // Do not apply any BCs to field
  GKYL_SPECIES_REFLECT, // perfect reflector
  GKYL_SPECIES_ABSORB, // Absorbing BCs
  GKYL_SPECIES_NO_SLIP, // no-slip boundary conditions
  GKYL_SPECIES_WEDGE, // specialized "wedge" BCs for RZ-theta
  GKYL_SPECIES_FUNC, // Function boundary conditions
  GKYL_SPECIES_FIXED_FUNC, // Fixed function, time-independent, boundary conditions
  GKYL_SPECIES_EMISSION, // Emission spectrum BCs
  GKYL_SPECIES_ZERO_FLUX, // Zero flux BCs; must be applied on both lower and upper BC
  GKYL_SPECIES_RECYCLE, // Recycling BCs
};

// Boundary conditions on fields
enum gkyl_field_bc_type {
  GKYL_FIELD_COPY = 0, // copy BCs
  GKYL_FIELD_SKIP, // Do not apply any BCs to field
  GKYL_FIELD_PEC_WALL, // Maxwell's perfect electrical conductor (zero normal B and zero tangent E)
  GKYL_FIELD_SYM_WALL, // Maxwell's symmetry BC (zero normal E and zero tangent B)
  GKYL_FIELD_RESERVOIR, // Reservoir Maxwell's BCs for heat flux problem
  GKYL_FIELD_WEDGE, // specialized "wedge" BCs for RZ-theta
  GKYL_FIELD_FUNC, // Function boundary conditions
  GKYL_FIELD_DIRICHLET, // Dirichlet boundary conditions
  GKYL_FIELD_NEUMANN, // Nemann boundary conditions
  GKYL_FIELD_NONE, // Do not apply any boundary conditions
};

// Type of file import for initial conditions
enum gkyl_ic_import_type {
  GKYL_IC_IMPORT_NONE = 0,
  GKYL_IC_IMPORT_F, // Import f only.
  GKYL_IC_IMPORT_AF, // Import f and scale by alpha(x).
};
// ended inlining gkyl_app.h 
// skipping file: gkyl_array.h 
// start inlining gkyl_array_ops.h 

// skipping file: gkyl_array.h 
// skipping file: gkyl_elem_type.h 
// start inlining gkyl_evalf_def.h 

#include <stddef.h>

// Wave equation object.
typedef struct gkyl_wv_eqn gkyl_wv_eqn;

/**
 * Type of function to project.
 *
 * @param t Time to evaluate function
 * @param xn Coordinates for evaluation
 * @param fout Output vector of 'num_ret_vals'
 * @param ctx Context for function evaluation. Can be NULL
 */
typedef void (*evalf_t)(double t, const double *xn, double *fout, void *ctx);

/**
 * Type of function to apply BC
 *
 * @param eqn Base equation object.
 * @param t Time at which BC is applied
 * @param ncomp Number of compontents (size of skin and ghost arrays)
 * @param skin Pointer to data in skin-cell
 * @param ghost Pointer to data in ghost-cell
 * @param ctx Context for function evaluation. Can be NULL
 */
typedef void (*wv_bc_func_t)(const struct gkyl_wv_eqn* eqn, double t, int ncomp, const double* skin, double* ghost, void* ctx);

/**
 * Type of function for use in array copy op.
 *
 * @param nc Number of elements in @a out and @a inp
 * @param out Output buffer
 * @param inp Input buffer
 * @param ctx Context for function evaluation. Can be NULL
 */
typedef void (*array_copy_func_t)(size_t nc, double *out, const double *inp, void *ctx);
// ended inlining gkyl_evalf_def.h 
// skipping file: gkyl_range.h 

GKYL_CU_DH
static inline void*
gkyl_flat_fetch(void *data, size_t loc)
{
  return ((char*) data) + loc;
}

// Struct used to pass function pointer and context to various buffer
// copy operators
struct gkyl_array_copy_func {
  array_copy_func_t func;
  void *ctx;
  uint32_t flags;

  void *ctx_on_dev; // pointer to on-device context (or itself)
  struct gkyl_array_copy_func *on_dev; // pointer to itself or device data
};

// To return diff of two arrays
struct gkyl_array_diff {
  bool is_compatible; // are arrays compatible

  // the following make sense only if is_compatible = true
  double max_abs_diff; // maximum absolute difference
  double min_abs_diff; // minmum absolute difference
  double max_rel_diff; // maximum relative difference
  double min_rel_diff; // minmum relative difference
};  

/**
 * Check if array_copy_func is on device.
 *
 * @param bc BC function to check
 * @return true if eqn on device, false otherwise
 */
bool
gkyl_array_copy_func_is_cu_dev(const struct gkyl_array_copy_func *bc);

/**
 * Clear out = val. Returns out.
 *
 * @param out Output array
 * @param val Factor to set 
 * @return out array
 */
struct gkyl_array* gkyl_array_clear(struct gkyl_array *out, double val);

/**
 * Compute out = out + a*inp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @return out array
 */
struct gkyl_array* gkyl_array_accumulate(struct gkyl_array *out,
  double a, const struct gkyl_array *inp);

/**
 * Compute out = out + a*inp[coff] where coff is a component-offset if
 * out->ncomp < inp->ncomp, or out[coff] = out[coff]+ a*inp if
 * out->ncomp > inp->ncomp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @param coff Component offset
 * @return out array
 */
struct gkyl_array* gkyl_array_accumulate_offset(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, int coff);

/**
 * Set out = a*inp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @return out array
 */
struct gkyl_array* gkyl_array_set(struct gkyl_array *out,
  double a, const struct gkyl_array *inp);

/**
 * Set out = a*inp[coff] where coff is a component-offset if
 * out->ncomp < inp->ncomp, or out[coff] = a*inp if
 * out->ncomp > inp->ncomp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @param coff Component offset
 * @return out array
 */
struct gkyl_array* gkyl_array_set_offset(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, int coff);

/**
 * Scale out = a*out. Returns out.
 *
 * @param out Output array
 * @param a Factor to scale
 * @return out array
 */
struct gkyl_array* gkyl_array_scale(struct gkyl_array *out, double a);

/**
 * Scale out = a*out. Returns out.
 *
 * @param out Output array
 * @param a Factor to scale that varies by cell
 * @return out array
 */
struct gkyl_array* gkyl_array_scale_by_cell(struct gkyl_array *out, const struct gkyl_array *a);

/**
 * Shift the k-th coefficient in every cell, out_k = a+out_k. Returns out.
 *
 * @param out Output array.
 * @param a Factor to shift k-th coefficient by.
 * @param k Coefficient to be shifted.
 * @return out array.
 */
struct gkyl_array* gkyl_array_shiftc(struct gkyl_array *out, double a, unsigned k);

/**
 * Clear out = val. Returns out.
 *
 * @param out Output array
 * @param val Factor to set 
 * @return out array
 */
struct gkyl_array* gkyl_array_clear_range(struct gkyl_array *out, double val,
  const struct gkyl_range *range);

/**
 * Compute out = out + a*inp over a range of indices.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @param range Range specifying region to accumulate
 * @return out array
 */
struct gkyl_array* gkyl_array_accumulate_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, const struct gkyl_range *range);

/**
 * Compute out = out + a*inp[coff] where coff is a component-offset if
 * out->ncomp < inp->ncomp, or out[coff] = out[coff]+ a*inp if
 * out->ncomp > inp->ncomp, over a range of indices. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @param coff Component offset
 * @return out array
 */
struct gkyl_array* gkyl_array_accumulate_offset_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, int coff, const struct gkyl_range *range);

/**
 * Set out = a*inp. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @return out array
 * @param range Range specifying region to set
 */
struct gkyl_array* gkyl_array_set_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, const struct gkyl_range *range);

/**
 * Set out = a*inp over specified ranges. Returns out.
 * input and output ranges must have the same volume.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @return out array
 * @param out_range Range specifying region of out to set
 * @param inp_range Range specifying region of inp to use
 */
struct gkyl_array* gkyl_array_set_range_to_range(struct gkyl_array *out, double a,
  const struct gkyl_array *inp, struct gkyl_range *out_range, struct gkyl_range *inp_range);

/**
 * Set out = a*inp[coff] where coff is a component-offset if
 * out->ncomp < inp->ncomp, or out[coff] = a*inp if
 * out->ncomp > inp->ncomp, over a range of indices. Returns out.
 *
 * @param out Output array
 * @param a Factor to multiply input array
 * @param inp Input array
 * @return out array
 * @param range Range specifying region to set
 */
struct gkyl_array* gkyl_array_set_offset_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, int coff, const struct gkyl_range *range);

/**
 * Scale out = a*ut. Returns out.
 *
 * @param out Output array
 * @param a Factor to scale by
 * @return out array
 * @param range Range specifying region to scale
 */
struct gkyl_array* gkyl_array_scale_range(struct gkyl_array *out,
  double a, const struct gkyl_range *range);

/**
 * Shift the k-th coefficient in every cell, out_k = a+out_k within
 * a given range. Returns out.
 *
 * @param out Output array.
 * @param a Factor to shift k-th coefficient by.
 * @param k Coefficient to be shifted.
 * @param range Range to shift coefficient k in.
 * @return out array.
 */
struct gkyl_array* gkyl_array_shiftc_range(struct gkyl_array *out, double a,
  unsigned k, const struct gkyl_range *range);

/**
 * Copy out inp. Returns out.
 *
 * @param out Output array
 * @param inp Input array
 * @param range Range specifying region to copy
 * @return out array
 */
struct gkyl_array* gkyl_array_copy_range(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range);

/**
 * Copy out inp over specified ranges. Returns out.
 * input and output ranges must have the same volume.
 *
 * @param out Output array
 * @param inp Input array
 * @param out_range Range specifying region to copy to in out array
 * @param inp_range Range specifying region to copy to from in inp array
 * @return out array
 */
struct gkyl_array* gkyl_array_copy_range_to_range(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *out_range, const struct gkyl_range *inp_range);

/**
 * Perform an "reduce" operation of data in the array.
 *
 * @param res On output, reduces values (ncomp size).
 * @param arr Array to perform reduction on.
 * @param op Reduction operators.
 */
void gkyl_array_reduce(double *res, const struct gkyl_array *arr, enum gkyl_array_op op);

/**
 * Perform an "reduce" operation of data in the array. Data is reduced
 * component-wise.
 *
 * @param res On output, reduced values (ncomp size).
 * @param arr Array to perform reduction on.
 * @param op Reduction operators.
 * @param range Range specifying region.
 */
void gkyl_array_reduce_range(double *res,
  const struct gkyl_array *arr, enum gkyl_array_op op, const struct gkyl_range *range);

/**
 * Copy region of array into a buffer. The buffer must be preallocated
 * and at least of size arr->size*arr->elemSz bytes.
 *
 * @param data Output data buffer.
 * @param arr Array to copy from
 * @param range Range specifying region to copy from
 */
void gkyl_array_copy_to_buffer(void *data, const struct gkyl_array *arr,
  const struct gkyl_range *range);

/**
 * Copy buffer into region of array. The array must be preallocated.
 *
 * @param arr Array to copy into
 * @param data Input data buffer.
 * @param range Range specifying region to copy into
 */
void gkyl_array_copy_from_buffer(struct gkyl_array *arr, const void *data,
  const struct gkyl_range *range);

/**
 * Copy region of array into a buffer, calling user-specified function
 * as the copying is done. The buffer must be preallocated and at
 * least of size arr->size*arr->elemSz bytes.
 *
 * @param data Output data buffer.
 * @param arr Array to copy from
 * @param range Range specifying region to copy from
 * @param cf Function pointer and context
 */
void gkyl_array_copy_to_buffer_fn(void *data, const struct gkyl_array *arr,
  const struct gkyl_range *range, struct gkyl_array_copy_func *cf);

/**
 * Copy region of array into a buffer, calling user-specified function
 * as the copying is done. While the copying is being performed the
 * index in @a dir is "flipped". (TODO: WHAT DOES THIS MEAN?). The
 * buffer must be preallocated and at least of size
 * arr->size*arr->elemSz bytes.
 *
 * @param data Output data buffer.
 * @param arr Array to copy from
 * @dir Direction to apply index flip
 * @param range Range specifying region to copy from
 * @param cf Function pointer and context
 */
void gkyl_array_flip_copy_to_buffer_fn(void *data, const struct gkyl_array *arr,
  int dir, const struct gkyl_range *range, struct gkyl_array_copy_func *cf);

/**
 * Return difference between two arrays. Mostly useful for testing.
 *
 * @param arr1 First array to compare
 * @param arr2 Second array to compare
 * @param range Range to compare over
 * @return diff between arrays
 */
struct gkyl_array_diff gkyl_array_diff(const struct gkyl_array *arr1,
  const struct gkyl_array *arr2, const struct gkyl_range *range);

/**
 * Host-side wrappers for array operations
 */
void gkyl_array_clear_cu(struct gkyl_array* out, double val);

void gkyl_array_accumulate_cu(struct gkyl_array* out, double a, const struct gkyl_array* inp);

void gkyl_array_accumulate_offset_cu(struct gkyl_array* out, double a, const struct gkyl_array* inp, int coff);

void gkyl_array_set_cu(struct gkyl_array* out, double a, const struct gkyl_array* inp);

void gkyl_array_set_offset_cu(struct gkyl_array* out, double a, const struct gkyl_array* inp, int coff);

void gkyl_array_scale_cu(struct gkyl_array* out, double a);

void gkyl_array_scale_by_cell_cu(struct gkyl_array* out, const struct gkyl_array* a);

void gkyl_array_shiftc_cu(struct gkyl_array* out, double a, unsigned k);

void gkyl_array_shiftc_range_cu(struct gkyl_array *out, double a, unsigned k, const struct gkyl_range *range);

/**
 * Host-side wrappers for range-based array operations
 */
void gkyl_array_clear_range_cu(struct gkyl_array *out, double val, const struct gkyl_range *range);

void gkyl_array_accumulate_range_cu(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, const struct gkyl_range *range);

void gkyl_array_accumulate_offset_range_cu(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, int coff, const struct gkyl_range *range);

void gkyl_array_set_range_cu(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, const struct gkyl_range *range);

void gkyl_array_set_range_to_range_cu(struct gkyl_array *out, double a,
  const struct gkyl_array *inp, struct gkyl_range *out_range, struct gkyl_range *inp_range);

void gkyl_array_set_offset_range_cu(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, int coff, const struct gkyl_range *range);

void gkyl_array_scale_range_cu(struct gkyl_array *out,
  double a, const struct gkyl_range *range);

void gkyl_array_copy_range_cu(struct gkyl_array *out, const struct gkyl_array* inp, 
  const struct gkyl_range *range);

void gkyl_array_copy_range_to_range_cu(struct gkyl_array *out, const struct gkyl_array* inp,
  const struct gkyl_range *out_range, const struct gkyl_range *inp_range);

void gkyl_array_copy_to_buffer_cu(void *data, const struct gkyl_array *arr, 
  const struct gkyl_range *range);

void gkyl_array_copy_from_buffer_cu(struct gkyl_array *arr, const void *data, 
  const struct gkyl_range *range);

void gkyl_array_copy_to_buffer_fn_cu(void *data, const struct gkyl_array *arr,
  const struct gkyl_range *range, struct gkyl_array_copy_func *cf);

void gkyl_array_flip_copy_to_buffer_fn_cu(void *data, const struct gkyl_array *arr,
  int dir, const struct gkyl_range *range, struct gkyl_array_copy_func *cf);
// ended inlining gkyl_array_ops.h 
// skipping file: gkyl_array_rio.h 
// start inlining gkyl_array_rio_format_desc.h 

/**

  Format description for raw Gkeyll output file.

  IF THIS FORMAT IF MODIFIED, PLEASE COPY AND THEN CHANGE THE
  DESCRIPTION SO WE HAVE THE OLDER VERSIONS DOCUMENTED HERE. UPDATE
  VERSION BY 1 EACH TIME YOU CHANGE THE FORMAT.

  The format of the gkyl binary output is as follows.

  ## Version 0: Jan 2021. Created by A.H. Note Version 0 has no header
     information

  Data      Type and meaning
  --------------------------
  ndim      uint64_t Dimension of field
  cells     uint64_t[ndim] number of cells in each direction
  lower     float64[ndim] Lower bounds of grid
  upper     float64[ndim] Upper bounds of grid
  esznc     uint64_t Element-size * number of components in field
  size      uint64_t Total number of cells in field
  DATA      size*esznc bytes of data

  ## Version 1: May 9th 2022. Created by A.H

  Data      Type and meaning
  --------------------------
  gkyl0     5 bytes
  version   uint64_t
  file_type uint64_t (See header gkyl_elem_type.h for file types)
  meta_size uint64_t Number of bytes of meta-data
  DATA      meta_size bytes of data. This is in msgpack format

  * For file_type = 1 (field) the above header is followed by

  real_type uint64_t. Indicates real type of data
  ndim      uint64_t Dimension of field
  cells     uint64_t[ndim] number of cells in each direction
  lower     float64[ndim] Lower bounds of grid
  upper     float64[ndim] Upper bounds of grid
  esznc     uint64_t Element-size * number of components in field
  size      uint64_t Total number of cells in field
  DATA      size*esznc bytes of data

  * For file_type = 2 (dynvec) the above header is followed by

  real_type uint64_t. Indicates real type of data
  esznc     uint64_t Element-size * number of components in field
  size      uint64_t Total number of cells in field
  TIME_DATA float64[size] bytes of data
  DATA      size*esznc bytes of data

  * For file_type = 3 (multi-range field) the above header is followed by

  real_type uint64_t. Indicates real type of data
  ndim      uint64_t Dimension of field
  cells     uint64_t[ndim] number of cells in each direction
  lower     float64[ndim] Lower bounds of grid
  upper     float64[ndim] Upper bounds of grid
  esznc     uint64_t Element-size * number of components in field
  size      uint64_t Total number of cells in field
  nrange    uint64_t Number of ranges stored in this file

  For each of the nrange ranges in the field the following data is
  present

  loidx     uint64_t[ndim] Index of lower-left corner of the range
  upidx     uint64_t[ndim] Index of upper-right corner of the range
  size      uint64_t Total number of cells in range
  DATA      size*esznc bytes of data

  Note: the global range in Gkeyll, of which each range is a part,
  is 1-indexed.

  * For file_type = 4 (block topology) there is no additional
    data. The topology is store in mpack format in the following way.

    {
      ndim = dimension of topology,
      num_blocks = number of blocks,
      connections = array of connection data
    }

    The array of connection data is a list of integers with bid, dir
    and edge (in this order) for each block, in each direction and
    lower and upper edges (in this order). See struct gkyl_block_topo
    and btopo_create_mpack in file block_topo.c file for details.

  * For file_type = 5 is written out for multi-block simulations and
    has no additional data. The information is all in the mpack format:

    {
      frame = fname number,
      stime = simulation time,
      topo_file_name = file in which topology is stored
    }

 */

#include <stdio.h>

// The following utility functions allow to determine the size in
// bytes of the headers for gkyl output files.

size_t gkyl_base_hdr_size(size_t meta_sz);
size_t gkyl_file_type_1_hrd_size(int ndim);
size_t gkyl_file_type_1_partial_hrd_size(int ndim);
size_t gkyl_file_type_2_hrd_size(void);
size_t gkyl_file_type_3_hrd_size(int ndim);
size_t gkyl_file_type_3_range_hrd_size(int ndim);

/**
 * Return .gkyl file type. Returns -1 if file does not exist or is not
 * a gkyl file.
 *
 * @param fname File name
 * @param file type.
 */
int gkyl_get_gkyl_file_type(const char *fname);
// ended inlining gkyl_array_rio_format_desc.h 
// start inlining gkyl_comm.h 

// skipping file: gkyl_array.h 
// skipping file: gkyl_array_ops.h 
// skipping file: gkyl_array_rio.h 
// skipping file: gkyl_elem_type.h 
// start inlining gkyl_rect_decomp.h 

// skipping file: gkyl_range.h 
// skipping file: gkyl_rect_grid.h 
// skipping file: gkyl_ref_count.h 

// Decomposition object
struct gkyl_rect_decomp {
  int ndim; // dimension of decomposition
  int ndecomp; // number of sub-domains
  struct gkyl_range parent_range; // range that was decomposed
  struct gkyl_range *ranges; // decomposed ranges

  struct gkyl_ref_count ref_count;
};

// List of neighbors 
struct gkyl_rect_decomp_neigh {
  int num_neigh; // number of neighbors
  const int *neigh; // list of neighbors

  // following information is not typically needed
  const int *dir; // direction in which neigh[i] range is located
  const int *edge; // edge on which neigh[i] range is located
};

/**
 * Create a new decomposition of @a range, given @a cuts in each
 * direction. The total number of decomposed ranges are product of all
 * cuts. The decomposed ranges index independent of @a range,
 * i.e. decomposed ranges are NOT sub-ranges of @a range.
 *
 * @param ndim Number of dimensions
 * @param cuts Cuts in each direction.
 * @param range Range to decompose
 * @return Decomposition of @a range
 */
struct gkyl_rect_decomp* gkyl_rect_decomp_new_from_cuts(int ndim, const int cuts[],
  const struct gkyl_range *range);

/**
 * Create a new decomposition given @a cuts and cells in each
 * direction. The total number of decomposed ranges are product of all
 * cuts.
 *
 * @param ndim Number of dimensions
 * @param cuts Cuts in each direction.
 * @param cells Number of cells in each direction
 * @return Decomposition of range based on cuts
 */
struct gkyl_rect_decomp *gkyl_rect_decomp_new_from_cuts_and_cells(int ndim,
  const int cuts[], const int cells[]);

/**
 * Create a new decomposition from a given decomposition. The new
 * decomposition extends each region by a tensor product with @a
 * arange.
 *
 * @param arange Range to extend by
 * @return New extended decomposition
 */
struct gkyl_rect_decomp *gkyl_rect_decomp_extended_new(const struct gkyl_range *arange,
  const struct gkyl_rect_decomp *decomp);

/**
 * Acquire a pointer to the decomposition.
 *
 * @param decomp Decom to acquire pointer to
 * @return New decomposition
 */
struct gkyl_rect_decomp* gkyl_rect_decomp_acquire(const struct gkyl_rect_decomp *decomp);

/**
 * Check if decomposition is  a valid covering of the range.
 *
 * NOTE: This function internally allocates memory over the complete
 * parent range. This can be a problem if the parent range is huge.
 *
 * @param decomp Demposition to check
 * @return true if this is a valid covering
 */
bool gkyl_rect_decomp_check_covering(const struct gkyl_rect_decomp *decomp);

/**
 * Compute the neighbor of range @a nidx. The returned object must be
 * freed using the gkyl_rect_decomp_neigh_release call.
 *
 * @param decomp Decomposition object
 * @param inc_corners If true, corner neighbors are also included
 * @param nidx Index of range for which neighbor data is needed
 * @return Neighbor list for range nidx
 */
struct gkyl_rect_decomp_neigh* gkyl_rect_decomp_calc_neigh(
  const struct gkyl_rect_decomp *decomp, bool inc_corners, int nidx);

/**
 * Compute the periodic neighbor of range @a nidx in the specified
 * direction. The returned object must be freed using the
 * gkyl_rect_decomp_neigh_release call.
 *
 * @param decomp Decomposition object
 * @param dir Direction to compute periodic neighbors
 * @param inc_corners If true, corner neighbors are also included
 * @param nidx Index of range for which neighbor data is needed
 * @return Periodic neighbor list for range nidx
 */
struct gkyl_rect_decomp_neigh* gkyl_rect_decomp_calc_periodic_neigh(
  const struct gkyl_rect_decomp *decomp, int dir, bool inc_corners, int nidx);

/**
 * Free neighbor memory
 *
 * @param ng Neighbor data to free
 */
void gkyl_rect_decomp_neigh_release(struct gkyl_rect_decomp_neigh *ng);

/**
 * Compute cumulative offet of @a nidx range in the decomp. The
 * cumulative offset is the global linear index of the first cell in
 * the local range.
 *
 * @param decomp Decomposition object
 * @param nidx Index of range for which offset is needed
 * @return Offest of first cell in range[nidx]
 */
long gkyl_rect_decomp_calc_offset(const struct gkyl_rect_decomp *decomp, int nidx);

/**
 * Free decomposition.
 *
 * @param decomp Decomposition to free
 */
void gkyl_rect_decomp_release(struct gkyl_rect_decomp *decomp);

// The functions below are utility functions to construct properly
// nested ranges that extend over the grid or over local ranges, given
// ghost cells.

/**
 * Create range over global region given cells in each direction.
 *
 * @param ndim Grid dimension
 * @param cells Number of cells in each direction
 * @param range On output, global range
 */
void gkyl_create_global_range(int ndim, const int *cells, struct gkyl_range *range);

/**
 * Create range and extended ranges from grid and ghost-cell data. The
 * range is a sub-range of the extended range.
 *
 * @param grid Grid to compute ranges for
 * @param nghost Number of ghost-cells in each direction
 * @param ext_range On output, extended range spanning grid+ghost-cells
 * @param range On output, range spanning grid. Sub-range of ext_range.
 */
void gkyl_create_grid_ranges(const struct gkyl_rect_grid *grid,
  const int *nghost, struct gkyl_range *ext_range,
  struct gkyl_range *range);

/**
 * Create range and extended ranges from given range and ghost-cell
 * data. The range is a sub-range of the extended range.
 *
 * @param inrange Input range to use
 * @param nghost Number of ghost-cells in each direction
 * @param ext_range On output, extended range spanning inrange+ghost-cells
 * @param range On output, range same as inrange, but sub-range of ext_range.
 */
void gkyl_create_ranges(const struct gkyl_range *inrange,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range);

/**
 * Return the cuts used to create the the decomposition object.
 * 
 * @param decomp Decomposition object.
 * @param cuts Output cuts in each direction.
 */
void gkyl_rect_decomp_get_cuts(struct gkyl_rect_decomp* decomp, int* cuts);
// ended inlining gkyl_rect_decomp.h 
// skipping file: gkyl_rect_grid.h 
// skipping file: gkyl_ref_count.h 

// Structure holding data and function pointers to communicate various
// Gkeyll objects across multi-region or multi-block domains
struct gkyl_comm {
  char id[128]; // string ID for communcator
  bool has_decomp; // flag to indicate if comm has an associated decomp

  struct gkyl_ref_count ref_count; // reference count
};

/**
 * Get rank of communicator.
 *
 * @param comm Communicator
 * @param rank On output, the rank
 * @return error code: 0 for success
 */
int gkyl_comm_get_rank(struct gkyl_comm *comm, int *rank);

/**
 * Get number of ranks in communicator
 *
 * @param comm Communicator
 * @param rank On output, the rank
 * @return error code: 0 for success
 */
int gkyl_comm_get_size(struct gkyl_comm *comm, int *sz);

/**
 * All reduce values across domains.
 *
 * @param comm Communicator
 * @param type Data-type of element
 * @param op Operator to use in reduction
 * @param nelem Number of elemets in inp and out
 * @param inp Local values on domain
 * @param out Reduced values
 * @return error code: 0 for success
 */
int gkyl_comm_allreduce(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp, void *out);

/**
 * All reduce values across domains on the host/MPI communicator.
 *
 * @param comm Communicator
 * @param type Data-type of element
 * @param op Operator to use in reduction
 * @param nelem Number of elemets in inp and out
 * @param inp Local values on domain
 * @param out Reduced values
 * @return error code: 0 for success
 */
int gkyl_comm_allreduce_host(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp, void *out);

/**
 * Gather all local data into a global array on each process.
 *
 * @param comm Communicator
 * @param local Local range for array
 * @param global Global range for array
 * @param array_local Local array
 * @param array_global Global array
 * @return error code: 0 for success
 */
int gkyl_comm_array_allgather(struct gkyl_comm *comm, 
  const struct gkyl_range *local, const struct gkyl_range *global,
  const struct gkyl_array *array_local, struct gkyl_array *array_global);

/**
 * Gather all local data on host into a global array on each process.
 *
 * @param comm Communicator
 * @param local Local range for array
 * @param global Global range for array
 * @param array_local Local array
 * @param array_global Global array
 * @return error code: 0 for success
 */
int gkyl_comm_array_allgather_host(struct gkyl_comm *comm, 
  const struct gkyl_range *local, const struct gkyl_range *global,
  const struct gkyl_array *array_local, struct gkyl_array *array_global);

/**
 * Broadcast an array to other processes.
 *
 * @param comm Communicator.
 * @param array_send Array to send (only in rank 'root').
 * @param array_recv Receive buffer array.
 * @param root Broadcasting process.
 * @return error code: 0 for success
 */
int gkyl_comm_array_bcast(struct gkyl_comm *comm, 
  const struct gkyl_array *array_send, struct gkyl_array *array_recv, int root);

/**
 * Broadcast a host side array to other processes.
 *
 * @param comm Communicator.
 * @param array_send Array to send (only in rank 'root').
 * @param array_recv Receive buffer array.
 * @param root Broadcasting process.
 * @return error code: 0 for success
 */
int gkyl_comm_array_bcast_host(struct gkyl_comm *comm, 
  const struct gkyl_array *array_send, struct gkyl_array *array_recv, int root);

/**
 * Synchronize array across domain.
 *
 * @param comm Communicator
 * @param local Local range for array: sub-range of local_ext
 * @param local_ext Extended range, i.e. range over which array is defined
 * @param array Array to synchronize
 * @return error code: 0 for success
 */
int gkyl_comm_array_sync(struct gkyl_comm *comm,
  const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  struct gkyl_array *array);

/**
 * Synchronize array across domain in periodic directions.
 *
 * @param comm Communicator
 * @param local Local range for array: sub-range of local_ext
 * @param local_ext Extended range, i.e. range over which array is defined
 * @param nper_dirs Number of periodic directions
 * @param per_dirs Directions that are periodic
 * @param array Array to synchronize
 * @return error code: 0 for success
 */
int gkyl_comm_array_per_sync(struct gkyl_comm *comm,
  const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  int nper_dirs, const int *per_dirs,
  struct gkyl_array *array);

/**
 * Barrier across domains
 *
 * @param comm Communcator
 * @return error code: 0 for success
 */
int gkyl_comm_barrier(struct gkyl_comm *comm);


/**
 * Start and end a group call
 * 
 * @param comm Communcator
 */
void gkyl_comm_group_call_start(struct gkyl_comm *comm);
void gkyl_comm_group_call_end(struct gkyl_comm *comm);

/**
 * Create a new communcator that extends the communcator to work on a
 * extended domain specified by erange. (Each range handled by the
 * communicator is extended by a tensor-product with erange). The
 * returned communicator must be freed by calling gkyl_comm_release.
 *
 * @param comm Communicator
 * @param erange Range to extend by
 * @return Newly created communicator
 */
struct gkyl_comm* gkyl_comm_extend_comm(const struct gkyl_comm *comm,
  const struct gkyl_range *erange);

/**
 * Split a communicator into a new communcator based on color. All
 * ranks with the same color will form the new communcator. In the input @a
 * new_decomp can be NULL.
 *
 * @param comm Communicator.
 * @param color All ranks of same color will share a communicator.
 * @param new_decomp Decomp object to associate new communicator. Can be NULL
 * @return Newly created communicator
 */
struct gkyl_comm* gkyl_comm_split_comm(const struct gkyl_comm *comm, int color,
  struct gkyl_rect_decomp *new_decomp);

/**
 * Create a new communicator that incudes a subset of ranks in @a
 * comm. This call can return a NULL if the communicator is not valid
 * on the parent calling rank. In this case the is_valid flag is also
 * set to false.
 *
 * @param comm Communicator.
 * @param nrank Number of ranks to include
 * @param ranks List of ranks to include
 * @param new_decomp Decomp object to associate new communicator. Can be NULL
 * @param is_valid On output, true if comm is usable, false otherwise
 * @return Newly created communicator
 */
struct gkyl_comm* gkyl_comm_create_comm_from_ranks(const struct gkyl_comm *comm, int nranks,
  const int *ranks, struct gkyl_rect_decomp *new_decomp,
  bool *is_valid);

/**
 * Acquire pointer to communicator
 *
 * @param comm Communicator to to get acquire
 * @return Acquired comm obj pointer
 */
struct gkyl_comm* gkyl_comm_acquire(const struct gkyl_comm *comm);

/**
 * Release communicator memory.
 *
 * @param comm Communicator to release
 */
void gkyl_comm_release(const struct gkyl_comm *comm);
// ended inlining gkyl_comm.h 
// start inlining gkyl_comm_io.h 
// skipping file: gkyl_comm.h 

/**
 * Write out grid and array data to file in .gkyl format so postgkyl
 * can understand it.
 *
 * @param comm Communicator
 * @param grid Grid object to write
 * @param range Range describing portion of the array to output.
 * @param meta Meta-data to write. Set to NULL or 0 if no metadata
 * @param arr Array object to write
 * @param fname Name of output file (include .gkyl extension)
 * @return Status flag: 0 if write succeeded, 'errno' otherwise
 */
int gkyl_comm_array_write(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_msgpack_data *meta,
  const struct gkyl_array *arr, const char *fname);

/**
 * Read array data from .gkyl format. The input grid must be
 * pre-computed and must match the grid in the array. An error is
 * returned if this is not the case.
 *
 * @param comm Communicator
 * @param grid Grid object for read
 * @param range Range describing portion of the array to read.
 * @param arr Array object to read
 * @param fname Name of output file (include .gkyl extension)
 * @return Status flag: 0 if write succeeded, 'errno' otherwise
 */
int gkyl_comm_array_read(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  struct gkyl_array *arr, const char *fname);
// ended inlining gkyl_comm_io.h 
// start inlining gkyl_const.h 

#define GKYL_PI (3.141592653589793238462643383279502884)
#define GKYL_E  (2.718281828459045235360287471352662497)
#define GKYL_SPEED_OF_LIGHT (299792458.0) // m/s
#define GKYL_PLANCKS_CONSTANT_H (6.62606896e-34) // joule*seconds
#define GKYL_ELECTRON_MASS (9.10938215e-31) // Kg
#define GKYL_PROTON_MASS (1.672621637e-27) // Kg
#define GKYL_MASS_UNIT (1.66053907e-27) // Kg
#define GKYL_ELEMENTARY_CHARGE (1.602176487e-19) // Coulomb
#define GKYL_BOLTZMANN_CONSTANT (1.3806488e-23)
#define GKYL_EPSILON0 (8.854187817620389850536563031710750260608e-12) // farad/meter
#define GKYL_MU0  (12.56637061435917295385057353311801153679e-7) // newtons/ampere/ampere
#define GKYL_EV2KELVIN (GKYL_ELEMENTARY_CHARGE/GKYL_BOLTZMANN_CONSTANT)
// ended inlining gkyl_const.h 
// start inlining gkyl_dynvec.h 

// skipping file: gkyl_array.h 
// skipping file: gkyl_elem_type.h 

#include <stdbool.h>
#include <stddef.h>

/** Dynamic vector to store time-dependent diagnostics */
typedef struct gkyl_dynvec_tag* gkyl_dynvec;

// Element type and number of components
struct gkyl_dynvec_etype_ncomp {
  enum gkyl_elem_type type; // type of data stored in vector
  size_t ncomp; // number of 'components'
};

/**
 * Create a new new dynvec. Delete using gkyl_dynvec_release method.
 * 
 * @param type Type of data in vector
 * @param ncomp Number of components
 * @return Newly allocated vector.
 */
gkyl_dynvec gkyl_dynvec_new(enum gkyl_elem_type type, size_t ncomp);

/**
 * Get element type stored in dynvec.
 *
 * @param vec Dynvec object
 * @return Element type
 */
int gkyl_dynvec_elem_type(gkyl_dynvec vec);

/**
 * Get number of components stored
 *
 * @param vec Dynvec object
 * @return Number of components
 */
int gkyl_dynvec_ncomp(gkyl_dynvec vec);

/**
 * Reserve @a rsize more elements so additional @a rsize append calls
 * do not require memory allocations.
 *
 * @param vec Vector to reserve data for
 * @param rsize Additional number of elements to reserve
 */
void gkyl_dynvec_reserve_more(gkyl_dynvec vec, size_t rsize);

/**
 * Append data to vector. You must ensure the data has the proper type
 * and ncomp number of elements.
 * 
 * @param vec Vector to append to
 * @param tm Time-stamp of data
 * @param data to append
 */
void gkyl_dynvec_append(gkyl_dynvec vec, double tm, const void *data);

/**
 * Get data at @a idx location. You must ensure the data has the
 * proper type and ncomp number of elements.
 * 
 * @param vec Vector
 * @param idx Index of data to fetch
 * @param data On return, data is copied in this buffer
 * @return If idx is not in range, false is returned
 */
bool gkyl_dynvec_get(const gkyl_dynvec vec, size_t idx, void *data);

/**
 * Get last appended data. You must ensure the data has the proper
 * type and ncomp number of elements.
 * 
 * @param vec Vector
 * @param data On return, data is copied in this buffer
 * @return If vector is empty, return false.
 */
bool gkyl_dynvec_getlast(const gkyl_dynvec vec, void *data);

/**
 * Get time-stamp for data at index idx.
 * 
 * @param vec Vector
 * @return Time-stamp of last data at idx
 */
double gkyl_dynvec_get_tm(const gkyl_dynvec vec, size_t idx);

/**
 * Get last appended data time-stamp.
 * 
 * @param vec Vector
 * @return Time-stamp of last appened data
 */
double gkyl_dynvec_getlast_tm(const gkyl_dynvec vec);

/**
 * Get size of dynvec.
 * 
 * @param vec Vector 
 * @return Number of elements in vector
 */
size_t gkyl_dynvec_size(const gkyl_dynvec vec);

/**
 * Get capacity of dynvec.
 *
 * @param vec Vector
 * @return Capacity (elements that can be stored witout reallocation)
 */
size_t gkyl_dynvec_capacity(const gkyl_dynvec vec);

/**
 * Get capacity of dynvec.
 * 
 * @param vec Vector 
 * @return Capacity of vector
 */
size_t gkyl_dynvec_capacity(const gkyl_dynvec vec);

/**
 * Clear contents of the vector
 *
 * @param vec Vector to clear
 */
void gkyl_dynvec_clear(gkyl_dynvec vec);

/**
 * Clear contents of the vector, but keep the last @a num inserted
 * elements, shifting them to the start of the vector. This is
 * typically useful when the vector has been written to file and the
 * data needs to be flushed, but the vector still is in use.
 *
 * If @a num is larger than the number of elements in the vector
 * then nothing is done.
 *
 * @param vec Vector to clear
 * @param num Number of final elements to keep.
 */
void gkyl_dynvec_clear_all_but(gkyl_dynvec vec, size_t num);

/**
 * Acquire a reference to the dynvec. Delete using gkyl_dynvec_release method.
 *
 * @param vec Vector to acquire reference from
 * @return Dynamic vector
 */
gkyl_dynvec gkyl_dynvec_acquire(const gkyl_dynvec vec);

/**
 * Write out dynvec to file. File is overwritten with new data, and
 * existing contents will lost.
 *
 * @param vec Vector to write
 * @param fname Name of output file.
 * @return 0 if succeeded.
 */
int gkyl_dynvec_write(const gkyl_dynvec vec, const char *fname);

/**
 * Write out dynvec to file. The dynvec is appened to the end of the
 * file if it already exists.
 *
 * @param vec Vector to write
 * @param fname Name of output file.
 * @return 0 if succeeded.
 */
int gkyl_dynvec_awrite(const gkyl_dynvec vec, const char *fname);

/**
 * Read number of components from the dynvec file.
 *
 * @param Name of input file
 * @return Element type and number of components
 */
struct gkyl_dynvec_etype_ncomp gkyl_dynvec_read_ncomp(const char *fname);

/**
 * Read dynvector from file, appending data to end of the
 * vector. Existing data in vector is retained, and read data is
 * appended.
 *
 * @param vec Vector to read into
 * @param fname Name of input file.
 */
bool gkyl_dynvec_read(gkyl_dynvec vec, const char *fname);

/**
 * Convert contents of dynvector to array. Time mesh is not copied to
 * the array but it returned as a seperate array. The input arrays
 * must be preallocated to be big enough to contain all the data.
 *
 * @param vec Dynvector to convert
 * @param tm_mesh On output, time-mesh of data
 * @param dyndata On output, data in dynamic array
 */
void gkyl_dynvec_to_array(const gkyl_dynvec vec, struct gkyl_array *tm_mesh,
  struct gkyl_array *dyndata);

/**
 * Release dynvec.
 *
 * @param vec Vector to release
 */
void gkyl_dynvec_release(gkyl_dynvec vec);

// ended inlining gkyl_dynvec.h 
// skipping file: gkyl_elem_type.h 
// start inlining gkyl_eval_offset_fd.h 

// skipping file: gkyl_array.h 
// skipping file: gkyl_evalf_def.h 
// skipping file: gkyl_range.h 
// skipping file: gkyl_rect_grid.h 

// Object type
typedef struct gkyl_eval_offset_fd gkyl_eval_offset_fd;

// Struct to describe offset with respect to the cell-center: 0.0 is
// cell center, -0.5 the left edge and 0.5 the right edge.
struct gkyl_offset_descr {
  double od_off[GKYL_MAX_DIM];
};

// input packaged as a struct
struct gkyl_eval_offset_fd_inp {
  const struct gkyl_rect_grid *grid; // grid on which to project

  int num_ret_vals; // number of return values in eval function
  struct gkyl_offset_descr *offsets; // size num_ret_vals
  
  evalf_t eval; // function to project
  void *ctx; // function context
};

/**
 * Create new updater to evaluate function on a finite-difference
 * grid. Each returned value is evaluated on an internal node that is
 * potentially offset from the cell-center. This updater is useful for
 * initializing finite-difference grids.
 *
 * @param inp Input parameters
 * @return New updater pointer.
 */
gkyl_eval_offset_fd* gkyl_eval_offset_fd_new(const struct gkyl_eval_offset_fd_inp *inp);

/**
 * Compute function of nodes. The update_rng MUST be a sub-range of
 * the range on which the array is defined. That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up Updater to run
 * @param tm Time at which projection must be computed
 * @param update_rng Range on which to run projection.
 * @param out Output array
 */
void gkyl_eval_offset_fd_advance(const gkyl_eval_offset_fd *up,
  double tm, const struct gkyl_range *update_rng, struct gkyl_array *out);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_eval_offset_fd_release(gkyl_eval_offset_fd *up);
// ended inlining gkyl_eval_offset_fd.h 
// skipping file: gkyl_evalf_def.h 
// start inlining gkyl_fv_proj.h 

// skipping file: gkyl_array.h 
// skipping file: gkyl_evalf_def.h 
// skipping file: gkyl_range.h 
// skipping file: gkyl_rect_decomp.h 

// Object type
typedef struct gkyl_fv_proj gkyl_fv_proj;

/**
 * Create new updater to compute cell-average of function on a
 * grid. Free using gkyl_fv_proj_release method.
 *
 * @param grid Grid object
 * @param num_quad Number of quadrature nodes
 * @param num_ret_vals Number of values 'eval' sets
 * @param eval Function to project (See gkyl_proj_on_basis for signature).
 * @param ctx Context for function evaluation. Can be NULL.
 * @return New updater pointer.
 */
gkyl_fv_proj* gkyl_fv_proj_new(const struct gkyl_rect_grid *grid,
  int num_quad, int num_ret_vals, evalf_t eval, void *ctx);

/**
 * Compute cell averages. The update_rng MUST be a sub-range of
 * the range on which the array is defined. That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param pob Project on basis updater to run
 * @param tm Time at which projection must be computed
 * @param update_rng Range on which to run projection.
 * @param out Output array
 */
void gkyl_fv_proj_advance(const gkyl_fv_proj *pob,
  double tm, const struct gkyl_range *update_rng, struct gkyl_array *out);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_fv_proj_release(gkyl_fv_proj* pob);
// ended inlining gkyl_fv_proj.h 
// start inlining gkyl_fvec.h 

// skipping file: gkyl_alloc.h 

struct gkyl_fvec_header { size_t len, capacity; };
#define _gkyl_fvec_header_(arr) ((struct gkyl_fvec_header*) (arr) -  1)
#define GKYL_FVEC_INIT_SIZE 256
#define _gkyl_fvec_resize_(hdr, arr)                                          \
    do {                                                                \
      if (hdr->len >= hdr->capacity) {                                  \
        size_t cap = hdr->capacity = 2*hdr->capacity;                   \
        hdr = gkyl_realloc(hdr, sizeof(struct gkyl_fvec_header) + cap*sizeof(*arr)); \
        arr = (void*) (hdr + 1);                                        \
      }                                                                 \
    } while(0)

/**
 * Push value @a val at the end of array @a arr.
 */
#define gkyl_fvec_push(arr, val)                                        \
    do {                                                                \
      if (0 == arr) {                                                   \
        size_t cap = GKYL_FVEC_INIT_SIZE;                               \
        struct gkyl_fvec_header *hdr = gkyl_malloc(sizeof(struct gkyl_fvec_header) + cap*sizeof(*arr)); \
        hdr->len = 0; hdr->capacity = cap;                              \
        arr = (void*) (hdr + 1);                                        \
      }                                                                 \
      struct gkyl_fvec_header *hdr = _gkyl_fvec_header_(arr);           \
      _gkyl_fvec_resize_(hdr, arr);                                     \
    (arr)[hdr->len++] = val;                                            \
    } while (0)

/**
 * Returns capacity of vector (amount of space allocated) in @a
 * arr. This is a multiple of GKYL_FVEC_INIT_SIZE
 */
#define gkyl_fvec_capacity(arr) (0 == arr ? 0 : _gkyl_fvec_header_(arr)->capacity)

/**
 * Returns number of elements in vector @a arr
 */
#define gkyl_fvec_size(arr) (0 == arr ? 0 :_gkyl_fvec_header_(arr)->len)

/**
 * Frees the array.
 */
#define gkyl_fvec_free(arr)                                             \
    do {                                                                \
      if (0 != arr) gkyl_free( _gkyl_fvec_header_(arr) );               \
    } while (0)
// ended inlining gkyl_fvec.h 
// start inlining gkyl_gauss_quad_data.h 

#include <math.h>
#include <stdlib.h>

// Maximum order ordinate/weight data
static const int gkyl_gauss_max = 8;

// Ordinates
static const double gkyl_gauss_ordinates_1[] =
{ 0.0 };
static const double gkyl_gauss_ordinates_2[] =
{ -0.5773502691896257645091, 0.5773502691896257645091 };
static const double gkyl_gauss_ordinates_3[] =
{ -0.7745966692414833770359, 0, 0.7745966692414833770359 };
static const double gkyl_gauss_ordinates_4[] =
{ -0.8611363115940525752239, -0.3399810435848562648027, 0.3399810435848562648027, 0.8611363115940525752239 };
static const double gkyl_gauss_ordinates_5[] =
{ -0.9061798459386639927976, -0.5384693101056830910363, 0, 0.5384693101056830910363, 0.9061798459386639927976 };
static const double gkyl_gauss_ordinates_6[] =
{ -0.9324695142031520278123, -0.6612093864662645136614, -0.2386191860831969086305, 0.2386191860831969086305, 0.6612093864662645136614, 0.9324695142031520278123 };
static const double gkyl_gauss_ordinates_7[] =
{ -0.9491079123427585245262, -0.7415311855993944398639, -0.4058451513773971669066, 0, 0.4058451513773971669066, 0.7415311855993944398639, 0.9491079123427585245262 };
static const double gkyl_gauss_ordinates_8[] =
{ -0.960289856497536231684, -0.7966664774136267395916, -0.5255324099163289858177, -0.1834346424956498049395, 0.1834346424956498049395, 0.525532409916328985818, 0.796666477413626739592, 0.9602898564975362316836 };

// Weights
static const double gkyl_gauss_weights_1[] =
{ 2.0 };
static const double gkyl_gauss_weights_2[] =
{ 1.0, 1.0 };
static const double gkyl_gauss_weights_3[] =
{ 0.5555555555555555555556, 0.888888888888888888889, 0.555555555555555555556 };
static const double gkyl_gauss_weights_4[] =
{ 0.3478548451374538573731, 0.6521451548625461426269, 0.652145154862546142627, 0.3478548451374538573731 };
static const double gkyl_gauss_weights_5[] =
{ 0.2369268850561890875143, 0.4786286704993664680413, 0.568888888888888888889, 0.478628670499366468041, 0.236926885056189087514 };
static const double gkyl_gauss_weights_6[] =
{ 0.1713244923791703450403, 0.36076157304813860757, 0.46791393457269104739, 0.46791393457269104739, 0.36076157304813860757, 0.1713244923791703450403 };
static const double gkyl_gauss_weights_7[] =
{ 0.1294849661688696932706, 0.279705391489276667901, 0.38183005050511894495, 0.4179591836734693877552, 0.38183005050511894495, 0.279705391489276667901, 0.1294849661688696932706 };
static const double gkyl_gauss_weights_8[] =
{ 0.1012285362903762591525, 0.222381034453374470544, 0.313706645877887287338, 0.3626837833783619829652, 0.3626837833783619829652, 0.31370664587788728734, 0.222381034453374470544, 0.1012285362903762591525 };

// gkyl_gauss_ordinates[N] are ordinates for N-point Guassian
// integration
static const double* gkyl_gauss_ordinates[] = {
  0, // N=0 makes no sense,
  gkyl_gauss_ordinates_1,
  gkyl_gauss_ordinates_2,
  gkyl_gauss_ordinates_3,
  gkyl_gauss_ordinates_4,
  gkyl_gauss_ordinates_5,
  gkyl_gauss_ordinates_6,
  gkyl_gauss_ordinates_7,
  gkyl_gauss_ordinates_8
};

// gkyl_gauss_weights[N] are weights for N-point Guassian integration
static const double *gkyl_gauss_weights[] = {
  0, // N=0 makes no sense,
  gkyl_gauss_weights_1,
  gkyl_gauss_weights_2,
  gkyl_gauss_weights_3,
  gkyl_gauss_weights_4,
  gkyl_gauss_weights_5,
  gkyl_gauss_weights_6,
  gkyl_gauss_weights_7,
  gkyl_gauss_weights_8
};


// Lobatto quadrature

// Ordinates
static const double gkyl_gauss_lobatto_ordinates_2[] = 
{ -1.0,1.0 }; 
static const double gkyl_gauss_lobatto_ordinates_3[] = 
{ -1.0,0.0,1.0 }; 
static const double gkyl_gauss_lobatto_ordinates_4[] = 
{ -1.0,-0.4472135954999579,0.4472135954999579,1.0 }; 
static const double gkyl_gauss_lobatto_ordinates_5[] = 
{ -1.0,-0.6546536707079771,0.0,0.6546536707079771,1.0 }; 
static const double gkyl_gauss_lobatto_ordinates_6[] = 
{ -1.0,-0.7650553239294646,-0.285231516480645,0.285231516480645,0.7650553239294646,1.0 }; 
static const double gkyl_gauss_lobatto_ordinates_7[] = 
{ -1.0,-0.8302238962785669,-0.4688487934707142,0.0,0.4688487934707142,0.8302238962785669,1.0 }; 
static const double gkyl_gauss_lobatto_ordinates_8[] = 
{ -1.0,-0.8717401485096066,-0.5917001814331423,-0.2092992179024789,0.2092992179024789,0.5917001814331423,0.8717401485096066,1.0 };

// Weights
static const double gkyl_gauss_lobatto_weights_2[] = 
{ 1.0,1.0 }; 
static const double gkyl_gauss_lobatto_weights_3[] = 
{ 0.3333333333333333,1.333333333333333,0.3333333333333333 }; 
static const double gkyl_gauss_lobatto_weights_4[] = 
{ 0.1666666666666667,0.8333333333333334,0.8333333333333334,0.1666666666666667 }; 
static const double gkyl_gauss_lobatto_weights_5[] = 
{ 0.1,0.5444444444444444,0.7111111111111111,0.5444444444444444,0.1 }; 
static const double gkyl_gauss_lobatto_weights_6[] = 
{ 0.06666666666666667,0.378474956297847,0.5548583770354863,0.5548583770354863,0.378474956297847,0.06666666666666667 }; 
static const double gkyl_gauss_lobatto_weights_7[] = 
{ 0.04761904761904762,0.276826047361566,0.4317453812098626,0.4876190476190476,0.4317453812098626,0.276826047361566,0.04761904761904762 }; 
static const double gkyl_gauss_lobatto_weights_8[] = 
{ 0.03571428571428571,0.210704227143506,0.3411226924835044,0.4124587946587039,0.4124587946587039,0.3411226924835044,0.210704227143506,0.03571428571428571 };

// gkyl_gauss_lobatto_ordinates[N] are ordinates for N-point
// Guass-Lobatto integration
static const double* gkyl_gauss_lobatto_ordinates[] = {
  0, // N=0 makes no sense,
  0, // N=1 makes no sense,
  gkyl_gauss_lobatto_ordinates_2,
  gkyl_gauss_lobatto_ordinates_3,
  gkyl_gauss_lobatto_ordinates_4,
  gkyl_gauss_lobatto_ordinates_5,
  gkyl_gauss_lobatto_ordinates_6,
  gkyl_gauss_lobatto_ordinates_7,
  gkyl_gauss_lobatto_ordinates_8
};

// gkyl_gauss_lobatto_weights[N] are weights for N-point Guass-Lobatto
// integration
static const double *gkyl_gauss_lobatto_weights[] = {
  0, // N=0 makes no sense,
  0, // N=1 makes no sense,  
  gkyl_gauss_lobatto_weights_2,
  gkyl_gauss_lobatto_weights_3,
  gkyl_gauss_lobatto_weights_4,
  gkyl_gauss_lobatto_weights_5,
  gkyl_gauss_lobatto_weights_6,
  gkyl_gauss_lobatto_weights_7,
  gkyl_gauss_lobatto_weights_8
};

/**
 * Compute ordinates and weights for use in Gaussian quadrature.
 *
 * @param x1 Left coordinate of domain.
 * @param x2 Right coordinate of domain.
 * @param x On output, ordinates.
 * @param w On output, weights.
 * @param n Order of the quadrature.
 */
void gkyl_gauleg(double x1, double x2,  double x[], double w[], int n);
// ended inlining gkyl_gauss_quad_data.h 
// start inlining gkyl_job_pool.h 

// skipping file: gkyl_ref_count.h 

// forward declare for use in function pointers
struct gkyl_job_pool;

// Function pointer sig for function that does the actual work 
typedef void (*jp_work_func)(void *ctx);

// Function sig that adds work to the pool
typedef bool (*jp_add_work)(const struct gkyl_job_pool *jp, jp_work_func func, void *ctx);

// Function sig for "wait" function that halts till all jobs have finished
typedef void (*jp_wait)(const struct gkyl_job_pool *jp);

struct gkyl_job_pool {
  int pool_size; // number of worked in pool
  jp_add_work add_work; // function to add work to pool
  jp_wait wait; // function to wait for jobs to finish

  struct gkyl_ref_count ref_count; // reference count  
};

/**
 * Add work to job pool. Work may start immediately on calling this method.
 *
 * @param jp Job-pool object.
 * @param func Pointer to function that does the work
 * @param ctx Context object to add to func
 * @return True if work was added, false otherwise
 */
bool gkyl_job_pool_add_work(const struct gkyl_job_pool *jp, jp_work_func func, void *ctx);

/**
 * Wait till all jobs are completed
 *
 * @param jp Job-pool object.
 */
void gkyl_job_pool_wait(const struct gkyl_job_pool *jp);

/**
 * Acquire pointer to job-pool. Delete using the release()
 * method
 *
 * @param jp Job-pool object.
 */
struct gkyl_job_pool* gkyl_job_pool_acquire(const struct gkyl_job_pool *jp);

/**
 * Delete job-pool object
 *
 * @param jp Object to delete.
 */
void gkyl_job_pool_release(const struct gkyl_job_pool* jp);
// ended inlining gkyl_job_pool.h 
// start inlining gkyl_null_comm.h 

// skipping file: gkyl_comm.h 

// skipping file: gkyl_rect_decomp.h 

// input to create new MPI communicator
struct gkyl_null_comm_inp {
  bool use_gpu; // flag to use if this communicator is on GPUs  
  const struct gkyl_rect_decomp *decomp; // pre-computed decomposition
  bool sync_corners; // should we sync corners?
};

/**
 * Return a new "null" communicator, i.e. a communicator for a single
 * core calculation.
 *
 * @param inp Input struct to use for initialization
 * @return New communicator
 */
struct gkyl_comm *gkyl_null_comm_inew(const struct gkyl_null_comm_inp *inp);

// ended inlining gkyl_null_comm.h 
// skipping file: gkyl_range.h 
// skipping file: gkyl_rect_decomp.h 
// skipping file: gkyl_rect_grid.h 
// skipping file: gkyl_ref_count.h 
// start inlining gkyl_thread_pool.h 

// skipping file: gkyl_job_pool.h 

/**
 * Create a new thread-pool object
 *
 * @param nthreads Number of threads to create
 * @return Pointer to new job-pool object
 */
struct gkyl_job_pool* gkyl_thread_pool_new(int nthreads);

// ended inlining gkyl_thread_pool.h 
// skipping file: gkyl_util.h 
// skipping file: gkyl_vargm.h 
