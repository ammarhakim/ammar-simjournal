
#define MPACK_HAS_CONFIG 0

// start inlining mpack.c 
/**
 * The MIT License (MIT)
 * 
 * Copyright (c) 2015-2021 Nicholas Fraser and the MPack authors
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 */

/*
 * This is the MPack 1.1.1 amalgamation package.
 *
 * http://github.com/ludocode/mpack
 */

#define MPACK_INTERNAL 1
#define MPACK_EMIT_INLINE_DEFS 1

// start inlining mpack.h 
/**
 * The MIT License (MIT)
 * 
 * Copyright (c) 2015-2021 Nicholas Fraser and the MPack authors
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 */

/*
 * This is the MPack 1.1.1 amalgamation package.
 *
 * http://github.com/ludocode/mpack
 */

#ifndef MPACK_H
#define MPACK_H 1

#define MPACK_AMALGAMATED 1
#define MPACK_RELEASE_VERSION 1

#if defined(MPACK_HAS_CONFIG) && MPACK_HAS_CONFIG
//#include "mpack-config.h"
#endif


/* mpack/mpack-platform.h.h */

/**
 * @file
 *
 * Abstracts all platform-specific code from MPack and handles configuration
 * options.
 *
 * This verifies the configuration and sets defaults based on the platform,
 * contains implementations of standard C functions when libc is not available,
 * and provides wrappers to all library functions.
 *
 * Documentation for configuration options is available here:
 *
 *     https://ludocode.github.io/mpack/group__config.html
 */

#ifndef MPACK_PLATFORM_H
#define MPACK_PLATFORM_H 1



/**
 * @defgroup config Configuration Options
 *
 * Defines the MPack configuration options.
 *
 * Custom configuration of MPack is not usually necessary. In almost all
 * cases you can ignore this and use the defaults.
 *
 * If you do want to configure MPack, you can define the below options as part
 * of your build system or project settings. This will override the below
 * defaults.
 *
 * If you'd like to use a file for configuration instead, define
 * @ref MPACK_HAS_CONFIG to 1 in your build system or project settings.
 * This will cause MPack to include a file you create called @c mpack-config.h
 * in which you can define your configuration. This is useful if you need to
 * include specific headers (such as a custom allocator) in order to configure
 * MPack to use it.
 *
 * @warning The value of all configuration options must be the same in
 * all translation units of your project, as well as in the mpack source
 * itself. These configuration options affect the layout of structs, among
 * other things, which cannot be different in source files that are linked
 * together.
 *
 * @note MPack does not contain defaults for building inside the Linux kernel.
 * There is a <a href="https://github.com/ludocode/mpack-linux-kernel">
 * configuration file for the Linux kernel</a> that can be used instead.
 *
 * @{
 */



/*
 * Pre-include checks
 *
 * These need to come before the user's mpack-config.h because they might be
 * including headers in it.
 */

/** @cond */
#if defined(_MSC_VER) && _MSC_VER < 1800 && !defined(__cplusplus)
    #error "In Visual Studio 2012 and earlier, MPack must be compiled as C++. Enable the /Tp compiler flag."
#endif

#if defined(_WIN32) && MPACK_INTERNAL
    #define _CRT_SECURE_NO_WARNINGS 1
#endif

#ifndef __STDC_LIMIT_MACROS
    #define __STDC_LIMIT_MACROS 1
#endif
#ifndef __STDC_FORMAT_MACROS
    #define __STDC_FORMAT_MACROS 1
#endif
#ifndef __STDC_CONSTANT_MACROS
    #define __STDC_CONSTANT_MACROS 1
#endif
/** @endcond */



/**
 * @name File Configuration
 * @{
 */

/**
 * @def MPACK_HAS_CONFIG
 *
 * Causes MPack to include a file you create called @c mpack-config.h .
 *
 * The file is included before MPack sets any defaults for undefined
 * configuration options. You can use it to configure MPack.
 *
 * This is off by default.
 */
#if defined(MPACK_HAS_CONFIG)
    #if MPACK_HAS_CONFIG
        #include "mpack-config.h"
    #endif
#else
    #define MPACK_HAS_CONFIG 0
#endif

/**
 * @}
 */

// this needs to come first since some stuff depends on it
/** @cond */
#ifndef MPACK_NO_BUILTINS
    #define MPACK_NO_BUILTINS 0
#endif
/** @endcond */



/**
 * @name Features
 * @{
 */

/**
 * @def MPACK_READER
 *
 * Enables compilation of the base Tag Reader.
 */
#ifndef MPACK_READER
#define MPACK_READER 1
#endif

/**
 * @def MPACK_EXPECT
 *
 * Enables compilation of the static Expect API.
 */
#ifndef MPACK_EXPECT
#define MPACK_EXPECT 1
#endif

/**
 * @def MPACK_NODE
 *
 * Enables compilation of the dynamic Node API.
 */
#ifndef MPACK_NODE
#define MPACK_NODE 1
#endif

/**
 * @def MPACK_WRITER
 *
 * Enables compilation of the Writer.
 */
#ifndef MPACK_WRITER
#define MPACK_WRITER 1
#endif

/**
 * @def MPACK_BUILDER
 *
 * Enables compilation of the Builder.
 *
 * The Builder API provides additional functions to the Writer for
 * automatically determining the element count of compound elements so you do
 * not have to specify them up-front.
 *
 * This requires a @c malloc(). It is enabled by default if MPACK_WRITER is
 * enabled and MPACK_MALLOC is defined.
 *
 * @see mpack_build_map()
 * @see mpack_build_array()
 * @see mpack_complete_map()
 * @see mpack_complete_array()
 */
// This is defined furthur below after we've resolved whether we have malloc().

/**
 * @def MPACK_COMPATIBILITY
 *
 * Enables compatibility features for reading and writing older
 * versions of MessagePack.
 *
 * This is disabled by default. When disabled, the behaviour is equivalent to
 * using the default version, @ref mpack_version_current.
 *
 * Enable this if you need to interoperate with applications or data that do
 * not support the new (v5) MessagePack spec. See the section on v4
 * compatibility in @ref docs/protocol.md for more information.
 */
#ifndef MPACK_COMPATIBILITY
#define MPACK_COMPATIBILITY 0
#endif

/**
 * @def MPACK_EXTENSIONS
 *
 * Enables the use of extension types.
 *
 * This is disabled by default. Define it to 1 to enable it. If disabled,
 * functions to read and write extensions will not exist, and any occurrence of
 * extension types in parsed messages will flag @ref mpack_error_invalid.
 *
 * MPack discourages the use of extension types. See the section on extension
 * types in @ref docs/protocol.md for more information.
 */
#ifndef MPACK_EXTENSIONS
#define MPACK_EXTENSIONS 0
#endif

/**
 * @}
 */



// workarounds for Doxygen
#if defined(MPACK_DOXYGEN)
#if MPACK_DOXYGEN
// We give these their default values of 0 here even though they are defined to
// 1 in the doxyfile. Doxygen will show this as the value in the docs, even
// though it ignores it when parsing the rest of the source. This is what we
// want, since we want the documentation to show these defaults but still
// generate documentation for the functions they add when they're on.
#define MPACK_COMPATIBILITY 0
#define MPACK_EXTENSIONS 0
#endif
#endif



/**
 * @name Dependencies
 * @{
 */

/**
 * @def MPACK_CONFORMING
 *
 * Enables the inclusion of basic C headers to define standard types and
 * macros.
 *
 * This causes MPack to include headers required for conforming implementations
 * of C99 even in freestanding, in particular <stddef.h>, <stdint.h>,
 * <stdbool.h> and <limits.h>. It also includes <inttypes.h>; this is
 * technically not required for freestanding but MPack needs it to detect
 * integer limits.
 *
 * You can disable this if these headers are unavailable or if they do not
 * define the standard types and macros (for example inside the Linux kernel.)
 * If this is disabled, MPack will include no headers and will assume a 32-bit
 * int. You will probably also want to define @ref MPACK_HAS_CONFIG to 1 and
 * include your own headers in the config file. You must provide definitions
 * for standard types such as @c size_t, @c bool, @c int32_t and so on.
 *
 * @see <a href="https://en.cppreference.com/w/c/language/conformance">
 * cppreference.com documentation on Conformance</a>
 */
#ifndef MPACK_CONFORMING
    #define MPACK_CONFORMING 1
#endif

/**
 * @def MPACK_STDLIB
 *
 * Enables the use of the C stdlib.
 *
 * This allows the library to use basic functions like @c memcmp() and @c
 * strlen(), as well as @c malloc() for debugging and in allocation helpers.
 *
 * If this is disabled, allocation helper functions will not be defined, and
 * MPack will attempt to detect compiler intrinsics for functions like @c
 * memcmp() (assuming @ref MPACK_NO_BUILTINS is not set.) It will fallback to
 * its own (slow) implementations if it cannot use builtins. You may want to
 * define @ref MPACK_MEMCMP and friends if you disable this.
 *
 * @see MPACK_MEMCMP
 * @see MPACK_MEMCPY
 * @see MPACK_MEMMOVE
 * @see MPACK_MEMSET
 * @see MPACK_STRLEN
 * @see MPACK_MALLOC
 * @see MPACK_REALLOC
 * @see MPACK_FREE
 */
#ifndef MPACK_STDLIB
    #if !MPACK_CONFORMING
        // If we don't even have a proper <limits.h> we assume we won't have
        // malloc() either.
        #define MPACK_STDLIB 0
    #else
        #define MPACK_STDLIB 1
    #endif
#endif

/**
 * @def MPACK_STDIO
 *
 * Enables the use of C stdio. This adds helpers for easily
 * reading/writing C files and makes debugging easier.
 */
#ifndef MPACK_STDIO
    #if !MPACK_STDLIB || defined(__AVR__)
        #define MPACK_STDIO 0
    #else
        #define MPACK_STDIO 1
    #endif
#endif

/**
 * Whether the 'float' type and floating point operations are supported.
 *
 * If @ref MPACK_FLOAT is disabled, floats are read and written as @c uint32_t
 * instead. This way messages with floats do not result in errors and you can
 * still perform manual float parsing yourself.
 */
#ifndef MPACK_FLOAT
    #define MPACK_FLOAT 1
#endif

/**
 * Whether the 'double' type is supported. This requires support for 'float'.
 *
 * If @ref MPACK_DOUBLE is disabled, doubles are read and written as @c
 * uint32_t instead. This way messages with doubles do not result in errors and
 * you can still perform manual doubles parsing yourself.
 *
 * If @ref MPACK_FLOAT is enabled but @ref MPACK_DOUBLE is not, doubles can be
 * read as floats using the shortening conversion functions, e.g. @ref
 * mpack_expect_float() or @ref mpack_node_float().
 */
#ifndef MPACK_DOUBLE
    #if !MPACK_FLOAT || defined(__AVR__)
        // AVR supports only float, not double.
        #define MPACK_DOUBLE 0
    #else
        #define MPACK_DOUBLE 1
    #endif
#endif

/**
 * @}
 */



/**
 * @name Allocation Functions
 * @{
 */

/**
 * @def MPACK_MALLOC
 *
 * Defines the memory allocation function used by MPack. This is used by
 * helpers for automatically allocating data the correct size, and for
 * debugging functions. If this macro is undefined, the allocation helpers
 * will not be compiled.
 *
 * Set this to use a custom @c malloc() function. Your function must have the
 * signature:
 *
 * @code
 * void* malloc(size_t size);
 * @endcode
 *
 * The default is @c malloc() if @ref MPACK_STDLIB is enabled.
 */
/**
 * @def MPACK_FREE
 *
 * Defines the memory free function used by MPack. This is used by helpers
 * for automatically allocating data the correct size. If this macro is
 * undefined, the allocation helpers will not be compiled.
 *
 * Set this to use a custom @c free() function. Your function must have the
 * signature:
 *
 * @code
 * void free(void* p);
 * @endcode
 *
 * The default is @c free() if @ref MPACK_MALLOC has not been customized and
 * @ref MPACK_STDLIB is enabled.
 */
/**
 * @def MPACK_REALLOC
 *
 * Defines the realloc function used by MPack. It is used by growable
 * buffers to resize more efficiently.
 *
 * The default is @c realloc() if @ref MPACK_MALLOC has not been customized and
 * @ref MPACK_STDLIB is enabled.
 *
 * Set this to use a custom @c realloc() function. Your function must have the
 * signature:
 *
 * @code
 * void* realloc(void* p, size_t new_size);
 * @endcode
 *
 * This is optional, even when @ref MPACK_MALLOC is used. If @ref MPACK_MALLOC is
 * set and @ref MPACK_REALLOC is not, @ref MPACK_MALLOC is used with a simple copy
 * to grow buffers.
 */

#if defined(MPACK_MALLOC) && !defined(MPACK_FREE)
    #error "MPACK_MALLOC requires MPACK_FREE."
#endif
#if !defined(MPACK_MALLOC) && defined(MPACK_FREE)
    #error "MPACK_FREE requires MPACK_MALLOC."
#endif

// These were never configurable in lowercase but we check anyway.
#ifdef mpack_malloc
    #error "Define MPACK_MALLOC, not mpack_malloc."
#endif
#ifdef mpack_realloc
    #error "Define MPACK_REALLOC, not mpack_realloc."
#endif
#ifdef mpack_free
    #error "Define MPACK_FREE, not mpack_free."
#endif

// We don't use calloc() at all.
#ifdef MPACK_CALLOC
    #error "Don't define MPACK_CALLOC. MPack does not use calloc()."
#endif
#ifdef mpack_calloc
    #error "Don't define mpack_calloc. MPack does not use calloc()."
#endif

// Use defaults in stdlib if we have them. Without it we don't use malloc.
#if defined(MPACK_STDLIB)
    #if MPACK_STDLIB && !defined(MPACK_MALLOC)
        #define MPACK_MALLOC malloc
        #define MPACK_REALLOC realloc
        #define MPACK_FREE free
    #endif
#endif

/**
 * @}
 */



// This needs to be defined after we've decided whether we have malloc().
#ifndef MPACK_BUILDER
    #if defined(MPACK_MALLOC) && MPACK_WRITER
        #define MPACK_BUILDER 1
    #else
        #define MPACK_BUILDER 0
    #endif
#endif



/**
 * @name System Functions
 * @{
 */

/**
 * @def MPACK_MEMCMP
 *
 * The function MPack will use to provide @c memcmp().
 *
 * Set this to use a custom @c memcmp() function. Your function must have the
 * signature:
 *
 * @code
 * int memcmp(const void* left, const void* right, size_t count);
 * @endcode
 */
/**
 * @def MPACK_MEMCPY
 *
 * The function MPack will use to provide @c memcpy().
 *
 * Set this to use a custom @c memcpy() function. Your function must have the
 * signature:
 *
 * @code
 * void* memcpy(void* restrict dest, const void* restrict src, size_t count);
 * @endcode
 */
/**
 * @def MPACK_MEMMOVE
 *
 * The function MPack will use to provide @c memmove().
 *
 * Set this to use a custom @c memmove() function. Your function must have the
 * signature:
 *
 * @code
 * void* memmove(void* dest, const void* src, size_t count);
 * @endcode
 */
/**
 * @def MPACK_MEMSET
 *
 * The function MPack will use to provide @c memset().
 *
 * Set this to use a custom @c memset() function. Your function must have the
 * signature:
 *
 * @code
 * void* memset(void* p, int c, size_t count);
 * @endcode
 */
/**
 * @def MPACK_STRLEN
 *
 * The function MPack will use to provide @c strlen().
 *
 * Set this to use a custom @c strlen() function. Your function must have the
 * signature:
 *
 * @code
 * size_t strlen(const char* str);
 * @endcode
 */

// These were briefly configurable in lowercase in an unreleased version. Just
// to make sure no one is doing this, we make sure these aren't already defined.
#ifdef mpack_memcmp
    #error "Define MPACK_MEMCMP, not mpack_memcmp."
#endif
#ifdef mpack_memcpy
    #error "Define MPACK_MEMCPY, not mpack_memcpy."
#endif
#ifdef mpack_memmove
    #error "Define MPACK_MEMMOVE, not mpack_memmove."
#endif
#ifdef mpack_memset
    #error "Define MPACK_MEMSET, not mpack_memset."
#endif
#ifdef mpack_strlen
    #error "Define MPACK_STRLEN, not mpack_strlen."
#endif

// If the standard library is available, we prefer to use its functions.
#if MPACK_STDLIB
    #ifndef MPACK_MEMCMP
        #define MPACK_MEMCMP memcmp
    #endif
    #ifndef MPACK_MEMCPY
        #define MPACK_MEMCPY memcpy
    #endif
    #ifndef MPACK_MEMMOVE
        #define MPACK_MEMMOVE memmove
    #endif
    #ifndef MPACK_MEMSET
        #define MPACK_MEMSET memset
    #endif
    #ifndef MPACK_STRLEN
        #define MPACK_STRLEN strlen
    #endif
#endif

#if !MPACK_NO_BUILTINS
    #ifdef __has_builtin
        #if !defined(MPACK_MEMCMP) && __has_builtin(__builtin_memcmp)
            #define MPACK_MEMCMP __builtin_memcmp
        #endif
        #if !defined(MPACK_MEMCPY) && __has_builtin(__builtin_memcpy)
            #define MPACK_MEMCPY __builtin_memcpy
        #endif
        #if !defined(MPACK_MEMMOVE) && __has_builtin(__builtin_memmove)
            #define MPACK_MEMMOVE __builtin_memmove
        #endif
        #if !defined(MPACK_MEMSET) && __has_builtin(__builtin_memset)
            #define MPACK_MEMSET __builtin_memset
        #endif
        #if !defined(MPACK_STRLEN) && __has_builtin(__builtin_strlen)
            #define MPACK_STRLEN __builtin_strlen
        #endif
    #elif defined(__GNUC__)
        #ifndef MPACK_MEMCMP
            #define MPACK_MEMCMP __builtin_memcmp
        #endif
        #ifndef MPACK_MEMCPY
            #define MPACK_MEMCPY __builtin_memcpy
        #endif
        // There's not always a builtin memmove for GCC. If we can't test for
        // it with __has_builtin above, we don't use it. It's been around for
        // much longer under clang, but then so has __has_builtin, so we let
        // the block above handle it.
        #ifndef MPACK_MEMSET
            #define MPACK_MEMSET __builtin_memset
        #endif
        #ifndef MPACK_STRLEN
            #define MPACK_STRLEN __builtin_strlen
        #endif
    #endif
#endif

/**
 * @}
 */



/**
 * @name Debugging Options
 * @{
 */

/**
 * @def MPACK_DEBUG
 *
 * Enables debug features. You may want to wrap this around your
 * own debug preprocs. By default, this is enabled if @c DEBUG or @c _DEBUG
 * are defined. (@c NDEBUG is not used since it is allowed to have
 * different values in different translation units.)
 */
#if !defined(MPACK_DEBUG)
    #if defined(DEBUG) || defined(_DEBUG)
        #define MPACK_DEBUG 1
    #else
        #define MPACK_DEBUG 0
    #endif
#endif

/**
 * @def MPACK_STRINGS
 *
 * Enables descriptive error and type strings.
 *
 * This can be turned off (by defining it to 0) to maximize space savings
 * on embedded devices. If this is disabled, string functions such as
 * mpack_error_to_string() and mpack_type_to_string() return an empty string.
 */
#ifndef MPACK_STRINGS
    #ifdef __AVR__
        #define MPACK_STRINGS 0
    #else
        #define MPACK_STRINGS 1
    #endif
#endif

/**
 * Set this to 1 to implement a custom @ref mpack_assert_fail() function.
 * See the documentation on @ref mpack_assert_fail() for details.
 *
 * Asserts are only used when @ref MPACK_DEBUG is enabled, and can be
 * triggered by bugs in MPack or bugs due to incorrect usage of MPack.
 */
#ifndef MPACK_CUSTOM_ASSERT
#define MPACK_CUSTOM_ASSERT 0
#endif

/**
 * @def MPACK_READ_TRACKING
 *
 * Enables compound type size tracking for readers. This ensures that the
 * correct number of elements or bytes are read from a compound type.
 *
 * This is enabled by default in debug builds (provided a @c malloc() is
 * available.)
 */
#if !defined(MPACK_READ_TRACKING)
    #if MPACK_DEBUG && MPACK_READER && defined(MPACK_MALLOC)
        #define MPACK_READ_TRACKING 1
    #else
        #define MPACK_READ_TRACKING 0
    #endif
#endif
#if MPACK_READ_TRACKING && !MPACK_READER
    #error "MPACK_READ_TRACKING requires MPACK_READER."
#endif

/**
 * @def MPACK_WRITE_TRACKING
 *
 * Enables compound type size tracking for writers. This ensures that the
 * correct number of elements or bytes are written in a compound type.
 *
 * Note that without write tracking enabled, it is possible for buggy code
 * to emit invalid MessagePack without flagging an error by writing the wrong
 * number of elements or bytes in a compound type. With tracking enabled,
 * MPack will catch such errors and break on the offending line of code.
 *
 * This is enabled by default in debug builds (provided a @c malloc() is
 * available.)
 */
#if !defined(MPACK_WRITE_TRACKING)
    #if MPACK_DEBUG && MPACK_WRITER && defined(MPACK_MALLOC)
        #define MPACK_WRITE_TRACKING 1
    #else
        #define MPACK_WRITE_TRACKING 0
    #endif
#endif
#if MPACK_WRITE_TRACKING && !MPACK_WRITER
    #error "MPACK_WRITE_TRACKING requires MPACK_WRITER."
#endif

/**
 * @}
 */




/**
 * @name Miscellaneous Options
 * @{
 */

/**
 * Whether to optimize for size or speed.
 *
 * Optimizing for size simplifies some parsing and encoding algorithms
 * at the expense of speed and saves a few kilobytes of space in the
 * resulting executable.
 *
 * This automatically detects -Os with GCC/Clang. Unfortunately there
 * doesn't seem to be a macro defined for /Os under MSVC.
 */
#ifndef MPACK_OPTIMIZE_FOR_SIZE
    #ifdef __OPTIMIZE_SIZE__
        #define MPACK_OPTIMIZE_FOR_SIZE 1
    #else
        #define MPACK_OPTIMIZE_FOR_SIZE 0
    #endif
#endif

/**
 * Stack space in bytes to use when initializing a reader or writer
 * with a stack-allocated buffer.
 *
 * @warning Make sure you have sufficient stack space. Some libc use relatively
 * small stacks even on desktop platforms, e.g. musl.
 */
#ifndef MPACK_STACK_SIZE
#define MPACK_STACK_SIZE 4096
#endif

/**
 * Buffer size to use for allocated buffers (such as for a file writer.)
 *
 * Starting with a single page and growing as needed seems to
 * provide the best performance with minimal memory waste.
 * Increasing this does not improve performance even when writing
 * huge messages.
 */
#ifndef MPACK_BUFFER_SIZE
#define MPACK_BUFFER_SIZE 4096
#endif

/**
 * Minimum size for paged allocations in bytes.
 *
 * This is the value used by default for MPACK_NODE_PAGE_SIZE and
 * MPACK_BUILDER_PAGE_SIZE.
 */
#ifndef MPACK_PAGE_SIZE
#define MPACK_PAGE_SIZE 4096
#endif

/**
 * Minimum size of an allocated node page in bytes.
 *
 * The children for a given compound element must be contiguous, so
 * larger pages than this may be allocated as needed. (Safety checks
 * exist to prevent malicious data from causing too large allocations.)
 *
 * See @ref mpack_node_data_t for the size of nodes.
 *
 * Using as many nodes fit in one memory page seems to provide the
 * best performance, and has very little waste when parsing small
 * messages.
 */
#ifndef MPACK_NODE_PAGE_SIZE
#define MPACK_NODE_PAGE_SIZE MPACK_PAGE_SIZE
#endif

/**
 * Minimum size of an allocated builder page in bytes.
 *
 * Builder writes are deferred to the allocated builder buffer which is
 * composed of a list of buffer pages. This defines the size of those pages.
 */
#ifndef MPACK_BUILDER_PAGE_SIZE
#define MPACK_BUILDER_PAGE_SIZE MPACK_PAGE_SIZE
#endif

/**
 * @def MPACK_BUILDER_INTERNAL_STORAGE
 *
 * Enables a small amount of internal storage within the writer to avoid some
 * allocations when using builders.
 *
 * This is disabled by default. Enable it to potentially improve performance at
 * the expense of a larger writer.
 *
 * @see MPACK_BUILDER_INTERNAL_STORAGE_SIZE to configure its size.
 */
#ifndef MPACK_BUILDER_INTERNAL_STORAGE
#define MPACK_BUILDER_INTERNAL_STORAGE 0
#endif

/**
 * Amount of space reserved inside @ref mpack_writer_t for the Builders. This
 * can allow small messages to be built with the Builder API without incurring
 * an allocation.
 *
 * Builder metadata is placed in this space in addition to the literal
 * MessagePack data. It needs to be big enough to be useful, but not so big as
 * to overflow the stack. If more space is needed, pages are allocated.
 *
 * This is only used if MPACK_BUILDER_INTERNAL_STORAGE is enabled.
 *
 * @see MPACK_BUILDER_PAGE_SIZE
 * @see MPACK_BUILDER_INTERNAL_STORAGE
 *
 * @warning Writers are typically placed on the stack so make sure you have
 * sufficient stack space. Some libc use relatively small stacks even on
 * desktop platforms, e.g. musl.
 */
#ifndef MPACK_BUILDER_INTERNAL_STORAGE_SIZE
#define MPACK_BUILDER_INTERNAL_STORAGE_SIZE 256
#endif

/**
 * The initial depth for the node parser. When MPACK_MALLOC is available,
 * the node parser has no practical depth limit, and it is not recursive
 * so there is no risk of overflowing the call stack.
 */
#ifndef MPACK_NODE_INITIAL_DEPTH
#define MPACK_NODE_INITIAL_DEPTH 8
#endif

/**
 * The maximum depth for the node parser if @ref MPACK_MALLOC is not available.
 */
#ifndef MPACK_NODE_MAX_DEPTH_WITHOUT_MALLOC
#define MPACK_NODE_MAX_DEPTH_WITHOUT_MALLOC 32
#endif

/**
 * @def MPACK_NO_BUILTINS
 *
 * Whether to disable compiler intrinsics and other built-in functions.
 *
 * If this is enabled, MPack won't use `__attribute__`, `__declspec`, any
 * function starting with `__builtin`, or pretty much anything else that isn't
 * standard C.
 */
#if defined(MPACK_DOXYGEN)
#if MPACK_DOXYGEN
    #define MPACK_NO_BUILTINS 0
#endif
#endif

/**
 * @}
 */



#if MPACK_DEBUG
/**
 * @name Debug Functions
 * @{
 */
/**
 * Implement this and define @ref MPACK_CUSTOM_ASSERT to use a custom
 * assertion function.
 *
 * This function should not return. If it does, MPack will @c abort().
 *
 * If you use C++, make sure you include @c mpack.h where you define
 * this to get the correct linkage (or define it <code>extern "C"</code>.)
 *
 * Asserts are only used when @ref MPACK_DEBUG is enabled, and can be
 * triggered by bugs in MPack or bugs due to incorrect usage of MPack.
 */
void mpack_assert_fail(const char* message);
/**
 * @}
 */
#endif



// The rest of this file shouldn't show up in Doxygen docs.
/** @cond */



/*
 * All remaining pseudo-configuration options that have not yet been set must
 * be defined here in order to support -Wundef.
 *
 * These aren't real configuration options; they are implementation details of
 * MPack.
 */
#ifndef MPACK_CUSTOM_BREAK
#define MPACK_CUSTOM_BREAK 0
#endif
#ifndef MPACK_EMIT_INLINE_DEFS
#define MPACK_EMIT_INLINE_DEFS 0
#endif
#ifndef MPACK_AMALGAMATED
#define MPACK_AMALGAMATED 0
#endif
#ifndef MPACK_RELEASE_VERSION
#define MPACK_RELEASE_VERSION 0
#endif
#ifndef MPACK_INTERNAL
#define MPACK_INTERNAL 0
#endif



/* System headers (based on configuration) */

#if MPACK_CONFORMING
    #include <stddef.h>
    #include <stdint.h>
    #include <stdbool.h>
    #include <inttypes.h>
    #include <limits.h>
#endif

#if MPACK_STDLIB
    #include <string.h>
    #include <stdlib.h>
#endif

#if MPACK_STDIO
    #include <stdio.h>
    #include <errno.h>
    #if MPACK_DEBUG
        #include <stdarg.h>
    #endif
#endif



/*
 * Integer Constants and Limits
 */

#if MPACK_CONFORMING
    #define MPACK_INT64_C INT64_C
    #define MPACK_UINT64_C UINT64_C

    #define MPACK_INT8_MIN INT8_MIN
    #define MPACK_INT16_MIN INT16_MIN
    #define MPACK_INT32_MIN INT32_MIN
    #define MPACK_INT64_MIN INT64_MIN
    #define MPACK_INT_MIN INT_MIN

    #define MPACK_INT8_MAX INT8_MAX
    #define MPACK_INT16_MAX INT16_MAX
    #define MPACK_INT32_MAX INT32_MAX
    #define MPACK_INT64_MAX INT64_MAX
    #define MPACK_INT_MAX INT_MAX

    #define MPACK_UINT8_MAX UINT8_MAX
    #define MPACK_UINT16_MAX UINT16_MAX
    #define MPACK_UINT32_MAX UINT32_MAX
    #define MPACK_UINT64_MAX UINT64_MAX
    #define MPACK_UINT_MAX UINT_MAX
#else
    // For a non-conforming implementation we assume int is 32 bits.

    #define MPACK_INT64_C(x) ((int64_t)(x##LL))
    #define MPACK_UINT64_C(x) ((uint64_t)(x##LLU))

    #define MPACK_INT8_MIN ((int8_t)(0x80))
    #define MPACK_INT16_MIN ((int16_t)(0x8000))
    #define MPACK_INT32_MIN ((int32_t)(0x80000000))
    #define MPACK_INT64_MIN MPACK_INT64_C(0x8000000000000000)
    #define MPACK_INT_MIN MPACK_INT32_MIN

    #define MPACK_INT8_MAX ((int8_t)(0x7f))
    #define MPACK_INT16_MAX ((int16_t)(0x7fff))
    #define MPACK_INT32_MAX ((int32_t)(0x7fffffff))
    #define MPACK_INT64_MAX MPACK_INT64_C(0x7fffffffffffffff)
    #define MPACK_INT_MAX MPACK_INT32_MAX

    #define MPACK_UINT8_MAX ((uint8_t)(0xffu))
    #define MPACK_UINT16_MAX ((uint16_t)(0xffffu))
    #define MPACK_UINT32_MAX ((uint32_t)(0xffffffffu))
    #define MPACK_UINT64_MAX MPACK_UINT64_C(0xffffffffffffffff)
    #define MPACK_UINT_MAX MPACK_UINT32_MAX
#endif



/*
 * Floating point support
 */

#if MPACK_DOUBLE && !MPACK_FLOAT
    #error "MPACK_DOUBLE requires MPACK_FLOAT."
#endif

// If we don't have support for float or double, we poison the identifiers to
// make sure we don't define anything related to them.
#if MPACK_INTERNAL
    #ifdef __GNUC__
        #if !MPACK_FLOAT
            #pragma GCC poison float
        #endif
        #if !MPACK_DOUBLE
            #pragma GCC poison double
        #endif
    #endif
#endif



/*
 * extern C
 */

#ifdef __cplusplus
    #define MPACK_EXTERN_C_BEGIN extern "C" {
    #define MPACK_EXTERN_C_END   }
#else
    #define MPACK_EXTERN_C_BEGIN /*nothing*/
    #define MPACK_EXTERN_C_END   /*nothing*/
#endif



/*
 * Warnings
 */

#if defined(__GNUC__)
    // Diagnostic push is not supported before GCC 4.6.
    #if defined(__clang__) || __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
        #define MPACK_SILENCE_WARNINGS_PUSH _Pragma ("GCC diagnostic push")
        #define MPACK_SILENCE_WARNINGS_POP _Pragma ("GCC diagnostic pop")
    #endif
#elif defined(_MSC_VER)
    // To support VS2017 and earlier we need to use __pragma and not _Pragma
    #define MPACK_SILENCE_WARNINGS_PUSH __pragma(warning(push))
    #define MPACK_SILENCE_WARNINGS_POP __pragma(warning(pop))
#endif

#if defined(_MSC_VER)
    // These are a bunch of mostly useless warnings emitted under MSVC /W4,
    // some as a result of the expansion of macros.
    #define MPACK_SILENCE_WARNINGS_MSVC_W4 \
            __pragma(warning(disable:4996)) /* _CRT_SECURE_NO_WARNINGS */ \
            __pragma(warning(disable:4127)) /* comparison is constant */ \
            __pragma(warning(disable:4702)) /* unreachable code */ \
            __pragma(warning(disable:4310)) /* cast truncates constant value */
#else
    #define MPACK_SILENCE_WARNINGS_MSVC_W4 /*nothing*/
#endif

/* GCC versions before 5.1 warn about defining a C99 non-static inline function
 * before declaring it (see issue #20). */
#if defined(__GNUC__) && !defined(__clang__)
    #if __GNUC__ < 5 || (__GNUC__ == 5 && __GNUC_MINOR__ < 1)
        #ifdef __cplusplus
            #define MPACK_SILENCE_WARNINGS_MISSING_PROTOTYPES \
                _Pragma ("GCC diagnostic ignored \"-Wmissing-declarations\"")
        #else
            #define MPACK_SILENCE_WARNINGS_MISSING_PROTOTYPES \
                _Pragma ("GCC diagnostic ignored \"-Wmissing-prototypes\"")
        #endif
    #endif
#endif
#ifndef MPACK_SILENCE_WARNINGS_MISSING_PROTOTYPES
    #define MPACK_SILENCE_WARNINGS_MISSING_PROTOTYPES /*nothing*/
#endif

/* GCC versions before 4.8 warn about shadowing a function with a variable that
 * isn't a function or function pointer (like "index"). */
#if defined(__GNUC__) && !defined(__clang__)
    #if __GNUC__ == 4 && __GNUC_MINOR__ < 8
        #define MPACK_SILENCE_WARNINGS_SHADOW \
            _Pragma ("GCC diagnostic ignored \"-Wshadow\"")
    #endif
#endif
#ifndef MPACK_SILENCE_WARNINGS_SHADOW
    #define MPACK_SILENCE_WARNINGS_SHADOW /*nothing*/
#endif

// On platforms with small size_t (e.g. AVR) we get type limits warnings where
// we compare a size_t to e.g. MPACK_UINT32_MAX.
#ifdef __AVR__
    #define MPACK_SILENCE_WARNINGS_TYPE_LIMITS \
        _Pragma ("GCC diagnostic ignored \"-Wtype-limits\"")
#else
    #define MPACK_SILENCE_WARNINGS_TYPE_LIMITS /*nothing*/
#endif

// MPack uses declarations after statements. This silences warnings about it
// (e.g. when using MPack in a Linux kernel module.)
#if defined(__GNUC__) && !defined(__cplusplus)
    #define MPACK_SILENCE_WARNINGS_DECLARATION_AFTER_STATEMENT \
        _Pragma ("GCC diagnostic ignored \"-Wdeclaration-after-statement\"")
#else
    #define MPACK_SILENCE_WARNINGS_DECLARATION_AFTER_STATEMENT /*nothing*/
#endif

#ifdef MPACK_SILENCE_WARNINGS_PUSH
    // We only silence warnings if push/pop is supported, that way we aren't
    // silencing warnings in code that uses MPack. If your compiler doesn't
    // support push/pop silencing of warnings, you'll have to turn off
    // conflicting warnings manually.

    #define MPACK_SILENCE_WARNINGS_BEGIN \
        MPACK_SILENCE_WARNINGS_PUSH \
        MPACK_SILENCE_WARNINGS_MSVC_W4 \
        MPACK_SILENCE_WARNINGS_MISSING_PROTOTYPES \
        MPACK_SILENCE_WARNINGS_SHADOW \
        MPACK_SILENCE_WARNINGS_TYPE_LIMITS \
        MPACK_SILENCE_WARNINGS_DECLARATION_AFTER_STATEMENT

    #define MPACK_SILENCE_WARNINGS_END \
        MPACK_SILENCE_WARNINGS_POP
#else
    #define MPACK_SILENCE_WARNINGS_BEGIN /*nothing*/
    #define MPACK_SILENCE_WARNINGS_END /*nothing*/
#endif

MPACK_SILENCE_WARNINGS_BEGIN
MPACK_EXTERN_C_BEGIN



/* Miscellaneous helper macros */

#define MPACK_UNUSED(var) ((void)(var))

#define MPACK_STRINGIFY_IMPL(arg) #arg
#define MPACK_STRINGIFY(arg) MPACK_STRINGIFY_IMPL(arg)

// This is a workaround for MSVC's incorrect expansion of __VA_ARGS__. It
// treats __VA_ARGS__ as a single preprocessor token when passed in the
// argument list of another macro unless we use an outer wrapper to expand it
// lexically first. (For safety/consistency we use this in all variadic macros
// that don't ignore the variadic arguments regardless of whether __VA_ARGS__
// is passed to another macro.)
//     https://stackoverflow.com/a/32400131
#define MPACK_EXPAND(x) x

// Extracts the first argument of a variadic macro list, where there might only
// be one argument.
#define MPACK_EXTRACT_ARG0_IMPL(first, ...) first
#define MPACK_EXTRACT_ARG0(...) MPACK_EXPAND(MPACK_EXTRACT_ARG0_IMPL( __VA_ARGS__ , ignored))

// Stringifies the first argument of a variadic macro list, where there might
// only be one argument.
#define MPACK_STRINGIFY_ARG0_impl(first, ...) #first
#define MPACK_STRINGIFY_ARG0(...) MPACK_EXPAND(MPACK_STRINGIFY_ARG0_impl( __VA_ARGS__ , ignored))



/*
 * Definition of inline macros.
 *
 * MPack does not use static inline in header files; only one non-inline definition
 * of each function should exist in the final build. This can reduce the binary size
 * in cases where the compiler cannot or chooses not to inline a function.
 * The addresses of functions should also compare equal across translation units
 * regardless of whether they are declared inline.
 *
 * The above requirements mean that the declaration and definition of non-trivial
 * inline functions must be separated so that the definitions will only
 * appear when necessary. In addition, three different linkage models need
 * to be supported:
 *
 *  - The C99 model, where a standalone function is emitted only if there is any
 *    `extern inline` or non-`inline` declaration (including the definition itself)
 *
 *  - The GNU model, where an `inline` definition emits a standalone function and an
 *    `extern inline` definition does not, regardless of other declarations
 *
 *  - The C++ model, where `inline` emits a standalone function with special
 *    (COMDAT) linkage
 *
 * The macros below wrap up everything above. All inline functions defined in header
 * files have a single non-inline definition emitted in the compilation of
 * mpack-platform.c. All inline declarations and definitions use the same MPACK_INLINE
 * specification to simplify the rules on when standalone functions are emitted.
 * Inline functions in source files are defined MPACK_STATIC_INLINE.
 *
 * Additional reading:
 *     http://www.greenend.org.uk/rjk/tech/inline.html
 */

#if defined(__cplusplus)
    // C++ rules
    // The linker will need COMDAT support to link C++ object files,
    // so we don't need to worry about emitting definitions from C++
    // translation units. If mpack-platform.c (or the amalgamation)
    // is compiled as C, its definition will be used, otherwise a
    // C++ definition will be used, and no other C files will emit
    // a definition.
    #define MPACK_INLINE inline

#elif defined(_MSC_VER)
    // MSVC 2013 always uses COMDAT linkage, but it doesn't treat 'inline' as a
    // keyword in C99 mode. (This appears to be fixed in a later version of
    // MSVC but we don't bother detecting it.)
    #define MPACK_INLINE __inline
    #define MPACK_STATIC_INLINE static __inline

#elif defined(__GNUC__) && (defined(__GNUC_GNU_INLINE__) || \
        (!defined(__GNUC_STDC_INLINE__) && !defined(__GNUC_GNU_INLINE__)))
    // GNU rules
    #if MPACK_EMIT_INLINE_DEFS
        #define MPACK_INLINE inline
    #else
        #define MPACK_INLINE extern inline
    #endif

#elif defined(__TINYC__)
    // tcc ignores the inline keyword, so we have to use static inline. We
    // issue a warning to make sure you are aware. You can define the below
    // macro to disable the warning. Hopefully this will be fixed soon:
    //     https://lists.nongnu.org/archive/html/tinycc-devel/2019-06/msg00000.html
    #ifndef MPACK_DISABLE_TINYC_INLINE_WARNING
        #warning "Single-definition inline is not supported by tcc. All inlines will be static. Define MPACK_DISABLE_TINYC_INLINE_WARNING to disable this warning."
    #endif
    #define MPACK_INLINE static inline

#else
    // C99 rules
    #if MPACK_EMIT_INLINE_DEFS
        #define MPACK_INLINE extern inline
    #else
        #define MPACK_INLINE inline
    #endif
#endif

#ifndef MPACK_STATIC_INLINE
#define MPACK_STATIC_INLINE static inline
#endif

#ifdef MPACK_OPTIMIZE_FOR_SPEED
    #error "You should define MPACK_OPTIMIZE_FOR_SIZE, not MPACK_OPTIMIZE_FOR_SPEED."
#endif



/*
 * Prevent inlining
 *
 * When a function is only used once, it is almost always inlined
 * automatically. This can cause poor instruction cache usage because a
 * function that should rarely be called (such as buffer exhaustion handling)
 * will get inlined into the middle of a hot code path.
 */

#if !MPACK_NO_BUILTINS
    #if defined(_MSC_VER)
        #define MPACK_NOINLINE __declspec(noinline)
    #elif defined(__GNUC__) || defined(__clang__)
        #define MPACK_NOINLINE __attribute__((__noinline__))
    #endif
#endif
#ifndef MPACK_NOINLINE
    #define MPACK_NOINLINE /* nothing */
#endif



/* restrict */

// We prefer the builtins even though e.g. MSVC's __restrict may not have
// exactly the same behaviour as the proper C99 restrict keyword because the
// builtins work in C++, so using the same keyword in both C and C++ prevents
// any incompatibilities when using MPack compiled as C in C++ code.
#if !MPACK_NO_BUILTINS
    #if defined(__GNUC__)
        #define MPACK_RESTRICT __restrict__
    #elif defined(_MSC_VER)
        #define MPACK_RESTRICT __restrict
    #endif
#endif

#ifndef MPACK_RESTRICT
    #ifdef __cplusplus
        #define MPACK_RESTRICT /* nothing, unavailable in C++ */
    #endif
#endif

#ifndef MPACK_RESTRICT
    #ifdef _MSC_VER
        // MSVC 2015 apparently doesn't properly support the restrict keyword
        // in C. We're using builtins above which do work on 2015, but when
        // MPACK_NO_BUILTINS is enabled we can't use it.
        #if _MSC_VER < 1910
            #define MPACK_RESTRICT /*nothing*/
        #endif
    #endif
#endif

#ifndef MPACK_RESTRICT
    #define MPACK_RESTRICT restrict /* required in C99 */
#endif



/* Some compiler-specific keywords and builtins */

#if !MPACK_NO_BUILTINS
    #if defined(__GNUC__) || defined(__clang__)
        #define MPACK_UNREACHABLE __builtin_unreachable()
        #define MPACK_NORETURN(fn) fn __attribute__((__noreturn__))
    #elif defined(_MSC_VER)
        #define MPACK_UNREACHABLE __assume(0)
        #define MPACK_NORETURN(fn) __declspec(noreturn) fn
    #endif
#endif

#ifndef MPACK_UNREACHABLE
#define MPACK_UNREACHABLE ((void)0)
#endif
#ifndef MPACK_NORETURN
#define MPACK_NORETURN(fn) fn
#endif



/*
 * Likely/unlikely
 *
 * These should only really be used when a branch is taken (or not taken) less
 * than 1/1000th of the time. Buffer flush checks when writing very small
 * elements are a good example.
 */

#if !MPACK_NO_BUILTINS
    #if defined(__GNUC__) || defined(__clang__)
        #ifndef MPACK_LIKELY
            #define MPACK_LIKELY(x) __builtin_expect((x),1)
        #endif
        #ifndef MPACK_UNLIKELY
            #define MPACK_UNLIKELY(x) __builtin_expect((x),0)
        #endif
    #endif
#endif

#ifndef MPACK_LIKELY
    #define MPACK_LIKELY(x) (x)
#endif
#ifndef MPACK_UNLIKELY
    #define MPACK_UNLIKELY(x) (x)
#endif



/* alignof */

#ifndef MPACK_ALIGNOF
    #if defined(__STDC_VERSION__)
        #if __STDC_VERSION__ >= 201112L
            #define MPACK_ALIGNOF(T) (_Alignof(T))
        #endif
    #endif
#endif

#ifndef MPACK_ALIGNOF
    #if defined(__cplusplus)
        #if __cplusplus >= 201103L
            #define MPACK_ALIGNOF(T) (alignof(T))
        #endif
    #endif
#endif

#ifndef MPACK_ALIGNOF
    #if defined(__GNUC__) && !defined(MPACK_NO_BUILTINS)
        #if defined(__clang__) || __GNUC__ >= 4
            #define MPACK_ALIGNOF(T) (__alignof__(T))
        #endif
    #endif
#endif

#ifndef MPACK_ALIGNOF
    #ifdef _MSC_VER
        #define MPACK_ALIGNOF(T) __alignof(T)
    #endif
#endif

// MPACK_ALIGNOF may not exist, in which case a workaround is used.



/* Static assert */

#ifndef MPACK_STATIC_ASSERT
    #if defined(__cplusplus)
        #if __cplusplus >= 201103L
            #define MPACK_STATIC_ASSERT static_assert
        #endif
    #elif defined(__STDC_VERSION__)
        #if __STDC_VERSION__ >= 201112L
            #define MPACK_STATIC_ASSERT _Static_assert
        #endif
    #endif
#endif

#if !MPACK_NO_BUILTINS
    #ifndef MPACK_STATIC_ASSERT
        #if defined(__has_feature)
            #if __has_feature(cxx_static_assert)
                #define MPACK_STATIC_ASSERT static_assert
            #elif __has_feature(c_static_assert)
                #define MPACK_STATIC_ASSERT _Static_assert
            #endif
        #endif
    #endif

    #ifndef MPACK_STATIC_ASSERT
        #if defined(__GNUC__)
            /* Diagnostic push is not supported before GCC 4.6. */
            #if defined(__clang__) || __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)
                #ifndef __cplusplus
                    #if defined(__clang__) || __GNUC__ >= 5
                    #define MPACK_IGNORE_PEDANTIC "GCC diagnostic ignored \"-Wpedantic\""
                    #else
                    #define MPACK_IGNORE_PEDANTIC "GCC diagnostic ignored \"-pedantic\""
                    #endif
                    #define MPACK_STATIC_ASSERT(expr, str) do { \
                        _Pragma ("GCC diagnostic push") \
                        _Pragma (MPACK_IGNORE_PEDANTIC) \
                        _Pragma ("GCC diagnostic ignored \"-Wc++-compat\"") \
                        _Static_assert(expr, str); \
                        _Pragma ("GCC diagnostic pop") \
                    } while (0)
                #endif
            #endif
        #endif
    #endif

    #ifndef MPACK_STATIC_ASSERT
        #ifdef _MSC_VER
            #if _MSC_VER >= 1600
                #define MPACK_STATIC_ASSERT static_assert
            #endif
        #endif
    #endif
#endif

#ifndef MPACK_STATIC_ASSERT
    #define MPACK_STATIC_ASSERT(expr, str) (MPACK_UNUSED(sizeof(char[1 - 2*!(expr)])))
#endif



/* _Generic */

#ifndef MPACK_HAS_GENERIC
    #if defined(__clang__) && defined(__has_feature)
        // With Clang we can test for _Generic support directly
        // and ignore C/C++ version
        #if __has_feature(c_generic_selections)
            #define MPACK_HAS_GENERIC 1
        #else
            #define MPACK_HAS_GENERIC 0
        #endif
    #endif
#endif

#ifndef MPACK_HAS_GENERIC
    #if defined(__STDC_VERSION__)
        #if __STDC_VERSION__ >= 201112L
            #if defined(__GNUC__) && !defined(__clang__)
                // GCC does not have full C11 support in GCC 4.7 and 4.8
                #if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 9)
                    #define MPACK_HAS_GENERIC 1
                #endif
            #else
                // We hope other compilers aren't lying about C11/_Generic support
                #define MPACK_HAS_GENERIC 1
            #endif
        #endif
    #endif
#endif

#ifndef MPACK_HAS_GENERIC
    #define MPACK_HAS_GENERIC 0
#endif



/*
 * Finite Math
 *
 * -ffinite-math-only, included in -ffast-math, breaks functions that
 * that check for non-finite real values such as isnan() and isinf().
 *
 * We should use this to trap errors when reading data that contains
 * non-finite reals. This isn't currently implemented.
 */

#ifndef MPACK_FINITE_MATH
#if defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__
#define MPACK_FINITE_MATH 1
#endif
#endif

#ifndef MPACK_FINITE_MATH
#define MPACK_FINITE_MATH 0
#endif



/*
 * Endianness checks
 *
 * These define MPACK_NHSWAP*() which swap network<->host byte
 * order when needed.
 *
 * We leave them undefined if we can't determine the endianness
 * at compile-time, in which case we fall back to bit-shifts.
 *
 * See the notes in mpack-common.h.
 */

#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && defined(__ORDER_BIG_ENDIAN__)
    #if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
        #define MPACK_NHSWAP16(x) (x)
        #define MPACK_NHSWAP32(x) (x)
        #define MPACK_NHSWAP64(x) (x)
    #elif __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__

        #if !MPACK_NO_BUILTINS
            #if defined(__clang__)
                #ifdef __has_builtin
                    // Unlike the GCC builtins, the bswap builtins in Clang
                    // significantly improve ARM performance.
                    #if __has_builtin(__builtin_bswap16)
                        #define MPACK_NHSWAP16(x) __builtin_bswap16(x)
                    #endif
                    #if __has_builtin(__builtin_bswap32)
                        #define MPACK_NHSWAP32(x) __builtin_bswap32(x)
                    #endif
                    #if __has_builtin(__builtin_bswap64)
                        #define MPACK_NHSWAP64(x) __builtin_bswap64(x)
                    #endif
                #endif

            #elif defined(__GNUC__)

                // The GCC bswap builtins are apparently poorly optimized on older
                // versions of GCC, so we set a minimum version here just in case.
                //     http://hardwarebug.org/2010/01/14/beware-the-builtins/

                #if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5)
                    #define MPACK_NHSWAP64(x) __builtin_bswap64(x)
                #endif

                // __builtin_bswap16() was not implemented on all platforms
                // until GCC 4.8.0:
                //     https://gcc.gnu.org/bugzilla/show_bug.cgi?id=52624
                //
                // The 16- and 32-bit versions in GCC significantly reduce performance
                // on ARM with little effect on code size so we don't use them.

            #endif
        #endif
    #endif

#elif defined(_MSC_VER) && defined(_WIN32) && MPACK_STDLIB && !MPACK_NO_BUILTINS

    // On Windows, we assume x86 and x86_64 are always little-endian.
    // We make no assumptions about ARM even though all current
    // Windows Phone devices are little-endian in case Microsoft's
    // compiler is ever used with a big-endian ARM device.

    // These are functions in <stdlib.h> so we depend on MPACK_STDLIB.
    // It's not clear if these are actually faster than just doing the
    // swap manually; maybe we shouldn't bother with this.

    #if defined(_M_IX86) || defined(_M_X64) || defined(_M_AMD64)
        #define MPACK_NHSWAP16(x) _byteswap_ushort(x)
        #define MPACK_NHSWAP32(x) _byteswap_ulong(x)
        #define MPACK_NHSWAP64(x) _byteswap_uint64(x)
    #endif

#endif

#if defined(__FLOAT_WORD_ORDER__) && defined(__BYTE_ORDER__)

    // We check where possible that the float byte order matches the
    // integer byte order. This is extremely unlikely to fail, but
    // we check anyway just in case.
    //
    // (The static assert is placed in float/double encoders instead
    // of here because our static assert fallback doesn't work at
    // file scope)

    #define MPACK_CHECK_FLOAT_ORDER() \
        MPACK_STATIC_ASSERT(__FLOAT_WORD_ORDER__ == __BYTE_ORDER__, \
            "float byte order does not match int byte order! float/double " \
            "encoding is not properly implemented on this platform.")

#endif

#ifndef MPACK_CHECK_FLOAT_ORDER
    #define MPACK_CHECK_FLOAT_ORDER() /* nothing */
#endif


/*
 * Here we define mpack_assert() and mpack_break(). They both work like a normal
 * assertion function in debug mode, causing a trap or abort. However, on some platforms
 * you can safely resume execution from mpack_break(), whereas mpack_assert() is
 * always fatal.
 *
 * In release mode, mpack_assert() is converted to an assurance to the compiler
 * that the expression cannot be false (via e.g. __assume() or __builtin_unreachable())
 * to improve optimization where supported. There is thus no point in "safely" handling
 * the case of this being false. Writing mpack_assert(0) rarely makes sense (except
 * possibly as a default handler in a switch) since the compiler will throw away any
 * code after it. If at any time an mpack_assert() is not true, the behaviour is
 * undefined. This also means the expression is evaluated even in release.
 *
 * mpack_break() on the other hand is compiled to nothing in release. It is
 * used in situations where we want to highlight a programming error as early as
 * possible (in the debugger), but we still handle the situation safely if it
 * happens in release to avoid producing incorrect results (such as in
 * MPACK_WRITE_TRACKING.) It does not take an expression to test because it
 * belongs in a safe-handling block after its failing condition has been tested.
 *
 * If stdio is available, we can add a format string describing the error, and
 * on some compilers we can declare it noreturn to get correct results from static
 * analysis tools. Note that the format string and arguments are not evaluated unless
 * the assertion is hit.
 *
 * Note that any arguments to mpack_assert() beyond the first are only evaluated
 * if the expression is false (and are never evaluated in release.)
 *
 * mpack_assert_fail() and mpack_break_hit() are defined separately
 * because assert is noreturn and break isn't. This distinction is very
 * important for static analysis tools to give correct results.
 */

#if MPACK_DEBUG
    MPACK_NORETURN(void mpack_assert_fail_wrapper(const char* message));
    #if MPACK_STDIO
        MPACK_NORETURN(void mpack_assert_fail_format(const char* format, ...));
        #define mpack_assert_fail_at(line, file, exprstr, format, ...) \
                MPACK_EXPAND(mpack_assert_fail_format("mpack assertion failed at " file ":" #line "\n%s\n" format, exprstr, __VA_ARGS__))
    #else
        #define mpack_assert_fail_at(line, file, exprstr, format, ...) \
                mpack_assert_fail_wrapper("mpack assertion failed at " file ":" #line "\n" exprstr "\n")
    #endif

    #define mpack_assert_fail_pos(line, file, exprstr, expr, ...) \
            MPACK_EXPAND(mpack_assert_fail_at(line, file, exprstr, __VA_ARGS__))

    // This contains a workaround to the pedantic C99 requirement of having at
    // least one argument to a variadic macro. The first argument is the
    // boolean expression, the optional second argument (if provided) must be a
    // literal format string, and any additional arguments are the format
    // argument list.
    //
    // Unfortunately this means macros are expanded in the expression before it
    // gets stringified. I haven't found a workaround to this.
    //
    // This adds two unused arguments to the format argument list when a
    // format string is provided, so this would complicate the use of
    // -Wformat and __attribute__((__format__)) on mpack_assert_fail_format()
    // if we ever bothered to implement it.
    #define mpack_assert(...) \
            MPACK_EXPAND(((!(MPACK_EXTRACT_ARG0(__VA_ARGS__))) ? \
                mpack_assert_fail_pos(__LINE__, __FILE__, MPACK_STRINGIFY_ARG0(__VA_ARGS__) , __VA_ARGS__ , "", NULL) : \
                (void)0))

    void mpack_break_hit(const char* message);
    #if MPACK_STDIO
        void mpack_break_hit_format(const char* format, ...);
        #define mpack_break_hit_at(line, file, ...) \
                MPACK_EXPAND(mpack_break_hit_format("mpack breakpoint hit at " file ":" #line "\n" __VA_ARGS__))
    #else
        #define mpack_break_hit_at(line, file, ...) \
                mpack_break_hit("mpack breakpoint hit at " file ":" #line )
    #endif
    #define mpack_break_hit_pos(line, file, ...) MPACK_EXPAND(mpack_break_hit_at(line, file, __VA_ARGS__))
    #define mpack_break(...) MPACK_EXPAND(mpack_break_hit_pos(__LINE__, __FILE__, __VA_ARGS__))
#else
    #define mpack_assert(...) \
            (MPACK_EXPAND((!(MPACK_EXTRACT_ARG0(__VA_ARGS__))) ? \
                (MPACK_UNREACHABLE, (void)0) : \
                (void)0))
    #define mpack_break(...) ((void)0)
#endif



// make sure we don't use the stdlib directly during development
#if MPACK_STDLIB && defined(MPACK_UNIT_TESTS) && MPACK_INTERNAL && defined(__GNUC__)
    #undef memcmp
    #undef memcpy
    #undef memmove
    #undef memset
    #undef strlen
    #undef malloc
    #undef calloc
    #undef realloc
    #undef free
    #pragma GCC poison memcmp
    #pragma GCC poison memcpy
    #pragma GCC poison memmove
    #pragma GCC poison memset
    #pragma GCC poison strlen
    #pragma GCC poison malloc
    #pragma GCC poison calloc
    #pragma GCC poison realloc
    #pragma GCC poison free
#endif



// If we don't have these stdlib functions, we need to define them ourselves.
// Either way we give them a lowercase name to make the code a bit nicer.

#ifdef MPACK_MEMCMP
    #define mpack_memcmp MPACK_MEMCMP
#else
    int mpack_memcmp(const void* s1, const void* s2, size_t n);
#endif

#ifdef MPACK_MEMCPY
    #define mpack_memcpy MPACK_MEMCPY
#else
    void* mpack_memcpy(void* MPACK_RESTRICT s1, const void* MPACK_RESTRICT s2, size_t n);
#endif

#ifdef MPACK_MEMMOVE
    #define mpack_memmove MPACK_MEMMOVE
#else
    void* mpack_memmove(void* s1, const void* s2, size_t n);
#endif

#ifdef MPACK_MEMSET
    #define mpack_memset MPACK_MEMSET
#else
    void* mpack_memset(void* s, int c, size_t n);
#endif

#ifdef MPACK_STRLEN
    #define mpack_strlen MPACK_STRLEN
#else
    size_t mpack_strlen(const char* s);
#endif



#if MPACK_STDIO
    #if defined(WIN32)
        #define mpack_snprintf _snprintf
    #else
        #define mpack_snprintf snprintf
    #endif
#endif



/* Debug logging */
#if 0
    #include <stdio.h>
    #define mpack_log(...) (MPACK_EXPAND(printf(__VA_ARGS__)), fflush(stdout))
#else
    #define mpack_log(...) ((void)0)
#endif



/* Make sure our configuration makes sense */
#ifndef MPACK_MALLOC
    #if MPACK_STDIO
        #error "MPACK_STDIO requires preprocessor definitions for MPACK_MALLOC and MPACK_FREE."
    #endif
    #if MPACK_READ_TRACKING
        #error "MPACK_READ_TRACKING requires preprocessor definitions for MPACK_MALLOC and MPACK_FREE."
    #endif
    #if MPACK_WRITE_TRACKING
        #error "MPACK_WRITE_TRACKING requires preprocessor definitions for MPACK_MALLOC and MPACK_FREE."
    #endif
#endif



/* Implement realloc if unavailable */
#ifdef MPACK_MALLOC
    #ifdef MPACK_REALLOC
        MPACK_INLINE void* mpack_realloc(void* old_ptr, size_t used_size, size_t new_size) {
            MPACK_UNUSED(used_size);
            return MPACK_REALLOC(old_ptr, new_size);
        }
    #else
        void* mpack_realloc(void* old_ptr, size_t used_size, size_t new_size);
    #endif
#endif



/** @endcond */
/**
 * @}
 */

MPACK_EXTERN_C_END
MPACK_SILENCE_WARNINGS_END

#endif

/* mpack/mpack-common.h.h */

/**
 * @file
 *
 * Defines types and functions shared by the MPack reader and writer.
 */

#ifndef MPACK_COMMON_H
#define MPACK_COMMON_H 1

/* #include "mpack-platform.h" */

#ifndef MPACK_PRINT_BYTE_COUNT
#define MPACK_PRINT_BYTE_COUNT 12
#endif

MPACK_SILENCE_WARNINGS_BEGIN
MPACK_EXTERN_C_BEGIN



/**
 * @defgroup common Tags and Common Elements
 *
 * Contains types, constants and functions shared by both the encoding
 * and decoding portions of MPack.
 *
 * @{
 */

/* Version information */

#define MPACK_VERSION_MAJOR 1  /**< The major version number of MPack. */
#define MPACK_VERSION_MINOR 1  /**< The minor version number of MPack. */
#define MPACK_VERSION_PATCH 1  /**< The patch version number of MPack. */

/** A number containing the version number of MPack for comparison purposes. */
#define MPACK_VERSION ((MPACK_VERSION_MAJOR * 10000) + \
        (MPACK_VERSION_MINOR * 100) + MPACK_VERSION_PATCH)

/** A macro to test for a minimum version of MPack. */
#define MPACK_VERSION_AT_LEAST(major, minor, patch) \
        (MPACK_VERSION >= (((major) * 10000) + ((minor) * 100) + (patch)))

/** @cond */
#if (MPACK_VERSION_PATCH > 0)
#define MPACK_VERSION_STRING_BASE \
        MPACK_STRINGIFY(MPACK_VERSION_MAJOR) "." \
        MPACK_STRINGIFY(MPACK_VERSION_MINOR) "." \
        MPACK_STRINGIFY(MPACK_VERSION_PATCH)
#else
#define MPACK_VERSION_STRING_BASE \
        MPACK_STRINGIFY(MPACK_VERSION_MAJOR) "." \
        MPACK_STRINGIFY(MPACK_VERSION_MINOR)
#endif
/** @endcond */

/**
 * @def MPACK_VERSION_STRING
 * @hideinitializer
 *
 * A string containing the MPack version.
 */
#if MPACK_RELEASE_VERSION
#define MPACK_VERSION_STRING MPACK_VERSION_STRING_BASE
#else
#define MPACK_VERSION_STRING MPACK_VERSION_STRING_BASE "dev"
#endif

/**
 * @def MPACK_LIBRARY_STRING
 * @hideinitializer
 *
 * A string describing MPack, containing the library name, version and debug mode.
 */
#if MPACK_DEBUG
#define MPACK_LIBRARY_STRING "MPack " MPACK_VERSION_STRING "-debug"
#else
#define MPACK_LIBRARY_STRING "MPack " MPACK_VERSION_STRING
#endif

/** @cond */
/**
 * @def MPACK_MAXIMUM_TAG_SIZE
 *
 * The maximum encoded size of a tag in bytes.
 */
#define MPACK_MAXIMUM_TAG_SIZE 9
/** @endcond */

#if MPACK_EXTENSIONS
/**
 * @def MPACK_TIMESTAMP_NANOSECONDS_MAX
 *
 * The maximum value of nanoseconds for a timestamp.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 */
#define MPACK_TIMESTAMP_NANOSECONDS_MAX 999999999
#endif



#if MPACK_COMPATIBILITY
/**
 * Versions of the MessagePack format.
 *
 * A reader, writer, or tree can be configured to serialize in an older
 * version of the MessagePack spec. This is necessary to interface with
 * older MessagePack libraries that do not support new MessagePack features.
 *
 * @note This requires @ref MPACK_COMPATIBILITY.
 */
typedef enum mpack_version_t {

    /**
     * Version 1.0/v4, supporting only the @c raw type without @c str8.
     */
    mpack_version_v4 = 4,

    /**
     * Version 2.0/v5, supporting the @c str8, @c bin and @c ext types.
     */
    mpack_version_v5 = 5,

    /**
     * The most recent supported version of MessagePack. This is the default.
     */
    mpack_version_current = mpack_version_v5,

} mpack_version_t;
#endif

/**
 * Error states for MPack objects.
 *
 * When a reader, writer, or tree is in an error state, all subsequent calls
 * are ignored and their return values are nil/zero. You should check whether
 * the source is in an error state before using such values.
 */
typedef enum mpack_error_t {
    mpack_ok = 0,        /**< No error. */
    mpack_error_io = 2,  /**< The reader or writer failed to fill or flush, or some other file or socket error occurred. */
    mpack_error_invalid, /**< The data read is not valid MessagePack. */
    mpack_error_unsupported, /**< The data read is not supported by this configuration of MPack. (See @ref MPACK_EXTENSIONS.) */
    mpack_error_type,    /**< The type or value range did not match what was expected by the caller. */
    mpack_error_too_big, /**< A read or write was bigger than the maximum size allowed for that operation. */
    mpack_error_memory,  /**< An allocation failure occurred. */
    mpack_error_bug,     /**< The MPack API was used incorrectly. (This will always assert in debug mode.) */
    mpack_error_data,    /**< The contained data is not valid. */
    mpack_error_eof,     /**< The reader failed to read because of file or socket EOF */
} mpack_error_t;

/**
 * Converts an MPack error to a string. This function returns an empty
 * string when MPACK_DEBUG is not set.
 */
const char* mpack_error_to_string(mpack_error_t error);

/**
 * Defines the type of a MessagePack tag.
 *
 * Note that extension types, both user defined and built-in, are represented
 * in tags as @ref mpack_type_ext. The value for an extension type is stored
 * separately.
 */
typedef enum mpack_type_t {
    mpack_type_missing = 0, /**< Special type indicating a missing optional value. */
    mpack_type_nil,         /**< A null value. */
    mpack_type_bool,        /**< A boolean (true or false.) */
    mpack_type_int,         /**< A 64-bit signed integer. */
    mpack_type_uint,        /**< A 64-bit unsigned integer. */
    mpack_type_float,       /**< A 32-bit IEEE 754 floating point number. */
    mpack_type_double,      /**< A 64-bit IEEE 754 floating point number. */
    mpack_type_str,         /**< A string. */
    mpack_type_bin,         /**< A chunk of binary data. */
    mpack_type_array,       /**< An array of MessagePack objects. */
    mpack_type_map,         /**< An ordered map of key/value pairs of MessagePack objects. */

    #if MPACK_EXTENSIONS
    /**
     * A typed MessagePack extension object containing a chunk of binary data.
     *
     * @note This requires @ref MPACK_EXTENSIONS.
     */
    mpack_type_ext,
    #endif
} mpack_type_t;

/**
 * Converts an MPack type to a string. This function returns an empty
 * string when MPACK_DEBUG is not set.
 */
const char* mpack_type_to_string(mpack_type_t type);

#if MPACK_EXTENSIONS
/**
 * A timestamp.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 */
typedef struct mpack_timestamp_t {
    int64_t seconds; /*< The number of seconds (signed) since 1970-01-01T00:00:00Z. */
    uint32_t nanoseconds; /*< The number of additional nanoseconds, between 0 and 999,999,999. */
} mpack_timestamp_t;
#endif

/**
 * An MPack tag is a MessagePack object header. It is a variant type
 * representing any kind of object, and includes the length of compound types
 * (e.g. map, array, string) or the value of non-compound types (e.g. boolean,
 * integer, float.)
 *
 * If the type is compound (str, bin, ext, array or map), the contained
 * elements or bytes are stored separately.
 *
 * This structure is opaque; its fields should not be accessed outside
 * of MPack.
 */
typedef struct mpack_tag_t mpack_tag_t;

/* Hide internals from documentation */
/** @cond */
struct mpack_tag_t {
    mpack_type_t type; /*< The type of value. */

    #if MPACK_EXTENSIONS
    int8_t exttype; /*< The extension type if the type is @ref mpack_type_ext. */
    #endif

    /* The value for non-compound types. */
    union {
        uint64_t u; /*< The value if the type is unsigned int. */
        int64_t  i; /*< The value if the type is signed int. */
        bool     b; /*< The value if the type is bool. */

        #if MPACK_FLOAT
        float    f; /*< The value if the type is float. */
        #else
        uint32_t f; /*< The raw value if the type is float. */
        #endif

        #if MPACK_DOUBLE
        double   d; /*< The value if the type is double. */
        #else
        uint64_t d; /*< The raw value if the type is double. */
        #endif

        /* The number of bytes if the type is str, bin or ext. */
        uint32_t l;

        /* The element count if the type is an array, or the number of
            key/value pairs if the type is map. */
        uint32_t n;
    } v;
};
/** @endcond */

/**
 * @name Tag Generators
 * @{
 */

/**
 * @def MPACK_TAG_ZERO
 *
 * An @ref mpack_tag_t initializer that zeroes the given tag.
 *
 * @warning This does not make the tag nil! The tag's type is invalid when
 * initialized this way. Use @ref mpack_tag_make_nil() to generate a nil tag.
 */
#if MPACK_EXTENSIONS
#define MPACK_TAG_ZERO {(mpack_type_t)0, 0, {0}}
#else
#define MPACK_TAG_ZERO {(mpack_type_t)0, {0}}
#endif

/** Generates a nil tag. */
MPACK_INLINE mpack_tag_t mpack_tag_make_nil(void) {
    mpack_tag_t ret = MPACK_TAG_ZERO;
    ret.type = mpack_type_nil;
    return ret;
}

/** Generates a bool tag. */
MPACK_INLINE mpack_tag_t mpack_tag_make_bool(bool value) {
    mpack_tag_t ret = MPACK_TAG_ZERO;
    ret.type = mpack_type_bool;
    ret.v.b = value;
    return ret;
}

/** Generates a bool tag with value true. */
MPACK_INLINE mpack_tag_t mpack_tag_make_true(void) {
    mpack_tag_t ret = MPACK_TAG_ZERO;
    ret.type = mpack_type_bool;
    ret.v.b = true;
    return ret;
}

/** Generates a bool tag with value false. */
MPACK_INLINE mpack_tag_t mpack_tag_make_false(void) {
    mpack_tag_t ret = MPACK_TAG_ZERO;
    ret.type = mpack_type_bool;
    ret.v.b = false;
    return ret;
}

/** Generates a signed int tag. */
MPACK_INLINE mpack_tag_t mpack_tag_make_int(int64_t value) {
    mpack_tag_t ret = MPACK_TAG_ZERO;
    ret.type = mpack_type_int;
    ret.v.i = value;
    return ret;
}

/** Generates an unsigned int tag. */
MPACK_INLINE mpack_tag_t mpack_tag_make_uint(uint64_t value) {
    mpack_tag_t ret = MPACK_TAG_ZERO;
    ret.type = mpack_type_uint;
    ret.v.u = value;
    return ret;
}

#if MPACK_FLOAT
/** Generates a float tag. */
MPACK_INLINE mpack_tag_t mpack_tag_make_float(float value)
#else
/** Generates a float tag from a raw uint32_t. */
MPACK_INLINE mpack_tag_t mpack_tag_make_raw_float(uint32_t value)
#endif
{
    mpack_tag_t ret = MPACK_TAG_ZERO;
    ret.type = mpack_type_float;
    ret.v.f = value;
    return ret;
}

#if MPACK_DOUBLE
/** Generates a double tag. */
MPACK_INLINE mpack_tag_t mpack_tag_make_double(double value)
#else
/** Generates a double tag from a raw uint64_t. */
MPACK_INLINE mpack_tag_t mpack_tag_make_raw_double(uint64_t value)
#endif
{
    mpack_tag_t ret = MPACK_TAG_ZERO;
    ret.type = mpack_type_double;
    ret.v.d = value;
    return ret;
}

/** Generates an array tag. */
MPACK_INLINE mpack_tag_t mpack_tag_make_array(uint32_t count) {
    mpack_tag_t ret = MPACK_TAG_ZERO;
    ret.type = mpack_type_array;
    ret.v.n = count;
    return ret;
}

/** Generates a map tag. */
MPACK_INLINE mpack_tag_t mpack_tag_make_map(uint32_t count) {
    mpack_tag_t ret = MPACK_TAG_ZERO;
    ret.type = mpack_type_map;
    ret.v.n = count;
    return ret;
}

/** Generates a str tag. */
MPACK_INLINE mpack_tag_t mpack_tag_make_str(uint32_t length) {
    mpack_tag_t ret = MPACK_TAG_ZERO;
    ret.type = mpack_type_str;
    ret.v.l = length;
    return ret;
}

/** Generates a bin tag. */
MPACK_INLINE mpack_tag_t mpack_tag_make_bin(uint32_t length) {
    mpack_tag_t ret = MPACK_TAG_ZERO;
    ret.type = mpack_type_bin;
    ret.v.l = length;
    return ret;
}

#if MPACK_EXTENSIONS
/**
 * Generates an ext tag.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 */
MPACK_INLINE mpack_tag_t mpack_tag_make_ext(int8_t exttype, uint32_t length) {
    mpack_tag_t ret = MPACK_TAG_ZERO;
    ret.type = mpack_type_ext;
    ret.exttype = exttype;
    ret.v.l = length;
    return ret;
}
#endif

/**
 * @}
 */

/**
 * @name Tag Querying Functions
 * @{
 */

/**
 * Gets the type of a tag.
 */
MPACK_INLINE mpack_type_t mpack_tag_type(mpack_tag_t* tag) {
    return tag->type;
}

/**
 * Gets the boolean value of a bool-type tag. The tag must be of type @ref
 * mpack_type_bool.
 *
 * This asserts that the type in the tag is @ref mpack_type_bool. (No check is
 * performed if MPACK_DEBUG is not set.)
 */
MPACK_INLINE bool mpack_tag_bool_value(mpack_tag_t* tag) {
    mpack_assert(tag->type == mpack_type_bool, "tag is not a bool!");
    return tag->v.b;
}

/**
 * Gets the signed integer value of an int-type tag.
 *
 * This asserts that the type in the tag is @ref mpack_type_int. (No check is
 * performed if MPACK_DEBUG is not set.)
 *
 * @warning This does not convert between signed and unsigned tags! A positive
 * integer may be stored in a tag as either @ref mpack_type_int or @ref
 * mpack_type_uint. You must check the type first; this can only be used if the
 * type is @ref mpack_type_int.
 *
 * @see mpack_type_int
 */
MPACK_INLINE int64_t mpack_tag_int_value(mpack_tag_t* tag) {
    mpack_assert(tag->type == mpack_type_int, "tag is not an int!");
    return tag->v.i;
}

/**
 * Gets the unsigned integer value of a uint-type tag.
 *
 * This asserts that the type in the tag is @ref mpack_type_uint. (No check is
 * performed if MPACK_DEBUG is not set.)
 *
 * @warning This does not convert between signed and unsigned tags! A positive
 * integer may be stored in a tag as either @ref mpack_type_int or @ref
 * mpack_type_uint. You must check the type first; this can only be used if the
 * type is @ref mpack_type_uint.
 *
 * @see mpack_type_uint
 */
MPACK_INLINE uint64_t mpack_tag_uint_value(mpack_tag_t* tag) {
    mpack_assert(tag->type == mpack_type_uint, "tag is not a uint!");
    return tag->v.u;
}

/**
 * Gets the float value of a float-type tag.
 *
 * This asserts that the type in the tag is @ref mpack_type_float. (No check is
 * performed if MPACK_DEBUG is not set.)
 *
 * @warning This does not convert between float and double tags! This can only
 * be used if the type is @ref mpack_type_float.
 *
 * @see mpack_type_float
 */
MPACK_INLINE
#if MPACK_FLOAT
float mpack_tag_float_value(mpack_tag_t* tag)
#else
uint32_t mpack_tag_raw_float_value(mpack_tag_t* tag)
#endif
{
    mpack_assert(tag->type == mpack_type_float, "tag is not a float!");
    return tag->v.f;
}

/**
 * Gets the double value of a double-type tag.
 *
 * This asserts that the type in the tag is @ref mpack_type_double. (No check
 * is performed if MPACK_DEBUG is not set.)
 *
 * @warning This does not convert between float and double tags! This can only
 * be used if the type is @ref mpack_type_double.
 *
 * @see mpack_type_double
 */
MPACK_INLINE
#if MPACK_DOUBLE
double mpack_tag_double_value(mpack_tag_t* tag)
#else
uint64_t mpack_tag_raw_double_value(mpack_tag_t* tag)
#endif
{
    mpack_assert(tag->type == mpack_type_double, "tag is not a double!");
    return tag->v.d;
}

/**
 * Gets the number of elements in an array tag.
 *
 * This asserts that the type in the tag is @ref mpack_type_array. (No check is
 * performed if MPACK_DEBUG is not set.)
 *
 * @see mpack_type_array
 */
MPACK_INLINE uint32_t mpack_tag_array_count(mpack_tag_t* tag) {
    mpack_assert(tag->type == mpack_type_array, "tag is not an array!");
    return tag->v.n;
}

/**
 * Gets the number of key-value pairs in a map tag.
 *
 * This asserts that the type in the tag is @ref mpack_type_map. (No check is
 * performed if MPACK_DEBUG is not set.)
 *
 * @see mpack_type_map
 */
MPACK_INLINE uint32_t mpack_tag_map_count(mpack_tag_t* tag) {
    mpack_assert(tag->type == mpack_type_map, "tag is not a map!");
    return tag->v.n;
}

/**
 * Gets the length in bytes of a str-type tag.
 *
 * This asserts that the type in the tag is @ref mpack_type_str. (No check is
 * performed if MPACK_DEBUG is not set.)
 *
 * @see mpack_type_str
 */
MPACK_INLINE uint32_t mpack_tag_str_length(mpack_tag_t* tag) {
    mpack_assert(tag->type == mpack_type_str, "tag is not a str!");
    return tag->v.l;
}

/**
 * Gets the length in bytes of a bin-type tag.
 *
 * This asserts that the type in the tag is @ref mpack_type_bin. (No check is
 * performed if MPACK_DEBUG is not set.)
 *
 * @see mpack_type_bin
 */
MPACK_INLINE uint32_t mpack_tag_bin_length(mpack_tag_t* tag) {
    mpack_assert(tag->type == mpack_type_bin, "tag is not a bin!");
    return tag->v.l;
}

#if MPACK_EXTENSIONS
/**
 * Gets the length in bytes of an ext-type tag.
 *
 * This asserts that the type in the tag is @ref mpack_type_ext. (No check is
 * performed if MPACK_DEBUG is not set.)
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 *
 * @see mpack_type_ext
 */
MPACK_INLINE uint32_t mpack_tag_ext_length(mpack_tag_t* tag) {
    mpack_assert(tag->type == mpack_type_ext, "tag is not an ext!");
    return tag->v.l;
}

/**
 * Gets the extension type (exttype) of an ext-type tag.
 *
 * This asserts that the type in the tag is @ref mpack_type_ext. (No check is
 * performed if MPACK_DEBUG is not set.)
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 *
 * @see mpack_type_ext
 */
MPACK_INLINE int8_t mpack_tag_ext_exttype(mpack_tag_t* tag) {
    mpack_assert(tag->type == mpack_type_ext, "tag is not an ext!");
    return tag->exttype;
}
#endif

/**
 * Gets the length in bytes of a str-, bin- or ext-type tag.
 *
 * This asserts that the type in the tag is @ref mpack_type_str, @ref
 * mpack_type_bin or @ref mpack_type_ext. (No check is performed if MPACK_DEBUG
 * is not set.)
 *
 * @see mpack_type_str
 * @see mpack_type_bin
 * @see mpack_type_ext
 */
MPACK_INLINE uint32_t mpack_tag_bytes(mpack_tag_t* tag) {
    #if MPACK_EXTENSIONS
    mpack_assert(tag->type == mpack_type_str || tag->type == mpack_type_bin
            || tag->type == mpack_type_ext, "tag is not a str, bin or ext!");
    #else
    mpack_assert(tag->type == mpack_type_str || tag->type == mpack_type_bin,
            "tag is not a str or bin!");
    #endif
    return tag->v.l;
}

/**
 * @}
 */

/**
 * @name Other tag functions
 * @{
 */

#if MPACK_EXTENSIONS
/**
 * The extension type for a timestamp.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 */
#define MPACK_EXTTYPE_TIMESTAMP ((int8_t)(-1))
#endif

/**
 * Compares two tags with an arbitrary fixed ordering. Returns 0 if the tags are
 * equal, a negative integer if left comes before right, or a positive integer
 * otherwise.
 *
 * \warning The ordering is not guaranteed to be preserved across MPack versions; do
 * not rely on it in persistent data.
 *
 * \warning Floating point numbers are compared bit-for-bit, not using the language's
 * operator==. This means that NaNs with matching representation will compare equal.
 * This behaviour is up for debate; see comments in the definition of mpack_tag_cmp().
 *
 * See mpack_tag_equal() for more information on when tags are considered equal.
 */
int mpack_tag_cmp(mpack_tag_t left, mpack_tag_t right);

/**
 * Compares two tags for equality. Tags are considered equal if the types are compatible
 * and the values (for non-compound types) are equal.
 *
 * The field width of variable-width fields is ignored (and in fact is not stored
 * in a tag), and positive numbers in signed integers are considered equal to their
 * unsigned counterparts. So for example the value 1 stored as a positive fixint
 * is equal to the value 1 stored in a 64-bit unsigned integer field.
 *
 * The "extension type" of an extension object is considered part of the value
 * and must match exactly.
 *
 * \warning Floating point numbers are compared bit-for-bit, not using the language's
 * operator==. This means that NaNs with matching representation will compare equal.
 * This behaviour is up for debate; see comments in the definition of mpack_tag_cmp().
 */
MPACK_INLINE bool mpack_tag_equal(mpack_tag_t left, mpack_tag_t right) {
    return mpack_tag_cmp(left, right) == 0;
}

#if MPACK_DEBUG && MPACK_STDIO
/**
 * Generates a json-like debug description of the given tag into the given buffer.
 *
 * This is only available in debug mode, and only if stdio is available (since
 * it uses snprintf().) It's strictly for debugging purposes.
 *
 * The prefix is used to print the first few hexadecimal bytes of a bin or ext
 * type. Pass NULL if not a bin or ext.
 */
void mpack_tag_debug_pseudo_json(mpack_tag_t tag, char* buffer, size_t buffer_size,
        const char* prefix, size_t prefix_size);

/**
 * Generates a debug string description of the given tag into the given buffer.
 *
 * This is only available in debug mode, and only if stdio is available (since
 * it uses snprintf().) It's strictly for debugging purposes.
 */
void mpack_tag_debug_describe(mpack_tag_t tag, char* buffer, size_t buffer_size);

/** @cond */

/*
 * A callback function for printing pseudo-JSON for debugging purposes.
 *
 * @see mpack_node_print_callback
 */
typedef void (*mpack_print_callback_t)(void* context, const char* data, size_t count);

// helpers for printing debug output
// i feel a bit like i'm re-implementing a buffered writer again...
typedef struct mpack_print_t {
    char* buffer;
    size_t size;
    size_t count;
    mpack_print_callback_t callback;
    void* context;
} mpack_print_t;

void mpack_print_append(mpack_print_t* print, const char* data, size_t count);

MPACK_INLINE void mpack_print_append_cstr(mpack_print_t* print, const char* cstr) {
    mpack_print_append(print, cstr, mpack_strlen(cstr));
}

void mpack_print_flush(mpack_print_t* print);

void mpack_print_file_callback(void* context, const char* data, size_t count);

/** @endcond */

#endif

/**
 * @}
 */

/**
 * @name Deprecated Tag Generators
 * @{
 */

/*
 * "make" has been added to their names to disambiguate them from the
 * value-fetching functions (e.g. mpack_tag_make_bool() vs
 * mpack_tag_bool_value().)
 *
 * The length and count for all compound types was the wrong sign (int32_t
 * instead of uint32_t.) These preserve the old behaviour; the new "make"
 * functions have the correct sign.
 */

/** \deprecated Renamed to mpack_tag_make_nil(). */
MPACK_INLINE mpack_tag_t mpack_tag_nil(void) {
    return mpack_tag_make_nil();
}

/** \deprecated Renamed to mpack_tag_make_bool(). */
MPACK_INLINE mpack_tag_t mpack_tag_bool(bool value) {
    return mpack_tag_make_bool(value);
}

/** \deprecated Renamed to mpack_tag_make_true(). */
MPACK_INLINE mpack_tag_t mpack_tag_true(void) {
    return mpack_tag_make_true();
}

/** \deprecated Renamed to mpack_tag_make_false(). */
MPACK_INLINE mpack_tag_t mpack_tag_false(void) {
    return mpack_tag_make_false();
}

/** \deprecated Renamed to mpack_tag_make_int(). */
MPACK_INLINE mpack_tag_t mpack_tag_int(int64_t value) {
    return mpack_tag_make_int(value);
}

/** \deprecated Renamed to mpack_tag_make_uint(). */
MPACK_INLINE mpack_tag_t mpack_tag_uint(uint64_t value) {
    return mpack_tag_make_uint(value);
}

#if MPACK_FLOAT
/** \deprecated Renamed to mpack_tag_make_float(). */
MPACK_INLINE mpack_tag_t mpack_tag_float(float value) {
    return mpack_tag_make_float(value);
}
#endif

#if MPACK_DOUBLE
/** \deprecated Renamed to mpack_tag_make_double(). */
MPACK_INLINE mpack_tag_t mpack_tag_double(double value) {
    return mpack_tag_make_double(value);
}
#endif

/** \deprecated Renamed to mpack_tag_make_array(). */
MPACK_INLINE mpack_tag_t mpack_tag_array(int32_t count) {
    return mpack_tag_make_array((uint32_t)count);
}

/** \deprecated Renamed to mpack_tag_make_map(). */
MPACK_INLINE mpack_tag_t mpack_tag_map(int32_t count) {
    return mpack_tag_make_map((uint32_t)count);
}

/** \deprecated Renamed to mpack_tag_make_str(). */
MPACK_INLINE mpack_tag_t mpack_tag_str(int32_t length) {
    return mpack_tag_make_str((uint32_t)length);
}

/** \deprecated Renamed to mpack_tag_make_bin(). */
MPACK_INLINE mpack_tag_t mpack_tag_bin(int32_t length) {
    return mpack_tag_make_bin((uint32_t)length);
}

#if MPACK_EXTENSIONS
/** \deprecated Renamed to mpack_tag_make_ext(). */
MPACK_INLINE mpack_tag_t mpack_tag_ext(int8_t exttype, int32_t length) {
    return mpack_tag_make_ext(exttype, (uint32_t)length);
}
#endif

/**
 * @}
 */

/** @cond */

/*
 * Helpers to perform unaligned network-endian loads and stores
 * at arbitrary addresses. Byte-swapping builtins are used if they
 * are available and if they improve performance.
 *
 * These will remain available in the public API so feel free to
 * use them for other purposes, but they are undocumented.
 */

MPACK_INLINE uint8_t mpack_load_u8(const char* p) {
    return (uint8_t)p[0];
}

MPACK_INLINE uint16_t mpack_load_u16(const char* p) {
    #ifdef MPACK_NHSWAP16
    uint16_t val;
    mpack_memcpy(&val, p, sizeof(val));
    return MPACK_NHSWAP16(val);
    #else
    return (uint16_t)((((uint16_t)(uint8_t)p[0]) << 8) |
           ((uint16_t)(uint8_t)p[1]));
    #endif
}

MPACK_INLINE uint32_t mpack_load_u32(const char* p) {
    #ifdef MPACK_NHSWAP32
    uint32_t val;
    mpack_memcpy(&val, p, sizeof(val));
    return MPACK_NHSWAP32(val);
    #else
    return (((uint32_t)(uint8_t)p[0]) << 24) |
           (((uint32_t)(uint8_t)p[1]) << 16) |
           (((uint32_t)(uint8_t)p[2]) <<  8) |
            ((uint32_t)(uint8_t)p[3]);
    #endif
}

MPACK_INLINE uint64_t mpack_load_u64(const char* p) {
    #ifdef MPACK_NHSWAP64
    uint64_t val;
    mpack_memcpy(&val, p, sizeof(val));
    return MPACK_NHSWAP64(val);
    #else
    return (((uint64_t)(uint8_t)p[0]) << 56) |
           (((uint64_t)(uint8_t)p[1]) << 48) |
           (((uint64_t)(uint8_t)p[2]) << 40) |
           (((uint64_t)(uint8_t)p[3]) << 32) |
           (((uint64_t)(uint8_t)p[4]) << 24) |
           (((uint64_t)(uint8_t)p[5]) << 16) |
           (((uint64_t)(uint8_t)p[6]) <<  8) |
            ((uint64_t)(uint8_t)p[7]);
    #endif
}

MPACK_INLINE void mpack_store_u8(char* p, uint8_t val) {
    uint8_t* u = (uint8_t*)p;
    u[0] = val;
}

MPACK_INLINE void mpack_store_u16(char* p, uint16_t val) {
    #ifdef MPACK_NHSWAP16
    val = MPACK_NHSWAP16(val);
    mpack_memcpy(p, &val, sizeof(val));
    #else
    uint8_t* u = (uint8_t*)p;
    u[0] = (uint8_t)((val >> 8) & 0xFF);
    u[1] = (uint8_t)( val       & 0xFF);
    #endif
}

MPACK_INLINE void mpack_store_u32(char* p, uint32_t val) {
    #ifdef MPACK_NHSWAP32
    val = MPACK_NHSWAP32(val);
    mpack_memcpy(p, &val, sizeof(val));
    #else
    uint8_t* u = (uint8_t*)p;
    u[0] = (uint8_t)((val >> 24) & 0xFF);
    u[1] = (uint8_t)((val >> 16) & 0xFF);
    u[2] = (uint8_t)((val >>  8) & 0xFF);
    u[3] = (uint8_t)( val        & 0xFF);
    #endif
}

MPACK_INLINE void mpack_store_u64(char* p, uint64_t val) {
    #ifdef MPACK_NHSWAP64
    val = MPACK_NHSWAP64(val);
    mpack_memcpy(p, &val, sizeof(val));
    #else
    uint8_t* u = (uint8_t*)p;
    u[0] = (uint8_t)((val >> 56) & 0xFF);
    u[1] = (uint8_t)((val >> 48) & 0xFF);
    u[2] = (uint8_t)((val >> 40) & 0xFF);
    u[3] = (uint8_t)((val >> 32) & 0xFF);
    u[4] = (uint8_t)((val >> 24) & 0xFF);
    u[5] = (uint8_t)((val >> 16) & 0xFF);
    u[6] = (uint8_t)((val >>  8) & 0xFF);
    u[7] = (uint8_t)( val        & 0xFF);
    #endif
}

MPACK_INLINE int8_t  mpack_load_i8 (const char* p) {return (int8_t) mpack_load_u8 (p);}
MPACK_INLINE int16_t mpack_load_i16(const char* p) {return (int16_t)mpack_load_u16(p);}
MPACK_INLINE int32_t mpack_load_i32(const char* p) {return (int32_t)mpack_load_u32(p);}
MPACK_INLINE int64_t mpack_load_i64(const char* p) {return (int64_t)mpack_load_u64(p);}
MPACK_INLINE void mpack_store_i8 (char* p, int8_t  val) {mpack_store_u8 (p, (uint8_t) val);}
MPACK_INLINE void mpack_store_i16(char* p, int16_t val) {mpack_store_u16(p, (uint16_t)val);}
MPACK_INLINE void mpack_store_i32(char* p, int32_t val) {mpack_store_u32(p, (uint32_t)val);}
MPACK_INLINE void mpack_store_i64(char* p, int64_t val) {mpack_store_u64(p, (uint64_t)val);}

#if MPACK_FLOAT
MPACK_INLINE float mpack_load_float(const char* p) {
    MPACK_CHECK_FLOAT_ORDER();
    MPACK_STATIC_ASSERT(sizeof(float) == sizeof(uint32_t), "float is wrong size??");
    union {
        float f;
        uint32_t u;
    } v;
    v.u = mpack_load_u32(p);
    return v.f;
}
#endif

#if MPACK_DOUBLE
MPACK_INLINE double mpack_load_double(const char* p) {
    MPACK_CHECK_FLOAT_ORDER();
    MPACK_STATIC_ASSERT(sizeof(double) == sizeof(uint64_t), "double is wrong size??");
    union {
        double d;
        uint64_t u;
    } v;
    v.u = mpack_load_u64(p);
    return v.d;
}
#endif

#if MPACK_FLOAT
MPACK_INLINE void mpack_store_float(char* p, float value) {
    MPACK_CHECK_FLOAT_ORDER();
    union {
        float f;
        uint32_t u;
    } v;
    v.f = value;
    mpack_store_u32(p, v.u);
}
#endif

#if MPACK_DOUBLE
MPACK_INLINE void mpack_store_double(char* p, double value) {
    MPACK_CHECK_FLOAT_ORDER();
    union {
        double d;
        uint64_t u;
    } v;
    v.d = value;
    mpack_store_u64(p, v.u);
}
#endif

#if MPACK_FLOAT && !MPACK_DOUBLE
/**
 * Performs a manual shortening conversion on the raw 64-bit representation of
 * a double. This is useful for parsing doubles on platforms that only support
 * floats (such as AVR.)
 *
 * The significand is truncated rather than rounded and subnormal numbers are
 * set to 0 so this may not be quite as accurate as a real double-to-float
 * conversion.
 */
MPACK_INLINE float mpack_shorten_raw_double_to_float(uint64_t d) {
    MPACK_CHECK_FLOAT_ORDER();
    union {
        float f;
        uint32_t u;
    } v;

    // float has  1 bit sign,  8 bits exponent, 23 bits significand
    // double has 1 bit sign, 11 bits exponent, 52 bits significand

    uint64_t d_sign = (uint64_t)(d >> 63);
    uint64_t d_exponent = (uint32_t)(d >> 52) & ((1 << 11) - 1);
    uint64_t d_significand = d & (((uint64_t)1 << 52) - 1);

    uint32_t f_sign = (uint32_t)d_sign;
    uint32_t f_exponent;
    uint32_t f_significand;

    if (MPACK_UNLIKELY(d_exponent == ((1 << 11) - 1))) {
        // infinity or NAN. shift down to preserve the top bit since it
        // indicates signaling NAN, but also set the low bit if any bits were
        // set (that way we can't shift NAN to infinity.)
        f_exponent = ((1 << 8) - 1);
        f_significand = (uint32_t)(d_significand >> 29) | (d_significand ? 1 : 0);

    } else {
        int fix_bias = (int)d_exponent - ((1 << 10) - 1) + ((1 << 7) - 1);
        if (MPACK_UNLIKELY(fix_bias <= 0)) {
            // we don't currently handle subnormal numbers. just set it to zero.
            f_exponent = 0;
            f_significand = 0;
        } else if (MPACK_UNLIKELY(fix_bias > 0xff)) {
            // exponent is too large; saturate to infinity
            f_exponent = 0xff;
            f_significand = 0;
        } else {
            // a normal number that fits in a float. this is the usual case.
            f_exponent = (uint32_t)fix_bias;
            f_significand = (uint32_t)(d_significand >> 29);
        }
    }

    #if 0
    printf("\n===============\n");
    for (size_t i = 0; i < 64; ++i)
        printf("%i%s",(int)((d>>(63-i))&1),((i%8)==7)?" ":"");
    printf("\n%lu %lu %lu\n", d_sign, d_exponent, d_significand);
    printf("%u %u %u\n", f_sign, f_exponent, f_significand);
    #endif

    v.u = (f_sign << 31) | (f_exponent << 23) | f_significand;
    return v.f;
}
#endif

/** @endcond */



/** @cond */

// Sizes in bytes for the various possible tags
#define MPACK_TAG_SIZE_FIXUINT  1
#define MPACK_TAG_SIZE_U8       2
#define MPACK_TAG_SIZE_U16      3
#define MPACK_TAG_SIZE_U32      5
#define MPACK_TAG_SIZE_U64      9
#define MPACK_TAG_SIZE_FIXINT   1
#define MPACK_TAG_SIZE_I8       2
#define MPACK_TAG_SIZE_I16      3
#define MPACK_TAG_SIZE_I32      5
#define MPACK_TAG_SIZE_I64      9
#define MPACK_TAG_SIZE_FLOAT    5
#define MPACK_TAG_SIZE_DOUBLE   9
#define MPACK_TAG_SIZE_FIXARRAY 1
#define MPACK_TAG_SIZE_ARRAY16  3
#define MPACK_TAG_SIZE_ARRAY32  5
#define MPACK_TAG_SIZE_FIXMAP   1
#define MPACK_TAG_SIZE_MAP16    3
#define MPACK_TAG_SIZE_MAP32    5
#define MPACK_TAG_SIZE_FIXSTR   1
#define MPACK_TAG_SIZE_STR8     2
#define MPACK_TAG_SIZE_STR16    3
#define MPACK_TAG_SIZE_STR32    5
#define MPACK_TAG_SIZE_BIN8     2
#define MPACK_TAG_SIZE_BIN16    3
#define MPACK_TAG_SIZE_BIN32    5
#define MPACK_TAG_SIZE_FIXEXT1  2
#define MPACK_TAG_SIZE_FIXEXT2  2
#define MPACK_TAG_SIZE_FIXEXT4  2
#define MPACK_TAG_SIZE_FIXEXT8  2
#define MPACK_TAG_SIZE_FIXEXT16 2
#define MPACK_TAG_SIZE_EXT8     3
#define MPACK_TAG_SIZE_EXT16    4
#define MPACK_TAG_SIZE_EXT32    6

// size in bytes for complete ext types
#define MPACK_EXT_SIZE_TIMESTAMP4 (MPACK_TAG_SIZE_FIXEXT4 + 4)
#define MPACK_EXT_SIZE_TIMESTAMP8 (MPACK_TAG_SIZE_FIXEXT8 + 8)
#define MPACK_EXT_SIZE_TIMESTAMP12 (MPACK_TAG_SIZE_EXT8 + 12)

/** @endcond */



#if MPACK_READ_TRACKING || MPACK_WRITE_TRACKING
/* Tracks the write state of compound elements (maps, arrays, */
/* strings, binary blobs and extension types) */
/** @cond */

typedef struct mpack_track_element_t {
    mpack_type_t type;
    uint32_t left;

    // indicates that a value still needs to be read/written for an already
    // read/written key. left is not decremented until both key and value are
    // read/written.
    bool key_needs_value;

    // tracks whether the map/array being written is using a builder. if true,
    // the number of elements is automatic, and left is 0.
    bool builder;
} mpack_track_element_t;

typedef struct mpack_track_t {
    size_t count;
    size_t capacity;
    mpack_track_element_t* elements;
} mpack_track_t;

#if MPACK_INTERNAL
mpack_error_t mpack_track_init(mpack_track_t* track);
mpack_error_t mpack_track_grow(mpack_track_t* track);
mpack_error_t mpack_track_push(mpack_track_t* track, mpack_type_t type, uint32_t count);
mpack_error_t mpack_track_push_builder(mpack_track_t* track, mpack_type_t type);
mpack_error_t mpack_track_pop(mpack_track_t* track, mpack_type_t type);
mpack_error_t mpack_track_pop_builder(mpack_track_t* track, mpack_type_t type);
mpack_error_t mpack_track_element(mpack_track_t* track, bool read);
mpack_error_t mpack_track_peek_element(mpack_track_t* track, bool read);
mpack_error_t mpack_track_bytes(mpack_track_t* track, bool read, size_t count);
mpack_error_t mpack_track_str_bytes_all(mpack_track_t* track, bool read, size_t count);
mpack_error_t mpack_track_check_empty(mpack_track_t* track);
mpack_error_t mpack_track_destroy(mpack_track_t* track, bool cancel);
#endif

/** @endcond */
#endif



#if MPACK_INTERNAL
/** @cond */



/* Miscellaneous string functions */

/**
 * Returns true if the given UTF-8 string is valid.
 */
bool mpack_utf8_check(const char* str, size_t bytes);

/**
 * Returns true if the given UTF-8 string is valid and contains no null characters.
 */
bool mpack_utf8_check_no_null(const char* str, size_t bytes);

/**
 * Returns true if the given string has no null bytes.
 */
bool mpack_str_check_no_null(const char* str, size_t bytes);



/** @endcond */
#endif



/**
 * @}
 */

MPACK_EXTERN_C_END
MPACK_SILENCE_WARNINGS_END

#endif


/* mpack/mpack-writer.h.h */

/**
 * @file
 *
 * Declares the MPack Writer.
 */

#ifndef MPACK_WRITER_H
#define MPACK_WRITER_H 1

/* #include "mpack-common.h" */

#if MPACK_WRITER

MPACK_SILENCE_WARNINGS_BEGIN
MPACK_EXTERN_C_BEGIN

#if MPACK_WRITE_TRACKING
struct mpack_track_t;
#endif

/**
 * @defgroup writer Write API
 *
 * The MPack Write API encodes structured data of a fixed (hardcoded) schema to MessagePack.
 *
 * @{
 */

/**
 * @def MPACK_WRITER_MINIMUM_BUFFER_SIZE
 *
 * The minimum buffer size for a writer with a flush function.
 */
#define MPACK_WRITER_MINIMUM_BUFFER_SIZE 32

/**
 * A buffered MessagePack encoder.
 *
 * The encoder wraps an existing buffer and, optionally, a flush function.
 * This allows efficiently encoding to an in-memory buffer or to a stream.
 *
 * All write operations are synchronous; they will block until the
 * data is fully written, or an error occurs.
 */
typedef struct mpack_writer_t mpack_writer_t;

/**
 * The MPack writer's flush function to flush the buffer to the output stream.
 * It should flag an appropriate error on the writer if flushing fails (usually
 * mpack_error_io or mpack_error_memory.)
 *
 * The specified context for callbacks is at writer->context.
 */
typedef void (*mpack_writer_flush_t)(mpack_writer_t* writer, const char* buffer, size_t count);

/**
 * An error handler function to be called when an error is flagged on
 * the writer.
 *
 * The error handler will only be called once on the first error flagged;
 * any subsequent writes and errors are ignored, and the writer is
 * permanently in that error state.
 *
 * MPack is safe against non-local jumps out of error handler callbacks.
 * This means you are allowed to longjmp or throw an exception (in C++,
 * Objective-C, or with SEH) out of this callback.
 *
 * Bear in mind when using longjmp that local non-volatile variables that
 * have changed are undefined when setjmp() returns, so you can't put the
 * writer on the stack in the same activation frame as the setjmp without
 * declaring it volatile.
 *
 * You must still eventually destroy the writer. It is not destroyed
 * automatically when an error is flagged. It is safe to destroy the
 * writer within this error callback, but you will either need to perform
 * a non-local jump, or store something in your context to identify
 * that the writer is destroyed since any future accesses to it cause
 * undefined behavior.
 */
typedef void (*mpack_writer_error_t)(mpack_writer_t* writer, mpack_error_t error);

/**
 * A teardown function to be called when the writer is destroyed.
 */
typedef void (*mpack_writer_teardown_t)(mpack_writer_t* writer);

/* Hide internals from documentation */
/** @cond */

#if MPACK_BUILDER
/**
 * Build buffer pages form a linked list.
 *
 * They don't always fill up. If there is not enough space within them to write
 * a tag or place an mpack_build_t, a new page is allocated. For this reason
 * they store the number of used bytes.
 */
typedef struct mpack_builder_page_t {
    struct mpack_builder_page_t* next;
    size_t bytes_used;
} mpack_builder_page_t;

/**
 * Builds form a linked list of mpack_build_t, interleaved with their encoded
 * contents directly in the paged builder buffer.
 */
typedef struct mpack_build_t {
    //mpack_builder_page_t* page;
    struct mpack_build_t* parent;
    //struct mpack_build_t* next;

    size_t bytes; // number of bytes between this build and the next one
    uint32_t count; // number of elements (or key/value pairs) in this map/array
    mpack_type_t type;

    // depth of nested non-build compound elements within this
    // build.
    uint32_t nested_compound_elements;

    // indicates that a value still needs to be written for an already
    // written key. count is not incremented until both key and value are
    // written.
    bool key_needs_value;
} mpack_build_t;

/**
 * The builder state. This is stored within mpack_writer_t.
 */
typedef struct mpack_builder_t {
    mpack_build_t* current_build; // build which is accumulating elements
    mpack_build_t* latest_build; // build which is accumulating bytes
    mpack_builder_page_t* current_page;
    mpack_builder_page_t* pages;
    char* stash_buffer;
    char* stash_position;
    char* stash_end;
    #if MPACK_BUILDER_INTERNAL_STORAGE
    char internal[MPACK_BUILDER_INTERNAL_STORAGE_SIZE];
    #endif
} mpack_builder_t;
#endif

struct mpack_writer_t {
    #if MPACK_COMPATIBILITY
    mpack_version_t version;          /* Version of the MessagePack spec to write */
    #endif
    mpack_writer_flush_t flush;       /* Function to write bytes to the output stream */
    mpack_writer_error_t error_fn;    /* Function to call on error */
    mpack_writer_teardown_t teardown; /* Function to teardown the context on destroy */
    void* context;                    /* Context for writer callbacks */

    char* buffer;         /* Byte buffer */
    char* position;       /* Current position within the buffer */
    char* end;            /* The end of the buffer */
    mpack_error_t error;  /* Error state */

    #if MPACK_WRITE_TRACKING
    mpack_track_t track; /* Stack of map/array/str/bin/ext writes */
    #endif

    #ifdef MPACK_MALLOC
    /* Reserved. You can use this space to allocate a custom
     * context in order to reduce heap allocations. */
    void* reserved[2];
    #endif

    #if MPACK_BUILDER
    mpack_builder_t builder;
    #endif
};


#if MPACK_WRITE_TRACKING
void mpack_writer_track_push(mpack_writer_t* writer, mpack_type_t type, uint32_t count);
void mpack_writer_track_push_builder(mpack_writer_t* writer, mpack_type_t type);
void mpack_writer_track_pop(mpack_writer_t* writer, mpack_type_t type);
void mpack_writer_track_pop_builder(mpack_writer_t* writer, mpack_type_t type);
void mpack_writer_track_bytes(mpack_writer_t* writer, size_t count);
#else
MPACK_INLINE void mpack_writer_track_push(mpack_writer_t* writer, mpack_type_t type, uint32_t count) {
    MPACK_UNUSED(writer);
    MPACK_UNUSED(type);
    MPACK_UNUSED(count);
}
MPACK_INLINE void mpack_writer_track_push_builder(mpack_writer_t* writer, mpack_type_t type) {
    MPACK_UNUSED(writer);
    MPACK_UNUSED(type);
}
MPACK_INLINE void mpack_writer_track_pop(mpack_writer_t* writer, mpack_type_t type) {
    MPACK_UNUSED(writer);
    MPACK_UNUSED(type);
}
MPACK_INLINE void mpack_writer_track_pop_builder(mpack_writer_t* writer, mpack_type_t type) {
    MPACK_UNUSED(writer);
    MPACK_UNUSED(type);
}
MPACK_INLINE void mpack_writer_track_bytes(mpack_writer_t* writer, size_t count) {
    MPACK_UNUSED(writer);
    MPACK_UNUSED(count);
}
#endif

/** @endcond */

/**
 * @name Lifecycle Functions
 * @{
 */

/**
 * Initializes an MPack writer with the given buffer. The writer
 * does not assume ownership of the buffer.
 *
 * Trying to write past the end of the buffer will result in mpack_error_too_big
 * unless a flush function is set with mpack_writer_set_flush(). To use the data
 * without flushing, call mpack_writer_buffer_used() to determine the number of
 * bytes written.
 *
 * @param writer The MPack writer.
 * @param buffer The buffer into which to write MessagePack data.
 * @param size The size of the buffer.
 */
void mpack_writer_init(mpack_writer_t* writer, char* buffer, size_t size);

#ifdef MPACK_MALLOC
/**
 * Initializes an MPack writer using a growable buffer.
 *
 * The data is placed in the given data pointer if and when the writer
 * is destroyed without error. The data pointer is NULL during writing,
 * and will remain NULL if an error occurs.
 *
 * The allocated data must be freed with MPACK_FREE() (or simply free()
 * if MPack's allocator hasn't been customized.)
 *
 * @throws mpack_error_memory if the buffer fails to grow when
 * flushing.
 *
 * @param writer The MPack writer.
 * @param data Where to place the allocated data.
 * @param size Where to write the size of the data.
 */
void mpack_writer_init_growable(mpack_writer_t* writer, char** data, size_t* size);
#endif

/**
 * Initializes an MPack writer directly into an error state. Use this if you
 * are writing a wrapper to mpack_writer_init() which can fail its setup.
 */
void mpack_writer_init_error(mpack_writer_t* writer, mpack_error_t error);

#if MPACK_STDIO
/**
 * Initializes an MPack writer that writes to a file.
 *
 * @throws mpack_error_memory if allocation fails
 * @throws mpack_error_io if the file cannot be opened
 */
void mpack_writer_init_filename(mpack_writer_t* writer, const char* filename);

/**
 * Deprecated.
 *
 * \deprecated Renamed to mpack_writer_init_filename().
 */
MPACK_INLINE void mpack_writer_init_file(mpack_writer_t* writer, const char* filename) {
    mpack_writer_init_filename(writer, filename);
}

/**
 * Initializes an MPack writer that writes to a libc FILE. This can be used to
 * write to stdout or stderr, or to a file opened separately.
 *
 * @param writer The MPack writer.
 * @param stdfile The FILE.
 * @param close_when_done If true, fclose() will be called on the FILE when it
 *         is no longer needed. If false, the file will not be flushed or
 *         closed when writing is done.
 *
 * @note The writer is buffered. If you want to write other data to the FILE in
 *         between messages, you must flush it first.
 *
 * @see mpack_writer_flush_message
 */
void mpack_writer_init_stdfile(mpack_writer_t* writer, FILE* stdfile, bool close_when_done);
#endif

/** @cond */

#define mpack_writer_init_stack_line_ex(line, writer) \
    char mpack_buf_##line[MPACK_STACK_SIZE]; \
    mpack_writer_init(writer, mpack_buf_##line, sizeof(mpack_buf_##line))

#define mpack_writer_init_stack_line(line, writer) \
    mpack_writer_init_stack_line_ex(line, writer)

/*
 * Initializes an MPack writer using stack space as a buffer. A flush function
 * should be added to the writer to flush the buffer.
 *
 * This is currently undocumented since it's not entirely useful on its own.
 */

#define mpack_writer_init_stack(writer) \
    mpack_writer_init_stack_line(__LINE__, (writer))

/** @endcond */

/**
 * Cleans up the MPack writer, flushing and closing the underlying stream,
 * if any. Returns the final error state of the writer.
 *
 * No flushing is performed if the writer is in an error state. The attached
 * teardown function is called whether or not the writer is in an error state.
 *
 * This will assert in tracking mode if the writer is not in an error
 * state and has any unclosed compound types. If you want to cancel
 * writing in the middle of a document, you need to flag an error on
 * the writer before destroying it (such as mpack_error_data).
 *
 * Note that a writer may raise an error and call your error handler during
 * the final flush. It is safe to longjmp or throw out of this error handler,
 * but if you do, the writer will not be destroyed, and the teardown function
 * will not be called. You can still get the writer's error state, and you
 * must call @ref mpack_writer_destroy() again. (The second call is guaranteed
 * not to call your error handler again since the writer is already in an error
 * state.)
 *
 * @see mpack_writer_set_error_handler
 * @see mpack_writer_set_flush
 * @see mpack_writer_set_teardown
 * @see mpack_writer_flag_error
 * @see mpack_error_data
 */
mpack_error_t mpack_writer_destroy(mpack_writer_t* writer);

/**
 * @}
 */

/**
 * @name Configuration
 * @{
 */

#if MPACK_COMPATIBILITY
/**
 * Sets the version of the MessagePack spec that will be generated.
 *
 * This can be used to interface with older libraries that do not support
 * the newest MessagePack features (such as the @c str8 type.)
 *
 * @note This requires @ref MPACK_COMPATIBILITY.
 */
MPACK_INLINE void mpack_writer_set_version(mpack_writer_t* writer, mpack_version_t version) {
    writer->version = version;
}
#endif

/**
 * Sets the custom pointer to pass to the writer callbacks, such as flush
 * or teardown.
 *
 * @param writer The MPack writer.
 * @param context User data to pass to the writer callbacks.
 *
 * @see mpack_writer_context()
 */
MPACK_INLINE void mpack_writer_set_context(mpack_writer_t* writer, void* context) {
    writer->context = context;
}

/**
 * Returns the custom context for writer callbacks.
 *
 * @see mpack_writer_set_context
 * @see mpack_writer_set_flush
 */
MPACK_INLINE void* mpack_writer_context(mpack_writer_t* writer) {
    return writer->context;
}

/**
 * Sets the flush function to write out the data when the buffer is full.
 *
 * If no flush function is used, trying to write past the end of the
 * buffer will result in mpack_error_too_big.
 *
 * This should normally be used with mpack_writer_set_context() to register
 * a custom pointer to pass to the flush function.
 *
 * @param writer The MPack writer.
 * @param flush The function to write out data from the buffer.
 *
 * @see mpack_writer_context()
 */
void mpack_writer_set_flush(mpack_writer_t* writer, mpack_writer_flush_t flush);

/**
 * Sets the error function to call when an error is flagged on the writer.
 *
 * This should normally be used with mpack_writer_set_context() to register
 * a custom pointer to pass to the error function.
 *
 * See the definition of mpack_writer_error_t for more information about
 * what you can do from an error callback.
 *
 * @see mpack_writer_error_t
 * @param writer The MPack writer.
 * @param error_fn The function to call when an error is flagged on the writer.
 */
MPACK_INLINE void mpack_writer_set_error_handler(mpack_writer_t* writer, mpack_writer_error_t error_fn) {
    writer->error_fn = error_fn;
}

/**
 * Sets the teardown function to call when the writer is destroyed.
 *
 * This should normally be used with mpack_writer_set_context() to register
 * a custom pointer to pass to the teardown function.
 *
 * @param writer The MPack writer.
 * @param teardown The function to call when the writer is destroyed.
 */
MPACK_INLINE void mpack_writer_set_teardown(mpack_writer_t* writer, mpack_writer_teardown_t teardown) {
    writer->teardown = teardown;
}

/**
 * @}
 */

/**
 * @name Core Writer Functions
 * @{
 */

/**
 * Flushes any buffered data to the underlying stream.
 *
 * If the writer is connected to a socket and you are keeping it open,
 * you will want to call this after writing a message (or set of
 * messages) so that the data is actually sent.
 *
 * It is not necessary to call this if you are not keeping the writer
 * open afterwards. You can just call `mpack_writer_destroy()` and it
 * will flush before cleaning up.
 *
 * This will assert if no flush function is assigned to the writer.
 *
 * If write tracking is enabled, this will break and flag @ref
 * mpack_error_bug if the writer has any open compound types, ensuring
 * that no compound types are still open. This prevents a "missing
 * finish" bug from causing a never-ending message.
 */
void mpack_writer_flush_message(mpack_writer_t* writer);

/**
 * Returns the number of bytes currently stored in the buffer. This
 * may be less than the total number of bytes written if bytes have
 * been flushed to an underlying stream.
 */
MPACK_INLINE size_t mpack_writer_buffer_used(mpack_writer_t* writer) {
    return (size_t)(writer->position - writer->buffer);
}

/**
 * Returns the amount of space left in the buffer. This may be reset
 * after a write if bytes are flushed to an underlying stream.
 */
MPACK_INLINE size_t mpack_writer_buffer_left(mpack_writer_t* writer) {
    return (size_t)(writer->end - writer->position);
}

/**
 * Returns the (current) size of the buffer. This may change after a write if
 * the flush callback changes the buffer.
 */
MPACK_INLINE size_t mpack_writer_buffer_size(mpack_writer_t* writer) {
    return (size_t)(writer->end - writer->buffer);
}

/**
 * Places the writer in the given error state, calling the error callback if one
 * is set.
 *
 * This allows you to externally flag errors, for example if you are validating
 * data as you write it, or if you want to cancel writing in the middle of a
 * document. (The writer will assert if you try to destroy it without error and
 * with unclosed compound types. In this case you should flag mpack_error_data
 * before destroying it.)
 *
 * If the writer is already in an error state, this call is ignored and no
 * error callback is called.
 *
 * @see mpack_writer_destroy
 * @see mpack_error_data
 */
void mpack_writer_flag_error(mpack_writer_t* writer, mpack_error_t error);

/**
 * Queries the error state of the MPack writer.
 *
 * If a writer is in an error state, you should discard all data since the
 * last time the error flag was checked. The error flag cannot be cleared.
 */
MPACK_INLINE mpack_error_t mpack_writer_error(mpack_writer_t* writer) {
    return writer->error;
}

/**
 * Writes a MessagePack object header (an MPack Tag.)
 *
 * If the value is a map, array, string, binary or extension type, the
 * containing elements or bytes must be written separately and the
 * appropriate finish function must be called (as though one of the
 * mpack_start_*() functions was called.)
 *
 * @see mpack_write_bytes()
 * @see mpack_finish_map()
 * @see mpack_finish_array()
 * @see mpack_finish_str()
 * @see mpack_finish_bin()
 * @see mpack_finish_ext()
 * @see mpack_finish_type()
 */
void mpack_write_tag(mpack_writer_t* writer, mpack_tag_t tag);

/**
 * @}
 */

/**
 * @name Integers
 * @{
 */

/** Writes an 8-bit integer in the most efficient packing available. */
void mpack_write_i8(mpack_writer_t* writer, int8_t value);

/** Writes a 16-bit integer in the most efficient packing available. */
void mpack_write_i16(mpack_writer_t* writer, int16_t value);

/** Writes a 32-bit integer in the most efficient packing available. */
void mpack_write_i32(mpack_writer_t* writer, int32_t value);

/** Writes a 64-bit integer in the most efficient packing available. */
void mpack_write_i64(mpack_writer_t* writer, int64_t value);

/** Writes an integer in the most efficient packing available. */
MPACK_INLINE void mpack_write_int(mpack_writer_t* writer, int64_t value) {
    mpack_write_i64(writer, value);
}

/** Writes an 8-bit unsigned integer in the most efficient packing available. */
void mpack_write_u8(mpack_writer_t* writer, uint8_t value);

/** Writes an 16-bit unsigned integer in the most efficient packing available. */
void mpack_write_u16(mpack_writer_t* writer, uint16_t value);

/** Writes an 32-bit unsigned integer in the most efficient packing available. */
void mpack_write_u32(mpack_writer_t* writer, uint32_t value);

/** Writes an 64-bit unsigned integer in the most efficient packing available. */
void mpack_write_u64(mpack_writer_t* writer, uint64_t value);

/** Writes an unsigned integer in the most efficient packing available. */
MPACK_INLINE void mpack_write_uint(mpack_writer_t* writer, uint64_t value) {
    mpack_write_u64(writer, value);
}

/**
 * @}
 */

/**
 * @name Other Basic Types
 * @{
 */

#if MPACK_FLOAT
/** Writes a float. */
void mpack_write_float(mpack_writer_t* writer, float value);
#else
/** Writes a float from a raw uint32_t. */
void mpack_write_raw_float(mpack_writer_t* writer, uint32_t raw_value);
#endif

#if MPACK_DOUBLE
/** Writes a double. */
void mpack_write_double(mpack_writer_t* writer, double value);
#else
/** Writes a double from a raw uint64_t. */
void mpack_write_raw_double(mpack_writer_t* writer, uint64_t raw_value);
#endif

/** Writes a boolean. */
void mpack_write_bool(mpack_writer_t* writer, bool value);

/** Writes a boolean with value true. */
void mpack_write_true(mpack_writer_t* writer);

/** Writes a boolean with value false. */
void mpack_write_false(mpack_writer_t* writer);

/** Writes a nil. */
void mpack_write_nil(mpack_writer_t* writer);

/** Write a pre-encoded messagepack object */
void mpack_write_object_bytes(mpack_writer_t* writer, const char* data, size_t bytes);

#if MPACK_EXTENSIONS
/**
 * Writes a timestamp.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 *
 * @param writer The writer
 * @param seconds The (signed) number of seconds since 1970-01-01T00:00:00Z.
 * @param nanoseconds The additional number of nanoseconds from 0 to 999,999,999 inclusive.
 */
void mpack_write_timestamp(mpack_writer_t* writer, int64_t seconds, uint32_t nanoseconds);

/**
 * Writes a timestamp with the given number of seconds (and zero nanoseconds).
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 *
 * @param writer The writer
 * @param seconds The (signed) number of seconds since 1970-01-01T00:00:00Z.
 */
MPACK_INLINE void mpack_write_timestamp_seconds(mpack_writer_t* writer, int64_t seconds) {
    mpack_write_timestamp(writer, seconds, 0);
}

/**
 * Writes a timestamp.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 */
MPACK_INLINE void mpack_write_timestamp_struct(mpack_writer_t* writer, mpack_timestamp_t timestamp) {
    mpack_write_timestamp(writer, timestamp.seconds, timestamp.nanoseconds);
}
#endif

/**
 * @}
 */

/**
 * @name Map and Array Functions
 * @{
 */

/**
 * Opens an array.
 *
 * `count` elements must follow, and mpack_finish_array() must be called
 * when done.
 *
 * If you do not know the number of elements to be written ahead of time, call
 * mpack_build_array() instead.
 *
 * @see mpack_finish_array()
 * @see mpack_build_array() to count the number of elements automatically
 */
void mpack_start_array(mpack_writer_t* writer, uint32_t count);

/**
 * Opens a map.
 *
 * `count * 2` elements must follow, and mpack_finish_map() must be called
 * when done.
 *
 * If you do not know the number of elements to be written ahead of time, call
 * mpack_build_map() instead.
 *
 * Remember that while map elements in MessagePack are implicitly ordered,
 * they are not ordered in JSON. If you need elements to be read back
 * in the order they are written, consider use an array instead.
 *
 * @see mpack_finish_map()
 * @see mpack_build_map() to count the number of key/value pairs automatically
 */
void mpack_start_map(mpack_writer_t* writer, uint32_t count);

MPACK_INLINE void mpack_builder_compound_push(mpack_writer_t* writer) {
    MPACK_UNUSED(writer);

    #if MPACK_BUILDER
    mpack_build_t* build = writer->builder.current_build;
    if (build != NULL) {
        ++build->nested_compound_elements;
    }
    #endif
}

MPACK_INLINE void mpack_builder_compound_pop(mpack_writer_t* writer) {
    MPACK_UNUSED(writer);

    #if MPACK_BUILDER
    mpack_build_t* build = writer->builder.current_build;
    if (build != NULL) {
        mpack_assert(build->nested_compound_elements > 0);
        --build->nested_compound_elements;
    }
    #endif
}

/**
 * Finishes writing an array.
 *
 * This should be called only after a corresponding call to mpack_start_array()
 * and after the array contents are written.
 *
 * In debug mode (or if MPACK_WRITE_TRACKING is not 0), this will track writes
 * to ensure that the correct number of elements are written.
 *
 * @see mpack_start_array()
 */
MPACK_INLINE void mpack_finish_array(mpack_writer_t* writer) {
    mpack_writer_track_pop(writer, mpack_type_array);
    mpack_builder_compound_pop(writer);
}

/**
 * Finishes writing a map.
 *
 * This should be called only after a corresponding call to mpack_start_map()
 * and after the map contents are written.
 *
 * In debug mode (or if MPACK_WRITE_TRACKING is not 0), this will track writes
 * to ensure that the correct number of elements are written.
 *
 * @see mpack_start_map()
 */
MPACK_INLINE void mpack_finish_map(mpack_writer_t* writer) {
    mpack_writer_track_pop(writer, mpack_type_map);
    mpack_builder_compound_pop(writer);
}

/**
 * Starts building an array.
 *
 * Elements must follow, and mpack_complete_array() must be called when done. The
 * number of elements is determined automatically.
 *
 * If you know ahead of time the number of elements in the array, it is more
 * efficient to call mpack_start_array() instead, even if you are already
 * within another open build.
 *
 * Builder containers can be nested within normal (known size) containers and
 * vice versa. You can call mpack_build_array(), then mpack_start_array()
 * inside it, then mpack_build_array() inside that, and so forth.
 *
 * @see mpack_complete_array() to complete this array
 * @see mpack_start_array() if you already know the size of the array
 * @see mpack_build_map() for implementation details
 */
void mpack_build_array(struct mpack_writer_t* writer);

/**
 * Starts building a map.
 *
 * An even number of elements must follow, and mpack_complete_map() must be
 * called when done. The number of elements is determined automatically.
 *
 * If you know ahead of time the number of elements in the map, it is more
 * efficient to call mpack_start_map() instead, even if you are already within
 * another open build.
 *
 * Builder containers can be nested within normal (known size) containers and
 * vice versa. You can call mpack_build_map(), then mpack_start_map() inside
 * it, then mpack_build_map() inside that, and so forth.
 *
 * A writer in build mode diverts writes to a builder buffer that allocates as
 * needed. Once the last map or array being built is completed, the deferred
 * message is composed with computed array and map sizes into the writer.
 * Builder maps and arrays are encoded exactly the same as ordinary maps and
 * arrays in the final message.
 *
 * This indirect encoding is costly, as it incurs at least an extra copy of all
 * data written within a builder (but not additional copies for nested
 * builders.) Expect a speed penalty of half or more.
 *
 * A good strategy is to use this during early development when your messages
 * are constantly changing, and then closer to release when your message
 * formats have stabilized, replace all your build calls with start calls with
 * pre-computed sizes. Or don't, if you find the builder has little impact on
 * performance, because even with builders MPack is extremely fast.
 *
 * @note When an array or map starts being built, nothing will be flushed
 *       until it is completed. If you are building a large message that
 *       does not fit in the output stream, you won't get an error about it
 *       until everything is written.
 *
 * @see mpack_complete_map() to complete this map
 * @see mpack_start_map() if you already know the size of the map
 */
void mpack_build_map(struct mpack_writer_t* writer);

/**
 * Completes an array being built.
 *
 * @see mpack_build_array()
 */
void mpack_complete_array(struct mpack_writer_t* writer);

/**
 * Completes a map being built.
 *
 * @see mpack_build_map()
 */
void mpack_complete_map(struct mpack_writer_t* writer);

/**
 * @}
 */

/**
 * @name Data Helpers
 * @{
 */

/**
 * Writes a string.
 *
 * To stream a string in chunks, use mpack_start_str() instead.
 *
 * MPack does not care about the underlying encoding, but UTF-8 is highly
 * recommended, especially for compatibility with JSON. You should consider
 * calling mpack_write_utf8() instead, especially if you will be reading
 * it back as UTF-8.
 *
 * You should not call mpack_finish_str() after calling this; this
 * performs both start and finish.
 */
void mpack_write_str(mpack_writer_t* writer, const char* str, uint32_t length);

/**
 * Writes a string, ensuring that it is valid UTF-8.
 *
 * This does not accept any UTF-8 variant such as Modified UTF-8, CESU-8 or
 * WTF-8. Only pure UTF-8 is allowed.
 *
 * You should not call mpack_finish_str() after calling this; this
 * performs both start and finish.
 *
 * @throws mpack_error_invalid if the string is not valid UTF-8
 */
void mpack_write_utf8(mpack_writer_t* writer, const char* str, uint32_t length);

/**
 * Writes a null-terminated string. (The null-terminator is not written.)
 *
 * MPack does not care about the underlying encoding, but UTF-8 is highly
 * recommended, especially for compatibility with JSON. You should consider
 * calling mpack_write_utf8_cstr() instead, especially if you will be reading
 * it back as UTF-8.
 *
 * You should not call mpack_finish_str() after calling this; this
 * performs both start and finish.
 */
void mpack_write_cstr(mpack_writer_t* writer, const char* cstr);

/**
 * Writes a null-terminated string, or a nil node if the given cstr pointer
 * is NULL. (The null-terminator is not written.)
 *
 * MPack does not care about the underlying encoding, but UTF-8 is highly
 * recommended, especially for compatibility with JSON. You should consider
 * calling mpack_write_utf8_cstr_or_nil() instead, especially if you will
 * be reading it back as UTF-8.
 *
 * You should not call mpack_finish_str() after calling this; this
 * performs both start and finish.
 */
void mpack_write_cstr_or_nil(mpack_writer_t* writer, const char* cstr);

/**
 * Writes a null-terminated string, ensuring that it is valid UTF-8. (The
 * null-terminator is not written.)
 *
 * This does not accept any UTF-8 variant such as Modified UTF-8, CESU-8 or
 * WTF-8. Only pure UTF-8 is allowed.
 *
 * You should not call mpack_finish_str() after calling this; this
 * performs both start and finish.
 *
 * @throws mpack_error_invalid if the string is not valid UTF-8
 */
void mpack_write_utf8_cstr(mpack_writer_t* writer, const char* cstr);

/**
 * Writes a null-terminated string ensuring that it is valid UTF-8, or
 * writes nil if the given cstr pointer is NULL. (The null-terminator
 * is not written.)
 *
 * This does not accept any UTF-8 variant such as Modified UTF-8, CESU-8 or
 * WTF-8. Only pure UTF-8 is allowed.
 *
 * You should not call mpack_finish_str() after calling this; this
 * performs both start and finish.
 *
 * @throws mpack_error_invalid if the string is not valid UTF-8
 */
void mpack_write_utf8_cstr_or_nil(mpack_writer_t* writer, const char* cstr);

/**
 * Writes a binary blob.
 *
 * To stream a binary blob in chunks, use mpack_start_bin() instead.
 *
 * You should not call mpack_finish_bin() after calling this; this
 * performs both start and finish.
 */
void mpack_write_bin(mpack_writer_t* writer, const char* data, uint32_t count);

#if MPACK_EXTENSIONS
/**
 * Writes an extension type.
 *
 * To stream an extension blob in chunks, use mpack_start_ext() instead.
 *
 * Extension types [0, 127] are available for application-specific types. Extension
 * types [-128, -1] are reserved for future extensions of MessagePack.
 *
 * You should not call mpack_finish_ext() after calling this; this
 * performs both start and finish.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 */
void mpack_write_ext(mpack_writer_t* writer, int8_t exttype, const char* data, uint32_t count);
#endif

/**
 * @}
 */

/**
 * @name Chunked Data Functions
 * @{
 */

/**
 * Opens a string. `count` bytes should be written with calls to
 * mpack_write_bytes(), and mpack_finish_str() should be called
 * when done.
 *
 * To write an entire string at once, use mpack_write_str() or
 * mpack_write_cstr() instead.
 *
 * MPack does not care about the underlying encoding, but UTF-8 is highly
 * recommended, especially for compatibility with JSON.
 */
void mpack_start_str(mpack_writer_t* writer, uint32_t count);

/**
 * Opens a binary blob. `count` bytes should be written with calls to
 * mpack_write_bytes(), and mpack_finish_bin() should be called
 * when done.
 */
void mpack_start_bin(mpack_writer_t* writer, uint32_t count);

#if MPACK_EXTENSIONS
/**
 * Opens an extension type. `count` bytes should be written with calls
 * to mpack_write_bytes(), and mpack_finish_ext() should be called
 * when done.
 *
 * Extension types [0, 127] are available for application-specific types. Extension
 * types [-128, -1] are reserved for future extensions of MessagePack.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 */
void mpack_start_ext(mpack_writer_t* writer, int8_t exttype, uint32_t count);
#endif

/**
 * Writes a portion of bytes for a string, binary blob or extension type which
 * was opened by mpack_write_tag() or one of the mpack_start_*() functions.
 *
 * This can be called multiple times to write the data in chunks, as long as
 * the total amount of bytes written matches the count given when the compound
 * type was started.
 *
 * The corresponding mpack_finish_*() function must be called when done.
 *
 * To write an entire string, binary blob or extension type at
 * once, use one of the mpack_write_*() functions instead.
 *
 * @see mpack_write_tag()
 * @see mpack_start_str()
 * @see mpack_start_bin()
 * @see mpack_start_ext()
 * @see mpack_finish_str()
 * @see mpack_finish_bin()
 * @see mpack_finish_ext()
 * @see mpack_finish_type()
 */
void mpack_write_bytes(mpack_writer_t* writer, const char* data, size_t count);

/**
 * Finishes writing a string.
 *
 * This should be called only after a corresponding call to mpack_start_str()
 * and after the string bytes are written with mpack_write_bytes().
 *
 * This will track writes to ensure that the correct number of elements are written.
 *
 * @see mpack_start_str()
 * @see mpack_write_bytes()
 */
MPACK_INLINE void mpack_finish_str(mpack_writer_t* writer) {
    mpack_writer_track_pop(writer, mpack_type_str);
}

/**
 * Finishes writing a binary blob.
 *
 * This should be called only after a corresponding call to mpack_start_bin()
 * and after the binary bytes are written with mpack_write_bytes().
 *
 * This will track writes to ensure that the correct number of bytes are written.
 *
 * @see mpack_start_bin()
 * @see mpack_write_bytes()
 */
MPACK_INLINE void mpack_finish_bin(mpack_writer_t* writer) {
    mpack_writer_track_pop(writer, mpack_type_bin);
}

#if MPACK_EXTENSIONS
/**
 * Finishes writing an extended type binary data blob.
 *
 * This should be called only after a corresponding call to mpack_start_bin()
 * and after the binary bytes are written with mpack_write_bytes().
 *
 * This will track writes to ensure that the correct number of bytes are written.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 *
 * @see mpack_start_ext()
 * @see mpack_write_bytes()
 */
MPACK_INLINE void mpack_finish_ext(mpack_writer_t* writer) {
    mpack_writer_track_pop(writer, mpack_type_ext);
}
#endif

/**
 * Finishes writing the given compound type.
 *
 * This will track writes to ensure that the correct number of elements
 * or bytes are written.
 *
 * This can be called with the appropriate type instead the corresponding
 * mpack_finish_*() function if you want to finish a dynamic type.
 */
MPACK_INLINE void mpack_finish_type(mpack_writer_t* writer, mpack_type_t type) {
    mpack_writer_track_pop(writer, type);
}

/**
 * @}
 */

#if MPACK_HAS_GENERIC && !defined(__cplusplus)

/**
 * @name Type-Generic Writers
 * @{
 */

/**
 * @def mpack_write(writer, value)
 *
 * Type-generic writer for primitive types.
 *
 * The compiler will dispatch to an appropriate write function based
 * on the type of the @a value parameter.
 *
 * @note This requires C11 `_Generic` support. (A set of inline overloads
 * are used in C++ to provide the same functionality.)
 *
 * @warning In C11, the indentifiers `true`, `false` and `NULL` are
 * all of type `int`, not `bool` or `void*`! They will emit unexpected
 * types when passed uncast, so be careful when using them.
 */
#if MPACK_FLOAT
    #define MPACK_WRITE_GENERIC_FLOAT float: mpack_write_float,
#else
    #define MPACK_WRITE_GENERIC_FLOAT /*nothing*/
#endif
#if MPACK_DOUBLE
    #define MPACK_WRITE_GENERIC_DOUBLE double: mpack_write_double,
#else
    #define MPACK_WRITE_GENERIC_DOUBLE /*nothing*/
#endif
#define mpack_write(writer, value) \
    _Generic(((void)0, value),                      \
              int8_t: mpack_write_i8,               \
             int16_t: mpack_write_i16,              \
             int32_t: mpack_write_i32,              \
             int64_t: mpack_write_i64,              \
             uint8_t: mpack_write_u8,               \
            uint16_t: mpack_write_u16,              \
            uint32_t: mpack_write_u32,              \
            uint64_t: mpack_write_u64,              \
                bool: mpack_write_bool,             \
            MPACK_WRITE_GENERIC_FLOAT               \
            MPACK_WRITE_GENERIC_DOUBLE              \
              char *: mpack_write_cstr_or_nil,      \
        const char *: mpack_write_cstr_or_nil       \
    )(writer, value)

/**
 * @def mpack_write_kv(writer, key, value)
 *
 * Type-generic writer for key-value pairs of null-terminated string
 * keys and primitive values.
 *
 * @warning @a writer may be evaluated multiple times.
 *
 * @warning In C11, the indentifiers `true`, `false` and `NULL` are
 * all of type `int`, not `bool` or `void*`! They will emit unexpected
 * types when passed uncast, so be careful when using them.
 *
 * @param writer The writer.
 * @param key A null-terminated C string.
 * @param value A primitive type supported by mpack_write().
 */
#define mpack_write_kv(writer, key, value) do {     \
    mpack_write_cstr(writer, key);                  \
    mpack_write(writer, value);                     \
} while (0)

/**
 * @}
 */

#endif // MPACK_HAS_GENERIC && !defined(__cplusplus)

// The rest of this file contains C++ overloads, so we end extern "C" here.
MPACK_EXTERN_C_END

#if defined(__cplusplus) || defined(MPACK_DOXYGEN)

/**
 * @name C++ write overloads
 * @{
 */

/*
 * C++ generic writers for primitive values
 */

#ifdef MPACK_DOXYGEN
#undef mpack_write
#undef mpack_write_kv
#endif

MPACK_INLINE void mpack_write(mpack_writer_t* writer, int8_t value) {
    mpack_write_i8(writer, value);
}

MPACK_INLINE void mpack_write(mpack_writer_t* writer, int16_t value) {
    mpack_write_i16(writer, value);
}

MPACK_INLINE void mpack_write(mpack_writer_t* writer, int32_t value) {
    mpack_write_i32(writer, value);
}

MPACK_INLINE void mpack_write(mpack_writer_t* writer, int64_t value) {
    mpack_write_i64(writer, value);
}

MPACK_INLINE void mpack_write(mpack_writer_t* writer, uint8_t value) {
    mpack_write_u8(writer, value);
}

MPACK_INLINE void mpack_write(mpack_writer_t* writer, uint16_t value) {
    mpack_write_u16(writer, value);
}

MPACK_INLINE void mpack_write(mpack_writer_t* writer, uint32_t value) {
    mpack_write_u32(writer, value);
}

MPACK_INLINE void mpack_write(mpack_writer_t* writer, uint64_t value) {
    mpack_write_u64(writer, value);
}

MPACK_INLINE void mpack_write(mpack_writer_t* writer, bool value) {
    mpack_write_bool(writer, value);
}

MPACK_INLINE void mpack_write(mpack_writer_t* writer, float value) {
    mpack_write_float(writer, value);
}

MPACK_INLINE void mpack_write(mpack_writer_t* writer, double value) {
    mpack_write_double(writer, value);
}

MPACK_INLINE void mpack_write(mpack_writer_t* writer, char *value) {
    mpack_write_cstr_or_nil(writer, value);
}

MPACK_INLINE void mpack_write(mpack_writer_t* writer, const char *value) {
    mpack_write_cstr_or_nil(writer, value);
}

/* C++ generic write for key-value pairs */

MPACK_INLINE void mpack_write_kv(mpack_writer_t* writer, const char *key, int8_t value) {
    mpack_write_cstr(writer, key);
    mpack_write_i8(writer, value);
}

MPACK_INLINE void mpack_write_kv(mpack_writer_t* writer, const char *key, int16_t value) {
    mpack_write_cstr(writer, key);
    mpack_write_i16(writer, value);
}

MPACK_INLINE void mpack_write_kv(mpack_writer_t* writer, const char *key, int32_t value) {
    mpack_write_cstr(writer, key);
    mpack_write_i32(writer, value);
}

MPACK_INLINE void mpack_write_kv(mpack_writer_t* writer, const char *key, int64_t value) {
    mpack_write_cstr(writer, key);
    mpack_write_i64(writer, value);
}

MPACK_INLINE void mpack_write_kv(mpack_writer_t* writer, const char *key, uint8_t value) {
    mpack_write_cstr(writer, key);
    mpack_write_u8(writer, value);
}

MPACK_INLINE void mpack_write_kv(mpack_writer_t* writer, const char *key, uint16_t value) {
    mpack_write_cstr(writer, key);
    mpack_write_u16(writer, value);
}

MPACK_INLINE void mpack_write_kv(mpack_writer_t* writer, const char *key, uint32_t value) {
    mpack_write_cstr(writer, key);
    mpack_write_u32(writer, value);
}

MPACK_INLINE void mpack_write_kv(mpack_writer_t* writer, const char *key, uint64_t value) {
    mpack_write_cstr(writer, key);
    mpack_write_u64(writer, value);
}

MPACK_INLINE void mpack_write_kv(mpack_writer_t* writer, const char *key, bool value) {
    mpack_write_cstr(writer, key);
    mpack_write_bool(writer, value);
}

MPACK_INLINE void mpack_write_kv(mpack_writer_t* writer, const char *key, float value) {
    mpack_write_cstr(writer, key);
    mpack_write_float(writer, value);
}

MPACK_INLINE void mpack_write_kv(mpack_writer_t* writer, const char *key, double value) {
    mpack_write_cstr(writer, key);
    mpack_write_double(writer, value);
}

MPACK_INLINE void mpack_write_kv(mpack_writer_t* writer, const char *key, char *value) {
    mpack_write_cstr(writer, key);
    mpack_write_cstr_or_nil(writer, value);
}

MPACK_INLINE void mpack_write_kv(mpack_writer_t* writer, const char *key, const char *value) {
    mpack_write_cstr(writer, key);
    mpack_write_cstr_or_nil(writer, value);
}

/**
 * @}
 */

#endif /* __cplusplus */

/**
 * @}
 */

MPACK_SILENCE_WARNINGS_END

#endif // MPACK_WRITER

#endif

/* mpack/mpack-reader.h.h */

/**
 * @file
 *
 * Declares the core MPack Tag Reader.
 */

#ifndef MPACK_READER_H
#define MPACK_READER_H 1

/* #include "mpack-common.h" */

MPACK_SILENCE_WARNINGS_BEGIN
MPACK_EXTERN_C_BEGIN

#if MPACK_READER

#if MPACK_READ_TRACKING
struct mpack_track_t;
#endif

// The denominator to determine whether a read is a small
// fraction of the buffer size.
#define MPACK_READER_SMALL_FRACTION_DENOMINATOR 32

/**
 * @defgroup reader Reader API
 *
 * The MPack Reader API contains functions for imperatively reading dynamically
 * typed data from a MessagePack stream.
 *
 * See @ref docs/reader.md for examples.
 *
 * @note If you are not writing code for an embedded device (or otherwise do
 * not need maximum performance with minimal memory usage), you should not use
 * this. You probably want to use the @link node Node API@endlink instead.
 *
 * This forms the basis of the @link expect Expect API@endlink, which can be
 * used to interpret the stream of elements in expected types and value ranges.
 *
 * @{
 */

/**
 * @def MPACK_READER_MINIMUM_BUFFER_SIZE
 *
 * The minimum buffer size for a reader with a fill function.
 */
#define MPACK_READER_MINIMUM_BUFFER_SIZE 32

/**
 * A buffered MessagePack decoder.
 *
 * The decoder wraps an existing buffer and, optionally, a fill function.
 * This allows efficiently decoding data from existing memory buffers, files,
 * streams, etc.
 *
 * All read operations are synchronous; they will block until the
 * requested data is fully read, or an error occurs.
 *
 * This structure is opaque; its fields should not be accessed outside
 * of MPack.
 */
typedef struct mpack_reader_t mpack_reader_t;

/**
 * The MPack reader's fill function. It should fill the buffer with at
 * least one byte and at most the given @c count, returning the number
 * of bytes written to the buffer.
 *
 * In case of error, it should flag an appropriate error on the reader
 * (usually @ref mpack_error_io), or simply return zero. If zero is
 * returned, mpack_error_io is raised.
 *
 * @note When reading from a stream, you should only copy and return
 * the bytes that are immediately available. It is always safe to return
 * less than the requested count as long as some non-zero number of bytes
 * are read; if more bytes are needed, the read function will simply be
 * called again.
 *
 * @see mpack_reader_context()
 */
typedef size_t (*mpack_reader_fill_t)(mpack_reader_t* reader, char* buffer, size_t count);

/**
 * The MPack reader's skip function. It should discard the given number
 * of bytes from the source (for example by seeking forward.)
 *
 * In case of error, it should flag an appropriate error on the reader.
 *
 * @see mpack_reader_context()
 */
typedef void (*mpack_reader_skip_t)(mpack_reader_t* reader, size_t count);

/**
 * An error handler function to be called when an error is flagged on
 * the reader.
 *
 * The error handler will only be called once on the first error flagged;
 * any subsequent reads and errors are ignored, and the reader is
 * permanently in that error state.
 *
 * MPack is safe against non-local jumps out of error handler callbacks.
 * This means you are allowed to longjmp or throw an exception (in C++,
 * Objective-C, or with SEH) out of this callback.
 *
 * Bear in mind when using longjmp that local non-volatile variables that
 * have changed are undefined when setjmp() returns, so you can't put the
 * reader on the stack in the same activation frame as the setjmp without
 * declaring it volatile.
 *
 * You must still eventually destroy the reader. It is not destroyed
 * automatically when an error is flagged. It is safe to destroy the
 * reader within this error callback, but you will either need to perform
 * a non-local jump, or store something in your context to identify
 * that the reader is destroyed since any future accesses to it cause
 * undefined behavior.
 */
typedef void (*mpack_reader_error_t)(mpack_reader_t* reader, mpack_error_t error);

/**
 * A teardown function to be called when the reader is destroyed.
 */
typedef void (*mpack_reader_teardown_t)(mpack_reader_t* reader);

/* Hide internals from documentation */
/** @cond */

struct mpack_reader_t {
    void* context;                    /* Context for reader callbacks */
    mpack_reader_fill_t fill;         /* Function to read bytes into the buffer */
    mpack_reader_error_t error_fn;    /* Function to call on error */
    mpack_reader_teardown_t teardown; /* Function to teardown the context on destroy */
    mpack_reader_skip_t skip;         /* Function to skip bytes from the source */

    char* buffer;       /* Writeable byte buffer */
    size_t size;        /* Size of the buffer */

    const char* data;   /* Current data pointer (in the buffer, if it is used) */
    const char* end;    /* The end of available data (in the buffer, if it is used) */

    mpack_error_t error;  /* Error state */

    #if MPACK_READ_TRACKING
    mpack_track_t track; /* Stack of map/array/str/bin/ext reads */
    #endif
};

/** @endcond */

/**
 * @name Lifecycle Functions
 * @{
 */

/**
 * Initializes an MPack reader with the given buffer. The reader does
 * not assume ownership of the buffer, but the buffer must be writeable
 * if a fill function will be used to refill it.
 *
 * @param reader The MPack reader.
 * @param buffer The buffer with which to read MessagePack data.
 * @param size The size of the buffer.
 * @param count The number of bytes already in the buffer.
 */
void mpack_reader_init(mpack_reader_t* reader, char* buffer, size_t size, size_t count);

/**
 * Initializes an MPack reader directly into an error state. Use this if you
 * are writing a wrapper to mpack_reader_init() which can fail its setup.
 */
void mpack_reader_init_error(mpack_reader_t* reader, mpack_error_t error);

/**
 * Initializes an MPack reader to parse a pre-loaded contiguous chunk of data. The
 * reader does not assume ownership of the data.
 *
 * @param reader The MPack reader.
 * @param data The data to parse.
 * @param count The number of bytes pointed to by data.
 */
void mpack_reader_init_data(mpack_reader_t* reader, const char* data, size_t count);

#if MPACK_STDIO
/**
 * Initializes an MPack reader that reads from a file.
 *
 * The file will be automatically opened and closed by the reader.
 */
void mpack_reader_init_filename(mpack_reader_t* reader, const char* filename);

/**
 * Deprecated.
 *
 * \deprecated Renamed to mpack_reader_init_filename().
 */
MPACK_INLINE void mpack_reader_init_file(mpack_reader_t* reader, const char* filename) {
    mpack_reader_init_filename(reader, filename);
}

/**
 * Initializes an MPack reader that reads from a libc FILE. This can be used to
 * read from stdin, or from a file opened separately.
 *
 * @param reader The MPack reader.
 * @param stdfile The FILE.
 * @param close_when_done If true, fclose() will be called on the FILE when it
 *         is no longer needed. If false, the file will not be closed when
 *         reading is done.
 *
 * @warning The reader is buffered. It will read data in advance of parsing it,
 * and it may read more data than it parsed. See mpack_reader_remaining() to
 * access the extra data.
 */
void mpack_reader_init_stdfile(mpack_reader_t* reader, FILE* stdfile, bool close_when_done);
#endif

/**
 * @def mpack_reader_init_stack(reader)
 * @hideinitializer
 *
 * Initializes an MPack reader using stack space as a buffer. A fill function
 * should be added to the reader to fill the buffer.
 *
 * @see mpack_reader_set_fill
 */

/** @cond */
#define mpack_reader_init_stack_line_ex(line, reader) \
    char mpack_buf_##line[MPACK_STACK_SIZE]; \
    mpack_reader_init((reader), mpack_buf_##line, sizeof(mpack_buf_##line), 0)

#define mpack_reader_init_stack_line(line, reader) \
    mpack_reader_init_stack_line_ex(line, reader)
/** @endcond */

#define mpack_reader_init_stack(reader) \
    mpack_reader_init_stack_line(__LINE__, (reader))

/**
 * Cleans up the MPack reader, ensuring that all compound elements
 * have been completely read. Returns the final error state of the
 * reader.
 *
 * This will assert in tracking mode if the reader is not in an error
 * state and has any incomplete reads. If you want to cancel reading
 * in the middle of a document, you need to flag an error on the reader
 * before destroying it (such as mpack_error_data).
 *
 * @see mpack_read_tag()
 * @see mpack_reader_flag_error()
 * @see mpack_error_data
 */
mpack_error_t mpack_reader_destroy(mpack_reader_t* reader);

/**
 * @}
 */

/**
 * @name Callbacks
 * @{
 */

/**
 * Sets the custom pointer to pass to the reader callbacks, such as fill
 * or teardown.
 *
 * @param reader The MPack reader.
 * @param context User data to pass to the reader callbacks.
 *
 * @see mpack_reader_context()
 */
MPACK_INLINE void mpack_reader_set_context(mpack_reader_t* reader, void* context) {
    reader->context = context;
}

/**
 * Returns the custom context for reader callbacks.
 *
 * @see mpack_reader_set_context
 * @see mpack_reader_set_fill
 * @see mpack_reader_set_skip
 */
MPACK_INLINE void* mpack_reader_context(mpack_reader_t* reader) {
    return reader->context;
}

/**
 * Sets the fill function to refill the data buffer when it runs out of data.
 *
 * If no fill function is used, truncated MessagePack data results in
 * mpack_error_invalid (since the buffer is assumed to contain a
 * complete MessagePack object.)
 *
 * If a fill function is used, truncated MessagePack data usually
 * results in mpack_error_io (since the fill function fails to get
 * the missing data.)
 *
 * This should normally be used with mpack_reader_set_context() to register
 * a custom pointer to pass to the fill function.
 *
 * @param reader The MPack reader.
 * @param fill The function to fetch additional data into the buffer.
 */
void mpack_reader_set_fill(mpack_reader_t* reader, mpack_reader_fill_t fill);

/**
 * Sets the skip function to discard bytes from the source stream.
 *
 * It's not necessary to implement this function. If the stream is not
 * seekable, don't set a skip callback. The reader will fall back to
 * using the fill function instead.
 *
 * This should normally be used with mpack_reader_set_context() to register
 * a custom pointer to pass to the skip function.
 *
 * The skip function is ignored in size-optimized builds to reduce code
 * size. Data will be skipped with the fill function when necessary.
 *
 * @param reader The MPack reader.
 * @param skip The function to discard bytes from the source stream.
 */
void mpack_reader_set_skip(mpack_reader_t* reader, mpack_reader_skip_t skip);

/**
 * Sets the error function to call when an error is flagged on the reader.
 *
 * This should normally be used with mpack_reader_set_context() to register
 * a custom pointer to pass to the error function.
 *
 * See the definition of mpack_reader_error_t for more information about
 * what you can do from an error callback.
 *
 * @see mpack_reader_error_t
 * @param reader The MPack reader.
 * @param error_fn The function to call when an error is flagged on the reader.
 */
MPACK_INLINE void mpack_reader_set_error_handler(mpack_reader_t* reader, mpack_reader_error_t error_fn) {
    reader->error_fn = error_fn;
}

/**
 * Sets the teardown function to call when the reader is destroyed.
 *
 * This should normally be used with mpack_reader_set_context() to register
 * a custom pointer to pass to the teardown function.
 *
 * @param reader The MPack reader.
 * @param teardown The function to call when the reader is destroyed.
 */
MPACK_INLINE void mpack_reader_set_teardown(mpack_reader_t* reader, mpack_reader_teardown_t teardown) {
    reader->teardown = teardown;
}

/**
 * @}
 */

/**
 * @name Core Reader Functions
 * @{
 */

/**
 * Queries the error state of the MPack reader.
 *
 * If a reader is in an error state, you should discard all data since the
 * last time the error flag was checked. The error flag cannot be cleared.
 */
MPACK_INLINE mpack_error_t mpack_reader_error(mpack_reader_t* reader) {
    return reader->error;
}

/**
 * Places the reader in the given error state, calling the error callback if one
 * is set.
 *
 * This allows you to externally flag errors, for example if you are validating
 * data as you read it.
 *
 * If the reader is already in an error state, this call is ignored and no
 * error callback is called.
 */
void mpack_reader_flag_error(mpack_reader_t* reader, mpack_error_t error);

/**
 * Places the reader in the given error state if the given error is not mpack_ok,
 * returning the resulting error state of the reader.
 *
 * This allows you to externally flag errors, for example if you are validating
 * data as you read it.
 *
 * If the given error is mpack_ok or if the reader is already in an error state,
 * this call is ignored and the actual error state of the reader is returned.
 */
MPACK_INLINE mpack_error_t mpack_reader_flag_if_error(mpack_reader_t* reader, mpack_error_t error) {
    if (error != mpack_ok)
        mpack_reader_flag_error(reader, error);
    return mpack_reader_error(reader);
}

/**
 * Returns bytes left in the reader's buffer.
 *
 * If you are done reading MessagePack data but there is other interesting data
 * following it, the reader may have buffered too much data. The number of bytes
 * remaining in the buffer and a pointer to the position of those bytes can be
 * queried here.
 *
 * If you know the length of the MPack chunk beforehand, it's better to instead
 * have your fill function limit the data it reads so that the reader does not
 * have extra data. In this case you can simply check that this returns zero.
 *
 * Returns 0 if the reader is in an error state.
 *
 * @param reader The MPack reader from which to query remaining data.
 * @param data [out] A pointer to the remaining data, or NULL.
 * @return The number of bytes remaining in the buffer.
 */
size_t mpack_reader_remaining(mpack_reader_t* reader, const char** data);

/**
 * Reads a MessagePack object header (an MPack tag.)
 *
 * If an error occurs, the reader is placed in an error state and a
 * nil tag is returned. If the reader is already in an error state,
 * a nil tag is returned.
 *
 * If the type is compound (i.e. is a map, array, string, binary or
 * extension type), additional reads are required to get the contained
 * data, and the corresponding done function must be called when done.
 *
 * @note Maps in JSON are unordered, so it is recommended not to expect
 * a specific ordering for your map values in case your data is converted
 * to/from JSON.
 *
 * @see mpack_read_bytes()
 * @see mpack_done_array()
 * @see mpack_done_map()
 * @see mpack_done_str()
 * @see mpack_done_bin()
 * @see mpack_done_ext()
 */
mpack_tag_t mpack_read_tag(mpack_reader_t* reader);

/**
 * Parses the next MessagePack object header (an MPack tag) without
 * advancing the reader.
 *
 * If an error occurs, the reader is placed in an error state and a
 * nil tag is returned. If the reader is already in an error state,
 * a nil tag is returned.
 *
 * @note Maps in JSON are unordered, so it is recommended not to expect
 * a specific ordering for your map values in case your data is converted
 * to/from JSON.
 *
 * @see mpack_read_tag()
 * @see mpack_discard()
 */
mpack_tag_t mpack_peek_tag(mpack_reader_t* reader);

/**
 * @}
 */

/**
 * @name String and Data Functions
 * @{
 */

/**
 * Skips bytes from the underlying stream. This is used only to
 * skip the contents of a string, binary blob or extension object.
 */
void mpack_skip_bytes(mpack_reader_t* reader, size_t count);

/**
 * Reads bytes from a string, binary blob or extension object, copying
 * them into the given buffer.
 *
 * A str, bin or ext must have been opened by a call to mpack_read_tag()
 * which yielded one of these types, or by a call to an expect function
 * such as mpack_expect_str() or mpack_expect_bin().
 *
 * If an error occurs, the buffer contents are undefined.
 *
 * This can be called multiple times for a single str, bin or ext
 * to read the data in chunks. The total data read must add up
 * to the size of the object.
 *
 * @param reader The MPack reader
 * @param p The buffer in which to copy the bytes
 * @param count The number of bytes to read
 */
void mpack_read_bytes(mpack_reader_t* reader, char* p, size_t count);

/**
 * Reads bytes from a string, ensures that the string is valid UTF-8,
 * and copies the bytes into the given buffer.
 *
 * A string must have been opened by a call to mpack_read_tag() which
 * yielded a string, or by a call to an expect function such as
 * mpack_expect_str().
 *
 * The given byte count must match the complete size of the string as
 * returned by the tag or expect function. You must ensure that the
 * buffer fits the data.
 *
 * This does not accept any UTF-8 variant such as Modified UTF-8, CESU-8 or
 * WTF-8. Only pure UTF-8 is allowed.
 *
 * If an error occurs, the buffer contents are undefined.
 *
 * Unlike mpack_read_bytes(), this cannot be used to read the data in
 * chunks (since this might split a character's UTF-8 bytes, and the
 * reader does not keep track of the UTF-8 decoding state between reads.)
 *
 * @throws mpack_error_type if the string contains invalid UTF-8.
 */
void mpack_read_utf8(mpack_reader_t* reader, char* p, size_t byte_count);

/**
 * Reads bytes from a string, ensures that the string contains no NUL
 * bytes, copies the bytes into the given buffer and adds a null-terminator.
 *
 * A string must have been opened by a call to mpack_read_tag() which
 * yielded a string, or by a call to an expect function such as
 * mpack_expect_str().
 *
 * The given byte count must match the size of the string as returned
 * by the tag or expect function. The string will only be copied if
 * the buffer is large enough to store it.
 *
 * If an error occurs, the buffer will contain an empty string.
 *
 * @note If you know the object will be a string before reading it,
 * it is highly recommended to use mpack_expect_cstr() instead.
 * Alternatively you could use mpack_peek_tag() and call
 * mpack_expect_cstr() if it's a string.
 *
 * @throws mpack_error_too_big if the string plus null-terminator is larger than the given buffer size
 * @throws mpack_error_type if the string contains a null byte.
 *
 * @see mpack_peek_tag()
 * @see mpack_expect_cstr()
 * @see mpack_expect_utf8_cstr()
 */
void mpack_read_cstr(mpack_reader_t* reader, char* buf, size_t buffer_size, size_t byte_count);

/**
 * Reads bytes from a string, ensures that the string is valid UTF-8
 * with no NUL bytes, copies the bytes into the given buffer and adds a
 * null-terminator.
 *
 * A string must have been opened by a call to mpack_read_tag() which
 * yielded a string, or by a call to an expect function such as
 * mpack_expect_str().
 *
 * The given byte count must match the size of the string as returned
 * by the tag or expect function. The string will only be copied if
 * the buffer is large enough to store it.
 *
 * This does not accept any UTF-8 variant such as Modified UTF-8, CESU-8 or
 * WTF-8. Only pure UTF-8 is allowed, but without the NUL character, since
 * it cannot be represented in a null-terminated string.
 *
 * If an error occurs, the buffer will contain an empty string.
 *
 * @note If you know the object will be a string before reading it,
 * it is highly recommended to use mpack_expect_utf8_cstr() instead.
 * Alternatively you could use mpack_peek_tag() and call
 * mpack_expect_utf8_cstr() if it's a string.
 *
 * @throws mpack_error_too_big if the string plus null-terminator is larger than the given buffer size
 * @throws mpack_error_type if the string contains invalid UTF-8 or a null byte.
 *
 * @see mpack_peek_tag()
 * @see mpack_expect_utf8_cstr()
 */
void mpack_read_utf8_cstr(mpack_reader_t* reader, char* buf, size_t buffer_size, size_t byte_count);

#ifdef MPACK_MALLOC
/** @cond */
// This can optionally add a null-terminator, but it does not check
// whether the data contains null bytes. This must be done separately
// in a cstring read function (possibly as part of a UTF-8 check.)
char* mpack_read_bytes_alloc_impl(mpack_reader_t* reader, size_t count, bool null_terminated);
/** @endcond */

/**
 * Reads bytes from a string, binary blob or extension object, allocating
 * storage for them and returning the allocated pointer.
 *
 * The allocated string must be freed with MPACK_FREE() (or simply free()
 * if MPack's allocator hasn't been customized.)
 *
 * Returns NULL if any error occurs, or if count is zero.
 */
MPACK_INLINE char* mpack_read_bytes_alloc(mpack_reader_t* reader, size_t count) {
    return mpack_read_bytes_alloc_impl(reader, count, false);
}
#endif

/**
 * Reads bytes from a string, binary blob or extension object in-place in
 * the buffer. This can be used to avoid copying the data.
 *
 * A str, bin or ext must have been opened by a call to mpack_read_tag()
 * which yielded one of these types, or by a call to an expect function
 * such as mpack_expect_str() or mpack_expect_bin().
 *
 * If the bytes are from a string, the string is not null-terminated! Use
 * mpack_read_cstr() to copy the string into a buffer and add a null-terminator.
 *
 * The returned pointer is invalidated on the next read, or when the buffer
 * is destroyed.
 *
 * The reader will move data around in the buffer if needed to ensure that
 * the pointer can always be returned, so this should only be used if
 * count is very small compared to the buffer size. If you need to check
 * whether a small size is reasonable (for example you intend to handle small and
 * large sizes differently), you can call mpack_should_read_bytes_inplace().
 *
 * This can be called multiple times for a single str, bin or ext
 * to read the data in chunks. The total data read must add up
 * to the size of the object.
 *
 * NULL is returned if the reader is in an error state.
 *
 * @throws mpack_error_too_big if the requested size is larger than the buffer size
 *
 * @see mpack_should_read_bytes_inplace()
 */
const char* mpack_read_bytes_inplace(mpack_reader_t* reader, size_t count);

/**
 * Reads bytes from a string in-place in the buffer and ensures they are
 * valid UTF-8. This can be used to avoid copying the data.
 *
 * A string must have been opened by a call to mpack_read_tag() which
 * yielded a string, or by a call to an expect function such as
 * mpack_expect_str().
 *
 * The string is not null-terminated! Use mpack_read_utf8_cstr() to
 * copy the string into a buffer and add a null-terminator.
 *
 * The returned pointer is invalidated on the next read, or when the buffer
 * is destroyed.
 *
 * The reader will move data around in the buffer if needed to ensure that
 * the pointer can always be returned, so this should only be used if
 * count is very small compared to the buffer size. If you need to check
 * whether a small size is reasonable (for example you intend to handle small and
 * large sizes differently), you can call mpack_should_read_bytes_inplace().
 *
 * This does not accept any UTF-8 variant such as Modified UTF-8, CESU-8 or
 * WTF-8. Only pure UTF-8 is allowed.
 *
 * Unlike mpack_read_bytes_inplace(), this cannot be used to read the data in
 * chunks (since this might split a character's UTF-8 bytes, and the
 * reader does not keep track of the UTF-8 decoding state between reads.)
 *
 * NULL is returned if the reader is in an error state.
 *
 * @throws mpack_error_type if the string contains invalid UTF-8
 * @throws mpack_error_too_big if the requested size is larger than the buffer size
 *
 * @see mpack_should_read_bytes_inplace()
 */
const char* mpack_read_utf8_inplace(mpack_reader_t* reader, size_t count);

/**
 * Returns true if it's a good idea to read the given number of bytes
 * in-place.
 *
 * If the read will be larger than some small fraction of the buffer size,
 * this will return false to avoid shuffling too much data back and forth
 * in the buffer.
 *
 * Use this if you're expecting arbitrary size data, and you want to read
 * in-place for the best performance when possible but will fall back to
 * a normal read if the data is too large.
 *
 * @see mpack_read_bytes_inplace()
 */
MPACK_INLINE bool mpack_should_read_bytes_inplace(mpack_reader_t* reader, size_t count) {
    return (reader->size == 0 || count <= reader->size / MPACK_READER_SMALL_FRACTION_DENOMINATOR);
}

#if MPACK_EXTENSIONS
/**
 * Reads a timestamp contained in an ext object of the given size, closing the
 * ext type.
 *
 * An ext object of exttype @ref MPACK_EXTTYPE_TIMESTAMP must have been opened
 * by a call to e.g. mpack_read_tag() or mpack_expect_ext().
 *
 * You must NOT call mpack_done_ext() after calling this. A timestamp ext
 * object can only contain a single timestamp value, so this calls
 * mpack_done_ext() automatically.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 *
 * @throws mpack_error_invalid if the size is not one of the supported
 * timestamp sizes, or if the nanoseconds are out of range.
 */
mpack_timestamp_t mpack_read_timestamp(mpack_reader_t* reader, size_t size);
#endif

/**
 * @}
 */

/**
 * @name Core Reader Functions
 * @{
 */

#if MPACK_READ_TRACKING
/**
 * Finishes reading the given type.
 *
 * This will track reads to ensure that the correct number of elements
 * or bytes are read.
 */
void mpack_done_type(mpack_reader_t* reader, mpack_type_t type);
#else
MPACK_INLINE void mpack_done_type(mpack_reader_t* reader, mpack_type_t type) {
    MPACK_UNUSED(reader);
    MPACK_UNUSED(type);
}
#endif

/**
 * Finishes reading an array.
 *
 * This will track reads to ensure that the correct number of elements are read.
 */
MPACK_INLINE void mpack_done_array(mpack_reader_t* reader) {
    mpack_done_type(reader, mpack_type_array);
}

/**
 * @fn mpack_done_map(mpack_reader_t* reader)
 *
 * Finishes reading a map.
 *
 * This will track reads to ensure that the correct number of elements are read.
 */
MPACK_INLINE void mpack_done_map(mpack_reader_t* reader) {
    mpack_done_type(reader, mpack_type_map);
}

/**
 * @fn mpack_done_str(mpack_reader_t* reader)
 *
 * Finishes reading a string.
 *
 * This will track reads to ensure that the correct number of bytes are read.
 */
MPACK_INLINE void mpack_done_str(mpack_reader_t* reader) {
    mpack_done_type(reader, mpack_type_str);
}

/**
 * @fn mpack_done_bin(mpack_reader_t* reader)
 *
 * Finishes reading a binary data blob.
 *
 * This will track reads to ensure that the correct number of bytes are read.
 */
MPACK_INLINE void mpack_done_bin(mpack_reader_t* reader) {
    mpack_done_type(reader, mpack_type_bin);
}

#if MPACK_EXTENSIONS
/**
 * @fn mpack_done_ext(mpack_reader_t* reader)
 *
 * Finishes reading an extended type binary data blob.
 *
 * This will track reads to ensure that the correct number of bytes are read.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 */
MPACK_INLINE void mpack_done_ext(mpack_reader_t* reader) {
    mpack_done_type(reader, mpack_type_ext);
}
#endif

/**
 * Reads and discards the next object. This will read and discard all
 * contained data as well if it is a compound type.
 */
void mpack_discard(mpack_reader_t* reader);

/**
 * @}
 */

/** @cond */

#if MPACK_DEBUG && MPACK_STDIO
/**
 * @name Debugging Functions
 * @{
 */
/*
 * Converts a blob of MessagePack to a pseudo-JSON string for debugging
 * purposes, placing the result in the given buffer with a null-terminator.
 *
 * If the buffer does not have enough space, the result will be truncated (but
 * it is guaranteed to be null-terminated.)
 *
 * This is only available in debug mode, and only if stdio is available (since
 * it uses snprintf().) It's strictly for debugging purposes.
 */
void mpack_print_data_to_buffer(const char* data, size_t data_size, char* buffer, size_t buffer_size);

/*
 * Converts a node to pseudo-JSON for debugging purposes, calling the given
 * callback as many times as is necessary to output the character data.
 *
 * No null-terminator or trailing newline will be written.
 *
 * This is only available in debug mode, and only if stdio is available (since
 * it uses snprintf().) It's strictly for debugging purposes.
 */
void mpack_print_data_to_callback(const char* data, size_t size, mpack_print_callback_t callback, void* context);

/*
 * Converts a blob of MessagePack to pseudo-JSON for debugging purposes
 * and pretty-prints it to the given file.
 */
void mpack_print_data_to_file(const char* data, size_t len, FILE* file);

/*
 * Converts a blob of MessagePack to pseudo-JSON for debugging purposes
 * and pretty-prints it to stdout.
 */
MPACK_INLINE void mpack_print_data_to_stdout(const char* data, size_t len) {
    mpack_print_data_to_file(data, len, stdout);
}

/*
 * Converts the MessagePack contained in the given `FILE*` to pseudo-JSON for
 * debugging purposes, calling the given callback as many times as is necessary
 * to output the character data.
 */
void mpack_print_stdfile_to_callback(FILE* file, mpack_print_callback_t callback, void* context);

/*
 * Deprecated.
 *
 * \deprecated Renamed to mpack_print_data_to_stdout().
 */
MPACK_INLINE void mpack_print(const char* data, size_t len) {
    mpack_print_data_to_stdout(data, len);
}

/**
 * @}
 */
#endif

/** @endcond */

/**
 * @}
 */



#if MPACK_INTERNAL

bool mpack_reader_ensure_straddle(mpack_reader_t* reader, size_t count);

/*
 * Ensures there are at least @c count bytes left in the
 * data, raising an error and returning false if more
 * data cannot be made available.
 */
MPACK_INLINE bool mpack_reader_ensure(mpack_reader_t* reader, size_t count) {
    mpack_assert(count != 0, "cannot ensure zero bytes!");
    mpack_assert(reader->error == mpack_ok, "reader cannot be in an error state!");

    if (count <= (size_t)(reader->end - reader->data))
        return true;
    return mpack_reader_ensure_straddle(reader, count);
}

void mpack_read_native_straddle(mpack_reader_t* reader, char* p, size_t count);

// Reads count bytes into p, deferring to mpack_read_native_straddle() if more
// bytes are needed than are available in the buffer.
MPACK_INLINE void mpack_read_native(mpack_reader_t* reader, char* p, size_t count) {
    mpack_assert(count == 0 || p != NULL, "data pointer for %i bytes is NULL", (int)count);

    if (count > (size_t)(reader->end - reader->data)) {
        mpack_read_native_straddle(reader, p, count);
    } else {
        mpack_memcpy(p, reader->data, count);
        reader->data += count;
    }
}

#if MPACK_READ_TRACKING
#define MPACK_READER_TRACK(reader, error_expr) \
    (((reader)->error == mpack_ok) ? mpack_reader_flag_if_error((reader), (error_expr)) : (reader)->error)
#else
#define MPACK_READER_TRACK(reader, error_expr) (MPACK_UNUSED(reader), mpack_ok)
#endif

MPACK_INLINE mpack_error_t mpack_reader_track_element(mpack_reader_t* reader) {
    return MPACK_READER_TRACK(reader, mpack_track_element(&reader->track, true));
}

MPACK_INLINE mpack_error_t mpack_reader_track_peek_element(mpack_reader_t* reader) {
    return MPACK_READER_TRACK(reader, mpack_track_peek_element(&reader->track, true));
}

MPACK_INLINE mpack_error_t mpack_reader_track_bytes(mpack_reader_t* reader, size_t count) {
    MPACK_UNUSED(count);
    return MPACK_READER_TRACK(reader, mpack_track_bytes(&reader->track, true, count));
}

MPACK_INLINE mpack_error_t mpack_reader_track_str_bytes_all(mpack_reader_t* reader, size_t count) {
    MPACK_UNUSED(count);
    return MPACK_READER_TRACK(reader, mpack_track_str_bytes_all(&reader->track, true, count));
}

#endif



#endif

MPACK_EXTERN_C_END
MPACK_SILENCE_WARNINGS_END

#endif


/* mpack/mpack-expect.h.h */

/**
 * @file
 *
 * Declares the MPack static Expect API.
 */

#ifndef MPACK_EXPECT_H
#define MPACK_EXPECT_H 1

/* #include "mpack-reader.h" */

MPACK_SILENCE_WARNINGS_BEGIN
MPACK_EXTERN_C_BEGIN

#if MPACK_EXPECT

#if !MPACK_READER
#error "MPACK_EXPECT requires MPACK_READER."
#endif

/**
 * @defgroup expect Expect API
 *
 * The MPack Expect API allows you to easily read MessagePack data when you
 * expect it to follow a predefined schema.
 *
 * @note If you are not writing code for an embedded device (or otherwise do
 * not need maximum performance with minimal memory usage), you should not use
 * this. You probably want to use the @link node Node API@endlink instead.
 *
 * See @ref docs/expect.md for examples.
 *
 * The main purpose of the Expect API is convenience, so the API is lax. It
 * automatically converts between similar types where there is no loss of
 * precision.
 *
 * When using any of the expect functions, if the type or value of what was
 * read does not match what is expected, @ref mpack_error_type is raised.
 *
 * @{
 */

/**
 * @name Basic Number Functions
 * @{
 */

/**
 * Reads an 8-bit unsigned integer.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in an 8-bit unsigned int.
 *
 * Returns zero if an error occurs.
 */
uint8_t mpack_expect_u8(mpack_reader_t* reader);

/**
 * Reads a 16-bit unsigned integer.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 16-bit unsigned int.
 *
 * Returns zero if an error occurs.
 */
uint16_t mpack_expect_u16(mpack_reader_t* reader);

/**
 * Reads a 32-bit unsigned integer.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 32-bit unsigned int.
 *
 * Returns zero if an error occurs.
 */
uint32_t mpack_expect_u32(mpack_reader_t* reader);

/**
 * Reads a 64-bit unsigned integer.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 64-bit unsigned int.
 *
 * Returns zero if an error occurs.
 */
uint64_t mpack_expect_u64(mpack_reader_t* reader);

/**
 * Reads an 8-bit signed integer.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in an 8-bit signed int.
 *
 * Returns zero if an error occurs.
 */
int8_t mpack_expect_i8(mpack_reader_t* reader);

/**
 * Reads a 16-bit signed integer.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 16-bit signed int.
 *
 * Returns zero if an error occurs.
 */
int16_t mpack_expect_i16(mpack_reader_t* reader);

/**
 * Reads a 32-bit signed integer.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 32-bit signed int.
 *
 * Returns zero if an error occurs.
 */
int32_t mpack_expect_i32(mpack_reader_t* reader);

/**
 * Reads a 64-bit signed integer.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 64-bit signed int.
 *
 * Returns zero if an error occurs.
 */
int64_t mpack_expect_i64(mpack_reader_t* reader);

#if MPACK_FLOAT
/**
 * Reads a number, returning the value as a float. The underlying value can be an
 * integer, float or double; the value is converted to a float.
 *
 * @note Reading a double or a large integer with this function can incur a
 * loss of precision.
 *
 * @throws mpack_error_type if the underlying value is not a float, double or integer.
 */
float mpack_expect_float(mpack_reader_t* reader);
#endif

#if MPACK_DOUBLE
/**
 * Reads a number, returning the value as a double. The underlying value can be an
 * integer, float or double; the value is converted to a double.
 *
 * @note Reading a very large integer with this function can incur a
 * loss of precision.
 *
 * @throws mpack_error_type if the underlying value is not a float, double or integer.
 */
double mpack_expect_double(mpack_reader_t* reader);
#endif

#if MPACK_FLOAT
/**
 * Reads a float. The underlying value must be a float, not a double or an integer.
 * This ensures no loss of precision can occur.
 *
 * @throws mpack_error_type if the underlying value is not a float.
 */
float mpack_expect_float_strict(mpack_reader_t* reader);
#endif

#if MPACK_DOUBLE
/**
 * Reads a double. The underlying value must be a float or double, not an integer.
 * This ensures no loss of precision can occur.
 *
 * @throws mpack_error_type if the underlying value is not a float or double.
 */
double mpack_expect_double_strict(mpack_reader_t* reader);
#endif

#if !MPACK_FLOAT
/**
 * Reads a float as a raw uint32_t. The underlying value must be a float, not a
 * double or an integer.
 *
 * @throws mpack_error_type if the underlying value is not a float.
 */
uint32_t mpack_expect_raw_float(mpack_reader_t* reader);
#endif

#if !MPACK_DOUBLE
/**
 * Reads a double as a raw uint64_t. The underlying value must be a double, not a
 * float or an integer.
 *
 * @throws mpack_error_type if the underlying value is not a double.
 */
uint64_t mpack_expect_raw_double(mpack_reader_t* reader);
#endif

/**
 * @}
 */

/**
 * @name Ranged Number Functions
 * @{
 */

/**
 * Reads an 8-bit unsigned integer, ensuring that it falls within the given range.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in an 8-bit unsigned int.
 *
 * Returns min_value if an error occurs.
 */
uint8_t mpack_expect_u8_range(mpack_reader_t* reader, uint8_t min_value, uint8_t max_value);

/**
 * Reads a 16-bit unsigned integer, ensuring that it falls within the given range.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 16-bit unsigned int.
 *
 * Returns min_value if an error occurs.
 */
uint16_t mpack_expect_u16_range(mpack_reader_t* reader, uint16_t min_value, uint16_t max_value);

/**
 * Reads a 32-bit unsigned integer, ensuring that it falls within the given range.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 32-bit unsigned int.
 *
 * Returns min_value if an error occurs.
 */
uint32_t mpack_expect_u32_range(mpack_reader_t* reader, uint32_t min_value, uint32_t max_value);

/**
 * Reads a 64-bit unsigned integer, ensuring that it falls within the given range.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 64-bit unsigned int.
 *
 * Returns min_value if an error occurs.
 */
uint64_t mpack_expect_u64_range(mpack_reader_t* reader, uint64_t min_value, uint64_t max_value);

/**
 * Reads an unsigned integer, ensuring that it falls within the given range.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in an unsigned int.
 *
 * Returns min_value if an error occurs.
 */
MPACK_INLINE unsigned int mpack_expect_uint_range(mpack_reader_t* reader, unsigned int min_value, unsigned int max_value) {
    // This should be true at compile-time, so this just wraps the 32-bit
    // function. We fallback to 64-bit if for some reason sizeof(int) isn't 4.
    if (sizeof(unsigned int) == 4)
        return (unsigned int)mpack_expect_u32_range(reader, (uint32_t)min_value, (uint32_t)max_value);
    return (unsigned int)mpack_expect_u64_range(reader, min_value, max_value);
}

/**
 * Reads an 8-bit unsigned integer, ensuring that it is at most @a max_value.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in an 8-bit unsigned int.
 *
 * Returns 0 if an error occurs.
 */
MPACK_INLINE uint8_t mpack_expect_u8_max(mpack_reader_t* reader, uint8_t max_value) {
    return mpack_expect_u8_range(reader, 0, max_value);
}

/**
 * Reads a 16-bit unsigned integer, ensuring that it is at most @a max_value.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 16-bit unsigned int.
 *
 * Returns 0 if an error occurs.
 */
MPACK_INLINE uint16_t mpack_expect_u16_max(mpack_reader_t* reader, uint16_t max_value) {
    return mpack_expect_u16_range(reader, 0, max_value);
}

/**
 * Reads a 32-bit unsigned integer, ensuring that it is at most @a max_value.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 32-bit unsigned int.
 *
 * Returns 0 if an error occurs.
 */
MPACK_INLINE uint32_t mpack_expect_u32_max(mpack_reader_t* reader, uint32_t max_value) {
    return mpack_expect_u32_range(reader, 0, max_value);
}

/**
 * Reads a 64-bit unsigned integer, ensuring that it is at most @a max_value.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 64-bit unsigned int.
 *
 * Returns 0 if an error occurs.
 */
MPACK_INLINE uint64_t mpack_expect_u64_max(mpack_reader_t* reader, uint64_t max_value) {
    return mpack_expect_u64_range(reader, 0, max_value);
}

/**
 * Reads an unsigned integer, ensuring that it is at most @a max_value.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in an unsigned int.
 *
 * Returns 0 if an error occurs.
 */
MPACK_INLINE unsigned int mpack_expect_uint_max(mpack_reader_t* reader, unsigned int max_value) {
    return mpack_expect_uint_range(reader, 0, max_value);
}

/**
 * Reads an 8-bit signed integer, ensuring that it falls within the given range.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in an 8-bit signed int.
 *
 * Returns min_value if an error occurs.
 */
int8_t mpack_expect_i8_range(mpack_reader_t* reader, int8_t min_value, int8_t max_value);

/**
 * Reads a 16-bit signed integer, ensuring that it falls within the given range.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 16-bit signed int.
 *
 * Returns min_value if an error occurs.
 */
int16_t mpack_expect_i16_range(mpack_reader_t* reader, int16_t min_value, int16_t max_value);

/**
 * Reads a 32-bit signed integer, ensuring that it falls within the given range.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 32-bit signed int.
 *
 * Returns min_value if an error occurs.
 */
int32_t mpack_expect_i32_range(mpack_reader_t* reader, int32_t min_value, int32_t max_value);

/**
 * Reads a 64-bit signed integer, ensuring that it falls within the given range.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 64-bit signed int.
 *
 * Returns min_value if an error occurs.
 */
int64_t mpack_expect_i64_range(mpack_reader_t* reader, int64_t min_value, int64_t max_value);

/**
 * Reads a signed integer, ensuring that it falls within the given range.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a signed int.
 *
 * Returns min_value if an error occurs.
 */
MPACK_INLINE int mpack_expect_int_range(mpack_reader_t* reader, int min_value, int max_value) {
    // This should be true at compile-time, so this just wraps the 32-bit
    // function. We fallback to 64-bit if for some reason sizeof(int) isn't 4.
    if (sizeof(int) == 4)
        return (int)mpack_expect_i32_range(reader, (int32_t)min_value, (int32_t)max_value);
    return (int)mpack_expect_i64_range(reader, min_value, max_value);
}

/**
 * Reads an 8-bit signed integer, ensuring that it is at least zero and at
 * most @a max_value.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in an 8-bit signed int.
 *
 * Returns 0 if an error occurs.
 */
MPACK_INLINE int8_t mpack_expect_i8_max(mpack_reader_t* reader, int8_t max_value) {
    return mpack_expect_i8_range(reader, 0, max_value);
}

/**
 * Reads a 16-bit signed integer, ensuring that it is at least zero and at
 * most @a max_value.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 16-bit signed int.
 *
 * Returns 0 if an error occurs.
 */
MPACK_INLINE int16_t mpack_expect_i16_max(mpack_reader_t* reader, int16_t max_value) {
    return mpack_expect_i16_range(reader, 0, max_value);
}

/**
 * Reads a 32-bit signed integer, ensuring that it is at least zero and at
 * most @a max_value.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 32-bit signed int.
 *
 * Returns 0 if an error occurs.
 */
MPACK_INLINE int32_t mpack_expect_i32_max(mpack_reader_t* reader, int32_t max_value) {
    return mpack_expect_i32_range(reader, 0, max_value);
}

/**
 * Reads a 64-bit signed integer, ensuring that it is at least zero and at
 * most @a max_value.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a 64-bit signed int.
 *
 * Returns 0 if an error occurs.
 */
MPACK_INLINE int64_t mpack_expect_i64_max(mpack_reader_t* reader, int64_t max_value) {
    return mpack_expect_i64_range(reader, 0, max_value);
}

/**
 * Reads an int, ensuring that it is at least zero and at most @a max_value.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a signed int.
 *
 * Returns 0 if an error occurs.
 */
MPACK_INLINE int mpack_expect_int_max(mpack_reader_t* reader, int max_value) {
    return mpack_expect_int_range(reader, 0, max_value);
}

#if MPACK_FLOAT
/**
 * Reads a number, ensuring that it falls within the given range and returning
 * the value as a float. The underlying value can be an integer, float or
 * double; the value is converted to a float.
 *
 * @note Reading a double or a large integer with this function can incur a
 * loss of precision.
 *
 * @throws mpack_error_type if the underlying value is not a float, double or integer.
 */
float mpack_expect_float_range(mpack_reader_t* reader, float min_value, float max_value);
#endif

#if MPACK_DOUBLE
/**
 * Reads a number, ensuring that it falls within the given range and returning
 * the value as a double. The underlying value can be an integer, float or
 * double; the value is converted to a double.
 *
 * @note Reading a very large integer with this function can incur a
 * loss of precision.
 *
 * @throws mpack_error_type if the underlying value is not a float, double or integer.
 */
double mpack_expect_double_range(mpack_reader_t* reader, double min_value, double max_value);
#endif

/**
 * @}
 */



// These are additional Basic Number functions that wrap inline range functions.

/**
 * @name Basic Number Functions
 * @{
 */

/**
 * Reads an unsigned int.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in an unsigned int.
 *
 * Returns zero if an error occurs.
 */
MPACK_INLINE unsigned int mpack_expect_uint(mpack_reader_t* reader) {

    // This should be true at compile-time, so this just wraps the 32-bit function.
    if (sizeof(unsigned int) == 4)
        return (unsigned int)mpack_expect_u32(reader);

    // Otherwise we wrap the max function to ensure it fits.
    return (unsigned int)mpack_expect_u64_max(reader, MPACK_UINT_MAX);

}

/**
 * Reads a signed int.
 *
 * The underlying type may be an integer type of any size and signedness,
 * as long as the value can be represented in a signed int.
 *
 * Returns zero if an error occurs.
 */
MPACK_INLINE int mpack_expect_int(mpack_reader_t* reader) {

    // This should be true at compile-time, so this just wraps the 32-bit function.
    if (sizeof(int) == 4)
        return (int)mpack_expect_i32(reader);

    // Otherwise we wrap the range function to ensure it fits.
    return (int)mpack_expect_i64_range(reader, MPACK_INT_MIN, MPACK_INT_MAX);

}

/**
 * @}
 */



/**
 * @name Matching Number Functions
 * @{
 */

/**
 * Reads an unsigned integer, ensuring that it exactly matches the given value.
 *
 * mpack_error_type is raised if the value is not representable as an unsigned
 * integer or if it does not exactly match the given value.
 */
void mpack_expect_uint_match(mpack_reader_t* reader, uint64_t value);

/**
 * Reads a signed integer, ensuring that it exactly matches the given value.
 *
 * mpack_error_type is raised if the value is not representable as a signed
 * integer or if it does not exactly match the given value.
 */
void mpack_expect_int_match(mpack_reader_t* reader, int64_t value);

/**
 * @}
 */

/**
 * @name Other Basic Types
 * @{
 */

/**
 * Reads a nil, raising @ref mpack_error_type if the value is not nil.
 */
void mpack_expect_nil(mpack_reader_t* reader);

/**
 * Reads a boolean.
 *
 * @note Integers will raise mpack_error_type; the value must be strictly a boolean.
 */
bool mpack_expect_bool(mpack_reader_t* reader);

/**
 * Reads a boolean, raising @ref mpack_error_type if its value is not @c true.
 */
void mpack_expect_true(mpack_reader_t* reader);

/**
 * Reads a boolean, raising @ref mpack_error_type if its value is not @c false.
 */
void mpack_expect_false(mpack_reader_t* reader);

/**
 * @}
 */

/**
 * @name Extension Functions
 * @{
 */

#if MPACK_EXTENSIONS
/**
 * Reads a timestamp.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 */
mpack_timestamp_t mpack_expect_timestamp(mpack_reader_t* reader);

/**
 * Reads a timestamp in seconds, truncating the nanoseconds (if any).
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 */
int64_t mpack_expect_timestamp_truncate(mpack_reader_t* reader);
#endif

/**
 * @}
 */

/**
 * @name Compound Types
 * @{
 */

/**
 * Reads the start of a map, returning its element count.
 *
 * A number of values follow equal to twice the element count of the map,
 * alternating between keys and values. @ref mpack_done_map() must be called
 * once all elements have been read.
 *
 * @note Maps in JSON are unordered, so it is recommended not to expect
 * a specific ordering for your map values in case your data is converted
 * to/from JSON.
 *
 * @warning This call is dangerous! It does not have a size limit, and it
 * does not have any way of checking whether there is enough data in the
 * message (since the data could be coming from a stream.) When looping
 * through the map's contents, you must check for errors on each iteration
 * of the loop. Otherwise an attacker could craft a message declaring a map
 * of a billion elements which would throw your parsing code into an
 * infinite loop! You should strongly consider using mpack_expect_map_max()
 * with a safe maximum size instead.
 *
 * @throws mpack_error_type if the value is not a map.
 */
uint32_t mpack_expect_map(mpack_reader_t* reader);

/**
 * Reads the start of a map with a number of elements in the given range, returning
 * its element count.
 *
 * A number of values follow equal to twice the element count of the map,
 * alternating between keys and values. @ref mpack_done_map() must be called
 * once all elements have been read.
 *
 * @note Maps in JSON are unordered, so it is recommended not to expect
 * a specific ordering for your map values in case your data is converted
 * to/from JSON.
 *
 * min_count is returned if an error occurs.
 *
 * @throws mpack_error_type if the value is not a map or if its size does
 * not fall within the given range.
 */
uint32_t mpack_expect_map_range(mpack_reader_t* reader, uint32_t min_count, uint32_t max_count);

/**
 * Reads the start of a map with a number of elements at most @a max_count,
 * returning its element count.
 *
 * A number of values follow equal to twice the element count of the map,
 * alternating between keys and values. @ref mpack_done_map() must be called
 * once all elements have been read.
 *
 * @note Maps in JSON are unordered, so it is recommended not to expect
 * a specific ordering for your map values in case your data is converted
 * to/from JSON.
 *
 * Zero is returned if an error occurs.
 *
 * @throws mpack_error_type if the value is not a map or if its size is
 * greater than max_count.
 */
MPACK_INLINE uint32_t mpack_expect_map_max(mpack_reader_t* reader, uint32_t max_count) {
    return mpack_expect_map_range(reader, 0, max_count);
}

/**
 * Reads the start of a map of the exact size given.
 *
 * A number of values follow equal to twice the element count of the map,
 * alternating between keys and values. @ref mpack_done_map() must be called
 * once all elements have been read.
 *
 * @note Maps in JSON are unordered, so it is recommended not to expect
 * a specific ordering for your map values in case your data is converted
 * to/from JSON.
 *
 * @throws mpack_error_type if the value is not a map or if its size
 * does not match the given count.
 */
void mpack_expect_map_match(mpack_reader_t* reader, uint32_t count);

/**
 * Reads a nil node or the start of a map, returning whether a map was
 * read and placing its number of key/value pairs in count.
 *
 * If a map was read, a number of values follow equal to twice the element count
 * of the map, alternating between keys and values. @ref mpack_done_map() should
 * also be called once all elements have been read (only if a map was read.)
 *
 * @note Maps in JSON are unordered, so it is recommended not to expect
 * a specific ordering for your map values in case your data is converted
 * to/from JSON.
 *
 * @warning This call is dangerous! It does not have a size limit, and it
 * does not have any way of checking whether there is enough data in the
 * message (since the data could be coming from a stream.) When looping
 * through the map's contents, you must check for errors on each iteration
 * of the loop. Otherwise an attacker could craft a message declaring a map
 * of a billion elements which would throw your parsing code into an
 * infinite loop! You should strongly consider using mpack_expect_map_max_or_nil()
 * with a safe maximum size instead.
 *
 * @returns @c true if a map was read successfully; @c false if nil was read
 *     or an error occurred.
 * @throws mpack_error_type if the value is not a nil or map.
 */
bool mpack_expect_map_or_nil(mpack_reader_t* reader, uint32_t* count);

/**
 * Reads a nil node or the start of a map with a number of elements at most
 * max_count, returning whether a map was read and placing its number of
 * key/value pairs in count.
 *
 * If a map was read, a number of values follow equal to twice the element count
 * of the map, alternating between keys and values. @ref mpack_done_map() should
 * anlso be called once all elements have been read (only if a map was read.)
 *
 * @note Maps in JSON are unordered, so it is recommended not to expect
 * a specific ordering for your map values in case your data is converted
 * to/from JSON. Consider using mpack_expect_key_cstr() or mpack_expect_key_uint()
 * to switch on the key; see @ref docs/expect.md for examples.
 *
 * @returns @c true if a map was read successfully; @c false if nil was read
 *     or an error occurred.
 * @throws mpack_error_type if the value is not a nil or map.
 */
bool mpack_expect_map_max_or_nil(mpack_reader_t* reader, uint32_t max_count, uint32_t* count);

/**
 * Reads the start of an array, returning its element count.
 *
 * A number of values follow equal to the element count of the array.
 * @ref mpack_done_array() must be called once all elements have been read.
 *
 * @warning This call is dangerous! It does not have a size limit, and it
 * does not have any way of checking whether there is enough data in the
 * message (since the data could be coming from a stream.) When looping
 * through the array's contents, you must check for errors on each iteration
 * of the loop. Otherwise an attacker could craft a message declaring an array
 * of a billion elements which would throw your parsing code into an
 * infinite loop! You should strongly consider using mpack_expect_array_max()
 * with a safe maximum size instead.
 */
uint32_t mpack_expect_array(mpack_reader_t* reader);

/**
 * Reads the start of an array with a number of elements in the given range,
 * returning its element count.
 *
 * A number of values follow equal to the element count of the array.
 * @ref mpack_done_array() must be called once all elements have been read.
 *
 * min_count is returned if an error occurs.
 *
 * @throws mpack_error_type if the value is not an array or if its size does
 * not fall within the given range.
 */
uint32_t mpack_expect_array_range(mpack_reader_t* reader, uint32_t min_count, uint32_t max_count);

/**
 * Reads the start of an array with a number of elements at most @a max_count,
 * returning its element count.
 *
 * A number of values follow equal to the element count of the array.
 * @ref mpack_done_array() must be called once all elements have been read.
 *
 * Zero is returned if an error occurs.
 *
 * @throws mpack_error_type if the value is not an array or if its size is
 * greater than max_count.
 */
MPACK_INLINE uint32_t mpack_expect_array_max(mpack_reader_t* reader, uint32_t max_count) {
    return mpack_expect_array_range(reader, 0, max_count);
}

/**
 * Reads the start of an array of the exact size given.
 *
 * A number of values follow equal to the element count of the array.
 * @ref mpack_done_array() must be called once all elements have been read.
 *
 * @throws mpack_error_type if the value is not an array or if its size does
 * not match the given count.
 */
void mpack_expect_array_match(mpack_reader_t* reader, uint32_t count);

/**
 * Reads a nil node or the start of an array, returning whether an array was
 * read and placing its number of elements in count.
 *
 * If an array was read, a number of values follow equal to the element count
 * of the array. @ref mpack_done_array() should also be called once all elements
 * have been read (only if an array was read.)
 *
 * @warning This call is dangerous! It does not have a size limit, and it
 * does not have any way of checking whether there is enough data in the
 * message (since the data could be coming from a stream.) When looping
 * through the array's contents, you must check for errors on each iteration
 * of the loop. Otherwise an attacker could craft a message declaring an array
 * of a billion elements which would throw your parsing code into an
 * infinite loop! You should strongly consider using mpack_expect_array_max_or_nil()
 * with a safe maximum size instead.
 *
 * @returns @c true if an array was read successfully; @c false if nil was read
 *     or an error occurred.
 * @throws mpack_error_type if the value is not a nil or array.
 */
bool mpack_expect_array_or_nil(mpack_reader_t* reader, uint32_t* count);

/**
 * Reads a nil node or the start of an array with a number of elements at most
 * max_count, returning whether an array was read and placing its number of
 * key/value pairs in count.
 *
 * If an array was read, a number of values follow equal to the element count
 * of the array. @ref mpack_done_array() should also be called once all elements
 * have been read (only if an array was read.)
 *
 * @returns @c true if an array was read successfully; @c false if nil was read
 *     or an error occurred.
 * @throws mpack_error_type if the value is not a nil or array.
 */
bool mpack_expect_array_max_or_nil(mpack_reader_t* reader, uint32_t max_count, uint32_t* count);

#ifdef MPACK_MALLOC
/**
 * @hideinitializer
 *
 * Reads the start of an array and allocates storage for it, placing its
 * size in out_count. A number of objects follow equal to the element count
 * of the array. You must call @ref mpack_done_array() when done (even
 * if the element count is zero.)
 *
 * If an error occurs, NULL is returned and the reader is placed in an
 * error state.
 *
 * If the count is zero, NULL is returned. This does not indicate error.
 * You should not check the return value for NULL to check for errors; only
 * check the reader's error state.
 *
 * The allocated array must be freed with MPACK_FREE() (or simply free()
 * if MPack's allocator hasn't been customized.)
 *
 * @throws mpack_error_type if the value is not an array or if its size is
 * greater than max_count.
 */
#define mpack_expect_array_alloc(reader, Type, max_count, out_count) \
    ((Type*)mpack_expect_array_alloc_impl(reader, sizeof(Type), max_count, out_count, false))

/**
 * @hideinitializer
 *
 * Reads a nil node or the start of an array and allocates storage for it,
 * placing its size in out_count. A number of objects follow equal to the element
 * count of the array if a non-empty array was read.
 *
 * If an error occurs, NULL is returned and the reader is placed in an
 * error state.
 *
 * If a nil node was read, NULL is returned. If an empty array was read,
 * mpack_done_array() is called automatically and NULL is returned. These
 * do not indicate error. You should not check the return value for NULL
 * to check for errors; only check the reader's error state.
 *
 * The allocated array must be freed with MPACK_FREE() (or simply free()
 * if MPack's allocator hasn't been customized.)
 *
 * @warning You must call @ref mpack_done_array() if and only if a non-zero
 * element count is read. This function does not differentiate between nil
 * and an empty array.
 *
 * @throws mpack_error_type if the value is not an array or if its size is
 * greater than max_count.
 */
#define mpack_expect_array_or_nil_alloc(reader, Type, max_count, out_count) \
    ((Type*)mpack_expect_array_alloc_impl(reader, sizeof(Type), max_count, out_count, true))
#endif

/**
 * @}
 */

/** @cond */
#ifdef MPACK_MALLOC
void* mpack_expect_array_alloc_impl(mpack_reader_t* reader,
        size_t element_size, uint32_t max_count, uint32_t* out_count, bool allow_nil);
#endif
/** @endcond */


/**
 * @name String Functions
 * @{
 */

/**
 * Reads the start of a string, returning its size in bytes.
 *
 * The bytes follow and must be read separately with mpack_read_bytes()
 * or mpack_read_bytes_inplace(). mpack_done_str() must be called
 * once all bytes have been read.
 *
 * NUL bytes are allowed in the string, and no encoding checks are done.
 *
 * mpack_error_type is raised if the value is not a string.
 */
uint32_t mpack_expect_str(mpack_reader_t* reader);

/**
 * Reads a string of at most the given size, writing it into the
 * given buffer and returning its size in bytes.
 *
 * This does not add a null-terminator! Use mpack_expect_cstr() to
 * add a null-terminator.
 *
 * NUL bytes are allowed in the string, and no encoding checks are done.
 */
size_t mpack_expect_str_buf(mpack_reader_t* reader, char* buf, size_t bufsize);

/**
 * Reads a string into the given buffer, ensuring it is a valid UTF-8 string
 * and returning its size in bytes.
 *
 * This does not add a null-terminator! Use mpack_expect_utf8_cstr() to
 * add a null-terminator.
 *
 * This does not accept any UTF-8 variant such as Modified UTF-8, CESU-8 or
 * WTF-8. Only pure UTF-8 is allowed.
 *
 * NUL bytes are allowed in the string (as they are in UTF-8.)
 *
 * Raises mpack_error_too_big if there is not enough room for the string.
 * Raises mpack_error_type if the value is not a string or is not a valid UTF-8 string.
 */
size_t mpack_expect_utf8(mpack_reader_t* reader, char* buf, size_t bufsize);

/**
 * Reads the start of a string, raising an error if its length is not
 * at most the given number of bytes (not including any null-terminator.)
 *
 * The bytes follow and must be read separately with mpack_read_bytes()
 * or mpack_read_bytes_inplace(). @ref mpack_done_str() must be called
 * once all bytes have been read.
 *
 * @throws mpack_error_type If the value is not a string.
 * @throws mpack_error_too_big If the string's length in bytes is larger than the given maximum size.
 */
MPACK_INLINE uint32_t mpack_expect_str_max(mpack_reader_t* reader, uint32_t maxsize) {
    uint32_t length = mpack_expect_str(reader);
    if (length > maxsize) {
        mpack_reader_flag_error(reader, mpack_error_too_big);
        return 0;
    }
    return length;
}

/**
 * Reads the start of a string, raising an error if its length is not
 * exactly the given number of bytes (not including any null-terminator.)
 *
 * The bytes follow and must be read separately with mpack_read_bytes()
 * or mpack_read_bytes_inplace(). @ref mpack_done_str() must be called
 * once all bytes have been read.
 *
 * mpack_error_type is raised if the value is not a string or if its
 * length does not match.
 */
MPACK_INLINE void mpack_expect_str_length(mpack_reader_t* reader, uint32_t count) {
    if (mpack_expect_str(reader) != count)
        mpack_reader_flag_error(reader, mpack_error_type);
}

/**
 * Reads a string, ensuring it exactly matches the given string.
 *
 * Remember that maps are unordered in JSON. Don't use this for map keys
 * unless the map has only a single key!
 */
void mpack_expect_str_match(mpack_reader_t* reader, const char* str, size_t length);

/**
 * Reads a string into the given buffer, ensures it has no null bytes,
 * and adds a null-terminator at the end.
 *
 * Raises mpack_error_too_big if there is not enough room for the string and null-terminator.
 * Raises mpack_error_type if the value is not a string or contains a null byte.
 */
void mpack_expect_cstr(mpack_reader_t* reader, char* buf, size_t size);

/**
 * Reads a string into the given buffer, ensures it is a valid UTF-8 string
 * without NUL characters, and adds a null-terminator at the end.
 *
 * This does not accept any UTF-8 variant such as Modified UTF-8, CESU-8 or
 * WTF-8. Only pure UTF-8 is allowed, but without the NUL character, since
 * it cannot be represented in a null-terminated string.
 *
 * Raises mpack_error_too_big if there is not enough room for the string and null-terminator.
 * Raises mpack_error_type if the value is not a string or is not a valid UTF-8 string.
 */
void mpack_expect_utf8_cstr(mpack_reader_t* reader, char* buf, size_t size);

#ifdef MPACK_MALLOC
/**
 * Reads a string with the given total maximum size (including space for a
 * null-terminator), allocates storage for it, ensures it has no null-bytes,
 * and adds a null-terminator at the end. You assume ownership of the
 * returned pointer if reading succeeds.
 *
 * The allocated string must be freed with MPACK_FREE() (or simply free()
 * if MPack's allocator hasn't been customized.)
 *
 * @throws mpack_error_too_big If the string plus null-terminator is larger than the given maxsize.
 * @throws mpack_error_type If the value is not a string or contains a null byte.
 */
char* mpack_expect_cstr_alloc(mpack_reader_t* reader, size_t maxsize);

/**
 * Reads a string with the given total maximum size (including space for a
 * null-terminator), allocates storage for it, ensures it is valid UTF-8
 * with no null-bytes, and adds a null-terminator at the end. You assume
 * ownership of the returned pointer if reading succeeds.
 *
 * The length in bytes of the string, not including the null-terminator,
 * will be written to size.
 *
 * This does not accept any UTF-8 variant such as Modified UTF-8, CESU-8 or
 * WTF-8. Only pure UTF-8 is allowed, but without the NUL character, since
 * it cannot be represented in a null-terminated string.
 *
 * The allocated string must be freed with MPACK_FREE() (or simply free()
 * if MPack's allocator hasn't been customized.)
 * if you want a null-terminator.
 *
 * @throws mpack_error_too_big If the string plus null-terminator is larger
 *     than the given maxsize.
 * @throws mpack_error_type If the value is not a string or contains
 *     invalid UTF-8 or a null byte.
 */
char* mpack_expect_utf8_cstr_alloc(mpack_reader_t* reader, size_t maxsize);
#endif

/**
 * Reads a string, ensuring it exactly matches the given null-terminated
 * string.
 *
 * Remember that maps are unordered in JSON. Don't use this for map keys
 * unless the map has only a single key!
 */
MPACK_INLINE void mpack_expect_cstr_match(mpack_reader_t* reader, const char* cstr) {
    mpack_assert(cstr != NULL, "cstr pointer is NULL");
    mpack_expect_str_match(reader, cstr, mpack_strlen(cstr));
}

/**
 * @}
 */

/**
 * @name Binary Data
 * @{
 */

/**
 * Reads the start of a binary blob, returning its size in bytes.
 *
 * The bytes follow and must be read separately with mpack_read_bytes()
 * or mpack_read_bytes_inplace(). @ref mpack_done_bin() must be called
 * once all bytes have been read.
 *
 * mpack_error_type is raised if the value is not a binary blob.
 */
uint32_t mpack_expect_bin(mpack_reader_t* reader);

/**
 * Reads the start of a binary blob, raising an error if its length is not
 * at most the given number of bytes.
 *
 * The bytes follow and must be read separately with mpack_read_bytes()
 * or mpack_read_bytes_inplace(). @ref mpack_done_bin() must be called
 * once all bytes have been read.
 *
 * mpack_error_type is raised if the value is not a binary blob or if its
 * length does not match.
 */
MPACK_INLINE uint32_t mpack_expect_bin_max(mpack_reader_t* reader, uint32_t maxsize) {
    uint32_t length = mpack_expect_bin(reader);
    if (length > maxsize) {
        mpack_reader_flag_error(reader, mpack_error_type);
        return 0;
    }
    return length;
}

/**
 * Reads the start of a binary blob, raising an error if its length is not
 * exactly the given number of bytes.
 *
 * The bytes follow and must be read separately with mpack_read_bytes()
 * or mpack_read_bytes_inplace(). @ref mpack_done_bin() must be called
 * once all bytes have been read.
 *
 * @throws mpack_error_type if the value is not a binary blob or if its size
 * does not match.
 */
MPACK_INLINE void mpack_expect_bin_size(mpack_reader_t* reader, uint32_t count) {
    if (mpack_expect_bin(reader) != count)
        mpack_reader_flag_error(reader, mpack_error_type);
}

/**
 * Reads a binary blob into the given buffer, returning its size in bytes.
 *
 * For compatibility, this will accept if the underlying type is string or
 * binary (since in MessagePack 1.0, strings and binary data were combined
 * under the "raw" type which became string in 1.1.)
 */
size_t mpack_expect_bin_buf(mpack_reader_t* reader, char* buf, size_t size);

/**
 * Reads a binary blob with the exact given size into the given buffer.
 *
 * For compatibility, this will accept if the underlying type is string or
 * binary (since in MessagePack 1.0, strings and binary data were combined
 * under the "raw" type which became string in 1.1.)
 *
 * @throws mpack_error_type if the value is not a binary blob or if its size
 * does not match.
 */
void mpack_expect_bin_size_buf(mpack_reader_t* reader, char* buf, uint32_t size);

/**
 * Reads a binary blob with the given total maximum size, allocating storage for it.
 */
char* mpack_expect_bin_alloc(mpack_reader_t* reader, size_t maxsize, size_t* size);

/**
 * @}
 */

/**
 * @name Extension Functions
 * @{
 */

#if MPACK_EXTENSIONS
/**
 * Reads the start of an extension blob, returning its size in bytes and
 * placing the type into @p type.
 *
 * The bytes follow and must be read separately with mpack_read_bytes()
 * or mpack_read_bytes_inplace(). @ref mpack_done_ext() must be called
 * once all bytes have been read.
 *
 * @p type will be a user-defined type in the range [0,127] or a reserved type
 * in the range [-128,-2].
 *
 * mpack_error_type is raised if the value is not an extension blob. The @p
 * type value is zero if an error occurs.
 *
 * @note This cannot be used to match a timestamp. @ref mpack_error_type will
 * be flagged if the value is a timestamp. Use mpack_expect_timestamp() or
 * mpack_expect_timestamp_truncate() instead.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 *
 * @warning Be careful when using reserved types. They may no longer be ext
 * types in the future, and previously valid data containing reserved types may
 * become invalid in the future.
 */
uint32_t mpack_expect_ext(mpack_reader_t* reader, int8_t* type);

/**
 * Reads the start of an extension blob, raising an error if its length is not
 * at most the given number of bytes and placing the type into @p type.
 *
 * The bytes follow and must be read separately with mpack_read_bytes()
 * or mpack_read_bytes_inplace(). @ref mpack_done_ext() must be called
 * once all bytes have been read.
 *
 * mpack_error_type is raised if the value is not an extension blob or if its
 * length does not match. The @p type value is zero if an error is raised.
 *
 * @p type will be a user-defined type in the range [0,127] or a reserved type
 * in the range [-128,-2].
 *
 * @note This cannot be used to match a timestamp. @ref mpack_error_type will
 * be flagged if the value is a timestamp. Use mpack_expect_timestamp() or
 * mpack_expect_timestamp_truncate() instead.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 *
 * @warning Be careful when using reserved types. They may no longer be ext
 * types in the future, and previously valid data containing reserved types may
 * become invalid in the future.
 *
 * @see mpack_expect_ext()
 */
MPACK_INLINE uint32_t mpack_expect_ext_max(mpack_reader_t* reader, int8_t* type, uint32_t maxsize) {
    uint32_t length = mpack_expect_ext(reader, type);
    if (length > maxsize) {
        mpack_reader_flag_error(reader, mpack_error_type);
        return 0;
    }
    return length;
}

/**
 * Reads the start of an extension blob, raising an error if its length is not
 * exactly the given number of bytes and placing the type into @p type.
 *
 * The bytes follow and must be read separately with mpack_read_bytes()
 * or mpack_read_bytes_inplace(). @ref mpack_done_ext() must be called
 * once all bytes have been read.
 *
 * mpack_error_type is raised if the value is not an extension blob or if its
 * length does not match. The @p type value is zero if an error is raised.
 *
 * @p type will be a user-defined type in the range [0,127] or a reserved type
 * in the range [-128,-2].
 *
 * @note This cannot be used to match a timestamp. @ref mpack_error_type will
 * be flagged if the value is a timestamp. Use mpack_expect_timestamp() or
 * mpack_expect_timestamp_truncate() instead.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 *
 * @warning Be careful when using reserved types. They may no longer be ext
 * types in the future, and previously valid data containing reserved types may
 * become invalid in the future.
 *
 * @see mpack_expect_ext()
 */
MPACK_INLINE void mpack_expect_ext_size(mpack_reader_t* reader, int8_t* type, uint32_t count) {
    if (mpack_expect_ext(reader, type) != count) {
        *type = 0;
        mpack_reader_flag_error(reader, mpack_error_type);
    }
}

/**
 * Reads an extension blob into the given buffer, returning its size in bytes
 * and placing the type into @p type.
 *
 * mpack_error_type is raised if the value is not an extension blob or if its
 * length does not match. The @p type value is zero if an error is raised.
 *
 * @p type will be a user-defined type in the range [0,127] or a reserved type
 * in the range [-128,-2].
 *
 * @note This cannot be used to match a timestamp. @ref mpack_error_type will
 * be flagged if the value is a timestamp. Use mpack_expect_timestamp() or
 * mpack_expect_timestamp_truncate() instead.
 *
 * @warning Be careful when using reserved types. They may no longer be ext
 * types in the future, and previously valid data containing reserved types may
 * become invalid in the future.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 *
 * @see mpack_expect_ext()
 */
size_t mpack_expect_ext_buf(mpack_reader_t* reader, int8_t* type, char* buf, size_t size);
#endif

#if MPACK_EXTENSIONS && defined(MPACK_MALLOC)
/**
 * Reads an extension blob with the given total maximum size, allocating
 * storage for it, and placing the type into @p type.
 *
 * mpack_error_type is raised if the value is not an extension blob or if its
 * length does not match. The @p type value is zero if an error is raised.
 *
 * @p type will be a user-defined type in the range [0,127] or a reserved type
 * in the range [-128,-2].
 *
 * @note This cannot be used to match a timestamp. @ref mpack_error_type will
 * be flagged if the value is a timestamp. Use mpack_expect_timestamp() or
 * mpack_expect_timestamp_truncate() instead.
 *
 * @warning Be careful when using reserved types. They may no longer be ext
 * types in the future, and previously valid data containing reserved types may
 * become invalid in the future.
 *
 * @note This requires @ref MPACK_EXTENSIONS and @ref MPACK_MALLOC.
 *
 * @see mpack_expect_ext()
 */
char* mpack_expect_ext_alloc(mpack_reader_t* reader, int8_t* type, size_t maxsize, size_t* size);
#endif

/**
 * @}
 */

/**
 * @name Special Functions
 * @{
 */

/**
 * Reads a MessagePack object header (an MPack tag), expecting it to exactly
 * match the given tag.
 *
 * If the type is compound (i.e. is a map, array, string, binary or
 * extension type), additional reads are required to get the contained
 * data, and the corresponding done function must be called when done.
 *
 * @throws mpack_error_type if the tag does not match
 *
 * @see mpack_read_bytes()
 * @see mpack_done_array()
 * @see mpack_done_map()
 * @see mpack_done_str()
 * @see mpack_done_bin()
 * @see mpack_done_ext()
 */
void mpack_expect_tag(mpack_reader_t* reader, mpack_tag_t tag);

/**
 * Expects a string matching one of the strings in the given array,
 * returning its array index.
 *
 * If the value does not match any of the given strings,
 * @ref mpack_error_type is flagged. Use mpack_expect_enum_optional()
 * if you want to allow other values than the given strings.
 *
 * If any error occurs or the reader is in an error state, @a count
 * is returned.
 *
 * This can be used to quickly parse a string into an enum when the
 * enum values range from 0 to @a count-1. If the last value in the
 * enum is a special "count" value, it can be passed as the count,
 * and the return value can be cast directly to the enum type.
 *
 * @code{.c}
 * typedef enum           { APPLE ,  BANANA ,  ORANGE , COUNT} fruit_t;
 * const char* fruits[] = {"apple", "banana", "orange"};
 *
 * fruit_t fruit = (fruit_t)mpack_expect_enum(reader, fruits, COUNT);
 * @endcode
 *
 * See @ref docs/expect.md for more examples.
 *
 * The maximum string length is the size of the buffer (strings are read in-place.)
 *
 * @param reader The reader
 * @param strings An array of expected strings of length count
 * @param count The number of strings
 * @return The index of the matched string, or @a count in case of error
 */
size_t mpack_expect_enum(mpack_reader_t* reader, const char* strings[], size_t count);

/**
 * Expects a string matching one of the strings in the given array
 * returning its array index, or @a count if no strings match.
 *
 * If the value is not a string, or it does not match any of the
 * given strings, @a count is returned and no error is flagged.
 *
 * If any error occurs or the reader is in an error state, @a count
 * is returned.
 *
 * This can be used to quickly parse a string into an enum when the
 * enum values range from 0 to @a count-1. If the last value in the
 * enum is a special "count" value, it can be passed as the count,
 * and the return value can be cast directly to the enum type.
 *
 * @code{.c}
 * typedef enum           { APPLE ,  BANANA ,  ORANGE , COUNT} fruit_t;
 * const char* fruits[] = {"apple", "banana", "orange"};
 *
 * fruit_t fruit = (fruit_t)mpack_expect_enum_optional(reader, fruits, COUNT);
 * @endcode
 *
 * See @ref docs/expect.md for more examples.
 *
 * The maximum string length is the size of the buffer (strings are read in-place.)
 *
 * @param reader The reader
 * @param strings An array of expected strings of length count
 * @param count The number of strings
 *
 * @return The index of the matched string, or @a count if it does not
 * match or an error occurs
 */
size_t mpack_expect_enum_optional(mpack_reader_t* reader, const char* strings[], size_t count);

/**
 * Expects an unsigned integer map key between 0 and count-1, marking it
 * as found in the given bool array and returning it.
 *
 * This is a helper for switching among int keys in a map. It is
 * typically used with an enum to define the key values. It should
 * be called in the expression of a switch() statement. See @ref
 * docs/expect.md for an example.
 *
 * The found array must be cleared before expecting the first key. If the
 * flag for a given key is already set when found (i.e. the map contains a
 * duplicate key), mpack_error_invalid is flagged.
 *
 * If the key is not a non-negative integer, or if the key is @a count or
 * larger, @a count is returned and no error is flagged. If you want an error
 * on unrecognized keys, flag an error in the default case in your switch;
 * otherwise you must call mpack_discard() to discard its content.
 *
 * @param reader The reader
 * @param found An array of bool flags of length count
 * @param count The number of values in the found array, and one more than the
 *              maximum allowed key
 *
 * @see @ref docs/expect.md
 */
size_t mpack_expect_key_uint(mpack_reader_t* reader, bool found[], size_t count);

/**
 * Expects a string map key matching one of the strings in the given key list,
 * marking it as found in the given bool array and returning its index.
 *
 * This is a helper for switching among string keys in a map. It is
 * typically used with an enum with names matching the strings in the
 * array to define the key indices. It should be called in the expression
 * of a switch() statement. See @ref docs/expect.md for an example.
 *
 * The found array must be cleared before expecting the first key. If the
 * flag for a given key is already set when found (i.e. the map contains a
 * duplicate key), mpack_error_invalid is flagged.
 *
 * If the key is unrecognized, count is returned and no error is flagged. If
 * you want an error on unrecognized keys, flag an error in the default case
 * in your switch; otherwise you must call mpack_discard() to discard its content.
 *
 * The maximum key length is the size of the buffer (keys are read in-place.)
 *
 * @param reader The reader
 * @param keys An array of expected string keys of length count
 * @param found An array of bool flags of length count
 * @param count The number of values in the keys and found arrays
 *
 * @see @ref docs/expect.md
 */
size_t mpack_expect_key_cstr(mpack_reader_t* reader, const char* keys[],
        bool found[], size_t count);

/**
 * @}
 */

/**
 * @}
 */

#endif

MPACK_EXTERN_C_END
MPACK_SILENCE_WARNINGS_END

#endif



/* mpack/mpack-node.h.h */

/**
 * @file
 *
 * Declares the MPack dynamic Node API.
 */

#ifndef MPACK_NODE_H
#define MPACK_NODE_H 1

/* #include "mpack-reader.h" */

MPACK_SILENCE_WARNINGS_BEGIN
MPACK_EXTERN_C_BEGIN

#if MPACK_NODE

/**
 * @defgroup node Node API
 *
 * The MPack Node API allows you to parse a chunk of MessagePack into a
 * dynamically typed data structure, providing random access to the parsed
 * data.
 *
 * See @ref docs/node.md for examples.
 *
 * @{
 */

/**
 * A handle to node data in a parsed MPack tree.
 *
 * Nodes represent either primitive values or compound types. If a
 * node is a compound type, it contains a pointer to its child nodes,
 * or a pointer to its underlying data.
 *
 * Nodes are immutable.
 *
 * @note @ref mpack_node_t is an opaque reference to the node data, not the
 * node data itself. (It contains pointers to both the node data and the tree.)
 * It is passed by value in the Node API.
 */
typedef struct mpack_node_t mpack_node_t;

/**
 * The storage for nodes in an MPack tree.
 *
 * You only need to use this if you intend to provide your own storage
 * for nodes instead of letting the tree allocate it.
 *
 * @ref mpack_node_data_t is 16 bytes on most common architectures (32-bit
 * and 64-bit.)
 */
typedef struct mpack_node_data_t mpack_node_data_t;

/**
 * An MPack tree parser to parse a blob or stream of MessagePack.
 *
 * When a message is parsed, the tree contains a single root node which
 * contains all parsed data. The tree and its nodes are immutable.
 */
typedef struct mpack_tree_t mpack_tree_t;

/**
 * An error handler function to be called when an error is flagged on
 * the tree.
 *
 * The error handler will only be called once on the first error flagged;
 * any subsequent node reads and errors are ignored, and the tree is
 * permanently in that error state.
 *
 * MPack is safe against non-local jumps out of error handler callbacks.
 * This means you are allowed to longjmp or throw an exception (in C++,
 * Objective-C, or with SEH) out of this callback.
 *
 * Bear in mind when using longjmp that local non-volatile variables that
 * have changed are undefined when setjmp() returns, so you can't put the
 * tree on the stack in the same activation frame as the setjmp without
 * declaring it volatile.
 *
 * You must still eventually destroy the tree. It is not destroyed
 * automatically when an error is flagged. It is safe to destroy the
 * tree within this error callback, but you will either need to perform
 * a non-local jump, or store something in your context to identify
 * that the tree is destroyed since any future accesses to it cause
 * undefined behavior.
 */
typedef void (*mpack_tree_error_t)(mpack_tree_t* tree, mpack_error_t error);

/**
 * The MPack tree's read function. It should fill the buffer with as many bytes
 * as are immediately available up to the given @c count, returning the number
 * of bytes written to the buffer.
 *
 * In case of error, it should flag an appropriate error on the reader
 * (usually @ref mpack_error_io.)
 *
 * The blocking or non-blocking behaviour of the read should match whether you
 * are using mpack_tree_parse() or mpack_tree_try_parse().
 *
 * If you are using mpack_tree_parse(), the read should block until at least
 * one byte is read. If you return 0, mpack_tree_parse() will raise @ref
 * mpack_error_io.
 *
 * If you are using mpack_tree_try_parse(), the read function can always
 * return 0, and must never block waiting for data (otherwise
 * mpack_tree_try_parse() would be equivalent to mpack_tree_parse().)
 * When you return 0, mpack_tree_try_parse() will return false without flagging
 * an error.
 */
typedef size_t (*mpack_tree_read_t)(mpack_tree_t* tree, char* buffer, size_t count);

/**
 * A teardown function to be called when the tree is destroyed.
 */
typedef void (*mpack_tree_teardown_t)(mpack_tree_t* tree);



/* Hide internals from documentation */
/** @cond */

struct mpack_node_t {
    mpack_node_data_t* data;
    mpack_tree_t* tree;
};

struct mpack_node_data_t {
    mpack_type_t type;

    /*
     * The element count if the type is an array;
     * the number of key/value pairs if the type is map;
     * or the number of bytes if the type is str, bin or ext.
     */
    uint32_t len;

    union {
        bool     b; /* The value if the type is bool. */

        #if MPACK_FLOAT
        float    f; /* The value if the type is float. */
        #else
        uint32_t f; /*< The raw value if the type is float. */
        #endif

        #if MPACK_DOUBLE
        double   d; /* The value if the type is double. */
        #else
        uint64_t d; /*< The raw value if the type is double. */
        #endif

        int64_t  i; /* The value if the type is signed int. */
        uint64_t u; /* The value if the type is unsigned int. */
        size_t offset; /* The byte offset for str, bin and ext */

        mpack_node_data_t* children; /* The children for map or array */
    } value;
};

typedef struct mpack_tree_page_t {
    struct mpack_tree_page_t* next;
    mpack_node_data_t nodes[1]; // variable size
} mpack_tree_page_t;

typedef enum mpack_tree_parse_state_t {
    mpack_tree_parse_state_not_started,
    mpack_tree_parse_state_in_progress,
    mpack_tree_parse_state_parsed,
} mpack_tree_parse_state_t;

typedef struct mpack_level_t {
    mpack_node_data_t* child;
    size_t left; // children left in level
} mpack_level_t;

typedef struct mpack_tree_parser_t {
    mpack_tree_parse_state_t state;

    // We keep track of the number of "possible nodes" left in the data rather
    // than the number of bytes.
    //
    // When a map or array is parsed, we ensure at least one byte for each child
    // exists and subtract them right away. This ensures that if ever a map or
    // array declares more elements than could possibly be contained in the data,
    // we will error out immediately rather than allocating storage for them.
    //
    // For example malicious data that repeats 0xDE 0xFF 0xFF (start of a map
    // with 65536 key-value pairs) would otherwise cause us to run out of
    // memory. With this, the parser can allocate at most as many nodes as
    // there are bytes in the data (plus the paging overhead, 12%.) An error
    // will be flagged immediately if and when there isn't enough data left to
    // fully read all children of all open compound types on the parsing stack.
    //
    // Once an entire message has been parsed (and there are no nodes left to
    // parse whose bytes have been subtracted), this matches the number of left
    // over bytes in the data.
    size_t possible_nodes_left;

    mpack_node_data_t* nodes; // next node in current page/pool
    size_t nodes_left; // nodes left in current page/pool

    size_t current_node_reserved;
    size_t level;

    #ifdef MPACK_MALLOC
    // It's much faster to allocate the initial parsing stack inline within the
    // parser. We replace it with a heap allocation if we need to grow it.
    mpack_level_t* stack;
    size_t stack_capacity;
    bool stack_owned;
    mpack_level_t stack_local[MPACK_NODE_INITIAL_DEPTH];
    #else
    // Without malloc(), we have to reserve a parsing stack the maximum allowed
    // parsing depth.
    mpack_level_t stack[MPACK_NODE_MAX_DEPTH_WITHOUT_MALLOC];
    #endif
} mpack_tree_parser_t;

struct mpack_tree_t {
    mpack_tree_error_t error_fn;    /* Function to call on error */
    mpack_tree_read_t read_fn;      /* Function to call to read more data */
    mpack_tree_teardown_t teardown; /* Function to teardown the context on destroy */
    void* context;                  /* Context for tree callbacks */

    mpack_node_data_t nil_node;     /* a nil node to be returned in case of error */
    mpack_node_data_t missing_node; /* a missing node to be returned in optional lookups */
    mpack_error_t error;

    #ifdef MPACK_MALLOC
    char* buffer;
    size_t buffer_capacity;
    #endif

    const char* data;
    size_t data_length; // length of data (and content of buffer, if used)

    size_t size; // size in bytes of tree (usually matches data_length, but not if tree has trailing data)
    size_t node_count; // total number of nodes in tree (across all pages)

    size_t max_size;  // maximum message size
    size_t max_nodes; // maximum nodes in a message

    mpack_tree_parser_t parser;
    mpack_node_data_t* root;

    mpack_node_data_t* pool; // pool, or NULL if no pool provided
    size_t pool_count;

    #ifdef MPACK_MALLOC
    mpack_tree_page_t* next;
    #endif
};

// internal functions

MPACK_INLINE mpack_node_t mpack_node(mpack_tree_t* tree, mpack_node_data_t* data) {
    mpack_node_t node;
    node.data = data;
    node.tree = tree;
    return node;
}

MPACK_INLINE mpack_node_data_t* mpack_node_child(mpack_node_t node, size_t child) {
    return node.data->value.children + child;
}

MPACK_INLINE mpack_node_t mpack_tree_nil_node(mpack_tree_t* tree) {
    return mpack_node(tree, &tree->nil_node);
}

MPACK_INLINE mpack_node_t mpack_tree_missing_node(mpack_tree_t* tree) {
    return mpack_node(tree, &tree->missing_node);
}

/** @endcond */



/**
 * @name Tree Initialization
 * @{
 */

#ifdef MPACK_MALLOC
/**
 * Initializes a tree parser with the given data.
 *
 * Configure the tree if desired, then call mpack_tree_parse() to parse it. The
 * tree will allocate pages of nodes as needed and will free them when
 * destroyed.
 *
 * The tree must be destroyed with mpack_tree_destroy().
 *
 * Any string or blob data types reference the original data, so the given data
 * pointer must remain valid until after the tree is destroyed.
 */
void mpack_tree_init_data(mpack_tree_t* tree, const char* data, size_t length);

/**
 * Deprecated.
 *
 * \deprecated Renamed to mpack_tree_init_data().
 */
MPACK_INLINE void mpack_tree_init(mpack_tree_t* tree, const char* data, size_t length) {
    mpack_tree_init_data(tree, data, length);
}

/**
 * Initializes a tree parser from an unbounded stream, or a stream of
 * unknown length.
 *
 * The parser can be used to read a single message from a stream of unknown
 * length, or multiple messages from an unbounded stream, allowing it to
 * be used for RPC communication. Call @ref mpack_tree_parse() to parse
 * a message from a blocking stream, or @ref mpack_tree_try_parse() for a
 * non-blocking stream.
 *
 * The stream will use a growable internal buffer to store the most recent
 * message, as well as allocated pages of nodes for the parse tree.
 *
 * Maximum allowances for message size and node count must be specified in this
 * function (since the stream is unbounded.) They can be changed later with
 * @ref mpack_tree_set_limits().
 *
 * @param tree The tree parser
 * @param read_fn The read function
 * @param context The context for the read function
 * @param max_message_size The maximum size of a message in bytes
 * @param max_message_nodes The maximum number of nodes per message. See
 *        @ref mpack_node_data_t for the size of nodes.
 *
 * @see mpack_tree_read_t
 * @see mpack_reader_context()
 */
void mpack_tree_init_stream(mpack_tree_t* tree, mpack_tree_read_t read_fn, void* context,
        size_t max_message_size, size_t max_message_nodes);
#endif

/**
 * Initializes a tree parser with the given data, using the given node data
 * pool to store the results.
 *
 * Configure the tree if desired, then call mpack_tree_parse() to parse it.
 *
 * If the data does not fit in the pool, @ref mpack_error_too_big will be flagged
 * on the tree.
 *
 * The tree must be destroyed with mpack_tree_destroy(), even if parsing fails.
 */
void mpack_tree_init_pool(mpack_tree_t* tree, const char* data, size_t length,
        mpack_node_data_t* node_pool, size_t node_pool_count);

/**
 * Initializes an MPack tree directly into an error state. Use this if you
 * are writing a wrapper to another <tt>mpack_tree_init*()</tt> function which
 * can fail its setup.
 */
void mpack_tree_init_error(mpack_tree_t* tree, mpack_error_t error);

#if MPACK_STDIO
/**
 * Initializes a tree to parse the given file. The tree must be destroyed with
 * mpack_tree_destroy(), even if parsing fails.
 *
 * The file is opened, loaded fully into memory, and closed before this call
 * returns.
 *
 * @param tree The tree to initialize
 * @param filename The filename passed to fopen() to read the file
 * @param max_bytes The maximum size of file to load, or 0 for unlimited size.
 */
void mpack_tree_init_filename(mpack_tree_t* tree, const char* filename, size_t max_bytes);

/**
 * Deprecated.
 *
 * \deprecated Renamed to mpack_tree_init_filename().
 */
MPACK_INLINE void mpack_tree_init_file(mpack_tree_t* tree, const char* filename, size_t max_bytes) {
    mpack_tree_init_filename(tree, filename, max_bytes);
}

/**
 * Initializes a tree to parse the given libc FILE. This can be used to
 * read from stdin, or from a file opened separately.
 *
 * The tree must be destroyed with mpack_tree_destroy(), even if parsing fails.
 *
 * The FILE is fully loaded fully into memory (and closed if requested) before
 * this call returns.
 *
 * @param tree The tree to initialize.
 * @param stdfile The FILE.
 * @param max_bytes The maximum size of file to load, or 0 for unlimited size.
 * @param close_when_done If true, fclose() will be called on the FILE when it
 *         is no longer needed. If false, the file will not be closed when
 *         reading is done.
 *
 * @warning The tree will read all data in the FILE before parsing it. If this
 *          is used on stdin, the parser will block until it is closed, even if
 *          a complete message has been written to it!
 */
void mpack_tree_init_stdfile(mpack_tree_t* tree, FILE* stdfile, size_t max_bytes, bool close_when_done);
#endif

/**
 * @}
 */

/**
 * @name Tree Functions
 * @{
 */

/**
 * Sets the maximum byte size and maximum number of nodes allowed per message.
 *
 * The default is SIZE_MAX (no limit) unless @ref mpack_tree_init_stream() is
 * called (where maximums are required.)
 *
 * If a pool of nodes is used, the node limit is the lesser of this limit and
 * the pool size.
 *
 * @param tree The tree parser
 * @param max_message_size The maximum size of a message in bytes
 * @param max_message_nodes The maximum number of nodes per message. See
 *        @ref mpack_node_data_t for the size of nodes.
 */
void mpack_tree_set_limits(mpack_tree_t* tree, size_t max_message_size,
        size_t max_message_nodes);

/**
 * Parses a MessagePack message into a tree of immutable nodes.
 *
 * If successful, the root node will be available under @ref mpack_tree_root().
 * If not, an appropriate error will be flagged.
 *
 * This can be called repeatedly to parse a series of messages from a data
 * source. When this is called, all previous nodes from this tree and their
 * contents (including the root node) are invalidated.
 *
 * If this is called with a stream (see @ref mpack_tree_init_stream()), the
 * stream must block until data is available. (Otherwise, if this is called on
 * a non-blocking stream, parsing will fail with @ref mpack_error_io when the
 * fill function returns 0.)
 *
 * There is no way to recover a tree in an error state. It must be destroyed.
 */
void mpack_tree_parse(mpack_tree_t* tree);

/**
 * Attempts to parse a MessagePack message from a non-blocking stream into a
 * tree of immutable nodes.
 *
 * A non-blocking read function must have been passed to the tree in
 * mpack_tree_init_stream().
 *
 * If this returns true, a message is available under
 * @ref mpack_tree_root(). The tree nodes and data will be valid until
 * the next time a parse is started.
 *
 * If this returns false, no message is available, because either not enough
 * data is available yet or an error has occurred. You must check the tree for
 * errors whenever this returns false. If there is no error, you should try
 * again later when more data is available. (You will want to select()/poll()
 * on the underlying socket or use some other asynchronous mechanism to
 * determine when it has data.)
 *
 * There is no way to recover a tree in an error state. It must be destroyed.
 *
 * @see mpack_tree_init_stream()
 */
bool mpack_tree_try_parse(mpack_tree_t* tree);

/**
 * Returns the root node of the tree, if the tree is not in an error state.
 * Returns a nil node otherwise.
 *
 * @warning You must call mpack_tree_parse() before calling this. If
 * @ref mpack_tree_parse() was never called, the tree will assert.
 */
mpack_node_t mpack_tree_root(mpack_tree_t* tree);

/**
 * Returns the error state of the tree.
 */
MPACK_INLINE mpack_error_t mpack_tree_error(mpack_tree_t* tree) {
    return tree->error;
}

/**
 * Returns the size in bytes of the current parsed message.
 *
 * If there is something in the buffer after the MessagePack object, this can
 * be used to find it.
 *
 * This is zero if an error occurred during tree parsing (since the
 * portion of the data that the first complete object occupies cannot
 * be determined if the data is invalid or corrupted.)
 */
MPACK_INLINE size_t mpack_tree_size(mpack_tree_t* tree) {
    return tree->size;
}

/**
 * Destroys the tree.
 */
mpack_error_t mpack_tree_destroy(mpack_tree_t* tree);

/**
 * Sets the custom pointer to pass to the tree callbacks, such as teardown.
 *
 * @param tree The MPack tree.
 * @param context User data to pass to the tree callbacks.
 *
 * @see mpack_reader_context()
 */
MPACK_INLINE void mpack_tree_set_context(mpack_tree_t* tree, void* context) {
    tree->context = context;
}

/**
 * Returns the custom context for tree callbacks.
 *
 * @see mpack_tree_set_context
 * @see mpack_tree_init_stream
 */
MPACK_INLINE void* mpack_tree_context(mpack_tree_t* tree) {
    return tree->context;
}

/**
 * Sets the error function to call when an error is flagged on the tree.
 *
 * This should normally be used with mpack_tree_set_context() to register
 * a custom pointer to pass to the error function.
 *
 * See the definition of mpack_tree_error_t for more information about
 * what you can do from an error callback.
 *
 * @see mpack_tree_error_t
 * @param tree The MPack tree.
 * @param error_fn The function to call when an error is flagged on the tree.
 */
MPACK_INLINE void mpack_tree_set_error_handler(mpack_tree_t* tree, mpack_tree_error_t error_fn) {
    tree->error_fn = error_fn;
}

/**
 * Sets the teardown function to call when the tree is destroyed.
 *
 * This should normally be used with mpack_tree_set_context() to register
 * a custom pointer to pass to the teardown function.
 *
 * @param tree The MPack tree.
 * @param teardown The function to call when the tree is destroyed.
 */
MPACK_INLINE void mpack_tree_set_teardown(mpack_tree_t* tree, mpack_tree_teardown_t teardown) {
    tree->teardown = teardown;
}

/**
 * Places the tree in the given error state, calling the error callback if one
 * is set.
 *
 * This allows you to externally flag errors, for example if you are validating
 * data as you read it.
 *
 * If the tree is already in an error state, this call is ignored and no
 * error callback is called.
 */
void mpack_tree_flag_error(mpack_tree_t* tree, mpack_error_t error);

/**
 * @}
 */

/**
 * @name Node Core Functions
 * @{
 */

/**
 * Places the node's tree in the given error state, calling the error callback
 * if one is set.
 *
 * This allows you to externally flag errors, for example if you are validating
 * data as you read it.
 *
 * If the tree is already in an error state, this call is ignored and no
 * error callback is called.
 */
void mpack_node_flag_error(mpack_node_t node, mpack_error_t error);

/**
 * Returns the error state of the node's tree.
 */
MPACK_INLINE mpack_error_t mpack_node_error(mpack_node_t node) {
    return mpack_tree_error(node.tree);
}

/**
 * Returns a tag describing the given node, or a nil tag if the
 * tree is in an error state.
 */
mpack_tag_t mpack_node_tag(mpack_node_t node);

/** @cond */

#if MPACK_DEBUG && MPACK_STDIO
/*
 * Converts a node to a pseudo-JSON string for debugging purposes, placing the
 * result in the given buffer with a null-terminator.
 *
 * If the buffer does not have enough space, the result will be truncated (but
 * it is guaranteed to be null-terminated.)
 *
 * This is only available in debug mode, and only if stdio is available (since
 * it uses snprintf().) It's strictly for debugging purposes.
 */
void mpack_node_print_to_buffer(mpack_node_t node, char* buffer, size_t buffer_size);

/*
 * Converts a node to pseudo-JSON for debugging purposes, calling the given
 * callback as many times as is necessary to output the character data.
 *
 * No null-terminator or trailing newline will be written.
 *
 * This is only available in debug mode, and only if stdio is available (since
 * it uses snprintf().) It's strictly for debugging purposes.
 */
void mpack_node_print_to_callback(mpack_node_t node, mpack_print_callback_t callback, void* context);

/*
 * Converts a node to pseudo-JSON for debugging purposes
 * and pretty-prints it to the given file.
 *
 * This is only available in debug mode, and only if stdio is available (since
 * it uses snprintf().) It's strictly for debugging purposes.
 */
void mpack_node_print_to_file(mpack_node_t node, FILE* file);

/*
 * Converts a node to pseudo-JSON for debugging purposes
 * and pretty-prints it to stdout.
 *
 * This is only available in debug mode, and only if stdio is available (since
 * it uses snprintf().) It's strictly for debugging purposes.
 */
MPACK_INLINE void mpack_node_print_to_stdout(mpack_node_t node) {
    mpack_node_print_to_file(node, stdout);
}

/*
 * Deprecated.
 *
 * \deprecated Renamed to mpack_node_print_to_stdout().
 */
MPACK_INLINE void mpack_node_print(mpack_node_t node) {
    mpack_node_print_to_stdout(node);
}
#endif

/** @endcond */

/**
 * @}
 */

/**
 * @name Node Primitive Value Functions
 * @{
 */

/**
 * Returns the type of the node.
 */
mpack_type_t mpack_node_type(mpack_node_t node);

/**
 * Returns true if the given node is a nil node; false otherwise.
 *
 * To ensure that a node is nil and flag an error otherwise, use
 * mpack_node_nil().
 */
bool mpack_node_is_nil(mpack_node_t node);

/**
 * Returns true if the given node handle indicates a missing node; false otherwise.
 *
 * To ensure that a node is missing and flag an error otherwise, use
 * mpack_node_missing().
 */
bool mpack_node_is_missing(mpack_node_t node);

/**
 * Checks that the given node is of nil type, raising @ref mpack_error_type
 * otherwise.
 *
 * Use mpack_node_is_nil() to return whether the node is nil.
 */
void mpack_node_nil(mpack_node_t node);

/**
 * Checks that the given node indicates a missing node, raising @ref
 * mpack_error_type otherwise.
 *
 * Use mpack_node_is_missing() to return whether the node is missing.
 */
void mpack_node_missing(mpack_node_t node);

/**
 * Returns the bool value of the node. If this node is not of the correct
 * type, false is returned and mpack_error_type is raised.
 */
bool mpack_node_bool(mpack_node_t node);

/**
 * Checks if the given node is of bool type with value true, raising
 * mpack_error_type otherwise.
 */
void mpack_node_true(mpack_node_t node);

/**
 * Checks if the given node is of bool type with value false, raising
 * mpack_error_type otherwise.
 */
void mpack_node_false(mpack_node_t node);

/**
 * Returns the 8-bit unsigned value of the node. If this node is not
 * of a compatible type, @ref mpack_error_type is raised and zero is returned.
 */
uint8_t mpack_node_u8(mpack_node_t node);

/**
 * Returns the 8-bit signed value of the node. If this node is not
 * of a compatible type, @ref mpack_error_type is raised and zero is returned.
 */
int8_t mpack_node_i8(mpack_node_t node);

/**
 * Returns the 16-bit unsigned value of the node. If this node is not
 * of a compatible type, @ref mpack_error_type is raised and zero is returned.
 */
uint16_t mpack_node_u16(mpack_node_t node);

/**
 * Returns the 16-bit signed value of the node. If this node is not
 * of a compatible type, @ref mpack_error_type is raised and zero is returned.
 */
int16_t mpack_node_i16(mpack_node_t node);

/**
 * Returns the 32-bit unsigned value of the node. If this node is not
 * of a compatible type, @ref mpack_error_type is raised and zero is returned.
 */
uint32_t mpack_node_u32(mpack_node_t node);

/**
 * Returns the 32-bit signed value of the node. If this node is not
 * of a compatible type, @ref mpack_error_type is raised and zero is returned.
 */
int32_t mpack_node_i32(mpack_node_t node);

/**
 * Returns the 64-bit unsigned value of the node. If this node is not
 * of a compatible type, @ref mpack_error_type is raised, and zero is returned.
 */
uint64_t mpack_node_u64(mpack_node_t node);

/**
 * Returns the 64-bit signed value of the node. If this node is not
 * of a compatible type, @ref mpack_error_type is raised and zero is returned.
 */
int64_t mpack_node_i64(mpack_node_t node);

/**
 * Returns the unsigned int value of the node.
 *
 * Returns zero if an error occurs.
 *
 * @throws mpack_error_type If the node is not an integer type or does not fit in the range of an unsigned int
 */
unsigned int mpack_node_uint(mpack_node_t node);

/**
 * Returns the int value of the node.
 *
 * Returns zero if an error occurs.
 *
 * @throws mpack_error_type If the node is not an integer type or does not fit in the range of an int
 */
int mpack_node_int(mpack_node_t node);

#if MPACK_FLOAT
/**
 * Returns the float value of the node. The underlying value can be an
 * integer, float or double; the value is converted to a float.
 *
 * @note Reading a double or a large integer with this function can incur a
 * loss of precision.
 *
 * @throws mpack_error_type if the underlying value is not a float, double or integer.
 */
float mpack_node_float(mpack_node_t node);
#endif

#if MPACK_DOUBLE
/**
 * Returns the double value of the node. The underlying value can be an
 * integer, float or double; the value is converted to a double.
 *
 * @note Reading a very large integer with this function can incur a
 * loss of precision.
 *
 * @throws mpack_error_type if the underlying value is not a float, double or integer.
 */
double mpack_node_double(mpack_node_t node);
#endif

#if MPACK_FLOAT
/**
 * Returns the float value of the node. The underlying value must be a float,
 * not a double or an integer. This ensures no loss of precision can occur.
 *
 * @throws mpack_error_type if the underlying value is not a float.
 */
float mpack_node_float_strict(mpack_node_t node);
#endif

#if MPACK_DOUBLE
/**
 * Returns the double value of the node. The underlying value must be a float
 * or double, not an integer. This ensures no loss of precision can occur.
 *
 * @throws mpack_error_type if the underlying value is not a float or double.
 */
double mpack_node_double_strict(mpack_node_t node);
#endif

#if !MPACK_FLOAT
/**
 * Returns the float value of the node as a raw uint32_t. The underlying value
 * must be a float, not a double or an integer.
 *
 * @throws mpack_error_type if the underlying value is not a float.
 */
uint32_t mpack_node_raw_float(mpack_node_t node);
#endif

#if !MPACK_DOUBLE
/**
 * Returns the double value of the node as a raw uint64_t. The underlying value
 * must be a double, not a float or an integer.
 *
 * @throws mpack_error_type if the underlying value is not a float or double.
 */
uint64_t mpack_node_raw_double(mpack_node_t node);
#endif


#if MPACK_EXTENSIONS
/**
 * Returns a timestamp.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 *
 * @throws mpack_error_type if the underlying value is not a timestamp.
 */
mpack_timestamp_t mpack_node_timestamp(mpack_node_t node);

/**
 * Returns a timestamp's (signed) seconds since 1970-01-01T00:00:00Z.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 *
 * @throws mpack_error_type if the underlying value is not a timestamp.
 */
int64_t mpack_node_timestamp_seconds(mpack_node_t node);

/**
 * Returns a timestamp's additional nanoseconds.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 *
 * @return A nanosecond count between 0 and 999,999,999 inclusive.
 * @throws mpack_error_type if the underlying value is not a timestamp.
 */
uint32_t mpack_node_timestamp_nanoseconds(mpack_node_t node);
#endif

/**
 * @}
 */

/**
 * @name Node String and Data Functions
 * @{
 */

/**
 * Checks that the given node contains a valid UTF-8 string.
 *
 * If the string is invalid, this flags an error, which would cause subsequent calls
 * to mpack_node_str() to return NULL and mpack_node_strlen() to return zero. So you
 * can check the node for error immediately after calling this, or you can call those
 * functions to use the data anyway and check for errors later.
 *
 * @throws mpack_error_type If this node is not a string or does not contain valid UTF-8.
 *
 * @param node The string node to test
 *
 * @see mpack_node_str()
 * @see mpack_node_strlen()
 */
void mpack_node_check_utf8(mpack_node_t node);

/**
 * Checks that the given node contains a valid UTF-8 string with no NUL bytes.
 *
 * This does not check that the string has a null-terminator! It only checks whether
 * the string could safely be represented as a C-string by appending a null-terminator.
 * (If the string does already contain a null-terminator, this will flag an error.)
 *
 * This is performed automatically by other UTF-8 cstr helper functions. Only
 * call this if you will do something else with the data directly, but you still
 * want to ensure it will be valid as a UTF-8 C-string.
 *
 * @throws mpack_error_type If this node is not a string, does not contain valid UTF-8,
 *     or contains a NUL byte.
 *
 * @param node The string node to test
 *
 * @see mpack_node_str()
 * @see mpack_node_strlen()
 * @see mpack_node_copy_utf8_cstr()
 * @see mpack_node_utf8_cstr_alloc()
 */
void mpack_node_check_utf8_cstr(mpack_node_t node);

#if MPACK_EXTENSIONS
/**
 * Returns the extension type of the given ext node.
 *
 * This returns zero if the tree is in an error state.
 *
 * @note This requires @ref MPACK_EXTENSIONS.
 */
int8_t mpack_node_exttype(mpack_node_t node);
#endif

/**
 * Returns the number of bytes in the given bin node.
 *
 * This returns zero if the tree is in an error state.
 *
 * If this node is not a bin, @ref mpack_error_type is raised and zero is returned.
 */
size_t mpack_node_bin_size(mpack_node_t node);

/**
 * Returns the length of the given str, bin or ext node.
 *
 * This returns zero if the tree is in an error state.
 *
 * If this node is not a str, bin or ext, @ref mpack_error_type is raised and zero
 * is returned.
 */
uint32_t mpack_node_data_len(mpack_node_t node);

/**
 * Returns the length in bytes of the given string node. This does not
 * include any null-terminator.
 *
 * This returns zero if the tree is in an error state.
 *
 * If this node is not a str, @ref mpack_error_type is raised and zero is returned.
 */
size_t mpack_node_strlen(mpack_node_t node);

/**
 * Returns a pointer to the data contained by this node, ensuring the node is a
 * string.
 *
 * @warning Strings are not null-terminated! Use one of the cstr functions
 * to get a null-terminated string.
 *
 * The pointer is valid as long as the data backing the tree is valid.
 *
 * If this node is not a string, @ref mpack_error_type is raised and @c NULL is returned.
 *
 * @see mpack_node_copy_cstr()
 * @see mpack_node_cstr_alloc()
 * @see mpack_node_utf8_cstr_alloc()
 */
const char* mpack_node_str(mpack_node_t node);

/**
 * Returns a pointer to the data contained by this node.
 *
 * @note Strings are not null-terminated! Use one of the cstr functions
 * to get a null-terminated string.
 *
 * The pointer is valid as long as the data backing the tree is valid.
 *
 * If this node is not of a str, bin or ext, @ref mpack_error_type is raised, and
 * @c NULL is returned.
 *
 * @see mpack_node_copy_cstr()
 * @see mpack_node_cstr_alloc()
 * @see mpack_node_utf8_cstr_alloc()
 */
const char* mpack_node_data(mpack_node_t node);

/**
 * Returns a pointer to the data contained by this bin node.
 *
 * The pointer is valid as long as the data backing the tree is valid.
 *
 * If this node is not a bin, @ref mpack_error_type is raised and @c NULL is
 * returned.
 */
const char* mpack_node_bin_data(mpack_node_t node);

/**
 * Copies the bytes contained by this node into the given buffer, returning the
 * number of bytes in the node.
 *
 * @throws mpack_error_type If this node is not a str, bin or ext type
 * @throws mpack_error_too_big If the string does not fit in the given buffer
 *
 * @param node The string node from which to copy data
 * @param buffer A buffer in which to copy the node's bytes
 * @param bufsize The size of the given buffer
 *
 * @return The number of bytes in the node, or zero if an error occurs.
 */
size_t mpack_node_copy_data(mpack_node_t node, char* buffer, size_t bufsize);

/**
 * Checks that the given node contains a valid UTF-8 string and copies the
 * string into the given buffer, returning the number of bytes in the string.
 *
 * @throws mpack_error_type If this node is not a string
 * @throws mpack_error_too_big If the string does not fit in the given buffer
 *
 * @param node The string node from which to copy data
 * @param buffer A buffer in which to copy the node's bytes
 * @param bufsize The size of the given buffer
 *
 * @return The number of bytes in the node, or zero if an error occurs.
 */
size_t mpack_node_copy_utf8(mpack_node_t node, char* buffer, size_t bufsize);

/**
 * Checks that the given node contains a string with no NUL bytes, copies the string
 * into the given buffer, and adds a null terminator.
 *
 * If this node is not of a string type, @ref mpack_error_type is raised. If the string
 * does not fit, @ref mpack_error_data is raised.
 *
 * If any error occurs, the buffer will contain an empty null-terminated string.
 *
 * @param node The string node from which to copy data
 * @param buffer A buffer in which to copy the node's string
 * @param size The size of the given buffer
 */
void mpack_node_copy_cstr(mpack_node_t node, char* buffer, size_t size);

/**
 * Checks that the given node contains a valid UTF-8 string with no NUL bytes,
 * copies the string into the given buffer, and adds a null terminator.
 *
 * If this node is not of a string type, @ref mpack_error_type is raised. If the string
 * does not fit, @ref mpack_error_data is raised.
 *
 * If any error occurs, the buffer will contain an empty null-terminated string.
 *
 * @param node The string node from which to copy data
 * @param buffer A buffer in which to copy the node's string
 * @param size The size of the given buffer
 */
void mpack_node_copy_utf8_cstr(mpack_node_t node, char* buffer, size_t size);

#ifdef MPACK_MALLOC
/**
 * Allocates a new chunk of data using MPACK_MALLOC() with the bytes
 * contained by this node.
 *
 * The allocated data must be freed with MPACK_FREE() (or simply free()
 * if MPack's allocator hasn't been customized.)
 *
 * @throws mpack_error_type If this node is not a str, bin or ext type
 * @throws mpack_error_too_big If the size of the data is larger than the
 *     given maximum size
 * @throws mpack_error_memory If an allocation failure occurs
 *
 * @param node The node from which to allocate and copy data
 * @param maxsize The maximum size to allocate
 *
 * @return The allocated data, or NULL if any error occurs.
 */
char* mpack_node_data_alloc(mpack_node_t node, size_t maxsize);

/**
 * Allocates a new null-terminated string using MPACK_MALLOC() with the string
 * contained by this node.
 *
 * The allocated string must be freed with MPACK_FREE() (or simply free()
 * if MPack's allocator hasn't been customized.)
 *
 * @throws mpack_error_type If this node is not a string or contains NUL bytes
 * @throws mpack_error_too_big If the size of the string plus null-terminator
 *     is larger than the given maximum size
 * @throws mpack_error_memory If an allocation failure occurs
 *
 * @param node The node from which to allocate and copy string data
 * @param maxsize The maximum size to allocate, including the null-terminator
 *
 * @return The allocated string, or NULL if any error occurs.
 */
char* mpack_node_cstr_alloc(mpack_node_t node, size_t maxsize);

/**
 * Allocates a new null-terminated string using MPACK_MALLOC() with the UTF-8
 * string contained by this node.
 *
 * The allocated string must be freed with MPACK_FREE() (or simply free()
 * if MPack's allocator hasn't been customized.)
 *
 * @throws mpack_error_type If this node is not a string, is not valid UTF-8,
 *     or contains NUL bytes
 * @throws mpack_error_too_big If the size of the string plus null-terminator
 *     is larger than the given maximum size
 * @throws mpack_error_memory If an allocation failure occurs
 *
 * @param node The node from which to allocate and copy string data
 * @param maxsize The maximum size to allocate, including the null-terminator
 *
 * @return The allocated string, or NULL if any error occurs.
 */
char* mpack_node_utf8_cstr_alloc(mpack_node_t node, size_t maxsize);
#endif

/**
 * Searches the given string array for a string matching the given
 * node and returns its index.
 *
 * If the node does not match any of the given strings,
 * @ref mpack_error_type is flagged. Use mpack_node_enum_optional()
 * if you want to allow values other than the given strings.
 *
 * If any error occurs or if the tree is in an error state, @a count
 * is returned.
 *
 * This can be used to quickly parse a string into an enum when the
 * enum values range from 0 to @a count-1. If the last value in the
 * enum is a special "count" value, it can be passed as the count,
 * and the return value can be cast directly to the enum type.
 *
 * @code{.c}
 * typedef enum           { APPLE ,  BANANA ,  ORANGE , COUNT} fruit_t;
 * const char* fruits[] = {"apple", "banana", "orange"};
 *
 * fruit_t fruit = (fruit_t)mpack_node_enum(node, fruits, COUNT);
 * @endcode
 *
 * @param node The node
 * @param strings An array of expected strings of length count
 * @param count The number of strings
 * @return The index of the matched string, or @a count in case of error
 */
size_t mpack_node_enum(mpack_node_t node, const char* strings[], size_t count);

/**
 * Searches the given string array for a string matching the given node,
 * returning its index or @a count if no strings match.
 *
 * If the value is not a string, or it does not match any of the
 * given strings, @a count is returned and no error is flagged.
 *
 * If any error occurs or if the tree is in an error state, @a count
 * is returned.
 *
 * This can be used to quickly parse a string into an enum when the
 * enum values range from 0 to @a count-1. If the last value in the
 * enum is a special "count" value, it can be passed as the count,
 * and the return value can be cast directly to the enum type.
 *
 * @code{.c}
 * typedef enum           { APPLE ,  BANANA ,  ORANGE , COUNT} fruit_t;
 * const char* fruits[] = {"apple", "banana", "orange"};
 *
 * fruit_t fruit = (fruit_t)mpack_node_enum_optional(node, fruits, COUNT);
 * @endcode
 *
 * @param node The node
 * @param strings An array of expected strings of length count
 * @param count The number of strings
 * @return The index of the matched string, or @a count in case of error
 */
size_t mpack_node_enum_optional(mpack_node_t node, const char* strings[], size_t count);

/**
 * @}
 */

/**
 * @name Compound Node Functions
 * @{
 */

/**
 * Returns the length of the given array node. Raises mpack_error_type
 * and returns 0 if the given node is not an array.
 */
size_t mpack_node_array_length(mpack_node_t node);

/**
 * Returns the node in the given array at the given index. If the node
 * is not an array, @ref mpack_error_type is raised and a nil node is returned.
 * If the given index is out of bounds, @ref mpack_error_data is raised and
 * a nil node is returned.
 */
mpack_node_t mpack_node_array_at(mpack_node_t node, size_t index);

/**
 * Returns the number of key/value pairs in the given map node. Raises
 * mpack_error_type and returns 0 if the given node is not a map.
 */
size_t mpack_node_map_count(mpack_node_t node);

/**
 * Returns the key node in the given map at the given index.
 *
 * A nil node is returned in case of error.
 *
 * @throws mpack_error_type if the node is not a map
 * @throws mpack_error_data if the given index is out of bounds
 */
mpack_node_t mpack_node_map_key_at(mpack_node_t node, size_t index);

/**
 * Returns the value node in the given map at the given index.
 *
 * A nil node is returned in case of error.
 *
 * @throws mpack_error_type if the node is not a map
 * @throws mpack_error_data if the given index is out of bounds
 */
mpack_node_t mpack_node_map_value_at(mpack_node_t node, size_t index);

/**
 * Returns the value node in the given map for the given integer key.
 *
 * The key must exist within the map. Use mpack_node_map_int_optional() to
 * check for optional keys.
 *
 * The key must be unique. An error is flagged if the node has multiple
 * entries with the given key.
 *
 * @throws mpack_error_type If the node is not a map
 * @throws mpack_error_data If the node does not contain exactly one entry with the given key
 *
 * @return The value node for the given key, or a nil node in case of error
 */
mpack_node_t mpack_node_map_int(mpack_node_t node, int64_t num);

/**
 * Returns the value node in the given map for the given integer key, or a
 * missing node if the map does not contain the given key.
 *
 * The key must be unique. An error is flagged if the node has multiple
 * entries with the given key.
 *
 * @throws mpack_error_type If the node is not a map
 * @throws mpack_error_data If the node contains more than one entry with the given key
 *
 * @return The value node for the given key, or a missing node if the key does
 *         not exist, or a nil node in case of error
 *
 * @see mpack_node_is_missing()
 */
mpack_node_t mpack_node_map_int_optional(mpack_node_t node, int64_t num);

/**
 * Returns the value node in the given map for the given unsigned integer key.
 *
 * The key must exist within the map. Use mpack_node_map_uint_optional() to
 * check for optional keys.
 *
 * The key must be unique. An error is flagged if the node has multiple
 * entries with the given key.
 *
 * @throws mpack_error_type If the node is not a map
 * @throws mpack_error_data If the node does not contain exactly one entry with the given key
 *
 * @return The value node for the given key, or a nil node in case of error
 */
mpack_node_t mpack_node_map_uint(mpack_node_t node, uint64_t num);

/**
 * Returns the value node in the given map for the given unsigned integer
 * key, or a missing node if the map does not contain the given key.
 *
 * The key must be unique. An error is flagged if the node has multiple
 * entries with the given key.
 *
 * @throws mpack_error_type If the node is not a map
 * @throws mpack_error_data If the node contains more than one entry with the given key
 *
 * @return The value node for the given key, or a missing node if the key does
 *         not exist, or a nil node in case of error
 *
 * @see mpack_node_is_missing()
 */
mpack_node_t mpack_node_map_uint_optional(mpack_node_t node, uint64_t num);

/**
 * Returns the value node in the given map for the given string key.
 *
 * The key must exist within the map. Use mpack_node_map_str_optional() to
 * check for optional keys.
 *
 * The key must be unique. An error is flagged if the node has multiple
 * entries with the given key.
 *
 * @throws mpack_error_type If the node is not a map
 * @throws mpack_error_data If the node does not contain exactly one entry with the given key
 *
 * @return The value node for the given key, or a nil node in case of error
 */
mpack_node_t mpack_node_map_str(mpack_node_t node, const char* str, size_t length);

/**
 * Returns the value node in the given map for the given string key, or a missing
 * node if the map does not contain the given key.
 *
 * The key must be unique. An error is flagged if the node has multiple
 * entries with the given key.
 *
 * @throws mpack_error_type If the node is not a map
 * @throws mpack_error_data If the node contains more than one entry with the given key
 *
 * @return The value node for the given key, or a missing node if the key does
 *         not exist, or a nil node in case of error
 *
 * @see mpack_node_is_missing()
 */
mpack_node_t mpack_node_map_str_optional(mpack_node_t node, const char* str, size_t length);

/**
 * Returns the value node in the given map for the given null-terminated
 * string key.
 *
 * The key must exist within the map. Use mpack_node_map_cstr_optional() to
 * check for optional keys.
 *
 * The key must be unique. An error is flagged if the node has multiple
 * entries with the given key.
 *
 * @throws mpack_error_type If the node is not a map
 * @throws mpack_error_data If the node does not contain exactly one entry with the given key
 *
 * @return The value node for the given key, or a nil node in case of error
 */
mpack_node_t mpack_node_map_cstr(mpack_node_t node, const char* cstr);

/**
 * Returns the value node in the given map for the given null-terminated
 * string key, or a missing node if the map does not contain the given key.
 *
 * The key must be unique. An error is flagged if the node has multiple
 * entries with the given key.
 *
 * @throws mpack_error_type If the node is not a map
 * @throws mpack_error_data If the node contains more than one entry with the given key
 *
 * @return The value node for the given key, or a missing node if the key does
 *         not exist, or a nil node in case of error
 *
 * @see mpack_node_is_missing()
 */
mpack_node_t mpack_node_map_cstr_optional(mpack_node_t node, const char* cstr);

/**
 * Returns true if the given node map contains exactly one entry with the
 * given integer key.
 *
 * The key must be unique. An error is flagged if the node has multiple
 * entries with the given key.
 *
 * @throws mpack_error_type If the node is not a map
 * @throws mpack_error_data If the node contains more than one entry with the given key
 */
bool mpack_node_map_contains_int(mpack_node_t node, int64_t num);

/**
 * Returns true if the given node map contains exactly one entry with the
 * given unsigned integer key.
 *
 * The key must be unique. An error is flagged if the node has multiple
 * entries with the given key.
 *
 * @throws mpack_error_type If the node is not a map
 * @throws mpack_error_data If the node contains more than one entry with the given key
 */
bool mpack_node_map_contains_uint(mpack_node_t node, uint64_t num);

/**
 * Returns true if the given node map contains exactly one entry with the
 * given string key.
 *
 * The key must be unique. An error is flagged if the node has multiple
 * entries with the given key.
 *
 * @throws mpack_error_type If the node is not a map
 * @throws mpack_error_data If the node contains more than one entry with the given key
 */
bool mpack_node_map_contains_str(mpack_node_t node, const char* str, size_t length);

/**
 * Returns true if the given node map contains exactly one entry with the
 * given null-terminated string key.
 *
 * The key must be unique. An error is flagged if the node has multiple
 * entries with the given key.
 *
 * @throws mpack_error_type If the node is not a map
 * @throws mpack_error_data If the node contains more than one entry with the given key
 */
bool mpack_node_map_contains_cstr(mpack_node_t node, const char* cstr);

/**
 * @}
 */

/**
 * @}
 */

#endif

MPACK_EXTERN_C_END
MPACK_SILENCE_WARNINGS_END

#endif


#endif

// ended inlining mpack.h 


/* mpack/mpack-platform.c.c */


// We define MPACK_EMIT_INLINE_DEFS and include mpack.h to emit
// standalone definitions of all (non-static) inline functions in MPack.

#define MPACK_INTERNAL 1
#define MPACK_EMIT_INLINE_DEFS 1

/* #include "mpack-platform.h" */
/* #include "mpack.h" */

MPACK_SILENCE_WARNINGS_BEGIN

#if MPACK_DEBUG

#if MPACK_STDIO
void mpack_assert_fail_format(const char* format, ...) {
    char buffer[512];
    va_list args;
    va_start(args, format);
    vsnprintf(buffer, sizeof(buffer), format, args);
    va_end(args);
    buffer[sizeof(buffer) - 1] = 0;
    mpack_assert_fail_wrapper(buffer);
}

void mpack_break_hit_format(const char* format, ...) {
    char buffer[512];
    va_list args;
    va_start(args, format);
    vsnprintf(buffer, sizeof(buffer), format, args);
    va_end(args);
    buffer[sizeof(buffer) - 1] = 0;
    mpack_break_hit(buffer);
}
#endif

#if !MPACK_CUSTOM_ASSERT
void mpack_assert_fail(const char* message) {
    MPACK_UNUSED(message);

    #if MPACK_STDIO
    fprintf(stderr, "%s\n", message);
    #endif
}
#endif

// We split the assert failure from the wrapper so that a
// custom assert function can return.
void mpack_assert_fail_wrapper(const char* message) {

    #ifdef MPACK_GCOV
    // gcov marks even __builtin_unreachable() as an uncovered line. this
    // silences it.
    (mpack_assert_fail(message), __builtin_unreachable());

    #else
    mpack_assert_fail(message);

    // mpack_assert_fail() is not supposed to return. in case it does, we
    // abort.

    #if !MPACK_NO_BUILTINS
    #if defined(__GNUC__) || defined(__clang__)
    __builtin_trap();
    #elif defined(WIN32)
    __debugbreak();
    #endif
    #endif

    #if (defined(__GNUC__) || defined(__clang__)) && !MPACK_NO_BUILTINS
    __builtin_abort();
    #elif MPACK_STDLIB
    abort();
    #endif

    MPACK_UNREACHABLE;
    #endif
}

#if !MPACK_CUSTOM_BREAK

// If we have a custom assert handler, break wraps it by default.
// This allows users of MPack to only implement mpack_assert_fail() without
// having to worry about the difference between assert and break.
//
// MPACK_CUSTOM_BREAK is available to define a separate break handler
// (which is needed by the unit test suite), but this is not offered in
// mpack-config.h for simplicity.

#if MPACK_CUSTOM_ASSERT
void mpack_break_hit(const char* message) {
    mpack_assert_fail_wrapper(message);
}
#else
void mpack_break_hit(const char* message) {
    MPACK_UNUSED(message);

    #if MPACK_STDIO
    fprintf(stderr, "%s\n", message);
    #endif

    #if defined(__GNUC__) || defined(__clang__) && !MPACK_NO_BUILTINS
    __builtin_trap();
    #elif defined(WIN32) && !MPACK_NO_BUILTINS
    __debugbreak();
    #elif MPACK_STDLIB
    abort();
    #endif
}
#endif

#endif

#endif



// The below are adapted from the C wikibook:
//     https://en.wikibooks.org/wiki/C_Programming/Strings

#ifndef mpack_memcmp
int mpack_memcmp(const void* s1, const void* s2, size_t n) {
     const unsigned char *us1 = (const unsigned char *) s1;
     const unsigned char *us2 = (const unsigned char *) s2;
     while (n-- != 0) {
         if (*us1 != *us2)
             return (*us1 < *us2) ? -1 : +1;
         us1++;
         us2++;
     }
     return 0;
}
#endif

#ifndef mpack_memcpy
void* mpack_memcpy(void* MPACK_RESTRICT s1, const void* MPACK_RESTRICT s2, size_t n) {
    char* MPACK_RESTRICT dst = (char *)s1;
    const char* MPACK_RESTRICT src = (const char *)s2;
    while (n-- != 0)
        *dst++ = *src++;
    return s1;
}
#endif

#ifndef mpack_memmove
void* mpack_memmove(void* s1, const void* s2, size_t n) {
    char *p1 = (char *)s1;
    const char *p2 = (const char *)s2;
    if (p2 < p1 && p1 < p2 + n) {
        p2 += n;
        p1 += n;
        while (n-- != 0)
            *--p1 = *--p2;
    } else
        while (n-- != 0)
            *p1++ = *p2++;
    return s1;
}
#endif

#ifndef mpack_memset
void* mpack_memset(void* s, int c, size_t n) {
    unsigned char *us = (unsigned char *)s;
    unsigned char uc = (unsigned char)c;
    while (n-- != 0)
        *us++ = uc;
    return s;
}
#endif

#ifndef mpack_strlen
size_t mpack_strlen(const char* s) {
    const char* p = s;
    while (*p != '\0')
        p++;
    return (size_t)(p - s);
}
#endif



#if defined(MPACK_MALLOC) && !defined(MPACK_REALLOC)
void* mpack_realloc(void* old_ptr, size_t used_size, size_t new_size) {
    if (new_size == 0) {
        if (old_ptr)
            MPACK_FREE(old_ptr);
        return NULL;
    }

    void* new_ptr = MPACK_MALLOC(new_size);
    if (new_ptr == NULL)
        return NULL;

    mpack_memcpy(new_ptr, old_ptr, used_size);
    MPACK_FREE(old_ptr);
    return new_ptr;
}
#endif

MPACK_SILENCE_WARNINGS_END

/* mpack/mpack-common.c.c */

#define MPACK_INTERNAL 1

/* #include "mpack-common.h" */

MPACK_SILENCE_WARNINGS_BEGIN

const char* mpack_error_to_string(mpack_error_t error) {
    #if MPACK_STRINGS
    switch (error) {
        #define MPACK_ERROR_STRING_CASE(e) case e: return #e
        MPACK_ERROR_STRING_CASE(mpack_ok);
        MPACK_ERROR_STRING_CASE(mpack_error_io);
        MPACK_ERROR_STRING_CASE(mpack_error_invalid);
        MPACK_ERROR_STRING_CASE(mpack_error_unsupported);
        MPACK_ERROR_STRING_CASE(mpack_error_type);
        MPACK_ERROR_STRING_CASE(mpack_error_too_big);
        MPACK_ERROR_STRING_CASE(mpack_error_memory);
        MPACK_ERROR_STRING_CASE(mpack_error_bug);
        MPACK_ERROR_STRING_CASE(mpack_error_data);
        MPACK_ERROR_STRING_CASE(mpack_error_eof);
        #undef MPACK_ERROR_STRING_CASE
    }
    mpack_assert(0, "unrecognized error %i", (int)error);
    return "(unknown mpack_error_t)";
    #else
    MPACK_UNUSED(error);
    return "";
    #endif
}

const char* mpack_type_to_string(mpack_type_t type) {
    #if MPACK_STRINGS
    switch (type) {
        #define MPACK_TYPE_STRING_CASE(e) case e: return #e
        MPACK_TYPE_STRING_CASE(mpack_type_missing);
        MPACK_TYPE_STRING_CASE(mpack_type_nil);
        MPACK_TYPE_STRING_CASE(mpack_type_bool);
        MPACK_TYPE_STRING_CASE(mpack_type_float);
        MPACK_TYPE_STRING_CASE(mpack_type_double);
        MPACK_TYPE_STRING_CASE(mpack_type_int);
        MPACK_TYPE_STRING_CASE(mpack_type_uint);
        MPACK_TYPE_STRING_CASE(mpack_type_str);
        MPACK_TYPE_STRING_CASE(mpack_type_bin);
        MPACK_TYPE_STRING_CASE(mpack_type_array);
        MPACK_TYPE_STRING_CASE(mpack_type_map);
        #if MPACK_EXTENSIONS
        MPACK_TYPE_STRING_CASE(mpack_type_ext);
        #endif
        #undef MPACK_TYPE_STRING_CASE
    }
    mpack_assert(0, "unrecognized type %i", (int)type);
    return "(unknown mpack_type_t)";
    #else
    MPACK_UNUSED(type);
    return "";
    #endif
}

int mpack_tag_cmp(mpack_tag_t left, mpack_tag_t right) {

    // positive numbers may be stored as int; convert to uint
    if (left.type == mpack_type_int && left.v.i >= 0) {
        left.type = mpack_type_uint;
        left.v.u = (uint64_t)left.v.i;
    }
    if (right.type == mpack_type_int && right.v.i >= 0) {
        right.type = mpack_type_uint;
        right.v.u = (uint64_t)right.v.i;
    }

    if (left.type != right.type)
        return ((int)left.type < (int)right.type) ? -1 : 1;

    switch (left.type) {
        case mpack_type_missing: // fallthrough
        case mpack_type_nil:
            return 0;

        case mpack_type_bool:
            return (int)left.v.b - (int)right.v.b;

        case mpack_type_int:
            if (left.v.i == right.v.i)
                return 0;
            return (left.v.i < right.v.i) ? -1 : 1;

        case mpack_type_uint:
            if (left.v.u == right.v.u)
                return 0;
            return (left.v.u < right.v.u) ? -1 : 1;

        case mpack_type_array:
        case mpack_type_map:
            if (left.v.n == right.v.n)
                return 0;
            return (left.v.n < right.v.n) ? -1 : 1;

        case mpack_type_str:
        case mpack_type_bin:
            if (left.v.l == right.v.l)
                return 0;
            return (left.v.l < right.v.l) ? -1 : 1;

        #if MPACK_EXTENSIONS
        case mpack_type_ext:
            if (left.exttype == right.exttype) {
                if (left.v.l == right.v.l)
                    return 0;
                return (left.v.l < right.v.l) ? -1 : 1;
            }
            return (int)left.exttype - (int)right.exttype;
        #endif

        // floats should not normally be compared for equality. we compare
        // with memcmp() to silence compiler warnings, but this will return
        // equal if both are NaNs with the same representation (though we may
        // want this, for instance if you are for some bizarre reason using
        // floats as map keys.) i'm not sure what the right thing to
        // do is here. check for NaN first? always return false if the type
        // is float? use operator== and pragmas to silence compiler warning?
        // please send me your suggestions.
        // note also that we don't convert floats to doubles, so when this is
        // used for ordering purposes, all floats are ordered before all
        // doubles.
        case mpack_type_float:
            return mpack_memcmp(&left.v.f, &right.v.f, sizeof(left.v.f));
        case mpack_type_double:
            return mpack_memcmp(&left.v.d, &right.v.d, sizeof(left.v.d));
    }

    mpack_assert(0, "unrecognized type %i", (int)left.type);
    return false;
}

#if MPACK_DEBUG && MPACK_STDIO
static char mpack_hex_char(uint8_t hex_value) {
    // Older compilers (e.g. GCC 4.4.7) promote the result of this ternary to
    // int and warn under -Wconversion, so we have to cast it back to char.
    return (char)((hex_value < 10) ? (char)('0' + hex_value) : (char)('a' + (hex_value - 10)));
}

static void mpack_tag_debug_complete_bin_ext(mpack_tag_t tag, size_t string_length, char* buffer, size_t buffer_size,
        const char* prefix, size_t prefix_size)
{
    // If at any point in this function we run out of space in the buffer, we
    // bail out. The outer tag print wrapper will make sure we have a
    // null-terminator.

    if (string_length == 0 || string_length >= buffer_size)
        return;
    buffer += string_length;
    buffer_size -= string_length;

    size_t total = mpack_tag_bytes(&tag);
    if (total == 0) {
        strncpy(buffer, ">", buffer_size);
        return;
    }

    strncpy(buffer, ": ", buffer_size);
    if (buffer_size < 2)
        return;
    buffer += 2;
    buffer_size -= 2;

    size_t hex_bytes = 0;
    size_t i;
    for (i = 0; i < MPACK_PRINT_BYTE_COUNT && i < prefix_size && buffer_size > 2; ++i) {
        uint8_t byte = (uint8_t)prefix[i];
        buffer[0] = mpack_hex_char((uint8_t)(byte >> 4));
        buffer[1] = mpack_hex_char((uint8_t)(byte & 0xfu));
        buffer += 2;
        buffer_size -= 2;
        ++hex_bytes;
    }

    if (buffer_size != 0)
        mpack_snprintf(buffer, buffer_size, "%s>", (total > hex_bytes) ? "..." : "");
}

static void mpack_tag_debug_pseudo_json_bin(mpack_tag_t tag, char* buffer, size_t buffer_size,
        const char* prefix, size_t prefix_size)
{
    mpack_assert(mpack_tag_type(&tag) == mpack_type_bin);
    size_t length = (size_t)mpack_snprintf(buffer, buffer_size, "<binary data of length %" PRIu32 "", tag.v.l);
    mpack_tag_debug_complete_bin_ext(tag, length, buffer, buffer_size, prefix, prefix_size);
}

#if MPACK_EXTENSIONS
static void mpack_tag_debug_pseudo_json_ext(mpack_tag_t tag, char* buffer, size_t buffer_size,
        const char* prefix, size_t prefix_size)
{
    mpack_assert(mpack_tag_type(&tag) == mpack_type_ext);
    size_t length = (size_t)mpack_snprintf(buffer, buffer_size, "<ext data of type %i and length %" PRIu32 "",
            mpack_tag_ext_exttype(&tag), mpack_tag_ext_length(&tag));
    mpack_tag_debug_complete_bin_ext(tag, length, buffer, buffer_size, prefix, prefix_size);
}
#endif

static void mpack_tag_debug_pseudo_json_impl(mpack_tag_t tag, char* buffer, size_t buffer_size,
        const char* prefix, size_t prefix_size)
{
    switch (tag.type) {
        case mpack_type_missing:
            mpack_snprintf(buffer, buffer_size, "<missing!>");
            return;
        case mpack_type_nil:
            mpack_snprintf(buffer, buffer_size, "null");
            return;
        case mpack_type_bool:
            mpack_snprintf(buffer, buffer_size, tag.v.b ? "true" : "false");
            return;
        case mpack_type_int:
            mpack_snprintf(buffer, buffer_size, "%" PRIi64, tag.v.i);
            return;
        case mpack_type_uint:
            mpack_snprintf(buffer, buffer_size, "%" PRIu64, tag.v.u);
            return;
        case mpack_type_float:
            #if MPACK_FLOAT
            mpack_snprintf(buffer, buffer_size, "%f", tag.v.f);
            #else
            mpack_snprintf(buffer, buffer_size, "<float>");
            #endif
            return;
        case mpack_type_double:
            #if MPACK_DOUBLE
            mpack_snprintf(buffer, buffer_size, "%f", tag.v.d);
            #else
            mpack_snprintf(buffer, buffer_size, "<double>");
            #endif
            return;

        case mpack_type_str:
            mpack_snprintf(buffer, buffer_size, "<string of %" PRIu32 " bytes>", tag.v.l);
            return;
        case mpack_type_bin:
            mpack_tag_debug_pseudo_json_bin(tag, buffer, buffer_size, prefix, prefix_size);
            return;
        #if MPACK_EXTENSIONS
        case mpack_type_ext:
            mpack_tag_debug_pseudo_json_ext(tag, buffer, buffer_size, prefix, prefix_size);
            return;
        #endif

        case mpack_type_array:
            mpack_snprintf(buffer, buffer_size, "<array of %" PRIu32 " elements>", tag.v.n);
            return;
        case mpack_type_map:
            mpack_snprintf(buffer, buffer_size, "<map of %" PRIu32 " key-value pairs>", tag.v.n);
            return;
    }

    mpack_snprintf(buffer, buffer_size, "<unknown!>");
}

void mpack_tag_debug_pseudo_json(mpack_tag_t tag, char* buffer, size_t buffer_size,
        const char* prefix, size_t prefix_size)
{
    mpack_assert(buffer_size > 0, "buffer size cannot be zero!");
    buffer[0] = 0;

    mpack_tag_debug_pseudo_json_impl(tag, buffer, buffer_size, prefix, prefix_size);

    // We always null-terminate the buffer manually just in case the snprintf()
    // function doesn't null-terminate when the string doesn't fit.
    buffer[buffer_size - 1] = 0;
}

static void mpack_tag_debug_describe_impl(mpack_tag_t tag, char* buffer, size_t buffer_size) {
    switch (tag.type) {
        case mpack_type_missing:
            mpack_snprintf(buffer, buffer_size, "missing");
            return;
        case mpack_type_nil:
            mpack_snprintf(buffer, buffer_size, "nil");
            return;
        case mpack_type_bool:
            mpack_snprintf(buffer, buffer_size, tag.v.b ? "true" : "false");
            return;
        case mpack_type_int:
            mpack_snprintf(buffer, buffer_size, "int %" PRIi64, tag.v.i);
            return;
        case mpack_type_uint:
            mpack_snprintf(buffer, buffer_size, "uint %" PRIu64, tag.v.u);
            return;
        case mpack_type_float:
            #if MPACK_FLOAT
            mpack_snprintf(buffer, buffer_size, "float %f", tag.v.f);
            #else
            mpack_snprintf(buffer, buffer_size, "float");
            #endif
            return;
        case mpack_type_double:
            #if MPACK_DOUBLE
            mpack_snprintf(buffer, buffer_size, "double %f", tag.v.d);
            #else
            mpack_snprintf(buffer, buffer_size, "double");
            #endif
            return;
        case mpack_type_str:
            mpack_snprintf(buffer, buffer_size, "str of %" PRIu32 " bytes", tag.v.l);
            return;
        case mpack_type_bin:
            mpack_snprintf(buffer, buffer_size, "bin of %" PRIu32 " bytes", tag.v.l);
            return;
        #if MPACK_EXTENSIONS
        case mpack_type_ext:
            mpack_snprintf(buffer, buffer_size, "ext of type %i, %" PRIu32 " bytes",
                    mpack_tag_ext_exttype(&tag), mpack_tag_ext_length(&tag));
            return;
        #endif
        case mpack_type_array:
            mpack_snprintf(buffer, buffer_size, "array of %" PRIu32 " elements", tag.v.n);
            return;
        case mpack_type_map:
            mpack_snprintf(buffer, buffer_size, "map of %" PRIu32 " key-value pairs", tag.v.n);
            return;
    }

    mpack_snprintf(buffer, buffer_size, "unknown!");
}

void mpack_tag_debug_describe(mpack_tag_t tag, char* buffer, size_t buffer_size) {
    mpack_assert(buffer_size > 0, "buffer size cannot be zero!");
    buffer[0] = 0;

    mpack_tag_debug_describe_impl(tag, buffer, buffer_size);

    // We always null-terminate the buffer manually just in case the snprintf()
    // function doesn't null-terminate when the string doesn't fit.
    buffer[buffer_size - 1] = 0;
}
#endif



#if MPACK_READ_TRACKING || MPACK_WRITE_TRACKING

#ifndef MPACK_TRACKING_INITIAL_CAPACITY
// seems like a reasonable number. we grow by doubling, and it only
// needs to be as long as the maximum depth of the message.
#define MPACK_TRACKING_INITIAL_CAPACITY 8
#endif

mpack_error_t mpack_track_init(mpack_track_t* track) {
    track->count = 0;
    track->capacity = MPACK_TRACKING_INITIAL_CAPACITY;
    track->elements = (mpack_track_element_t*)MPACK_MALLOC(sizeof(mpack_track_element_t) * track->capacity);
    if (track->elements == NULL)
        return mpack_error_memory;
    return mpack_ok;
}

mpack_error_t mpack_track_grow(mpack_track_t* track) {
    mpack_assert(track->elements, "null track elements!");
    mpack_assert(track->count == track->capacity, "incorrect growing?");

    size_t new_capacity = track->capacity * 2;

    mpack_track_element_t* new_elements = (mpack_track_element_t*)mpack_realloc(track->elements,
            sizeof(mpack_track_element_t) * track->count, sizeof(mpack_track_element_t) * new_capacity);
    if (new_elements == NULL)
        return mpack_error_memory;

    track->elements = new_elements;
    track->capacity = new_capacity;
    return mpack_ok;
}

mpack_error_t mpack_track_push(mpack_track_t* track, mpack_type_t type, uint32_t count) {
    mpack_assert(track->elements, "null track elements!");
    mpack_log("track pushing %s count %i\n", mpack_type_to_string(type), (int)count);

    // grow if needed
    if (track->count == track->capacity) {
        mpack_error_t error = mpack_track_grow(track);
        if (error != mpack_ok)
            return error;
    }

    // insert new track
    track->elements[track->count].type = type;
    track->elements[track->count].left = count;
    track->elements[track->count].builder = false;
    track->elements[track->count].key_needs_value = false;
    ++track->count;
    return mpack_ok;
}

// TODO dedupe this
mpack_error_t mpack_track_push_builder(mpack_track_t* track, mpack_type_t type) {
    mpack_assert(track->elements, "null track elements!");
    mpack_log("track pushing %s builder\n", mpack_type_to_string(type));

    // grow if needed
    if (track->count == track->capacity) {
        mpack_error_t error = mpack_track_grow(track);
        if (error != mpack_ok)
            return error;
    }

    // insert new track
    track->elements[track->count].type = type;
    track->elements[track->count].left = 0;
    track->elements[track->count].builder = true;
    track->elements[track->count].key_needs_value = false;
    ++track->count;
    return mpack_ok;
}

static mpack_error_t mpack_track_pop_impl(mpack_track_t* track, mpack_type_t type, bool builder) {
    mpack_assert(track->elements, "null track elements!");
    mpack_log("track popping %s\n", mpack_type_to_string(type));

    if (track->count == 0) {
        mpack_break("attempting to close a %s but nothing was opened!", mpack_type_to_string(type));
        return mpack_error_bug;
    }

    mpack_track_element_t* element = &track->elements[track->count - 1];

    if (element->type != type) {
        mpack_break("attempting to close a %s but the open element is a %s!",
                mpack_type_to_string(type), mpack_type_to_string(element->type));
        return mpack_error_bug;
    }

    if (element->key_needs_value) {
        mpack_assert(type == mpack_type_map, "key_needs_value can only be true for maps!");
        mpack_break("attempting to close a %s but an odd number of elements were written",
                mpack_type_to_string(type));
        return mpack_error_bug;
    }

    if (element->left != 0) {
        mpack_break("attempting to close a %s but there are %i %s left",
                mpack_type_to_string(type), element->left,
                (type == mpack_type_map || type == mpack_type_array) ? "elements" : "bytes");
        return mpack_error_bug;
    }

    if (element->builder != builder) {
        mpack_break("attempting to pop a %sbuilder but the open element is %sa builder",
                builder ? "" : "non-",
                element->builder ? "" : "not ");
        return mpack_error_bug;
    }

    --track->count;
    return mpack_ok;
}

mpack_error_t mpack_track_pop(mpack_track_t* track, mpack_type_t type) {
    return mpack_track_pop_impl(track, type, false);
}

mpack_error_t mpack_track_pop_builder(mpack_track_t* track, mpack_type_t type) {
    return mpack_track_pop_impl(track, type, true);
}

mpack_error_t mpack_track_peek_element(mpack_track_t* track, bool read) {
    MPACK_UNUSED(read);
    mpack_assert(track->elements, "null track elements!");

    // if there are no open elements, that's fine, we can read/write elements at will
    if (track->count == 0)
        return mpack_ok;

    mpack_track_element_t* element = &track->elements[track->count - 1];

    if (element->type != mpack_type_map && element->type != mpack_type_array) {
        mpack_break("elements cannot be %s within an %s", read ? "read" : "written",
                mpack_type_to_string(element->type));
        return mpack_error_bug;
    }

    if (!element->builder && element->left == 0 && !element->key_needs_value) {
        mpack_break("too many elements %s for %s", read ? "read" : "written",
                mpack_type_to_string(element->type));
        return mpack_error_bug;
    }

    return mpack_ok;
}

mpack_error_t mpack_track_element(mpack_track_t* track, bool read) {
    mpack_error_t error = mpack_track_peek_element(track, read);
    if (track->count == 0 || error != mpack_ok)
        return error;

    mpack_track_element_t* element = &track->elements[track->count - 1];

    if (element->type == mpack_type_map) {
        if (!element->key_needs_value) {
            element->key_needs_value = true;
            return mpack_ok; // don't decrement
        }
        element->key_needs_value = false;
    }

    if (!element->builder)
        --element->left;
    return mpack_ok;
}

mpack_error_t mpack_track_bytes(mpack_track_t* track, bool read, size_t count) {
    MPACK_UNUSED(read);
    mpack_assert(track->elements, "null track elements!");

    if (count > MPACK_UINT32_MAX) {
        mpack_break("%s more bytes than could possibly fit in a str/bin/ext!",
                read ? "reading" : "writing");
        return mpack_error_bug;
    }

    if (track->count == 0) {
        mpack_break("bytes cannot be %s with no open bin, str or ext", read ? "read" : "written");
        return mpack_error_bug;
    }

    mpack_track_element_t* element = &track->elements[track->count - 1];

    if (element->type == mpack_type_map || element->type == mpack_type_array) {
        mpack_break("bytes cannot be %s within an %s", read ? "read" : "written",
                mpack_type_to_string(element->type));
        return mpack_error_bug;
    }

    if (element->left < count) {
        mpack_break("too many bytes %s for %s", read ? "read" : "written",
                mpack_type_to_string(element->type));
        return mpack_error_bug;
    }

    element->left -= (uint32_t)count;
    return mpack_ok;
}

mpack_error_t mpack_track_str_bytes_all(mpack_track_t* track, bool read, size_t count) {
    mpack_error_t error = mpack_track_bytes(track, read, count);
    if (error != mpack_ok)
        return error;

    mpack_track_element_t* element = &track->elements[track->count - 1];

    if (element->type != mpack_type_str) {
        mpack_break("the open type must be a string, not a %s", mpack_type_to_string(element->type));
        return mpack_error_bug;
    }

    if (element->left != 0) {
        mpack_break("not all bytes were read; the wrong byte count was requested for a string read.");
        return mpack_error_bug;
    }

    return mpack_ok;
}

mpack_error_t mpack_track_check_empty(mpack_track_t* track) {
    if (track->count != 0) {
        mpack_break("unclosed %s", mpack_type_to_string(track->elements[0].type));
        return mpack_error_bug;
    }
    return mpack_ok;
}

mpack_error_t mpack_track_destroy(mpack_track_t* track, bool cancel) {
    mpack_error_t error = cancel ? mpack_ok : mpack_track_check_empty(track);
    if (track->elements) {
        MPACK_FREE(track->elements);
        track->elements = NULL;
    }
    return error;
}
#endif



static bool mpack_utf8_check_impl(const uint8_t* str, size_t count, bool allow_null) {
    while (count > 0) {
        uint8_t lead = str[0];

        // NUL
        if (!allow_null && lead == '\0') // we don't allow NUL bytes in MPack C-strings
            return false;

        // ASCII
        if (lead <= 0x7F) {
            ++str;
            --count;

        // 2-byte sequence
        } else if ((lead & 0xE0) == 0xC0) {
            if (count < 2) // truncated sequence
                return false;

            uint8_t cont = str[1];
            if ((cont & 0xC0) != 0x80) // not a continuation byte
                return false;

            str += 2;
            count -= 2;

            uint32_t z = ((uint32_t)(lead & ~0xE0) << 6) |
                          (uint32_t)(cont & ~0xC0);

            if (z < 0x80) // overlong sequence
                return false;

        // 3-byte sequence
        } else if ((lead & 0xF0) == 0xE0) {
            if (count < 3) // truncated sequence
                return false;

            uint8_t cont1 = str[1];
            if ((cont1 & 0xC0) != 0x80) // not a continuation byte
                return false;
            uint8_t cont2 = str[2];
            if ((cont2 & 0xC0) != 0x80) // not a continuation byte
                return false;

            str += 3;
            count -= 3;

            uint32_t z = ((uint32_t)(lead  & ~0xF0) << 12) |
                         ((uint32_t)(cont1 & ~0xC0) <<  6) |
                          (uint32_t)(cont2 & ~0xC0);

            if (z < 0x800) // overlong sequence
                return false;
            if (z >= 0xD800 && z <= 0xDFFF) // surrogate
                return false;

        // 4-byte sequence
        } else if ((lead & 0xF8) == 0xF0) {
            if (count < 4) // truncated sequence
                return false;

            uint8_t cont1 = str[1];
            if ((cont1 & 0xC0) != 0x80) // not a continuation byte
                return false;
            uint8_t cont2 = str[2];
            if ((cont2 & 0xC0) != 0x80) // not a continuation byte
                return false;
            uint8_t cont3 = str[3];
            if ((cont3 & 0xC0) != 0x80) // not a continuation byte
                return false;

            str += 4;
            count -= 4;

            uint32_t z = ((uint32_t)(lead  & ~0xF8) << 18) |
                         ((uint32_t)(cont1 & ~0xC0) << 12) |
                         ((uint32_t)(cont2 & ~0xC0) <<  6) |
                          (uint32_t)(cont3 & ~0xC0);

            if (z < 0x10000) // overlong sequence
                return false;
            if (z > 0x10FFFF) // codepoint limit
                return false;

        } else {
            return false; // continuation byte without a lead, or lead for a 5-byte sequence or longer
        }
    }
    return true;
}

bool mpack_utf8_check(const char* str, size_t bytes) {
    return mpack_utf8_check_impl((const uint8_t*)str, bytes, true);
}

bool mpack_utf8_check_no_null(const char* str, size_t bytes) {
    return mpack_utf8_check_impl((const uint8_t*)str, bytes, false);
}

bool mpack_str_check_no_null(const char* str, size_t bytes) {
    size_t i;
    for (i = 0; i < bytes; ++i)
        if (str[i] == '\0')
            return false;
    return true;
}

#if MPACK_DEBUG && MPACK_STDIO
void mpack_print_append(mpack_print_t* print, const char* data, size_t count) {

    // copy whatever fits into the buffer
    size_t copy = print->size - print->count;
    if (copy > count)
        copy = count;
    mpack_memcpy(print->buffer + print->count, data, copy);
    print->count += copy;
    data += copy;
    count -= copy;

    // if we don't need to flush or can't flush there's nothing else to do
    if (count == 0 || print->callback == NULL)
        return;

    // flush the buffer
    print->callback(print->context, print->buffer, print->count);

    if (count > print->size / 2) {
        // flush the rest of the data
        print->count = 0;
        print->callback(print->context, data, count);
    } else {
        // copy the rest of the data into the buffer
        mpack_memcpy(print->buffer, data, count);
        print->count = count;
    }

}

void mpack_print_flush(mpack_print_t* print) {
    if (print->count > 0 && print->callback != NULL) {
        print->callback(print->context, print->buffer, print->count);
        print->count = 0;
    }
}

void mpack_print_file_callback(void* context, const char* data, size_t count) {
    FILE* file = (FILE*)context;
    fwrite(data, 1, count, file);
}
#endif

MPACK_SILENCE_WARNINGS_END

/* mpack/mpack-writer.c.c */

#define MPACK_INTERNAL 1

/* #include "mpack-writer.h" */

MPACK_SILENCE_WARNINGS_BEGIN

#if MPACK_WRITER

#if MPACK_BUILDER
static void mpack_builder_flush(mpack_writer_t* writer);
#endif

#if MPACK_WRITE_TRACKING
static void mpack_writer_flag_if_error(mpack_writer_t* writer, mpack_error_t error) {
    if (error != mpack_ok)
        mpack_writer_flag_error(writer, error);
}

void mpack_writer_track_push(mpack_writer_t* writer, mpack_type_t type, uint32_t count) {
    if (writer->error == mpack_ok)
        mpack_writer_flag_if_error(writer, mpack_track_push(&writer->track, type, count));
}

void mpack_writer_track_push_builder(mpack_writer_t* writer, mpack_type_t type) {
    if (writer->error == mpack_ok)
        mpack_writer_flag_if_error(writer, mpack_track_push_builder(&writer->track, type));
}

void mpack_writer_track_pop(mpack_writer_t* writer, mpack_type_t type) {
    if (writer->error == mpack_ok)
        mpack_writer_flag_if_error(writer, mpack_track_pop(&writer->track, type));
}

void mpack_writer_track_pop_builder(mpack_writer_t* writer, mpack_type_t type) {
    if (writer->error == mpack_ok)
        mpack_writer_flag_if_error(writer, mpack_track_pop_builder(&writer->track, type));
}

void mpack_writer_track_bytes(mpack_writer_t* writer, size_t count) {
    if (writer->error == mpack_ok)
        mpack_writer_flag_if_error(writer, mpack_track_bytes(&writer->track, false, count));
}
#endif

// This should probably be renamed. It's not solely used for tracking.
static inline void mpack_writer_track_element(mpack_writer_t* writer) {
    (void)writer;

    #if MPACK_WRITE_TRACKING
    if (writer->error == mpack_ok)
        mpack_writer_flag_if_error(writer, mpack_track_element(&writer->track, false));
    #endif

    #if MPACK_BUILDER
    if (writer->builder.current_build != NULL) {
        mpack_build_t* build = writer->builder.current_build;
        // We only track this write if it's not nested within another non-build
        // map or array.
        if (build->nested_compound_elements == 0) {
            if (build->type != mpack_type_map) {
                ++build->count;
                mpack_log("adding element to build %p, now %" PRIu32 " elements\n", (void*)build, build->count);
            } else if (build->key_needs_value) {
                build->key_needs_value = false;
                ++build->count;
            } else {
                build->key_needs_value = true;
            }
        }
    }
    #endif
}

static void mpack_writer_clear(mpack_writer_t* writer) {
    #if MPACK_COMPATIBILITY
    writer->version = mpack_version_current;
    #endif
    writer->flush = NULL;
    writer->error_fn = NULL;
    writer->teardown = NULL;
    writer->context = NULL;

    writer->buffer = NULL;
    writer->position = NULL;
    writer->end = NULL;
    writer->error = mpack_ok;

    #if MPACK_WRITE_TRACKING
    mpack_memset(&writer->track, 0, sizeof(writer->track));
    #endif

    #if MPACK_BUILDER
    writer->builder.current_build = NULL;
    writer->builder.latest_build = NULL;
    writer->builder.current_page = NULL;
    writer->builder.pages = NULL;
    writer->builder.stash_buffer = NULL;
    writer->builder.stash_position = NULL;
    writer->builder.stash_end = NULL;
    #endif
}

void mpack_writer_init(mpack_writer_t* writer, char* buffer, size_t size) {
    mpack_assert(buffer != NULL, "cannot initialize writer with empty buffer");
    mpack_writer_clear(writer);
    writer->buffer = buffer;
    writer->position = buffer;
    writer->end = writer->buffer + size;

    #if MPACK_WRITE_TRACKING
    mpack_writer_flag_if_error(writer, mpack_track_init(&writer->track));
    #endif

    mpack_log("===========================\n");
    mpack_log("initializing writer with buffer size %i\n", (int)size);
}

void mpack_writer_init_error(mpack_writer_t* writer, mpack_error_t error) {
    mpack_writer_clear(writer);
    writer->error = error;

    mpack_log("===========================\n");
    mpack_log("initializing writer in error state %i\n", (int)error);
}

void mpack_writer_set_flush(mpack_writer_t* writer, mpack_writer_flush_t flush) {
    MPACK_STATIC_ASSERT(MPACK_WRITER_MINIMUM_BUFFER_SIZE >= MPACK_MAXIMUM_TAG_SIZE,
            "minimum buffer size must fit any tag!");
    MPACK_STATIC_ASSERT(31 + MPACK_TAG_SIZE_FIXSTR >= MPACK_WRITER_MINIMUM_BUFFER_SIZE,
            "minimum buffer size must fit the largest possible fixstr!");

    if (mpack_writer_buffer_size(writer) < MPACK_WRITER_MINIMUM_BUFFER_SIZE) {
        mpack_break("buffer size is %i, but minimum buffer size for flush is %i",
                (int)mpack_writer_buffer_size(writer), MPACK_WRITER_MINIMUM_BUFFER_SIZE);
        mpack_writer_flag_error(writer, mpack_error_bug);
        return;
    }

    writer->flush = flush;
}

#ifdef MPACK_MALLOC
typedef struct mpack_growable_writer_t {
    char** target_data;
    size_t* target_size;
} mpack_growable_writer_t;

static char* mpack_writer_get_reserved(mpack_writer_t* writer) {
    // This is in a separate function in order to avoid false strict aliasing
    // warnings. We aren't actually violating strict aliasing (the reserved
    // space is only ever dereferenced as an mpack_growable_writer_t.)
    return (char*)writer->reserved;
}

static void mpack_growable_writer_flush(mpack_writer_t* writer, const char* data, size_t count) {

    // This is an intrusive flush function which modifies the writer's buffer
    // in response to a flush instead of emptying it in order to add more
    // capacity for data. This removes the need to copy data from a fixed buffer
    // into a growable one, improving performance.
    //
    // There are three ways flush can be called:
    //   - flushing the buffer during writing (used is zero, count is all data, data is buffer)
    //   - flushing extra data during writing (used is all flushed data, count is extra data, data is not buffer)
    //   - flushing during teardown (used and count are both all flushed data, data is buffer)
    //
    // In the first two cases, we grow the buffer by at least double, enough
    // to ensure that new data will fit. We ignore the teardown flush.

    if (data == writer->buffer) {

        // teardown, do nothing
        if (mpack_writer_buffer_used(writer) == count)
            return;

        // otherwise leave the data in the buffer and just grow
        writer->position = writer->buffer + count;
        count = 0;
    }

    size_t used = mpack_writer_buffer_used(writer);
    size_t size = mpack_writer_buffer_size(writer);

    mpack_log("flush size %i used %i data %p buffer %p\n",
            (int)count, (int)used, data, writer->buffer);

    mpack_assert(data == writer->buffer || used + count > size,
            "extra flush for %i but there is %i space left in the buffer! (%i/%i)",
            (int)count, (int)mpack_writer_buffer_left(writer), (int)used, (int)size);

    // grow to fit the data
    // TODO: this really needs to correctly test for overflow
    size_t new_size = size * 2;
    while (new_size < used + count)
        new_size *= 2;

    mpack_log("flush growing buffer size from %i to %i\n", (int)size, (int)new_size);

    // grow the buffer
    char* new_buffer = (char*)mpack_realloc(writer->buffer, used, new_size);
    if (new_buffer == NULL) {
        mpack_writer_flag_error(writer, mpack_error_memory);
        return;
    }
    writer->position = new_buffer + used;
    writer->buffer = new_buffer;
    writer->end = writer->buffer + new_size;

    // append the extra data
    if (count > 0) {
        mpack_memcpy(writer->position, data, count);
        writer->position += count;
    }

    mpack_log("new buffer %p, used %i\n", new_buffer, (int)mpack_writer_buffer_used(writer));
}

static void mpack_growable_writer_teardown(mpack_writer_t* writer) {
    mpack_growable_writer_t* growable_writer = (mpack_growable_writer_t*)mpack_writer_get_reserved(writer);

    if (mpack_writer_error(writer) == mpack_ok) {

        // shrink the buffer to an appropriate size if the data is
        // much smaller than the buffer
        if (mpack_writer_buffer_used(writer) < mpack_writer_buffer_size(writer) / 2) {
            size_t used = mpack_writer_buffer_used(writer);

            // We always return a non-null pointer that must be freed, even if
            // nothing was written. malloc() and realloc() do not necessarily
            // do this so we enforce it ourselves.
            size_t size = (used != 0) ? used : 1;

            char* buffer = (char*)mpack_realloc(writer->buffer, used, size);
            if (!buffer) {
                MPACK_FREE(writer->buffer);
                mpack_writer_flag_error(writer, mpack_error_memory);
                return;
            }
            writer->buffer = buffer;
            writer->end = (writer->position = writer->buffer + used);
        }

        *growable_writer->target_data = writer->buffer;
        *growable_writer->target_size = mpack_writer_buffer_used(writer);
        writer->buffer = NULL;

    } else if (writer->buffer) {
        MPACK_FREE(writer->buffer);
        writer->buffer = NULL;
    }

    writer->context = NULL;
}

void mpack_writer_init_growable(mpack_writer_t* writer, char** target_data, size_t* target_size) {
    mpack_assert(target_data != NULL, "cannot initialize writer without a destination for the data");
    mpack_assert(target_size != NULL, "cannot initialize writer without a destination for the size");

    *target_data = NULL;
    *target_size = 0;

    MPACK_STATIC_ASSERT(sizeof(mpack_growable_writer_t) <= sizeof(writer->reserved),
            "not enough reserved space for growable writer!");
    mpack_growable_writer_t* growable_writer = (mpack_growable_writer_t*)mpack_writer_get_reserved(writer);

    growable_writer->target_data = target_data;
    growable_writer->target_size = target_size;

    size_t capacity = MPACK_BUFFER_SIZE;
    char* buffer = (char*)MPACK_MALLOC(capacity);
    if (buffer == NULL) {
        mpack_writer_init_error(writer, mpack_error_memory);
        return;
    }

    mpack_writer_init(writer, buffer, capacity);
    mpack_writer_set_flush(writer, mpack_growable_writer_flush);
    mpack_writer_set_teardown(writer, mpack_growable_writer_teardown);
}
#endif

#if MPACK_STDIO
static void mpack_file_writer_flush(mpack_writer_t* writer, const char* buffer, size_t count) {
    FILE* file = (FILE*)writer->context;
    size_t written = fwrite((const void*)buffer, 1, count, file);
    if (written != count)
        mpack_writer_flag_error(writer, mpack_error_io);
}

static void mpack_file_writer_teardown(mpack_writer_t* writer) {
    MPACK_FREE(writer->buffer);
    writer->buffer = NULL;
    writer->context = NULL;
}

static void mpack_file_writer_teardown_close(mpack_writer_t* writer) {
    FILE* file = (FILE*)writer->context;

    if (file) {
        int ret = fclose(file);
        if (ret != 0)
            mpack_writer_flag_error(writer, mpack_error_io);
    }

    mpack_file_writer_teardown(writer);
}

void mpack_writer_init_stdfile(mpack_writer_t* writer, FILE* file, bool close_when_done) {
    mpack_assert(file != NULL, "file is NULL");

    size_t capacity = MPACK_BUFFER_SIZE;
    char* buffer = (char*)MPACK_MALLOC(capacity);
    if (buffer == NULL) {
        mpack_writer_init_error(writer, mpack_error_memory);
        if (close_when_done) {
            fclose(file);
        }
        return;
    }

    mpack_writer_init(writer, buffer, capacity);
    mpack_writer_set_context(writer, file);
    mpack_writer_set_flush(writer, mpack_file_writer_flush);
    mpack_writer_set_teardown(writer, close_when_done ?
            mpack_file_writer_teardown_close :
            mpack_file_writer_teardown);
}

void mpack_writer_init_filename(mpack_writer_t* writer, const char* filename) {
    mpack_assert(filename != NULL, "filename is NULL");

    FILE* file = fopen(filename, "wb");
    if (file == NULL) {
        mpack_writer_init_error(writer, mpack_error_io);
        return;
    }

    mpack_writer_init_stdfile(writer, file, true);
}
#endif

void mpack_writer_flag_error(mpack_writer_t* writer, mpack_error_t error) {
    mpack_log("writer %p setting error %i: %s\n", (void*)writer, (int)error, mpack_error_to_string(error));

    if (writer->error == mpack_ok) {
        writer->error = error;
        if (writer->error_fn)
            writer->error_fn(writer, writer->error);
    }
}

MPACK_STATIC_INLINE void mpack_writer_flush_unchecked(mpack_writer_t* writer) {
    // This is a bit ugly; we reset used before calling flush so that
    // a flush function can distinguish between flushing the buffer
    // versus flushing external data. see mpack_growable_writer_flush()
    size_t used = mpack_writer_buffer_used(writer);
    writer->position = writer->buffer;
    writer->flush(writer, writer->buffer, used);
}

void mpack_writer_flush_message(mpack_writer_t* writer) {
    if (writer->error != mpack_ok)
        return;

    #if MPACK_WRITE_TRACKING
    // You cannot flush while there are elements open.
    mpack_writer_flag_if_error(writer, mpack_track_check_empty(&writer->track));
    if (writer->error != mpack_ok)
        return;
    #endif

    #if MPACK_BUILDER
    if (writer->builder.current_build != NULL) {
        mpack_break("cannot call mpack_writer_flush_message() while there are elements open!");
        mpack_writer_flag_error(writer, mpack_error_bug);
        return;
    }
    #endif

    if (writer->flush == NULL) {
        mpack_break("cannot call mpack_writer_flush_message() without a flush function!");
        mpack_writer_flag_error(writer, mpack_error_bug);
        return;
    }

    if (mpack_writer_buffer_used(writer) > 0)
        mpack_writer_flush_unchecked(writer);
}

// Ensures there are at least count bytes free in the buffer. This
// will flag an error if the flush function fails to make enough
// room in the buffer.
MPACK_NOINLINE static bool mpack_writer_ensure(mpack_writer_t* writer, size_t count) {
    mpack_assert(count != 0, "cannot ensure zero bytes!");
    mpack_assert(count <= MPACK_WRITER_MINIMUM_BUFFER_SIZE,
            "cannot ensure %i bytes, this is more than the minimum buffer size %i!",
            (int)count, (int)MPACK_WRITER_MINIMUM_BUFFER_SIZE);
    mpack_assert(count > mpack_writer_buffer_left(writer),
            "request to ensure %i bytes but there are already %i left in the buffer!",
            (int)count, (int)mpack_writer_buffer_left(writer));

    mpack_log("ensuring %i bytes, %i left\n", (int)count, (int)mpack_writer_buffer_left(writer));

    if (mpack_writer_error(writer) != mpack_ok)
        return false;

    #if MPACK_BUILDER
    // if we have a build in progress, we just ask the builder for a page.
    // either it will have space for a tag, or it will flag a memory error.
    if (writer->builder.current_build != NULL) {
        mpack_builder_flush(writer);
        return mpack_writer_error(writer) == mpack_ok;
    }
    #endif

    if (writer->flush == NULL) {
        mpack_writer_flag_error(writer, mpack_error_too_big);
        return false;
    }

    mpack_writer_flush_unchecked(writer);
    if (mpack_writer_error(writer) != mpack_ok)
        return false;

    if (mpack_writer_buffer_left(writer) >= count)
        return true;

    mpack_writer_flag_error(writer, mpack_error_io);
    return false;
}

// Writes encoded bytes to the buffer when we already know the data
// does not fit in the buffer (i.e. it straddles the edge of the
// buffer.) If there is a flush function, it is guaranteed to be
// called; otherwise mpack_error_too_big is raised.
MPACK_NOINLINE static void mpack_write_native_straddle(mpack_writer_t* writer, const char* p, size_t count) {
    mpack_assert(count == 0 || p != NULL, "data pointer for %i bytes is NULL", (int)count);

    if (mpack_writer_error(writer) != mpack_ok)
        return;
    mpack_log("big write for %i bytes from %p, %i space left in buffer\n",
            (int)count, p, (int)mpack_writer_buffer_left(writer));
    mpack_assert(count > mpack_writer_buffer_left(writer),
            "big write requested for %i bytes, but there is %i available "
            "space in buffer. should have called mpack_write_native() instead",
            (int)count, (int)(mpack_writer_buffer_left(writer)));

    #if MPACK_BUILDER
    // if we have a build in progress, we can't flush. we need to copy all
    // bytes into as many build buffer pages as it takes.
    if (writer->builder.current_build != NULL) {
        while (true) {
            size_t step = (size_t)(writer->end - writer->position);
            if (step > count)
                step = count;
            mpack_memcpy(writer->position, p, step);
            writer->position += step;
            p += step;
            count -= step;

            if (count == 0)
                return;

            mpack_builder_flush(writer);
            if (mpack_writer_error(writer) != mpack_ok)
                return;
            mpack_assert(writer->position != writer->end);
        }
    }
    #endif

    // we'll need a flush function
    if (!writer->flush) {
        mpack_writer_flag_error(writer, mpack_error_too_big);
        return;
    }

    // flush the buffer
    mpack_writer_flush_unchecked(writer);
    if (mpack_writer_error(writer) != mpack_ok)
        return;

    // note that an intrusive flush function (such as mpack_growable_writer_flush())
    // may have changed size and/or reset used to a non-zero value. we treat both as
    // though they may have changed, and there may still be data in the buffer.

    // flush the extra data directly if it doesn't fit in the buffer
    if (count > mpack_writer_buffer_left(writer)) {
        writer->flush(writer, p, count);
        if (mpack_writer_error(writer) != mpack_ok)
            return;
    } else {
        mpack_memcpy(writer->position, p, count);
        writer->position += count;
    }
}

// Writes encoded bytes to the buffer, flushing if necessary.
MPACK_STATIC_INLINE void mpack_write_native(mpack_writer_t* writer, const char* p, size_t count) {
    mpack_assert(count == 0 || p != NULL, "data pointer for %i bytes is NULL", (int)count);

    if (mpack_writer_buffer_left(writer) < count) {
        mpack_write_native_straddle(writer, p, count);
    } else {
        mpack_memcpy(writer->position, p, count);
        writer->position += count;
    }
}

mpack_error_t mpack_writer_destroy(mpack_writer_t* writer) {

    // clean up tracking, asserting if we're not already in an error state
    #if MPACK_WRITE_TRACKING
    mpack_track_destroy(&writer->track, writer->error != mpack_ok);
    #endif

    #if MPACK_BUILDER
    mpack_builder_t* builder = &writer->builder;
    if (builder->current_build != NULL) {
        // A builder is open!

        // Flag an error, if there's not already an error. You can only skip
        // closing any open compound types if a write error occurred. If there
        // wasn't already an error, it's a bug, which will assert in debug.
        if (mpack_writer_error(writer) == mpack_ok) {
            mpack_break("writer cannot be destroyed with an incomplete builder unless "
                    "an error was flagged!");
            mpack_writer_flag_error(writer, mpack_error_bug);
        }

        // Free any remaining builder pages
        mpack_builder_page_t* page = builder->pages;
        #if MPACK_BUILDER_INTERNAL_STORAGE
        mpack_assert(page == (mpack_builder_page_t*)builder->internal);
        page = page->next;
        #endif
        while (page != NULL) {
            mpack_builder_page_t* next = page->next;
            MPACK_FREE(page);
            page = next;
        }

        // Restore the stashed pointers. The teardown function may need to free
        // them (e.g. mpack_growable_writer_teardown().)
        writer->buffer = builder->stash_buffer;
        writer->position = builder->stash_position;
        writer->end = builder->stash_end;

        // Note: It's not necessary to clean up the current_build or other
        // pointers at this point because we're guaranteed to be in an error
        // state already so a user error callback can't longjmp out. This
        // destroy function will complete no matter what so it doesn't matter
        // what junk is left in the writer.
    }
    #endif

    // flush any outstanding data
    if (mpack_writer_error(writer) == mpack_ok && mpack_writer_buffer_used(writer) != 0 && writer->flush != NULL) {
        writer->flush(writer, writer->buffer, mpack_writer_buffer_used(writer));
        writer->flush = NULL;
    }

    if (writer->teardown) {
        writer->teardown(writer);
        writer->teardown = NULL;
    }

    return writer->error;
}

void mpack_write_tag(mpack_writer_t* writer, mpack_tag_t value) {
    switch (value.type) {
        case mpack_type_missing:
            mpack_break("cannot write a missing value!");
            mpack_writer_flag_error(writer, mpack_error_bug);
            return;

        case mpack_type_nil:    mpack_write_nil   (writer);            return;
        case mpack_type_bool:   mpack_write_bool  (writer, value.v.b); return;
        case mpack_type_int:    mpack_write_int   (writer, value.v.i); return;
        case mpack_type_uint:   mpack_write_uint  (writer, value.v.u); return;

        case mpack_type_float:
            #if MPACK_FLOAT
            mpack_write_float
            #else
            mpack_write_raw_float
            #endif
                (writer, value.v.f);
            return;
        case mpack_type_double:
            #if MPACK_DOUBLE
            mpack_write_double
            #else
            mpack_write_raw_double
            #endif
                (writer, value.v.d);
            return;

        case mpack_type_str: mpack_start_str(writer, value.v.l); return;
        case mpack_type_bin: mpack_start_bin(writer, value.v.l); return;

        #if MPACK_EXTENSIONS
        case mpack_type_ext:
            mpack_start_ext(writer, mpack_tag_ext_exttype(&value), mpack_tag_ext_length(&value));
            return;
        #endif

        case mpack_type_array: mpack_start_array(writer, value.v.n); return;
        case mpack_type_map:   mpack_start_map(writer, value.v.n);   return;
    }

    mpack_break("unrecognized type %i", (int)value.type);
    mpack_writer_flag_error(writer, mpack_error_bug);
}

MPACK_STATIC_INLINE void mpack_write_byte_element(mpack_writer_t* writer, char value) {
    mpack_writer_track_element(writer);
    if (MPACK_LIKELY(mpack_writer_buffer_left(writer) >= 1) || mpack_writer_ensure(writer, 1))
        *(writer->position++) = value;
}

void mpack_write_nil(mpack_writer_t* writer) {
    mpack_write_byte_element(writer, (char)0xc0);
}

void mpack_write_bool(mpack_writer_t* writer, bool value) {
    mpack_write_byte_element(writer, (char)(0xc2 | (value ? 1 : 0)));
}

void mpack_write_true(mpack_writer_t* writer) {
    mpack_write_byte_element(writer, (char)0xc3);
}

void mpack_write_false(mpack_writer_t* writer) {
    mpack_write_byte_element(writer, (char)0xc2);
}

void mpack_write_object_bytes(mpack_writer_t* writer, const char* data, size_t bytes) {
    mpack_writer_track_element(writer);
    mpack_write_native(writer, data, bytes);
}

/*
 * Encode functions
 */

MPACK_STATIC_INLINE void mpack_encode_fixuint(char* p, uint8_t value) {
    mpack_assert(value <= 127);
    mpack_store_u8(p, value);
}

MPACK_STATIC_INLINE void mpack_encode_u8(char* p, uint8_t value) {
    mpack_assert(value > 127);
    mpack_store_u8(p, 0xcc);
    mpack_store_u8(p + 1, value);
}

MPACK_STATIC_INLINE void mpack_encode_u16(char* p, uint16_t value) {
    mpack_assert(value > MPACK_UINT8_MAX);
    mpack_store_u8(p, 0xcd);
    mpack_store_u16(p + 1, value);
}

MPACK_STATIC_INLINE void mpack_encode_u32(char* p, uint32_t value) {
    mpack_assert(value > MPACK_UINT16_MAX);
    mpack_store_u8(p, 0xce);
    mpack_store_u32(p + 1, value);
}

MPACK_STATIC_INLINE void mpack_encode_u64(char* p, uint64_t value) {
    mpack_assert(value > MPACK_UINT32_MAX);
    mpack_store_u8(p, 0xcf);
    mpack_store_u64(p + 1, value);
}

MPACK_STATIC_INLINE void mpack_encode_fixint(char* p, int8_t value) {
    // this can encode positive or negative fixints
    mpack_assert(value >= -32);
    mpack_store_i8(p, value);
}

MPACK_STATIC_INLINE void mpack_encode_i8(char* p, int8_t value) {
    mpack_assert(value < -32);
    mpack_store_u8(p, 0xd0);
    mpack_store_i8(p + 1, value);
}

MPACK_STATIC_INLINE void mpack_encode_i16(char* p, int16_t value) {
    mpack_assert(value < MPACK_INT8_MIN);
    mpack_store_u8(p, 0xd1);
    mpack_store_i16(p + 1, value);
}

MPACK_STATIC_INLINE void mpack_encode_i32(char* p, int32_t value) {
    mpack_assert(value < MPACK_INT16_MIN);
    mpack_store_u8(p, 0xd2);
    mpack_store_i32(p + 1, value);
}

MPACK_STATIC_INLINE void mpack_encode_i64(char* p, int64_t value) {
    mpack_assert(value < MPACK_INT32_MIN);
    mpack_store_u8(p, 0xd3);
    mpack_store_i64(p + 1, value);
}

#if MPACK_FLOAT
MPACK_STATIC_INLINE void mpack_encode_float(char* p, float value) {
    mpack_store_u8(p, 0xca);
    mpack_store_float(p + 1, value);
}
#else
MPACK_STATIC_INLINE void mpack_encode_raw_float(char* p, uint32_t value) {
    mpack_store_u8(p, 0xca);
    mpack_store_u32(p + 1, value);
}
#endif

#if MPACK_DOUBLE
MPACK_STATIC_INLINE void mpack_encode_double(char* p, double value) {
    mpack_store_u8(p, 0xcb);
    mpack_store_double(p + 1, value);
}
#else
MPACK_STATIC_INLINE void mpack_encode_raw_double(char* p, uint64_t value) {
    mpack_store_u8(p, 0xcb);
    mpack_store_u64(p + 1, value);
}
#endif

MPACK_STATIC_INLINE void mpack_encode_fixarray(char* p, uint8_t count) {
    mpack_assert(count <= 15);
    mpack_store_u8(p, (uint8_t)(0x90 | count));
}

MPACK_STATIC_INLINE void mpack_encode_array16(char* p, uint16_t count) {
    mpack_assert(count > 15);
    mpack_store_u8(p, 0xdc);
    mpack_store_u16(p + 1, count);
}

MPACK_STATIC_INLINE void mpack_encode_array32(char* p, uint32_t count) {
    mpack_assert(count > MPACK_UINT16_MAX);
    mpack_store_u8(p, 0xdd);
    mpack_store_u32(p + 1, count);
}

MPACK_STATIC_INLINE void mpack_encode_fixmap(char* p, uint8_t count) {
    mpack_assert(count <= 15);
    mpack_store_u8(p, (uint8_t)(0x80 | count));
}

MPACK_STATIC_INLINE void mpack_encode_map16(char* p, uint16_t count) {
    mpack_assert(count > 15);
    mpack_store_u8(p, 0xde);
    mpack_store_u16(p + 1, count);
}

MPACK_STATIC_INLINE void mpack_encode_map32(char* p, uint32_t count) {
    mpack_assert(count > MPACK_UINT16_MAX);
    mpack_store_u8(p, 0xdf);
    mpack_store_u32(p + 1, count);
}

MPACK_STATIC_INLINE void mpack_encode_fixstr(char* p, uint8_t count) {
    mpack_assert(count <= 31);
    mpack_store_u8(p, (uint8_t)(0xa0 | count));
}

MPACK_STATIC_INLINE void mpack_encode_str8(char* p, uint8_t count) {
    mpack_assert(count > 31);
    mpack_store_u8(p, 0xd9);
    mpack_store_u8(p + 1, count);
}

MPACK_STATIC_INLINE void mpack_encode_str16(char* p, uint16_t count) {
    // we might be encoding a raw in compatibility mode, so we
    // allow count to be in the range [32, MPACK_UINT8_MAX].
    mpack_assert(count > 31);
    mpack_store_u8(p, 0xda);
    mpack_store_u16(p + 1, count);
}

MPACK_STATIC_INLINE void mpack_encode_str32(char* p, uint32_t count) {
    mpack_assert(count > MPACK_UINT16_MAX);
    mpack_store_u8(p, 0xdb);
    mpack_store_u32(p + 1, count);
}

MPACK_STATIC_INLINE void mpack_encode_bin8(char* p, uint8_t count) {
    mpack_store_u8(p, 0xc4);
    mpack_store_u8(p + 1, count);
}

MPACK_STATIC_INLINE void mpack_encode_bin16(char* p, uint16_t count) {
    mpack_assert(count > MPACK_UINT8_MAX);
    mpack_store_u8(p, 0xc5);
    mpack_store_u16(p + 1, count);
}

MPACK_STATIC_INLINE void mpack_encode_bin32(char* p, uint32_t count) {
    mpack_assert(count > MPACK_UINT16_MAX);
    mpack_store_u8(p, 0xc6);
    mpack_store_u32(p + 1, count);
}

#if MPACK_EXTENSIONS
MPACK_STATIC_INLINE void mpack_encode_fixext1(char* p, int8_t exttype) {
    mpack_store_u8(p, 0xd4);
    mpack_store_i8(p + 1, exttype);
}

MPACK_STATIC_INLINE void mpack_encode_fixext2(char* p, int8_t exttype) {
    mpack_store_u8(p, 0xd5);
    mpack_store_i8(p + 1, exttype);
}

MPACK_STATIC_INLINE void mpack_encode_fixext4(char* p, int8_t exttype) {
    mpack_store_u8(p, 0xd6);
    mpack_store_i8(p + 1, exttype);
}

MPACK_STATIC_INLINE void mpack_encode_fixext8(char* p, int8_t exttype) {
    mpack_store_u8(p, 0xd7);
    mpack_store_i8(p + 1, exttype);
}

MPACK_STATIC_INLINE void mpack_encode_fixext16(char* p, int8_t exttype) {
    mpack_store_u8(p, 0xd8);
    mpack_store_i8(p + 1, exttype);
}

MPACK_STATIC_INLINE void mpack_encode_ext8(char* p, int8_t exttype, uint8_t count) {
    mpack_assert(count != 1 && count != 2 && count != 4 && count != 8 && count != 16);
    mpack_store_u8(p, 0xc7);
    mpack_store_u8(p + 1, count);
    mpack_store_i8(p + 2, exttype);
}

MPACK_STATIC_INLINE void mpack_encode_ext16(char* p, int8_t exttype, uint16_t count) {
    mpack_assert(count > MPACK_UINT8_MAX);
    mpack_store_u8(p, 0xc8);
    mpack_store_u16(p + 1, count);
    mpack_store_i8(p + 3, exttype);
}

MPACK_STATIC_INLINE void mpack_encode_ext32(char* p, int8_t exttype, uint32_t count) {
    mpack_assert(count > MPACK_UINT16_MAX);
    mpack_store_u8(p, 0xc9);
    mpack_store_u32(p + 1, count);
    mpack_store_i8(p + 5, exttype);
}

MPACK_STATIC_INLINE void mpack_encode_timestamp_4(char* p, uint32_t seconds) {
    mpack_encode_fixext4(p, MPACK_EXTTYPE_TIMESTAMP);
    mpack_store_u32(p + MPACK_TAG_SIZE_FIXEXT4, seconds);
}

MPACK_STATIC_INLINE void mpack_encode_timestamp_8(char* p, int64_t seconds, uint32_t nanoseconds) {
    mpack_assert(nanoseconds <= MPACK_TIMESTAMP_NANOSECONDS_MAX);
    mpack_encode_fixext8(p, MPACK_EXTTYPE_TIMESTAMP);
    uint64_t encoded = ((uint64_t)nanoseconds << 34) | (uint64_t)seconds;
    mpack_store_u64(p + MPACK_TAG_SIZE_FIXEXT8, encoded);
}

MPACK_STATIC_INLINE void mpack_encode_timestamp_12(char* p, int64_t seconds, uint32_t nanoseconds) {
    mpack_assert(nanoseconds <= MPACK_TIMESTAMP_NANOSECONDS_MAX);
    mpack_encode_ext8(p, MPACK_EXTTYPE_TIMESTAMP, 12);
    mpack_store_u32(p + MPACK_TAG_SIZE_EXT8, nanoseconds);
    mpack_store_i64(p + MPACK_TAG_SIZE_EXT8 + 4, seconds);
}
#endif



/*
 * Write functions
 */

// This is a macro wrapper to the encode functions to encode
// directly into the buffer. If mpack_writer_ensure() fails
// it will flag an error so we don't have to do anything.
#define MPACK_WRITE_ENCODED(encode_fn, size, ...) do {                                                 \
    if (MPACK_LIKELY(mpack_writer_buffer_left(writer) >= size) || mpack_writer_ensure(writer, size)) { \
        MPACK_EXPAND(encode_fn(writer->position, __VA_ARGS__));                                        \
        writer->position += size;                                                                      \
    }                                                                                                  \
} while (0)

void mpack_write_u8(mpack_writer_t* writer, uint8_t value) {
    #if MPACK_OPTIMIZE_FOR_SIZE
    mpack_write_u64(writer, value);
    #else
    mpack_writer_track_element(writer);
    if (value <= 127) {
        MPACK_WRITE_ENCODED(mpack_encode_fixuint, MPACK_TAG_SIZE_FIXUINT, value);
    } else {
        MPACK_WRITE_ENCODED(mpack_encode_u8, MPACK_TAG_SIZE_U8, value);
    }
    #endif
}

void mpack_write_u16(mpack_writer_t* writer, uint16_t value) {
    #if MPACK_OPTIMIZE_FOR_SIZE
    mpack_write_u64(writer, value);
    #else
    mpack_writer_track_element(writer);
    if (value <= 127) {
        MPACK_WRITE_ENCODED(mpack_encode_fixuint, MPACK_TAG_SIZE_FIXUINT, (uint8_t)value);
    } else if (value <= MPACK_UINT8_MAX) {
        MPACK_WRITE_ENCODED(mpack_encode_u8, MPACK_TAG_SIZE_U8, (uint8_t)value);
    } else {
        MPACK_WRITE_ENCODED(mpack_encode_u16, MPACK_TAG_SIZE_U16, value);
    }
    #endif
}

void mpack_write_u32(mpack_writer_t* writer, uint32_t value) {
    #if MPACK_OPTIMIZE_FOR_SIZE
    mpack_write_u64(writer, value);
    #else
    mpack_writer_track_element(writer);
    if (value <= 127) {
        MPACK_WRITE_ENCODED(mpack_encode_fixuint, MPACK_TAG_SIZE_FIXUINT, (uint8_t)value);
    } else if (value <= MPACK_UINT8_MAX) {
        MPACK_WRITE_ENCODED(mpack_encode_u8, MPACK_TAG_SIZE_U8, (uint8_t)value);
    } else if (value <= MPACK_UINT16_MAX) {
        MPACK_WRITE_ENCODED(mpack_encode_u16, MPACK_TAG_SIZE_U16, (uint16_t)value);
    } else {
        MPACK_WRITE_ENCODED(mpack_encode_u32, MPACK_TAG_SIZE_U32, value);
    }
    #endif
}

void mpack_write_u64(mpack_writer_t* writer, uint64_t value) {
    mpack_writer_track_element(writer);

    if (value <= 127) {
        MPACK_WRITE_ENCODED(mpack_encode_fixuint, MPACK_TAG_SIZE_FIXUINT, (uint8_t)value);
    } else if (value <= MPACK_UINT8_MAX) {
        MPACK_WRITE_ENCODED(mpack_encode_u8, MPACK_TAG_SIZE_U8, (uint8_t)value);
    } else if (value <= MPACK_UINT16_MAX) {
        MPACK_WRITE_ENCODED(mpack_encode_u16, MPACK_TAG_SIZE_U16, (uint16_t)value);
    } else if (value <= MPACK_UINT32_MAX) {
        MPACK_WRITE_ENCODED(mpack_encode_u32, MPACK_TAG_SIZE_U32, (uint32_t)value);
    } else {
        MPACK_WRITE_ENCODED(mpack_encode_u64, MPACK_TAG_SIZE_U64, value);
    }
}

void mpack_write_i8(mpack_writer_t* writer, int8_t value) {
    #if MPACK_OPTIMIZE_FOR_SIZE
    mpack_write_i64(writer, value);
    #else
    mpack_writer_track_element(writer);
    if (value >= -32) {
        // we encode positive and negative fixints together
        MPACK_WRITE_ENCODED(mpack_encode_fixint, MPACK_TAG_SIZE_FIXINT, (int8_t)value);
    } else {
        MPACK_WRITE_ENCODED(mpack_encode_i8, MPACK_TAG_SIZE_I8, (int8_t)value);
    }
    #endif
}

void mpack_write_i16(mpack_writer_t* writer, int16_t value) {
    #if MPACK_OPTIMIZE_FOR_SIZE
    mpack_write_i64(writer, value);
    #else
    mpack_writer_track_element(writer);
    if (value >= -32) {
        if (value <= 127) {
            // we encode positive and negative fixints together
            MPACK_WRITE_ENCODED(mpack_encode_fixint, MPACK_TAG_SIZE_FIXINT, (int8_t)value);
        } else if (value <= MPACK_UINT8_MAX) {
            MPACK_WRITE_ENCODED(mpack_encode_u8, MPACK_TAG_SIZE_U8, (uint8_t)value);
        } else {
            MPACK_WRITE_ENCODED(mpack_encode_u16, MPACK_TAG_SIZE_U16, (uint16_t)value);
        }
    } else if (value >= MPACK_INT8_MIN) {
        MPACK_WRITE_ENCODED(mpack_encode_i8, MPACK_TAG_SIZE_I8, (int8_t)value);
    } else {
        MPACK_WRITE_ENCODED(mpack_encode_i16, MPACK_TAG_SIZE_I16, (int16_t)value);
    }
    #endif
}

void mpack_write_i32(mpack_writer_t* writer, int32_t value) {
    #if MPACK_OPTIMIZE_FOR_SIZE
    mpack_write_i64(writer, value);
    #else
    mpack_writer_track_element(writer);
    if (value >= -32) {
        if (value <= 127) {
            // we encode positive and negative fixints together
            MPACK_WRITE_ENCODED(mpack_encode_fixint, MPACK_TAG_SIZE_FIXINT, (int8_t)value);
        } else if (value <= MPACK_UINT8_MAX) {
            MPACK_WRITE_ENCODED(mpack_encode_u8, MPACK_TAG_SIZE_U8, (uint8_t)value);
        } else if (value <= MPACK_UINT16_MAX) {
            MPACK_WRITE_ENCODED(mpack_encode_u16, MPACK_TAG_SIZE_U16, (uint16_t)value);
        } else {
            MPACK_WRITE_ENCODED(mpack_encode_u32, MPACK_TAG_SIZE_U32, (uint32_t)value);
        }
    } else if (value >= MPACK_INT8_MIN) {
        MPACK_WRITE_ENCODED(mpack_encode_i8, MPACK_TAG_SIZE_I8, (int8_t)value);
    } else if (value >= MPACK_INT16_MIN) {
        MPACK_WRITE_ENCODED(mpack_encode_i16, MPACK_TAG_SIZE_I16, (int16_t)value);
    } else {
        MPACK_WRITE_ENCODED(mpack_encode_i32, MPACK_TAG_SIZE_I32, value);
    }
    #endif
}

void mpack_write_i64(mpack_writer_t* writer, int64_t value) {
    #if MPACK_OPTIMIZE_FOR_SIZE
    if (value > 127) {
        // for non-fix positive ints we call the u64 writer to save space
        mpack_write_u64(writer, (uint64_t)value);
        return;
    }
    #endif

    mpack_writer_track_element(writer);
    if (value >= -32) {
        #if MPACK_OPTIMIZE_FOR_SIZE
        MPACK_WRITE_ENCODED(mpack_encode_fixint, MPACK_TAG_SIZE_FIXINT, (int8_t)value);
        #else
        if (value <= 127) {
            MPACK_WRITE_ENCODED(mpack_encode_fixint, MPACK_TAG_SIZE_FIXINT, (int8_t)value);
        } else if (value <= MPACK_UINT8_MAX) {
            MPACK_WRITE_ENCODED(mpack_encode_u8, MPACK_TAG_SIZE_U8, (uint8_t)value);
        } else if (value <= MPACK_UINT16_MAX) {
            MPACK_WRITE_ENCODED(mpack_encode_u16, MPACK_TAG_SIZE_U16, (uint16_t)value);
        } else if (value <= MPACK_UINT32_MAX) {
            MPACK_WRITE_ENCODED(mpack_encode_u32, MPACK_TAG_SIZE_U32, (uint32_t)value);
        } else {
            MPACK_WRITE_ENCODED(mpack_encode_u64, MPACK_TAG_SIZE_U64, (uint64_t)value);
        }
        #endif
    } else if (value >= MPACK_INT8_MIN) {
        MPACK_WRITE_ENCODED(mpack_encode_i8, MPACK_TAG_SIZE_I8, (int8_t)value);
    } else if (value >= MPACK_INT16_MIN) {
        MPACK_WRITE_ENCODED(mpack_encode_i16, MPACK_TAG_SIZE_I16, (int16_t)value);
    } else if (value >= MPACK_INT32_MIN) {
        MPACK_WRITE_ENCODED(mpack_encode_i32, MPACK_TAG_SIZE_I32, (int32_t)value);
    } else {
        MPACK_WRITE_ENCODED(mpack_encode_i64, MPACK_TAG_SIZE_I64, value);
    }
}

#if MPACK_FLOAT
void mpack_write_float(mpack_writer_t* writer, float value) {
    mpack_writer_track_element(writer);
    MPACK_WRITE_ENCODED(mpack_encode_float, MPACK_TAG_SIZE_FLOAT, value);
}
#else
void mpack_write_raw_float(mpack_writer_t* writer, uint32_t value) {
    mpack_writer_track_element(writer);
    MPACK_WRITE_ENCODED(mpack_encode_raw_float, MPACK_TAG_SIZE_FLOAT, value);
}
#endif

#if MPACK_DOUBLE
void mpack_write_double(mpack_writer_t* writer, double value) {
    mpack_writer_track_element(writer);
    MPACK_WRITE_ENCODED(mpack_encode_double, MPACK_TAG_SIZE_DOUBLE, value);
}
#else
void mpack_write_raw_double(mpack_writer_t* writer, uint64_t value) {
    mpack_writer_track_element(writer);
    MPACK_WRITE_ENCODED(mpack_encode_raw_double, MPACK_TAG_SIZE_DOUBLE, value);
}
#endif

#if MPACK_EXTENSIONS
void mpack_write_timestamp(mpack_writer_t* writer, int64_t seconds, uint32_t nanoseconds) {
    #if MPACK_COMPATIBILITY
    if (writer->version <= mpack_version_v4) {
        mpack_break("Timestamps require spec version v5 or later. This writer is in v%i mode.", (int)writer->version);
        mpack_writer_flag_error(writer, mpack_error_bug);
        return;
    }
    #endif

    if (nanoseconds > MPACK_TIMESTAMP_NANOSECONDS_MAX) {
        mpack_break("timestamp nanoseconds out of bounds: %" PRIu32 , nanoseconds);
        mpack_writer_flag_error(writer, mpack_error_bug);
        return;
    }

    mpack_writer_track_element(writer);

    if (seconds < 0 || seconds >= (MPACK_INT64_C(1) << 34)) {
        MPACK_WRITE_ENCODED(mpack_encode_timestamp_12, MPACK_EXT_SIZE_TIMESTAMP12, seconds, nanoseconds);
    } else if (seconds > MPACK_UINT32_MAX || nanoseconds > 0) {
        MPACK_WRITE_ENCODED(mpack_encode_timestamp_8, MPACK_EXT_SIZE_TIMESTAMP8, seconds, nanoseconds);
    } else {
        MPACK_WRITE_ENCODED(mpack_encode_timestamp_4, MPACK_EXT_SIZE_TIMESTAMP4, (uint32_t)seconds);
    }
}
#endif

static void mpack_write_array_notrack(mpack_writer_t* writer, uint32_t count) {
    if (count <= 15) {
        MPACK_WRITE_ENCODED(mpack_encode_fixarray, MPACK_TAG_SIZE_FIXARRAY, (uint8_t)count);
    } else if (count <= MPACK_UINT16_MAX) {
        MPACK_WRITE_ENCODED(mpack_encode_array16, MPACK_TAG_SIZE_ARRAY16, (uint16_t)count);
    } else {
        MPACK_WRITE_ENCODED(mpack_encode_array32, MPACK_TAG_SIZE_ARRAY32, (uint32_t)count);
    }
}

static void mpack_write_map_notrack(mpack_writer_t* writer, uint32_t count) {
    if (count <= 15) {
        MPACK_WRITE_ENCODED(mpack_encode_fixmap, MPACK_TAG_SIZE_FIXMAP, (uint8_t)count);
    } else if (count <= MPACK_UINT16_MAX) {
        MPACK_WRITE_ENCODED(mpack_encode_map16, MPACK_TAG_SIZE_MAP16, (uint16_t)count);
    } else {
        MPACK_WRITE_ENCODED(mpack_encode_map32, MPACK_TAG_SIZE_MAP32, (uint32_t)count);
    }
}

void mpack_start_array(mpack_writer_t* writer, uint32_t count) {
    mpack_writer_track_element(writer);
    mpack_write_array_notrack(writer, count);
    mpack_writer_track_push(writer, mpack_type_array, count);
    mpack_builder_compound_push(writer);
}

void mpack_start_map(mpack_writer_t* writer, uint32_t count) {
    mpack_writer_track_element(writer);
    mpack_write_map_notrack(writer, count);
    mpack_writer_track_push(writer, mpack_type_map, count);
    mpack_builder_compound_push(writer);
}

static void mpack_start_str_notrack(mpack_writer_t* writer, uint32_t count) {
    if (count <= 31) {
        MPACK_WRITE_ENCODED(mpack_encode_fixstr, MPACK_TAG_SIZE_FIXSTR, (uint8_t)count);

    // str8 is only supported in v5 or later.
    } else if (count <= MPACK_UINT8_MAX
            #if MPACK_COMPATIBILITY
            && writer->version >= mpack_version_v5
            #endif
            ) {
        MPACK_WRITE_ENCODED(mpack_encode_str8, MPACK_TAG_SIZE_STR8, (uint8_t)count);

    } else if (count <= MPACK_UINT16_MAX) {
        MPACK_WRITE_ENCODED(mpack_encode_str16, MPACK_TAG_SIZE_STR16, (uint16_t)count);
    } else {
        MPACK_WRITE_ENCODED(mpack_encode_str32, MPACK_TAG_SIZE_STR32, (uint32_t)count);
    }
}

static void mpack_start_bin_notrack(mpack_writer_t* writer, uint32_t count) {
    #if MPACK_COMPATIBILITY
    // In the v4 spec, there was only the raw type for any kind of
    // variable-length data. In v4 mode, we support the bin functions,
    // but we produce an old-style raw.
    if (writer->version <= mpack_version_v4) {
        mpack_start_str_notrack(writer, count);
        return;
    }
    #endif

    if (count <= MPACK_UINT8_MAX) {
        MPACK_WRITE_ENCODED(mpack_encode_bin8, MPACK_TAG_SIZE_BIN8, (uint8_t)count);
    } else if (count <= MPACK_UINT16_MAX) {
        MPACK_WRITE_ENCODED(mpack_encode_bin16, MPACK_TAG_SIZE_BIN16, (uint16_t)count);
    } else {
        MPACK_WRITE_ENCODED(mpack_encode_bin32, MPACK_TAG_SIZE_BIN32, (uint32_t)count);
    }
}

void mpack_start_str(mpack_writer_t* writer, uint32_t count) {
    mpack_writer_track_element(writer);
    mpack_start_str_notrack(writer, count);
    mpack_writer_track_push(writer, mpack_type_str, count);
}

void mpack_start_bin(mpack_writer_t* writer, uint32_t count) {
    mpack_writer_track_element(writer);
    mpack_start_bin_notrack(writer, count);
    mpack_writer_track_push(writer, mpack_type_bin, count);
}

#if MPACK_EXTENSIONS
void mpack_start_ext(mpack_writer_t* writer, int8_t exttype, uint32_t count) {
    #if MPACK_COMPATIBILITY
    if (writer->version <= mpack_version_v4) {
        mpack_break("Ext types require spec version v5 or later. This writer is in v%i mode.", (int)writer->version);
        mpack_writer_flag_error(writer, mpack_error_bug);
        return;
    }
    #endif

    mpack_writer_track_element(writer);

    if (count == 1) {
        MPACK_WRITE_ENCODED(mpack_encode_fixext1, MPACK_TAG_SIZE_FIXEXT1, exttype);
    } else if (count == 2) {
        MPACK_WRITE_ENCODED(mpack_encode_fixext2, MPACK_TAG_SIZE_FIXEXT2, exttype);
    } else if (count == 4) {
        MPACK_WRITE_ENCODED(mpack_encode_fixext4, MPACK_TAG_SIZE_FIXEXT4, exttype);
    } else if (count == 8) {
        MPACK_WRITE_ENCODED(mpack_encode_fixext8, MPACK_TAG_SIZE_FIXEXT8, exttype);
    } else if (count == 16) {
        MPACK_WRITE_ENCODED(mpack_encode_fixext16, MPACK_TAG_SIZE_FIXEXT16, exttype);
    } else if (count <= MPACK_UINT8_MAX) {
        MPACK_WRITE_ENCODED(mpack_encode_ext8, MPACK_TAG_SIZE_EXT8, exttype, (uint8_t)count);
    } else if (count <= MPACK_UINT16_MAX) {
        MPACK_WRITE_ENCODED(mpack_encode_ext16, MPACK_TAG_SIZE_EXT16, exttype, (uint16_t)count);
    } else {
        MPACK_WRITE_ENCODED(mpack_encode_ext32, MPACK_TAG_SIZE_EXT32, exttype, (uint32_t)count);
    }

    mpack_writer_track_push(writer, mpack_type_ext, count);
}
#endif



/*
 * Compound helpers and other functions
 */

void mpack_write_str(mpack_writer_t* writer, const char* data, uint32_t count) {
    mpack_assert(count == 0 || data != NULL, "data for string of length %i is NULL", (int)count);

    #if MPACK_OPTIMIZE_FOR_SIZE
    mpack_writer_track_element(writer);
    mpack_start_str_notrack(writer, count);
    mpack_write_native(writer, data, count);
    #else

    mpack_writer_track_element(writer);

    if (count <= 31) {
        // The minimum buffer size when using a flush function is guaranteed to
        // fit the largest possible fixstr.
        size_t size = count + MPACK_TAG_SIZE_FIXSTR;
        if (MPACK_LIKELY(mpack_writer_buffer_left(writer) >= size) || mpack_writer_ensure(writer, size)) {
            char* MPACK_RESTRICT p = writer->position;
            mpack_encode_fixstr(p, (uint8_t)count);
            mpack_memcpy(p + MPACK_TAG_SIZE_FIXSTR, data, count);
            writer->position += count + MPACK_TAG_SIZE_FIXSTR;
        }
        return;
    }

    if (count <= MPACK_UINT8_MAX
            #if MPACK_COMPATIBILITY
            && writer->version >= mpack_version_v5
            #endif
            ) {
        if (count + MPACK_TAG_SIZE_STR8 <= mpack_writer_buffer_left(writer)) {
            char* MPACK_RESTRICT p = writer->position;
            mpack_encode_str8(p, (uint8_t)count);
            mpack_memcpy(p + MPACK_TAG_SIZE_STR8, data, count);
            writer->position += count + MPACK_TAG_SIZE_STR8;
        } else {
            MPACK_WRITE_ENCODED(mpack_encode_str8, MPACK_TAG_SIZE_STR8, (uint8_t)count);
            mpack_write_native(writer, data, count);
        }
        return;
    }

    // str16 and str32 are likely to be a significant fraction of the buffer
    // size, so we don't bother with a combined space check in order to
    // minimize code size.
    if (count <= MPACK_UINT16_MAX) {
        MPACK_WRITE_ENCODED(mpack_encode_str16, MPACK_TAG_SIZE_STR16, (uint16_t)count);
        mpack_write_native(writer, data, count);
    } else {
        MPACK_WRITE_ENCODED(mpack_encode_str32, MPACK_TAG_SIZE_STR32, (uint32_t)count);
        mpack_write_native(writer, data, count);
    }

    #endif
}

void mpack_write_bin(mpack_writer_t* writer, const char* data, uint32_t count) {
    mpack_assert(count == 0 || data != NULL, "data pointer for bin of %i bytes is NULL", (int)count);
    mpack_start_bin(writer, count);
    mpack_write_bytes(writer, data, count);
    mpack_finish_bin(writer);
}

#if MPACK_EXTENSIONS
void mpack_write_ext(mpack_writer_t* writer, int8_t exttype, const char* data, uint32_t count) {
    mpack_assert(count == 0 || data != NULL, "data pointer for ext of type %i and %i bytes is NULL", exttype, (int)count);
    mpack_start_ext(writer, exttype, count);
    mpack_write_bytes(writer, data, count);
    mpack_finish_ext(writer);
}
#endif

void mpack_write_bytes(mpack_writer_t* writer, const char* data, size_t count) {
    mpack_assert(count == 0 || data != NULL, "data pointer for %i bytes is NULL", (int)count);
    mpack_writer_track_bytes(writer, count);
    mpack_write_native(writer, data, count);
}

void mpack_write_cstr(mpack_writer_t* writer, const char* cstr) {
    mpack_assert(cstr != NULL, "cstr pointer is NULL");
    size_t length = mpack_strlen(cstr);
    if (length > MPACK_UINT32_MAX)
        mpack_writer_flag_error(writer, mpack_error_invalid);
    mpack_write_str(writer, cstr, (uint32_t)length);
}

void mpack_write_cstr_or_nil(mpack_writer_t* writer, const char* cstr) {
    if (cstr)
        mpack_write_cstr(writer, cstr);
    else
        mpack_write_nil(writer);
}

void mpack_write_utf8(mpack_writer_t* writer, const char* str, uint32_t length) {
    mpack_assert(length == 0 || str != NULL, "data for string of length %i is NULL", (int)length);
    if (!mpack_utf8_check(str, length)) {
        mpack_writer_flag_error(writer, mpack_error_invalid);
        return;
    }
    mpack_write_str(writer, str, length);
}

void mpack_write_utf8_cstr(mpack_writer_t* writer, const char* cstr) {
    mpack_assert(cstr != NULL, "cstr pointer is NULL");
    size_t length = mpack_strlen(cstr);
    if (length > MPACK_UINT32_MAX) {
        mpack_writer_flag_error(writer, mpack_error_invalid);
        return;
    }
    mpack_write_utf8(writer, cstr, (uint32_t)length);
}

void mpack_write_utf8_cstr_or_nil(mpack_writer_t* writer, const char* cstr) {
    if (cstr)
        mpack_write_utf8_cstr(writer, cstr);
    else
        mpack_write_nil(writer);
}

/*
 * Builder implementation
 *
 * When a writer is in build mode, it diverts writes to an internal growable
 * buffer. All elements other than builder start tags are encoded as normal
 * into the builder buffer (even nested maps and arrays of known size, e.g.
 * `mpack_start_array()`.) But for compound elements of unknown size, an
 * mpack_build_t is written to the buffer instead.
 *
 * The mpack_build_t tracks everything needed to re-constitute the final
 * message once all sizes are known. When the last build element is completed,
 * the builder resolves the build by walking through the builds, outputting the
 * final encoded tag, and copying everything in between to the writer's true
 * buffer.
 *
 * To make things extra complicated, the builder buffer is not contiguous. It's
 * allocated in pages, where the first page may be an internal page in the
 * writer. But, each mpack_build_t must itself be contiguous and aligned
 * properly within the buffer. This means bytes can be skipped (and wasted)
 * before the builds or at the end of pages.
 *
 * To keep track of this, builds store both their element count and the number
 * of encoded bytes that follow, and pages store the number of bytes used. As
 * elements are written, each element adds to the count in the current open
 * build, and the number of bytes written adds to the current page and the byte
 * count in the last started build (whether or not it is completed.)
 */

#if MPACK_BUILDER

#ifdef MPACK_ALIGNOF
    #define MPACK_BUILD_ALIGNMENT MPACK_ALIGNOF(mpack_build_t)
#else
    // without alignof, we just align to the greater of size_t, void* and uint64_t.
    // (we do this even though we don't have uint64_t in it in case we add it later.)
    #define MPACK_BUILD_ALIGNMENT_MAX(x, y) ((x) > (y) ? (x) : (y))
    #define MPACK_BUILD_ALIGNMENT (MPACK_BUILD_ALIGNMENT_MAX(sizeof(void*), \
                MPACK_BUILD_ALIGNMENT_MAX(sizeof(size_t), sizeof(uint64_t))))
#endif

static inline void mpack_builder_check_sizes(mpack_writer_t* writer) {

    // We check internal and page sizes here so that we don't have to check
    // them again. A new page with a build in it will have a page header,
    // build, and minimum space for a tag. This will perform horribly and waste
    // tons of memory if the page size is small, so you're best off just
    // sticking with the defaults.
    //
    // These are all known at compile time, so if they are large
    // enough this function should trivially optimize to a no-op.

    #if MPACK_BUILDER_INTERNAL_STORAGE
    // make sure the internal storage is big enough to be useful
    MPACK_STATIC_ASSERT(MPACK_BUILDER_INTERNAL_STORAGE_SIZE >= (sizeof(mpack_builder_page_t) +
            sizeof(mpack_build_t) + MPACK_WRITER_MINIMUM_BUFFER_SIZE),
            "MPACK_BUILDER_INTERNAL_STORAGE_SIZE is too small to be useful!");
    if (MPACK_BUILDER_INTERNAL_STORAGE_SIZE < (sizeof(mpack_builder_page_t) +
            sizeof(mpack_build_t) + MPACK_WRITER_MINIMUM_BUFFER_SIZE))
    {
        mpack_break("MPACK_BUILDER_INTERNAL_STORAGE_SIZE is too small to be useful!");
        mpack_writer_flag_error(writer, mpack_error_bug);
    }
    #endif

    // make sure the builder page size is big enough to be useful
    MPACK_STATIC_ASSERT(MPACK_BUILDER_PAGE_SIZE >= (sizeof(mpack_builder_page_t) +
            sizeof(mpack_build_t) + MPACK_WRITER_MINIMUM_BUFFER_SIZE),
            "MPACK_BUILDER_PAGE_SIZE is too small to be useful!");
    if (MPACK_BUILDER_PAGE_SIZE < (sizeof(mpack_builder_page_t) +
            sizeof(mpack_build_t) + MPACK_WRITER_MINIMUM_BUFFER_SIZE))
    {
        mpack_break("MPACK_BUILDER_PAGE_SIZE is too small to be useful!");
        mpack_writer_flag_error(writer, mpack_error_bug);
    }
}

static inline size_t mpack_builder_page_size(mpack_writer_t* writer, mpack_builder_page_t* page) {
    #if MPACK_BUILDER_INTERNAL_STORAGE
    if ((char*)page == writer->builder.internal)
        return sizeof(writer->builder.internal);
    #else
    (void)writer;
    (void)page;
    #endif
    return MPACK_BUILDER_PAGE_SIZE;
}

static inline size_t mpack_builder_align_build(size_t bytes_used) {
    size_t offset = bytes_used;
    offset += MPACK_BUILD_ALIGNMENT - 1;
    offset -= offset % MPACK_BUILD_ALIGNMENT;
    mpack_log("aligned %zi to %zi\n", bytes_used, offset);
    return offset;
}

static inline void mpack_builder_free_page(mpack_writer_t* writer, mpack_builder_page_t* page) {
    mpack_log("freeing page %p\n", (void*)page);
    #if MPACK_BUILDER_INTERNAL_STORAGE
    if ((char*)page == writer->builder.internal)
        return;
    #else
    (void)writer;
    #endif
    MPACK_FREE(page);
}

static inline size_t mpack_builder_page_remaining(mpack_writer_t* writer, mpack_builder_page_t* page) {
    return mpack_builder_page_size(writer, page) - page->bytes_used;
}

static void mpack_builder_configure_buffer(mpack_writer_t* writer) {
    if (mpack_writer_error(writer) != mpack_ok)
        return;
    mpack_builder_t* builder = &writer->builder;

    mpack_builder_page_t* page = builder->current_page;
    mpack_assert(page != NULL, "page is null??");

    // This diverts the writer into the remainder of the current page of our
    // build buffer.
    writer->buffer = (char*)page + page->bytes_used;
    writer->position = (char*)page + page->bytes_used;
    writer->end = (char*)page + mpack_builder_page_size(writer, page);
    mpack_log("configuring buffer from %p to %p\n", (void*)writer->position, (void*)writer->end);
}

static void mpack_builder_add_page(mpack_writer_t* writer) {
    mpack_builder_t* builder = &writer->builder;
    mpack_assert(writer->error == mpack_ok);

    mpack_log("adding a page.\n");
    mpack_builder_page_t* page = (mpack_builder_page_t*)MPACK_MALLOC(MPACK_BUILDER_PAGE_SIZE);
    if (page == NULL) {
        mpack_writer_flag_error(writer, mpack_error_memory);
        return;
    }

    page->next = NULL;
    page->bytes_used = sizeof(mpack_builder_page_t);
    builder->current_page->next = page;
    builder->current_page = page;
}

// Checks how many bytes the writer wrote to the page, adding it to the page's
// bytes_used. This must be followed up with mpack_builder_configure_buffer()
// (after adding a new page, build, etc) to reset the writer's buffer pointers.
static void mpack_builder_apply_writes(mpack_writer_t* writer) {
    mpack_assert(writer->error == mpack_ok);
    mpack_builder_t* builder = &writer->builder;
    mpack_log("latest build is %p\n", (void*)builder->latest_build);

    // The difference between buffer and current is the number of bytes that
    // were written to the page.
    size_t bytes_written = (size_t)(writer->position - writer->buffer);
    mpack_log("applying write of %zi bytes to build %p\n", bytes_written, (void*)builder->latest_build);

    mpack_assert(builder->current_page != NULL);
    mpack_assert(builder->latest_build != NULL);
    builder->current_page->bytes_used += bytes_written;
    builder->latest_build->bytes += bytes_written;
    mpack_log("latest build %p now has %zi bytes\n", (void*)builder->latest_build, builder->latest_build->bytes);
}

static void mpack_builder_flush(mpack_writer_t* writer) {
    mpack_assert(writer->error == mpack_ok);
    mpack_builder_apply_writes(writer);
    mpack_builder_add_page(writer);
    mpack_builder_configure_buffer(writer);
}

MPACK_NOINLINE static void mpack_builder_begin(mpack_writer_t* writer) {
    mpack_builder_t* builder = &writer->builder;
    mpack_assert(writer->error == mpack_ok);
    mpack_assert(builder->current_build == NULL);
    mpack_assert(builder->latest_build == NULL);
    mpack_assert(builder->pages == NULL);

    // If this is the first build, we need to stash the real buffer backing our
    // writer. We'll be diverting the writer to our build buffer.
    builder->stash_buffer = writer->buffer;
    builder->stash_position = writer->position;
    builder->stash_end = writer->end;

    mpack_builder_page_t* page;

    // we've checked that both these sizes are large enough above.
    #if MPACK_BUILDER_INTERNAL_STORAGE
    page = (mpack_builder_page_t*)builder->internal;
    mpack_log("beginning builder with internal storage %p\n", (void*)page);
    #else
    page = (mpack_builder_page_t*)MPACK_MALLOC(MPACK_BUILDER_PAGE_SIZE);
    if (page == NULL) {
        mpack_writer_flag_error(writer, mpack_error_memory);
        return;
    }
    mpack_log("beginning builder with allocated page %p\n", (void*)page);
    #endif

    page->next = NULL;
    page->bytes_used = sizeof(mpack_builder_page_t);
    builder->pages = page;
    builder->current_page = page;
}

static void mpack_builder_build(mpack_writer_t* writer, mpack_type_t type) {
    mpack_builder_check_sizes(writer);
    if (mpack_writer_error(writer) != mpack_ok)
        return;

    mpack_writer_track_element(writer);
    mpack_writer_track_push_builder(writer, type);

    mpack_builder_t* builder = &writer->builder;

    if (builder->current_build == NULL) {
        mpack_builder_begin(writer);
    } else {
        mpack_builder_apply_writes(writer);
    }
    if (mpack_writer_error(writer) != mpack_ok)
        return;

    // find aligned space for a new build. if there isn't enough space in the
    // current page, we discard the remaining space in it and allocate a new
    // page.
    size_t offset = mpack_builder_align_build(builder->current_page->bytes_used);
    if (offset + sizeof(mpack_build_t) > mpack_builder_page_size(writer, builder->current_page)) {
        mpack_log("not enough space for a build. %zi bytes used of %zi in this page\n",
                builder->current_page->bytes_used, mpack_builder_page_size(writer, builder->current_page));
        mpack_builder_add_page(writer);
        // there is always enough space in a fresh page.
        offset = mpack_builder_align_build(builder->current_page->bytes_used);
    }

    // allocate the build within the page. note that we don't keep track of the
    // space wasted due to the offset. instead the previous build has stored
    // how many bytes follow it, and we'll redo this offset calculation to find
    // this build after it.
    mpack_builder_page_t* page = builder->current_page;
    page->bytes_used = offset + sizeof(mpack_build_t);
    mpack_assert(page->bytes_used <= mpack_builder_page_size(writer, page));
    mpack_build_t* build = (mpack_build_t*)((char*)page + offset);
    mpack_log("created new build %p within page %p, which now has %zi bytes used\n",
            (void*)build, (void*)page, page->bytes_used);

    // configure the new build
    build->parent = builder->current_build;
    build->bytes = 0;
    build->count = 0;
    build->type = type;
    build->key_needs_value = false;
    build->nested_compound_elements = 0;

    mpack_log("setting current and latest build to new build %p\n", (void*)build);
    builder->current_build = build;
    builder->latest_build = build;

    // we always need to provide a buffer that meets the minimum buffer size.
    // if there isn't enough space, we discard the remaining space in the
    // current page and allocate a new one.
    if (mpack_builder_page_remaining(writer, page) < MPACK_WRITER_MINIMUM_BUFFER_SIZE) {
        mpack_log("less than minimum buffer size in current page. %zi bytes used of %zi in this page\n",
                builder->current_page->bytes_used, mpack_builder_page_size(writer, builder->current_page));
        mpack_builder_add_page(writer);
        if (mpack_writer_error(writer) != mpack_ok)
            return;
    }
    mpack_assert(mpack_builder_page_remaining(writer, builder->current_page) >= MPACK_WRITER_MINIMUM_BUFFER_SIZE);
    mpack_builder_configure_buffer(writer);
}

MPACK_NOINLINE
static void mpack_builder_resolve(mpack_writer_t* writer) {
    mpack_builder_t* builder = &writer->builder;

    // We should not have gotten here if we are in an error state. If an error
    // occurs with an open builder, the writer will free the open builder pages
    // when destroyed.
    mpack_assert(mpack_writer_error(writer) == mpack_ok, "can't resolve in error state!");

    // We don't want the user to longjmp out of any I/O errors while we are
    // walking the page list, so defer error callbacks to after we're done.
    mpack_writer_error_t error_fn = writer->error_fn;
    writer->error_fn = NULL;

    // The starting page is the internal storage (if we have it), otherwise
    // it's the first page in the array
    mpack_builder_page_t* page =
        #if MPACK_BUILDER_INTERNAL_STORAGE
        (mpack_builder_page_t*)builder->internal
        #else
        builder->pages
        #endif
        ;

    // We start by restoring the writer's original buffer so we can write the
    // data for real.
    writer->buffer = builder->stash_buffer;
    writer->position = builder->stash_position;
    writer->end = builder->stash_end;

    // We can also close out the build now.
    builder->current_build = NULL;
    builder->latest_build = NULL;
    builder->current_page = NULL;
    builder->pages = NULL;

    // the starting page always starts with the first build
    size_t offset = mpack_builder_align_build(sizeof(mpack_builder_page_t));
    mpack_build_t* build = (mpack_build_t*)((char*)page + offset);
    mpack_log("starting resolve with build %p in page %p\n", (void*)build, (void*)page);

    // encoded data immediately follows the build
    offset += sizeof(mpack_build_t);

    // Walk the list of builds, writing everything out in the buffer. Note that
    // we don't check for errors anywhere. The lower-level write functions will
    // all check for errors and do nothing after an error occurs. We need to
    // walk all pages anyway to free them, so there's not much point in
    // optimizing an error path at the expense of the normal path.
    while (true) {

        // write out the container tag
        mpack_log("writing out an %s with count %" PRIu32 " followed by %zi bytes\n",
                mpack_type_to_string(build->type), build->count, build->bytes);
        switch (build->type) {
            case mpack_type_map:
                mpack_write_map_notrack(writer, build->count);
                break;
            case mpack_type_array:
                mpack_write_array_notrack(writer, build->count);
                break;
            default:
                mpack_break("invalid type in builder?");
                mpack_writer_flag_error(writer, mpack_error_bug);
                return;
        }

        // figure out how many bytes follow this container. we're going to be
        // freeing pages as we write, so we need to be done with this build.
        size_t left = build->bytes;
        build = NULL;

        // write out all bytes following this container
        while (left > 0) {
            size_t bytes_used = page->bytes_used;
            if (offset < bytes_used) {
                size_t step = bytes_used - offset;
                if (step > left)
                    step = left;
                mpack_log("writing out %zi bytes starting at %p in page %p\n",
                        step, (void*)((char*)page + offset), (void*)page);
                mpack_write_native(writer, (char*)page + offset, step);
                offset += step;
                left -= step;
            }

            if (left == 0) {
                mpack_log("done writing bytes for this build\n");
                break;
            }

            // still need to write more bytes. free this page and jump to the
            // next one.
            mpack_builder_page_t* next_page = page->next;
            mpack_builder_free_page(writer, page);
            page = next_page;
            // bytes on the next page immediately follow the header.
            offset = sizeof(mpack_builder_page_t);
        }

        // now see if we can find another build.
        offset = mpack_builder_align_build(offset);
        if (offset + sizeof(mpack_build_t) > mpack_builder_page_size(writer, page)) {
            mpack_log("not enough room in this page for another build\n");
            mpack_builder_page_t* next_page = page->next;
            mpack_builder_free_page(writer, page);
            page = next_page;
            if (page == NULL) {
                mpack_log("no more pages\n");
                // there are no more pages. we're done.
                break;
            }
            offset = mpack_builder_align_build(sizeof(mpack_builder_page_t));
        }
        if (offset + sizeof(mpack_build_t) > page->bytes_used) {
            // there is no more data. we're done.
            mpack_log("no more data\n");
            mpack_builder_free_page(writer, page);
            break;
        }

        // we've found another build. loop around!
        build = (mpack_build_t*)((char*)page + offset);
        offset += sizeof(mpack_build_t);
        mpack_log("found build %p\n", (void*)build);
    }

    mpack_log("done resolve.\n");

    // We can now restore the error handler and call it if an error occurred.
    writer->error_fn = error_fn;
    if (writer->error_fn && mpack_writer_error(writer) != mpack_ok)
        writer->error_fn(writer, writer->error);
}

static void mpack_builder_complete(mpack_writer_t* writer, mpack_type_t type) {
    mpack_writer_track_pop_builder(writer, type);
    if (mpack_writer_error(writer) != mpack_ok)
        return;

    mpack_builder_t* builder = &writer->builder;
    mpack_assert(builder->current_build != NULL, "no build in progress!");
    mpack_assert(builder->latest_build != NULL, "missing latest build!");
    mpack_assert(builder->current_build->type == type, "completing wrong type!");
    mpack_log("completing build %p\n", (void*)builder->current_build);

    if (builder->current_build->key_needs_value) {
        mpack_break("an odd number of elements were written in a map!");
        mpack_writer_flag_error(writer, mpack_error_bug);
        return;
    }

    if (builder->current_build->nested_compound_elements != 0) {
        mpack_break("there is a nested unfinished non-build map or array in this build.");
        mpack_writer_flag_error(writer, mpack_error_bug);
        return;
    }

    // We need to apply whatever writes have been made to the current build
    // before popping it.
    mpack_builder_apply_writes(writer);

    // For a nested build, we just switch the current build back to its parent.
    if (builder->current_build->parent != NULL) {
        mpack_log("setting current build to parent build %p. latest is still %p.\n",
                (void*)builder->current_build->parent, (void*)builder->latest_build);
        builder->current_build = builder->current_build->parent;
        mpack_builder_configure_buffer(writer);
    } else {
        // We're completing the final build.
        mpack_builder_resolve(writer);
    }
}

void mpack_build_map(mpack_writer_t* writer) {
    mpack_builder_build(writer, mpack_type_map);
}

void mpack_build_array(mpack_writer_t* writer) {
    mpack_builder_build(writer, mpack_type_array);
}

void mpack_complete_map(mpack_writer_t* writer) {
    mpack_builder_complete(writer, mpack_type_map);
}

void mpack_complete_array(mpack_writer_t* writer) {
    mpack_builder_complete(writer, mpack_type_array);
}

#endif // MPACK_BUILDER
#endif // MPACK_WRITER

MPACK_SILENCE_WARNINGS_END

/* mpack/mpack-reader.c.c */

#define MPACK_INTERNAL 1

/* #include "mpack-reader.h" */

MPACK_SILENCE_WARNINGS_BEGIN

#if MPACK_READER

static void mpack_reader_skip_using_fill(mpack_reader_t* reader, size_t count);

void mpack_reader_init(mpack_reader_t* reader, char* buffer, size_t size, size_t count) {
    mpack_assert(buffer != NULL, "buffer is NULL");

    mpack_memset(reader, 0, sizeof(*reader));
    reader->buffer = buffer;
    reader->size = size;
    reader->data = buffer;
    reader->end = buffer + count;

    #if MPACK_READ_TRACKING
    mpack_reader_flag_if_error(reader, mpack_track_init(&reader->track));
    #endif

    mpack_log("===========================\n");
    mpack_log("initializing reader with buffer size %i\n", (int)size);
}

void mpack_reader_init_error(mpack_reader_t* reader, mpack_error_t error) {
    mpack_memset(reader, 0, sizeof(*reader));
    reader->error = error;

    mpack_log("===========================\n");
    mpack_log("initializing reader error state %i\n", (int)error);
}

void mpack_reader_init_data(mpack_reader_t* reader, const char* data, size_t count) {
    mpack_assert(data != NULL, "data is NULL");

    mpack_memset(reader, 0, sizeof(*reader));
    reader->data = data;
    reader->end = data + count;

    #if MPACK_READ_TRACKING
    mpack_reader_flag_if_error(reader, mpack_track_init(&reader->track));
    #endif

    mpack_log("===========================\n");
    mpack_log("initializing reader with data size %i\n", (int)count);
}

void mpack_reader_set_fill(mpack_reader_t* reader, mpack_reader_fill_t fill) {
    MPACK_STATIC_ASSERT(MPACK_READER_MINIMUM_BUFFER_SIZE >= MPACK_MAXIMUM_TAG_SIZE,
            "minimum buffer size must fit any tag!");

    if (reader->size == 0) {
        mpack_break("cannot use fill function without a writeable buffer!");
        mpack_reader_flag_error(reader, mpack_error_bug);
        return;
    }

    if (reader->size < MPACK_READER_MINIMUM_BUFFER_SIZE) {
        mpack_break("buffer size is %i, but minimum buffer size for fill is %i",
                (int)reader->size, MPACK_READER_MINIMUM_BUFFER_SIZE);
        mpack_reader_flag_error(reader, mpack_error_bug);
        return;
    }

    reader->fill = fill;
}

void mpack_reader_set_skip(mpack_reader_t* reader, mpack_reader_skip_t skip) {
    mpack_assert(reader->size != 0, "cannot use skip function without a writeable buffer!");
    reader->skip = skip;
}

#if MPACK_STDIO
static size_t mpack_file_reader_fill(mpack_reader_t* reader, char* buffer, size_t count) {
    if (feof((FILE *)reader->context)) {
       mpack_reader_flag_error(reader, mpack_error_eof);
       return 0;
    }
    return fread((void*)buffer, 1, count, (FILE*)reader->context);
}

static void mpack_file_reader_skip(mpack_reader_t* reader, size_t count) {
    if (mpack_reader_error(reader) != mpack_ok)
        return;
    FILE* file = (FILE*)reader->context;

    // We call ftell() to test whether the stream is seekable
    // without causing a file error.
    if (ftell(file) >= 0) {
        mpack_log("seeking forward %i bytes\n", (int)count);
        if (fseek(file, (long int)count, SEEK_CUR) == 0)
            return;
        mpack_log("fseek() didn't return zero!\n");
        if (ferror(file)) {
            mpack_reader_flag_error(reader, mpack_error_io);
            return;
        }
    }

    // If the stream is not seekable, fall back to the fill function.
    mpack_reader_skip_using_fill(reader, count);
}

static void mpack_file_reader_teardown(mpack_reader_t* reader) {
    MPACK_FREE(reader->buffer);
    reader->buffer = NULL;
    reader->context = NULL;
    reader->size = 0;
    reader->fill = NULL;
    reader->skip = NULL;
    reader->teardown = NULL;
}

static void mpack_file_reader_teardown_close(mpack_reader_t* reader) {
    FILE* file = (FILE*)reader->context;

    if (file) {
        int ret = fclose(file);
        if (ret != 0)
            mpack_reader_flag_error(reader, mpack_error_io);
    }

    mpack_file_reader_teardown(reader);
}

void mpack_reader_init_stdfile(mpack_reader_t* reader, FILE* file, bool close_when_done) {
    mpack_assert(file != NULL, "file is NULL");

    size_t capacity = MPACK_BUFFER_SIZE;
    char* buffer = (char*)MPACK_MALLOC(capacity);
    if (buffer == NULL) {
        mpack_reader_init_error(reader, mpack_error_memory);
        if (close_when_done) {
            fclose(file);
        }
        return;
    }

    mpack_reader_init(reader, buffer, capacity, 0);
    mpack_reader_set_context(reader, file);
    mpack_reader_set_fill(reader, mpack_file_reader_fill);
    mpack_reader_set_skip(reader, mpack_file_reader_skip);
    mpack_reader_set_teardown(reader, close_when_done ?
            mpack_file_reader_teardown_close :
            mpack_file_reader_teardown);
}

void mpack_reader_init_filename(mpack_reader_t* reader, const char* filename) {
    mpack_assert(filename != NULL, "filename is NULL");

    FILE* file = fopen(filename, "rb");
    if (file == NULL) {
        mpack_reader_init_error(reader, mpack_error_io);
        return;
    }

    mpack_reader_init_stdfile(reader, file, true);
}
#endif

mpack_error_t mpack_reader_destroy(mpack_reader_t* reader) {

    // clean up tracking, asserting if we're not already in an error state
    #if MPACK_READ_TRACKING
    mpack_reader_flag_if_error(reader, mpack_track_destroy(&reader->track, mpack_reader_error(reader) != mpack_ok));
    #endif

    if (reader->teardown)
        reader->teardown(reader);
    reader->teardown = NULL;

    return reader->error;
}

size_t mpack_reader_remaining(mpack_reader_t* reader, const char** data) {
    if (mpack_reader_error(reader) != mpack_ok)
        return 0;

    #if MPACK_READ_TRACKING
    if (mpack_reader_flag_if_error(reader, mpack_track_check_empty(&reader->track)) != mpack_ok)
        return 0;
    #endif

    if (data)
        *data = reader->data;
    return (size_t)(reader->end - reader->data);
}

void mpack_reader_flag_error(mpack_reader_t* reader, mpack_error_t error) {
    mpack_log("reader %p setting error %i: %s\n", (void*)reader, (int)error, mpack_error_to_string(error));

    if (reader->error == mpack_ok) {
        reader->error = error;
        reader->end = reader->data;
        if (reader->error_fn)
            reader->error_fn(reader, error);
    }
}

// Loops on the fill function, reading between the minimum and
// maximum number of bytes and flagging an error if it fails.
MPACK_NOINLINE static size_t mpack_fill_range(mpack_reader_t* reader, char* p, size_t min_bytes, size_t max_bytes) {
    mpack_assert(reader->fill != NULL, "mpack_fill_range() called with no fill function?");
    mpack_assert(min_bytes > 0, "cannot fill zero bytes!");
    mpack_assert(max_bytes >= min_bytes, "min_bytes %i cannot be larger than max_bytes %i!",
            (int)min_bytes, (int)max_bytes);

    size_t count = 0;
    while (count < min_bytes) {
        size_t read = reader->fill(reader, p + count, max_bytes - count);

        // Reader fill functions can flag an error or return 0 on failure. We
        // also guard against functions that return -1 just in case.
        if (mpack_reader_error(reader) != mpack_ok)
            return 0;
        if (read == 0 || read == ((size_t)(-1))) {
            mpack_reader_flag_error(reader, mpack_error_io);
            return 0;
        }

        count += read;
    }
    return count;
}

MPACK_NOINLINE bool mpack_reader_ensure_straddle(mpack_reader_t* reader, size_t count) {
    mpack_assert(count != 0, "cannot ensure zero bytes!");
    mpack_assert(reader->error == mpack_ok, "reader cannot be in an error state!");

    mpack_assert(count > (size_t)(reader->end - reader->data),
            "straddling ensure requested for %i bytes, but there are %i bytes "
            "left in buffer. call mpack_reader_ensure() instead",
            (int)count, (int)(reader->end - reader->data));

    // we'll need a fill function to get more data. if there's no
    // fill function, the buffer should contain an entire MessagePack
    // object, so we raise mpack_error_invalid instead of mpack_error_io
    // on truncated data.
    if (reader->fill == NULL) {
        mpack_reader_flag_error(reader, mpack_error_invalid);
        return false;
    }

    // we need enough space in the buffer. if the buffer is not
    // big enough, we return mpack_error_too_big (since this is
    // for an in-place read larger than the buffer size.)
    if (count > reader->size) {
        mpack_reader_flag_error(reader, mpack_error_too_big);
        return false;
    }

    // move the existing data to the start of the buffer
    size_t left = (size_t)(reader->end - reader->data);
    mpack_memmove(reader->buffer, reader->data, left);
    reader->end -= reader->data - reader->buffer;
    reader->data = reader->buffer;

    // read at least the necessary number of bytes, accepting up to the
    // buffer size
    size_t read = mpack_fill_range(reader, reader->buffer + left,
            count - left, reader->size - left);
    if (mpack_reader_error(reader) != mpack_ok)
        return false;
    reader->end += read;
    return true;
}

// Reads count bytes into p. Used when there are not enough bytes
// left in the buffer to satisfy a read.
MPACK_NOINLINE void mpack_read_native_straddle(mpack_reader_t* reader, char* p, size_t count) {
    mpack_assert(count == 0 || p != NULL, "data pointer for %i bytes is NULL", (int)count);

    if (mpack_reader_error(reader) != mpack_ok) {
        mpack_memset(p, 0, count);
        return;
    }

    size_t left = (size_t)(reader->end - reader->data);
    mpack_log("big read for %i bytes into %p, %i left in buffer, buffer size %i\n",
            (int)count, p, (int)left, (int)reader->size);

    if (count <= left) {
        mpack_assert(0,
                "big read requested for %i bytes, but there are %i bytes "
                "left in buffer. call mpack_read_native() instead",
                (int)count, (int)left);
        mpack_reader_flag_error(reader, mpack_error_bug);
        mpack_memset(p, 0, count);
        return;
    }

    // we'll need a fill function to get more data. if there's no
    // fill function, the buffer should contain an entire MessagePack
    // object, so we raise mpack_error_invalid instead of mpack_error_io
    // on truncated data.
    if (reader->fill == NULL) {
        mpack_reader_flag_error(reader, mpack_error_invalid);
        mpack_memset(p, 0, count);
        return;
    }

    if (reader->size == 0) {
        // somewhat debatable what error should be returned here. when
        // initializing a reader with an in-memory buffer it's not
        // necessarily a bug if the data is blank; it might just have
        // been truncated to zero. for this reason we return the same
        // error as if the data was truncated.
        mpack_reader_flag_error(reader, mpack_error_io);
        mpack_memset(p, 0, count);
        return;
    }

    // flush what's left of the buffer
    if (left > 0) {
        mpack_log("flushing %i bytes remaining in buffer\n", (int)left);
        mpack_memcpy(p, reader->data, left);
        count -= left;
        p += left;
        reader->data += left;
    }

    // if the remaining data needed is some small fraction of the
    // buffer size, we'll try to fill the buffer as much as possible
    // and copy the needed data out.
    if (count <= reader->size / MPACK_READER_SMALL_FRACTION_DENOMINATOR) {
        size_t read = mpack_fill_range(reader, reader->buffer, count, reader->size);
        if (mpack_reader_error(reader) != mpack_ok)
            return;
        mpack_memcpy(p, reader->buffer, count);
        reader->data = reader->buffer + count;
        reader->end = reader->buffer + read;

    // otherwise we read the remaining data directly into the target.
    } else {
        mpack_log("reading %i additional bytes\n", (int)count);
        mpack_fill_range(reader, p, count, count);
    }
}

MPACK_NOINLINE static void mpack_skip_bytes_straddle(mpack_reader_t* reader, size_t count) {

    // we'll need at least a fill function to skip more data. if there's
    // no fill function, the buffer should contain an entire MessagePack
    // object, so we raise mpack_error_invalid instead of mpack_error_io
    // on truncated data. (see mpack_read_native_straddle())
    if (reader->fill == NULL) {
        mpack_log("reader has no fill function!\n");
        mpack_reader_flag_error(reader, mpack_error_invalid);
        return;
    }

    // discard whatever's left in the buffer
    size_t left = (size_t)(reader->end - reader->data);
    mpack_log("discarding %i bytes still in buffer\n", (int)left);
    count -= left;
    reader->data = reader->end;

    // use the skip function if we've got one, and if we're trying
    // to skip a lot of data. if we only need to skip some tiny
    // fraction of the buffer size, it's probably better to just
    // fill the buffer and skip from it instead of trying to seek.
    if (reader->skip && count > reader->size / 16) {
        mpack_log("calling skip function for %i bytes\n", (int)count);
        reader->skip(reader, count);
        return;
    }

    mpack_reader_skip_using_fill(reader, count);
}

void mpack_skip_bytes(mpack_reader_t* reader, size_t count) {
    if (mpack_reader_error(reader) != mpack_ok)
        return;
    mpack_log("skip requested for %i bytes\n", (int)count);

    mpack_reader_track_bytes(reader, count);

    // check if we have enough in the buffer already
    size_t left = (size_t)(reader->end - reader->data);
    if (left >= count) {
        mpack_log("skipping %" PRIu32 " bytes still in buffer\n", (uint32_t)count);
        reader->data += count;
        return;
    }

    mpack_skip_bytes_straddle(reader, count);
}

MPACK_NOINLINE static void mpack_reader_skip_using_fill(mpack_reader_t* reader, size_t count) {
    mpack_assert(reader->fill != NULL, "missing fill function!");
    mpack_assert(reader->data == reader->end, "there are bytes left in the buffer!");
    mpack_assert(reader->error == mpack_ok, "should not have called this in an error state (%i)", reader->error);
    mpack_log("skip using fill for %i bytes\n", (int)count);

    // fill and discard multiples of the buffer size
    while (count > reader->size) {
        mpack_log("filling and discarding buffer of %i bytes\n", (int)reader->size);
        if (mpack_fill_range(reader, reader->buffer, reader->size, reader->size) < reader->size) {
            mpack_reader_flag_error(reader, mpack_error_io);
            return;
        }
        count -= reader->size;
    }

    // fill the buffer as much as possible
    reader->data = reader->buffer;
    size_t read = mpack_fill_range(reader, reader->buffer, count, reader->size);
    if (read < count) {
        mpack_reader_flag_error(reader, mpack_error_io);
        return;
    }
    reader->end = reader->data + read;
    mpack_log("filled %i bytes into buffer; discarding %i bytes\n", (int)read, (int)count);
    reader->data += count;
}

void mpack_read_bytes(mpack_reader_t* reader, char* p, size_t count) {
    mpack_assert(p != NULL, "destination for read of %i bytes is NULL", (int)count);
    mpack_reader_track_bytes(reader, count);
    mpack_read_native(reader, p, count);
}

void mpack_read_utf8(mpack_reader_t* reader, char* p, size_t byte_count) {
    mpack_assert(p != NULL, "destination for read of %i bytes is NULL", (int)byte_count);
    mpack_reader_track_str_bytes_all(reader, byte_count);
    mpack_read_native(reader, p, byte_count);

    if (mpack_reader_error(reader) == mpack_ok && !mpack_utf8_check(p, byte_count))
        mpack_reader_flag_error(reader, mpack_error_type);
}

static void mpack_read_cstr_unchecked(mpack_reader_t* reader, char* buf, size_t buffer_size, size_t byte_count) {
    mpack_assert(buf != NULL, "destination for read of %i bytes is NULL", (int)byte_count);
    mpack_assert(buffer_size >= 1, "buffer size is zero; you must have room for at least a null-terminator");

    if (mpack_reader_error(reader)) {
        buf[0] = 0;
        return;
    }

    if (byte_count > buffer_size - 1) {
        mpack_reader_flag_error(reader, mpack_error_too_big);
        buf[0] = 0;
        return;
    }

    mpack_reader_track_str_bytes_all(reader, byte_count);
    mpack_read_native(reader, buf, byte_count);
    buf[byte_count] = 0;
}

void mpack_read_cstr(mpack_reader_t* reader, char* buf, size_t buffer_size, size_t byte_count) {
    mpack_read_cstr_unchecked(reader, buf, buffer_size, byte_count);

    // check for null bytes
    if (mpack_reader_error(reader) == mpack_ok && !mpack_str_check_no_null(buf, byte_count)) {
        buf[0] = 0;
        mpack_reader_flag_error(reader, mpack_error_type);
    }
}

void mpack_read_utf8_cstr(mpack_reader_t* reader, char* buf, size_t buffer_size, size_t byte_count) {
    mpack_read_cstr_unchecked(reader, buf, buffer_size, byte_count);

    // check encoding
    if (mpack_reader_error(reader) == mpack_ok && !mpack_utf8_check_no_null(buf, byte_count)) {
        buf[0] = 0;
        mpack_reader_flag_error(reader, mpack_error_type);
    }
}

#ifdef MPACK_MALLOC
// Reads native bytes with error callback disabled. This allows MPack reader functions
// to hold an allocated buffer and read native data into it without leaking it in
// case of a non-local jump (longjmp, throw) out of an error handler.
static void mpack_read_native_noerrorfn(mpack_reader_t* reader, char* p, size_t count) {
    mpack_assert(reader->error == mpack_ok, "cannot call if an error is already flagged!");
    mpack_reader_error_t error_fn = reader->error_fn;
    reader->error_fn = NULL;
    mpack_read_native(reader, p, count);
    reader->error_fn = error_fn;
}

char* mpack_read_bytes_alloc_impl(mpack_reader_t* reader, size_t count, bool null_terminated) {

    // track the bytes first in case it jumps
    mpack_reader_track_bytes(reader, count);
    if (mpack_reader_error(reader) != mpack_ok)
        return NULL;

    // cannot allocate zero bytes. this is not an error.
    if (count == 0 && null_terminated == false)
        return NULL;

    // allocate data
    char* data = (char*)MPACK_MALLOC(count + (null_terminated ? 1 : 0)); // TODO: can this overflow?
    if (data == NULL) {
        mpack_reader_flag_error(reader, mpack_error_memory);
        return NULL;
    }

    // read with error callback disabled so we don't leak our buffer
    mpack_read_native_noerrorfn(reader, data, count);

    // report flagged errors
    if (mpack_reader_error(reader) != mpack_ok) {
        MPACK_FREE(data);
        if (reader->error_fn)
            reader->error_fn(reader, mpack_reader_error(reader));
        return NULL;
    }

    if (null_terminated)
        data[count] = '\0';
    return data;
}
#endif

// read inplace without tracking (since there are different
// tracking modes for different inplace readers)
static const char* mpack_read_bytes_inplace_notrack(mpack_reader_t* reader, size_t count) {
    if (mpack_reader_error(reader) != mpack_ok)
        return NULL;

    // if we have enough bytes already in the buffer, we can return it directly.
    if ((size_t)(reader->end - reader->data) >= count) {
        const char* bytes = reader->data;
        reader->data += count;
        return bytes;
    }

    if (!mpack_reader_ensure(reader, count))
        return NULL;

    const char* bytes = reader->data;
    reader->data += count;
    return bytes;
}

const char* mpack_read_bytes_inplace(mpack_reader_t* reader, size_t count) {
    mpack_reader_track_bytes(reader, count);
    return mpack_read_bytes_inplace_notrack(reader, count);
}

const char* mpack_read_utf8_inplace(mpack_reader_t* reader, size_t count) {
    mpack_reader_track_str_bytes_all(reader, count);
    const char* str = mpack_read_bytes_inplace_notrack(reader, count);

    if (mpack_reader_error(reader) == mpack_ok && !mpack_utf8_check(str, count)) {
        mpack_reader_flag_error(reader, mpack_error_type);
        return NULL;
    }

    return str;
}

static size_t mpack_parse_tag(mpack_reader_t* reader, mpack_tag_t* tag) {
    mpack_assert(reader->error == mpack_ok, "reader cannot be in an error state!");

    if (!mpack_reader_ensure(reader, 1))
        return 0;
    uint8_t type = mpack_load_u8(reader->data);

    // unfortunately, by far the fastest way to parse a tag is to switch
    // on the first byte, and to explicitly list every possible byte. so for
    // infix types, the list of cases is quite large.
    //
    // in size-optimized builds, we switch on the top four bits first to
    // handle most infix types with a smaller jump table to save space.

    #if MPACK_OPTIMIZE_FOR_SIZE
    switch (type >> 4) {

        // positive fixnum
        case 0x0: case 0x1: case 0x2: case 0x3:
        case 0x4: case 0x5: case 0x6: case 0x7:
            *tag = mpack_tag_make_uint(type);
            return 1;

        // negative fixnum
        case 0xe: case 0xf:
            *tag = mpack_tag_make_int((int8_t)type);
            return 1;

        // fixmap
        case 0x8:
            *tag = mpack_tag_make_map(type & ~0xf0u);
            return 1;

        // fixarray
        case 0x9:
            *tag = mpack_tag_make_array(type & ~0xf0u);
            return 1;

        // fixstr
        case 0xa: case 0xb:
            *tag = mpack_tag_make_str(type & ~0xe0u);
            return 1;

        // not one of the common infix types
        default:
            break;

    }
    #endif

    // handle individual type tags
    switch (type) {

        #if !MPACK_OPTIMIZE_FOR_SIZE
        // positive fixnum
        case 0x00: case 0x01: case 0x02: case 0x03: case 0x04: case 0x05: case 0x06: case 0x07:
        case 0x08: case 0x09: case 0x0a: case 0x0b: case 0x0c: case 0x0d: case 0x0e: case 0x0f:
        case 0x10: case 0x11: case 0x12: case 0x13: case 0x14: case 0x15: case 0x16: case 0x17:
        case 0x18: case 0x19: case 0x1a: case 0x1b: case 0x1c: case 0x1d: case 0x1e: case 0x1f:
        case 0x20: case 0x21: case 0x22: case 0x23: case 0x24: case 0x25: case 0x26: case 0x27:
        case 0x28: case 0x29: case 0x2a: case 0x2b: case 0x2c: case 0x2d: case 0x2e: case 0x2f:
        case 0x30: case 0x31: case 0x32: case 0x33: case 0x34: case 0x35: case 0x36: case 0x37:
        case 0x38: case 0x39: case 0x3a: case 0x3b: case 0x3c: case 0x3d: case 0x3e: case 0x3f:
        case 0x40: case 0x41: case 0x42: case 0x43: case 0x44: case 0x45: case 0x46: case 0x47:
        case 0x48: case 0x49: case 0x4a: case 0x4b: case 0x4c: case 0x4d: case 0x4e: case 0x4f:
        case 0x50: case 0x51: case 0x52: case 0x53: case 0x54: case 0x55: case 0x56: case 0x57:
        case 0x58: case 0x59: case 0x5a: case 0x5b: case 0x5c: case 0x5d: case 0x5e: case 0x5f:
        case 0x60: case 0x61: case 0x62: case 0x63: case 0x64: case 0x65: case 0x66: case 0x67:
        case 0x68: case 0x69: case 0x6a: case 0x6b: case 0x6c: case 0x6d: case 0x6e: case 0x6f:
        case 0x70: case 0x71: case 0x72: case 0x73: case 0x74: case 0x75: case 0x76: case 0x77:
        case 0x78: case 0x79: case 0x7a: case 0x7b: case 0x7c: case 0x7d: case 0x7e: case 0x7f:
            *tag = mpack_tag_make_uint(type);
            return 1;

        // negative fixnum
        case 0xe0: case 0xe1: case 0xe2: case 0xe3: case 0xe4: case 0xe5: case 0xe6: case 0xe7:
        case 0xe8: case 0xe9: case 0xea: case 0xeb: case 0xec: case 0xed: case 0xee: case 0xef:
        case 0xf0: case 0xf1: case 0xf2: case 0xf3: case 0xf4: case 0xf5: case 0xf6: case 0xf7:
        case 0xf8: case 0xf9: case 0xfa: case 0xfb: case 0xfc: case 0xfd: case 0xfe: case 0xff:
            *tag = mpack_tag_make_int((int8_t)type);
            return 1;

        // fixmap
        case 0x80: case 0x81: case 0x82: case 0x83: case 0x84: case 0x85: case 0x86: case 0x87:
        case 0x88: case 0x89: case 0x8a: case 0x8b: case 0x8c: case 0x8d: case 0x8e: case 0x8f:
            *tag = mpack_tag_make_map(type & ~0xf0u);
            return 1;

        // fixarray
        case 0x90: case 0x91: case 0x92: case 0x93: case 0x94: case 0x95: case 0x96: case 0x97:
        case 0x98: case 0x99: case 0x9a: case 0x9b: case 0x9c: case 0x9d: case 0x9e: case 0x9f:
            *tag = mpack_tag_make_array(type & ~0xf0u);
            return 1;

        // fixstr
        case 0xa0: case 0xa1: case 0xa2: case 0xa3: case 0xa4: case 0xa5: case 0xa6: case 0xa7:
        case 0xa8: case 0xa9: case 0xaa: case 0xab: case 0xac: case 0xad: case 0xae: case 0xaf:
        case 0xb0: case 0xb1: case 0xb2: case 0xb3: case 0xb4: case 0xb5: case 0xb6: case 0xb7:
        case 0xb8: case 0xb9: case 0xba: case 0xbb: case 0xbc: case 0xbd: case 0xbe: case 0xbf:
            *tag = mpack_tag_make_str(type & ~0xe0u);
            return 1;
        #endif

        // nil
        case 0xc0:
            *tag = mpack_tag_make_nil();
            return 1;

        // bool
        case 0xc2: case 0xc3:
            *tag = mpack_tag_make_bool((bool)(type & 1));
            return 1;

        // bin8
        case 0xc4:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_BIN8))
                return 0;
            *tag = mpack_tag_make_bin(mpack_load_u8(reader->data + 1));
            return MPACK_TAG_SIZE_BIN8;

        // bin16
        case 0xc5:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_BIN16))
                return 0;
            *tag = mpack_tag_make_bin(mpack_load_u16(reader->data + 1));
            return MPACK_TAG_SIZE_BIN16;

        // bin32
        case 0xc6:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_BIN32))
                return 0;
            *tag = mpack_tag_make_bin(mpack_load_u32(reader->data + 1));
            return MPACK_TAG_SIZE_BIN32;

        #if MPACK_EXTENSIONS
        // ext8
        case 0xc7:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_EXT8))
                return 0;
            *tag = mpack_tag_make_ext(mpack_load_i8(reader->data + 2), mpack_load_u8(reader->data + 1));
            return MPACK_TAG_SIZE_EXT8;

        // ext16
        case 0xc8:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_EXT16))
                return 0;
            *tag = mpack_tag_make_ext(mpack_load_i8(reader->data + 3), mpack_load_u16(reader->data + 1));
            return MPACK_TAG_SIZE_EXT16;

        // ext32
        case 0xc9:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_EXT32))
                return 0;
            *tag = mpack_tag_make_ext(mpack_load_i8(reader->data + 5), mpack_load_u32(reader->data + 1));
            return MPACK_TAG_SIZE_EXT32;
        #endif

        // float
        case 0xca:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_FLOAT))
                return 0;
            #if MPACK_FLOAT
            *tag = mpack_tag_make_float(mpack_load_float(reader->data + 1));
            #else
            *tag = mpack_tag_make_raw_float(mpack_load_u32(reader->data + 1));
            #endif
            return MPACK_TAG_SIZE_FLOAT;

        // double
        case 0xcb:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_DOUBLE))
                return 0;
            #if MPACK_DOUBLE
            *tag = mpack_tag_make_double(mpack_load_double(reader->data + 1));
            #else
            *tag = mpack_tag_make_raw_double(mpack_load_u64(reader->data + 1));
            #endif
            return MPACK_TAG_SIZE_DOUBLE;

        // uint8
        case 0xcc:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_U8))
                return 0;
            *tag = mpack_tag_make_uint(mpack_load_u8(reader->data + 1));
            return MPACK_TAG_SIZE_U8;

        // uint16
        case 0xcd:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_U16))
                return 0;
            *tag = mpack_tag_make_uint(mpack_load_u16(reader->data + 1));
            return MPACK_TAG_SIZE_U16;

        // uint32
        case 0xce:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_U32))
                return 0;
            *tag = mpack_tag_make_uint(mpack_load_u32(reader->data + 1));
            return MPACK_TAG_SIZE_U32;

        // uint64
        case 0xcf:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_U64))
                return 0;
            *tag = mpack_tag_make_uint(mpack_load_u64(reader->data + 1));
            return MPACK_TAG_SIZE_U64;

        // int8
        case 0xd0:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_I8))
                return 0;
            *tag = mpack_tag_make_int(mpack_load_i8(reader->data + 1));
            return MPACK_TAG_SIZE_I8;

        // int16
        case 0xd1:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_I16))
                return 0;
            *tag = mpack_tag_make_int(mpack_load_i16(reader->data + 1));
            return MPACK_TAG_SIZE_I16;

        // int32
        case 0xd2:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_I32))
                return 0;
            *tag = mpack_tag_make_int(mpack_load_i32(reader->data + 1));
            return MPACK_TAG_SIZE_I32;

        // int64
        case 0xd3:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_I64))
                return 0;
            *tag = mpack_tag_make_int(mpack_load_i64(reader->data + 1));
            return MPACK_TAG_SIZE_I64;

        #if MPACK_EXTENSIONS
        // fixext1
        case 0xd4:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_FIXEXT1))
                return 0;
            *tag = mpack_tag_make_ext(mpack_load_i8(reader->data + 1), 1);
            return MPACK_TAG_SIZE_FIXEXT1;

        // fixext2
        case 0xd5:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_FIXEXT2))
                return 0;
            *tag = mpack_tag_make_ext(mpack_load_i8(reader->data + 1), 2);
            return MPACK_TAG_SIZE_FIXEXT2;

        // fixext4
        case 0xd6:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_FIXEXT4))
                return 0;
            *tag = mpack_tag_make_ext(mpack_load_i8(reader->data + 1), 4);
            return 2;

        // fixext8
        case 0xd7:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_FIXEXT8))
                return 0;
            *tag = mpack_tag_make_ext(mpack_load_i8(reader->data + 1), 8);
            return MPACK_TAG_SIZE_FIXEXT8;

        // fixext16
        case 0xd8:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_FIXEXT16))
                return 0;
            *tag = mpack_tag_make_ext(mpack_load_i8(reader->data + 1), 16);
            return MPACK_TAG_SIZE_FIXEXT16;
        #endif

        // str8
        case 0xd9:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_STR8))
                return 0;
            *tag = mpack_tag_make_str(mpack_load_u8(reader->data + 1));
            return MPACK_TAG_SIZE_STR8;

        // str16
        case 0xda:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_STR16))
                return 0;
            *tag = mpack_tag_make_str(mpack_load_u16(reader->data + 1));
            return MPACK_TAG_SIZE_STR16;

        // str32
        case 0xdb:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_STR32))
                return 0;
            *tag = mpack_tag_make_str(mpack_load_u32(reader->data + 1));
            return MPACK_TAG_SIZE_STR32;

        // array16
        case 0xdc:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_ARRAY16))
                return 0;
            *tag = mpack_tag_make_array(mpack_load_u16(reader->data + 1));
            return MPACK_TAG_SIZE_ARRAY16;

        // array32
        case 0xdd:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_ARRAY32))
                return 0;
            *tag = mpack_tag_make_array(mpack_load_u32(reader->data + 1));
            return MPACK_TAG_SIZE_ARRAY32;

        // map16
        case 0xde:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_MAP16))
                return 0;
            *tag = mpack_tag_make_map(mpack_load_u16(reader->data + 1));
            return MPACK_TAG_SIZE_MAP16;

        // map32
        case 0xdf:
            if (!mpack_reader_ensure(reader, MPACK_TAG_SIZE_MAP32))
                return 0;
            *tag = mpack_tag_make_map(mpack_load_u32(reader->data + 1));
            return MPACK_TAG_SIZE_MAP32;

        // reserved
        case 0xc1:
            mpack_reader_flag_error(reader, mpack_error_invalid);
            return 0;

        #if !MPACK_EXTENSIONS
        // ext
        case 0xc7: // fallthrough
        case 0xc8: // fallthrough
        case 0xc9: // fallthrough
        // fixext
        case 0xd4: // fallthrough
        case 0xd5: // fallthrough
        case 0xd6: // fallthrough
        case 0xd7: // fallthrough
        case 0xd8:
            mpack_reader_flag_error(reader, mpack_error_unsupported);
            return 0;
        #endif

        #if MPACK_OPTIMIZE_FOR_SIZE
        // any other bytes should have been handled by the infix switch
        default:
            break;
        #endif
    }

    mpack_assert(0, "unreachable");
    return 0;
}

mpack_tag_t mpack_read_tag(mpack_reader_t* reader) {
    mpack_log("reading tag\n");

    // make sure we can read a tag
    if (mpack_reader_error(reader) != mpack_ok)
        return mpack_tag_nil();
    if (mpack_reader_track_element(reader) != mpack_ok)
        return mpack_tag_nil();

    mpack_tag_t tag = MPACK_TAG_ZERO;
    size_t count = mpack_parse_tag(reader, &tag);
    if (count == 0)
        return mpack_tag_nil();

    #if MPACK_READ_TRACKING
    mpack_error_t track_error = mpack_ok;

    switch (tag.type) {
        case mpack_type_map:
        case mpack_type_array:
            track_error = mpack_track_push(&reader->track, tag.type, tag.v.n);
            break;
        #if MPACK_EXTENSIONS
        case mpack_type_ext:
        #endif
        case mpack_type_str:
        case mpack_type_bin:
            track_error = mpack_track_push(&reader->track, tag.type, tag.v.l);
            break;
        default:
            break;
    }

    if (track_error != mpack_ok) {
        mpack_reader_flag_error(reader, track_error);
        return mpack_tag_nil();
    }
    #endif

    reader->data += count;
    return tag;
}

mpack_tag_t mpack_peek_tag(mpack_reader_t* reader) {
    mpack_log("peeking tag\n");

    // make sure we can peek a tag
    if (mpack_reader_error(reader) != mpack_ok)
        return mpack_tag_nil();
    if (mpack_reader_track_peek_element(reader) != mpack_ok)
        return mpack_tag_nil();

    mpack_tag_t tag = MPACK_TAG_ZERO;
    if (mpack_parse_tag(reader, &tag) == 0)
        return mpack_tag_nil();
    return tag;
}

void mpack_discard(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (mpack_reader_error(reader))
        return;
    switch (var.type) {
        case mpack_type_str:
            mpack_skip_bytes(reader, var.v.l);
            mpack_done_str(reader);
            break;
        case mpack_type_bin:
            mpack_skip_bytes(reader, var.v.l);
            mpack_done_bin(reader);
            break;
        #if MPACK_EXTENSIONS
        case mpack_type_ext:
            mpack_skip_bytes(reader, var.v.l);
            mpack_done_ext(reader);
            break;
        #endif
        case mpack_type_array: {
            for (; var.v.n > 0; --var.v.n) {
                mpack_discard(reader);
                if (mpack_reader_error(reader))
                    break;
            }
            mpack_done_array(reader);
            break;
        }
        case mpack_type_map: {
            for (; var.v.n > 0; --var.v.n) {
                mpack_discard(reader);
                mpack_discard(reader);
                if (mpack_reader_error(reader))
                    break;
            }
            mpack_done_map(reader);
            break;
        }
        default:
            break;
    }
}

#if MPACK_EXTENSIONS
mpack_timestamp_t mpack_read_timestamp(mpack_reader_t* reader, size_t size) {
    mpack_timestamp_t timestamp = {0, 0};

    if (size != 4 && size != 8 && size != 12) {
        mpack_reader_flag_error(reader, mpack_error_invalid);
        return timestamp;
    }

    char buf[12];
    mpack_read_bytes(reader, buf, size);
    mpack_done_ext(reader);
    if (mpack_reader_error(reader) != mpack_ok)
        return timestamp;

    switch (size) {
        case 4:
            timestamp.seconds = (int64_t)(uint64_t)mpack_load_u32(buf);
            break;

        case 8: {
            uint64_t packed = mpack_load_u64(buf);
            timestamp.seconds = (int64_t)(packed & ((MPACK_UINT64_C(1) << 34) - 1));
            timestamp.nanoseconds = (uint32_t)(packed >> 34);
            break;
        }

        case 12:
            timestamp.nanoseconds = mpack_load_u32(buf);
            timestamp.seconds = mpack_load_i64(buf + 4);
            break;

        default:
            mpack_assert(false, "unreachable");
            break;
    }

    if (timestamp.nanoseconds > MPACK_TIMESTAMP_NANOSECONDS_MAX) {
        mpack_reader_flag_error(reader, mpack_error_invalid);
        mpack_timestamp_t zero = {0, 0};
        return zero;
    }

    return timestamp;
}
#endif

#if MPACK_READ_TRACKING
void mpack_done_type(mpack_reader_t* reader, mpack_type_t type) {
    if (mpack_reader_error(reader) == mpack_ok)
        mpack_reader_flag_if_error(reader, mpack_track_pop(&reader->track, type));
}
#endif

#if MPACK_DEBUG && MPACK_STDIO
static size_t mpack_print_read_prefix(mpack_reader_t* reader, size_t length, char* buffer, size_t buffer_size) {
    if (length == 0)
        return 0;

    size_t read = (length < buffer_size) ? length : buffer_size;
    mpack_read_bytes(reader, buffer, read);
    if (mpack_reader_error(reader) != mpack_ok)
        return 0;

    mpack_skip_bytes(reader, length - read);
    return read;
}

static void mpack_print_element(mpack_reader_t* reader, mpack_print_t* print, size_t depth) {
    mpack_tag_t val = mpack_read_tag(reader);
    if (mpack_reader_error(reader) != mpack_ok)
        return;

    // We read some bytes from bin and ext so we can print its prefix in hex.
    char buffer[MPACK_PRINT_BYTE_COUNT];
    size_t count = 0;
    size_t i, j;

    switch (val.type) {
        case mpack_type_str:
            mpack_print_append_cstr(print, "\"");
            for (i = 0; i < val.v.l; ++i) {
                char c;
                mpack_read_bytes(reader, &c, 1);
                if (mpack_reader_error(reader) != mpack_ok)
                    return;
                switch (c) {
                    case '\n': mpack_print_append_cstr(print, "\\n"); break;
                    case '\\': mpack_print_append_cstr(print, "\\\\"); break;
                    case '"': mpack_print_append_cstr(print, "\\\""); break;
                    default: mpack_print_append(print, &c, 1); break;
                }
            }
            mpack_print_append_cstr(print, "\"");
            mpack_done_str(reader);
            return;

        case mpack_type_array:
            mpack_print_append_cstr(print, "[\n");
            for (i = 0; i < val.v.n; ++i) {
                for (j = 0; j < depth + 1; ++j)
                    mpack_print_append_cstr(print, "    ");
                mpack_print_element(reader, print, depth + 1);
                if (mpack_reader_error(reader) != mpack_ok)
                    return;
                if (i != val.v.n - 1)
                    mpack_print_append_cstr(print, ",");
                mpack_print_append_cstr(print, "\n");
            }
            for (i = 0; i < depth; ++i)
                mpack_print_append_cstr(print, "    ");
            mpack_print_append_cstr(print, "]");
            mpack_done_array(reader);
            return;

        case mpack_type_map:
            mpack_print_append_cstr(print, "{\n");
            for (i = 0; i < val.v.n; ++i) {
                for (j = 0; j < depth + 1; ++j)
                    mpack_print_append_cstr(print, "    ");
                mpack_print_element(reader, print, depth + 1);
                if (mpack_reader_error(reader) != mpack_ok)
                    return;
                mpack_print_append_cstr(print, ": ");
                mpack_print_element(reader, print, depth + 1);
                if (mpack_reader_error(reader) != mpack_ok)
                    return;
                if (i != val.v.n - 1)
                    mpack_print_append_cstr(print, ",");
                mpack_print_append_cstr(print, "\n");
            }
            for (i = 0; i < depth; ++i)
                mpack_print_append_cstr(print, "    ");
            mpack_print_append_cstr(print, "}");
            mpack_done_map(reader);
            return;

        // The above cases return so as not to print a pseudo-json value. The
        // below cases break and print pseudo-json.

        case mpack_type_bin:
            count = mpack_print_read_prefix(reader, mpack_tag_bin_length(&val), buffer, sizeof(buffer));
            mpack_done_bin(reader);
            break;

        #if MPACK_EXTENSIONS
        case mpack_type_ext:
            count = mpack_print_read_prefix(reader, mpack_tag_ext_length(&val), buffer, sizeof(buffer));
            mpack_done_ext(reader);
            break;
        #endif

        default:
            break;
    }

    char buf[256];
    mpack_tag_debug_pseudo_json(val, buf, sizeof(buf), buffer, count);
    mpack_print_append_cstr(print, buf);
}

static void mpack_print_and_destroy(mpack_reader_t* reader, mpack_print_t* print, size_t depth) {
    size_t i;
    for (i = 0; i < depth; ++i)
        mpack_print_append_cstr(print, "    ");
    mpack_print_element(reader, print, depth);

    size_t remaining = mpack_reader_remaining(reader, NULL);

    char buf[256];
    if (mpack_reader_destroy(reader) != mpack_ok) {
        mpack_snprintf(buf, sizeof(buf), "\n<mpack parsing error %s>", mpack_error_to_string(mpack_reader_error(reader)));
        buf[sizeof(buf) - 1] = '\0';
        mpack_print_append_cstr(print, buf);
    } else if (remaining > 0) {
        mpack_snprintf(buf, sizeof(buf), "\n<%i extra bytes at end of message>", (int)remaining);
        buf[sizeof(buf) - 1] = '\0';
        mpack_print_append_cstr(print, buf);
    }
}

static void mpack_print_data(const char* data, size_t len, mpack_print_t* print, size_t depth) {
    mpack_reader_t reader;
    mpack_reader_init_data(&reader, data, len);
    mpack_print_and_destroy(&reader, print, depth);
}

void mpack_print_data_to_buffer(const char* data, size_t data_size, char* buffer, size_t buffer_size) {
    if (buffer_size == 0) {
        mpack_assert(false, "buffer size is zero!");
        return;
    }

    mpack_print_t print;
    mpack_memset(&print, 0, sizeof(print));
    print.buffer = buffer;
    print.size = buffer_size;
    mpack_print_data(data, data_size, &print, 0);
    mpack_print_append(&print, "",  1); // null-terminator
    mpack_print_flush(&print);

    // we always make sure there's a null-terminator at the end of the buffer
    // in case we ran out of space.
    print.buffer[print.size - 1] = '\0';
}

void mpack_print_data_to_callback(const char* data, size_t size, mpack_print_callback_t callback, void* context) {
    char buffer[1024];
    mpack_print_t print;
    mpack_memset(&print, 0, sizeof(print));
    print.buffer = buffer;
    print.size = sizeof(buffer);
    print.callback = callback;
    print.context = context;
    mpack_print_data(data, size, &print, 0);
    mpack_print_flush(&print);
}

void mpack_print_data_to_file(const char* data, size_t len, FILE* file) {
    mpack_assert(data != NULL, "data is NULL");
    mpack_assert(file != NULL, "file is NULL");

    char buffer[1024];
    mpack_print_t print;
    mpack_memset(&print, 0, sizeof(print));
    print.buffer = buffer;
    print.size = sizeof(buffer);
    print.callback = &mpack_print_file_callback;
    print.context = file;

    mpack_print_data(data, len, &print, 2);
    mpack_print_append_cstr(&print, "\n");
    mpack_print_flush(&print);
}

void mpack_print_stdfile_to_callback(FILE* file, mpack_print_callback_t callback, void* context) {
    char buffer[1024];
    mpack_print_t print;
    mpack_memset(&print, 0, sizeof(print));
    print.buffer = buffer;
    print.size = sizeof(buffer);
    print.callback = callback;
    print.context = context;

    mpack_reader_t reader;
    mpack_reader_init_stdfile(&reader, file, false);
    mpack_print_and_destroy(&reader, &print, 0);
    mpack_print_flush(&print);
}
#endif

#endif

MPACK_SILENCE_WARNINGS_END

/* mpack/mpack-expect.c.c */

#define MPACK_INTERNAL 1

/* #include "mpack-expect.h" */

MPACK_SILENCE_WARNINGS_BEGIN

#if MPACK_EXPECT


// Helpers

MPACK_STATIC_INLINE uint8_t mpack_expect_native_u8(mpack_reader_t* reader) {
    if (mpack_reader_error(reader) != mpack_ok)
        return 0;
    uint8_t type;
    if (!mpack_reader_ensure(reader, sizeof(type)))
        return 0;
    type = mpack_load_u8(reader->data);
    reader->data += sizeof(type);
    return type;
}

#if !MPACK_OPTIMIZE_FOR_SIZE
MPACK_STATIC_INLINE uint16_t mpack_expect_native_u16(mpack_reader_t* reader) {
    if (mpack_reader_error(reader) != mpack_ok)
        return 0;
    uint16_t type;
    if (!mpack_reader_ensure(reader, sizeof(type)))
        return 0;
    type = mpack_load_u16(reader->data);
    reader->data += sizeof(type);
    return type;
}

MPACK_STATIC_INLINE uint32_t mpack_expect_native_u32(mpack_reader_t* reader) {
    if (mpack_reader_error(reader) != mpack_ok)
        return 0;
    uint32_t type;
    if (!mpack_reader_ensure(reader, sizeof(type)))
        return 0;
    type = mpack_load_u32(reader->data);
    reader->data += sizeof(type);
    return type;
}
#endif

MPACK_STATIC_INLINE uint8_t mpack_expect_type_byte(mpack_reader_t* reader) {
    mpack_reader_track_element(reader);
    return mpack_expect_native_u8(reader);
}


// Basic Number Functions

uint8_t mpack_expect_u8(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_uint) {
        if (var.v.u <= MPACK_UINT8_MAX)
            return (uint8_t)var.v.u;
    } else if (var.type == mpack_type_int) {
        if (var.v.i >= 0 && var.v.i <= MPACK_UINT8_MAX)
            return (uint8_t)var.v.i;
    }
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0;
}

uint16_t mpack_expect_u16(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_uint) {
        if (var.v.u <= MPACK_UINT16_MAX)
            return (uint16_t)var.v.u;
    } else if (var.type == mpack_type_int) {
        if (var.v.i >= 0 && var.v.i <= MPACK_UINT16_MAX)
            return (uint16_t)var.v.i;
    }
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0;
}

uint32_t mpack_expect_u32(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_uint) {
        if (var.v.u <= MPACK_UINT32_MAX)
            return (uint32_t)var.v.u;
    } else if (var.type == mpack_type_int) {
        if (var.v.i >= 0 && var.v.i <= MPACK_UINT32_MAX)
            return (uint32_t)var.v.i;
    }
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0;
}

uint64_t mpack_expect_u64(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_uint) {
        return var.v.u;
    } else if (var.type == mpack_type_int) {
        if (var.v.i >= 0)
            return (uint64_t)var.v.i;
    }
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0;
}

int8_t mpack_expect_i8(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_uint) {
        if (var.v.u <= MPACK_INT8_MAX)
            return (int8_t)var.v.u;
    } else if (var.type == mpack_type_int) {
        if (var.v.i >= MPACK_INT8_MIN && var.v.i <= MPACK_INT8_MAX)
            return (int8_t)var.v.i;
    }
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0;
}

int16_t mpack_expect_i16(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_uint) {
        if (var.v.u <= MPACK_INT16_MAX)
            return (int16_t)var.v.u;
    } else if (var.type == mpack_type_int) {
        if (var.v.i >= MPACK_INT16_MIN && var.v.i <= MPACK_INT16_MAX)
            return (int16_t)var.v.i;
    }
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0;
}

int32_t mpack_expect_i32(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_uint) {
        if (var.v.u <= MPACK_INT32_MAX)
            return (int32_t)var.v.u;
    } else if (var.type == mpack_type_int) {
        if (var.v.i >= MPACK_INT32_MIN && var.v.i <= MPACK_INT32_MAX)
            return (int32_t)var.v.i;
    }
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0;
}

int64_t mpack_expect_i64(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_uint) {
        if (var.v.u <= MPACK_INT64_MAX)
            return (int64_t)var.v.u;
    } else if (var.type == mpack_type_int) {
        return var.v.i;
    }
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0;
}

#if MPACK_FLOAT
float mpack_expect_float(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_uint)
        return (float)var.v.u;
    if (var.type == mpack_type_int)
        return (float)var.v.i;
    if (var.type == mpack_type_float)
        return var.v.f;

    if (var.type == mpack_type_double) {
        #if MPACK_DOUBLE
        return (float)var.v.d;
        #else
        return mpack_shorten_raw_double_to_float(var.v.d);
        #endif
    }

    mpack_reader_flag_error(reader, mpack_error_type);
    return 0.0f;
}
#endif

#if MPACK_DOUBLE
double mpack_expect_double(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_uint)
        return (double)var.v.u;
    else if (var.type == mpack_type_int)
        return (double)var.v.i;
    else if (var.type == mpack_type_float)
        return (double)var.v.f;
    else if (var.type == mpack_type_double)
        return var.v.d;
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0.0;
}
#endif

#if MPACK_FLOAT
float mpack_expect_float_strict(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_float)
        return var.v.f;
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0.0f;
}
#endif

#if MPACK_DOUBLE
double mpack_expect_double_strict(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_float)
        return (double)var.v.f;
    else if (var.type == mpack_type_double)
        return var.v.d;
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0.0;
}
#endif

#if !MPACK_FLOAT
uint32_t mpack_expect_raw_float(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_float)
        return var.v.f;
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0;
}
#endif

#if !MPACK_DOUBLE
uint64_t mpack_expect_raw_double(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_double)
        return var.v.d;
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0;
}
#endif


// Ranged Number Functions
//
// All ranged functions are identical other than the type, so we
// define their content with a macro. The prototypes are still written
// out in full to support ctags/IDE tools.

#define MPACK_EXPECT_RANGE_IMPL(name, type_t)                           \
                                                                        \
    /* make sure the range is sensible */                               \
    mpack_assert(min_value <= max_value,                                \
            "min_value %i must be less than or equal to max_value %i",  \
            min_value, max_value);                                      \
                                                                        \
    /* read the value */                                                \
    type_t val = mpack_expect_##name(reader);                           \
    if (mpack_reader_error(reader) != mpack_ok)                         \
        return min_value;                                               \
                                                                        \
    /* make sure it fits */                                             \
    if (val < min_value || val > max_value) {                           \
        mpack_reader_flag_error(reader, mpack_error_type);              \
        return min_value;                                               \
    }                                                                   \
                                                                        \
    return val;

uint8_t mpack_expect_u8_range(mpack_reader_t* reader, uint8_t min_value, uint8_t max_value) {MPACK_EXPECT_RANGE_IMPL(u8, uint8_t)}
uint16_t mpack_expect_u16_range(mpack_reader_t* reader, uint16_t min_value, uint16_t max_value) {MPACK_EXPECT_RANGE_IMPL(u16, uint16_t)}
uint32_t mpack_expect_u32_range(mpack_reader_t* reader, uint32_t min_value, uint32_t max_value) {MPACK_EXPECT_RANGE_IMPL(u32, uint32_t)}
uint64_t mpack_expect_u64_range(mpack_reader_t* reader, uint64_t min_value, uint64_t max_value) {MPACK_EXPECT_RANGE_IMPL(u64, uint64_t)}

int8_t mpack_expect_i8_range(mpack_reader_t* reader, int8_t min_value, int8_t max_value) {MPACK_EXPECT_RANGE_IMPL(i8, int8_t)}
int16_t mpack_expect_i16_range(mpack_reader_t* reader, int16_t min_value, int16_t max_value) {MPACK_EXPECT_RANGE_IMPL(i16, int16_t)}
int32_t mpack_expect_i32_range(mpack_reader_t* reader, int32_t min_value, int32_t max_value) {MPACK_EXPECT_RANGE_IMPL(i32, int32_t)}
int64_t mpack_expect_i64_range(mpack_reader_t* reader, int64_t min_value, int64_t max_value) {MPACK_EXPECT_RANGE_IMPL(i64, int64_t)}

#if MPACK_FLOAT
float mpack_expect_float_range(mpack_reader_t* reader, float min_value, float max_value) {MPACK_EXPECT_RANGE_IMPL(float, float)}
#endif
#if MPACK_DOUBLE
double mpack_expect_double_range(mpack_reader_t* reader, double min_value, double max_value) {MPACK_EXPECT_RANGE_IMPL(double, double)}
#endif

uint32_t mpack_expect_map_range(mpack_reader_t* reader, uint32_t min_value, uint32_t max_value) {MPACK_EXPECT_RANGE_IMPL(map, uint32_t)}
uint32_t mpack_expect_array_range(mpack_reader_t* reader, uint32_t min_value, uint32_t max_value) {MPACK_EXPECT_RANGE_IMPL(array, uint32_t)}


// Matching Number Functions

void mpack_expect_uint_match(mpack_reader_t* reader, uint64_t value) {
    if (mpack_expect_u64(reader) != value)
        mpack_reader_flag_error(reader, mpack_error_type);
}

void mpack_expect_int_match(mpack_reader_t* reader, int64_t value) {
    if (mpack_expect_i64(reader) != value)
        mpack_reader_flag_error(reader, mpack_error_type);
}


// Other Basic Types

void mpack_expect_nil(mpack_reader_t* reader) {
    if (mpack_expect_type_byte(reader) != 0xc0)
        mpack_reader_flag_error(reader, mpack_error_type);
}

bool mpack_expect_bool(mpack_reader_t* reader) {
    uint8_t type = mpack_expect_type_byte(reader);
    if ((type & ~1) != 0xc2)
        mpack_reader_flag_error(reader, mpack_error_type);
    return (bool)(type & 1);
}

void mpack_expect_true(mpack_reader_t* reader) {
    if (mpack_expect_bool(reader) != true)
        mpack_reader_flag_error(reader, mpack_error_type);
}

void mpack_expect_false(mpack_reader_t* reader) {
    if (mpack_expect_bool(reader) != false)
        mpack_reader_flag_error(reader, mpack_error_type);
}

#if MPACK_EXTENSIONS
mpack_timestamp_t mpack_expect_timestamp(mpack_reader_t* reader) {
    mpack_timestamp_t zero = {0, 0};

    mpack_tag_t tag = mpack_read_tag(reader);
    if (tag.type != mpack_type_ext) {
        mpack_reader_flag_error(reader, mpack_error_type);
        return zero;
    }
    if (mpack_tag_ext_exttype(&tag) != MPACK_EXTTYPE_TIMESTAMP) {
        mpack_reader_flag_error(reader, mpack_error_type);
        return zero;
    }

    return mpack_read_timestamp(reader, mpack_tag_ext_length(&tag));
}

int64_t mpack_expect_timestamp_truncate(mpack_reader_t* reader) {
    return mpack_expect_timestamp(reader).seconds;
}
#endif


// Compound Types

uint32_t mpack_expect_map(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_map)
        return var.v.n;
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0;
}

void mpack_expect_map_match(mpack_reader_t* reader, uint32_t count) {
    if (mpack_expect_map(reader) != count)
        mpack_reader_flag_error(reader, mpack_error_type);
}

bool mpack_expect_map_or_nil(mpack_reader_t* reader, uint32_t* count) {
    mpack_assert(count != NULL, "count cannot be NULL");

    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_nil) {
        *count = 0;
        return false;
    }
    if (var.type == mpack_type_map) {
        *count = var.v.n;
        return true;
    }
    mpack_reader_flag_error(reader, mpack_error_type);
    *count = 0;
    return false;
}

bool mpack_expect_map_max_or_nil(mpack_reader_t* reader, uint32_t max_count, uint32_t* count) {
    mpack_assert(count != NULL, "count cannot be NULL");

    bool has_map = mpack_expect_map_or_nil(reader, count);
    if (has_map && *count > max_count) {
        *count = 0;
        mpack_reader_flag_error(reader, mpack_error_type);
        return false;
    }
    return has_map;
}

uint32_t mpack_expect_array(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_array)
        return var.v.n;
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0;
}

void mpack_expect_array_match(mpack_reader_t* reader, uint32_t count) {
    if (mpack_expect_array(reader) != count)
        mpack_reader_flag_error(reader, mpack_error_type);
}

bool mpack_expect_array_or_nil(mpack_reader_t* reader, uint32_t* count) {
    mpack_assert(count != NULL, "count cannot be NULL");

    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_nil) {
        *count = 0;
        return false;
    }
    if (var.type == mpack_type_array) {
        *count = var.v.n;
        return true;
    }
    mpack_reader_flag_error(reader, mpack_error_type);
    *count = 0;
    return false;
}

bool mpack_expect_array_max_or_nil(mpack_reader_t* reader, uint32_t max_count, uint32_t* count) {
    mpack_assert(count != NULL, "count cannot be NULL");

    bool has_array = mpack_expect_array_or_nil(reader, count);
    if (has_array && *count > max_count) {
        *count = 0;
        mpack_reader_flag_error(reader, mpack_error_type);
        return false;
    }
    return has_array;
}

#ifdef MPACK_MALLOC
void* mpack_expect_array_alloc_impl(mpack_reader_t* reader, size_t element_size, uint32_t max_count, uint32_t* out_count, bool allow_nil) {
    mpack_assert(out_count != NULL, "out_count cannot be NULL");
    *out_count = 0;

    uint32_t count;
    bool has_array = true;
    if (allow_nil)
        has_array = mpack_expect_array_max_or_nil(reader, max_count, &count);
    else
        count = mpack_expect_array_max(reader, max_count);
    if (mpack_reader_error(reader))
        return NULL;

    // size 0 is not an error; we return NULL for no elements.
    if (count == 0) {
        // we call mpack_done_array() automatically ONLY if we are using
        // the _or_nil variant. this is the only way to allow nil and empty
        // to work the same way.
        if (allow_nil && has_array)
            mpack_done_array(reader);
        return NULL;
    }

    void* p = MPACK_MALLOC(element_size * count);
    if (p == NULL) {
        mpack_reader_flag_error(reader, mpack_error_memory);
        return NULL;
    }

    *out_count = count;
    return p;
}
#endif


// Str, Bin and Ext Functions

uint32_t mpack_expect_str(mpack_reader_t* reader) {
    #if MPACK_OPTIMIZE_FOR_SIZE
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_str)
        return var.v.l;
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0;
    #else
    uint8_t type = mpack_expect_type_byte(reader);
    uint32_t count;

    if ((type >> 5) == 5) {
        count = type & (uint8_t)~0xe0;
    } else if (type == 0xd9) {
        count = mpack_expect_native_u8(reader);
    } else if (type == 0xda) {
        count = mpack_expect_native_u16(reader);
    } else if (type == 0xdb) {
        count = mpack_expect_native_u32(reader);
    } else {
        mpack_reader_flag_error(reader, mpack_error_type);
        return 0;
    }

    #if MPACK_READ_TRACKING
    mpack_reader_flag_if_error(reader, mpack_track_push(&reader->track, mpack_type_str, count));
    #endif
    return count;
    #endif
}

size_t mpack_expect_str_buf(mpack_reader_t* reader, char* buf, size_t bufsize) {
    mpack_assert(buf != NULL, "buf cannot be NULL");

    size_t length = mpack_expect_str(reader);
    if (mpack_reader_error(reader))
        return 0;

    if (length > bufsize) {
        mpack_reader_flag_error(reader, mpack_error_too_big);
        return 0;
    }

    mpack_read_bytes(reader, buf, length);
    if (mpack_reader_error(reader))
        return 0;

    mpack_done_str(reader);
    return length;
}

size_t mpack_expect_utf8(mpack_reader_t* reader, char* buf, size_t size) {
    mpack_assert(buf != NULL, "buf cannot be NULL");

    size_t length = mpack_expect_str_buf(reader, buf, size);

    if (!mpack_utf8_check(buf, length)) {
        mpack_reader_flag_error(reader, mpack_error_type);
        return 0;
    }

    return length;
}

uint32_t mpack_expect_bin(mpack_reader_t* reader) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_bin)
        return var.v.l;
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0;
}

size_t mpack_expect_bin_buf(mpack_reader_t* reader, char* buf, size_t bufsize) {
    mpack_assert(buf != NULL, "buf cannot be NULL");

    size_t binsize = mpack_expect_bin(reader);
    if (mpack_reader_error(reader))
        return 0;
    if (binsize > bufsize) {
        mpack_reader_flag_error(reader, mpack_error_too_big);
        return 0;
    }
    mpack_read_bytes(reader, buf, binsize);
    if (mpack_reader_error(reader))
        return 0;
    mpack_done_bin(reader);
    return binsize;
}

void mpack_expect_bin_size_buf(mpack_reader_t* reader, char* buf, uint32_t size) {
    mpack_assert(buf != NULL, "buf cannot be NULL");
    mpack_expect_bin_size(reader, size);
    mpack_read_bytes(reader, buf, size);
    mpack_done_bin(reader);
}

#if MPACK_EXTENSIONS
uint32_t mpack_expect_ext(mpack_reader_t* reader, int8_t* type) {
    mpack_tag_t var = mpack_read_tag(reader);
    if (var.type == mpack_type_ext) {
        *type = mpack_tag_ext_exttype(&var);
        return mpack_tag_ext_length(&var);
    }
    *type = 0;
    mpack_reader_flag_error(reader, mpack_error_type);
    return 0;
}

size_t mpack_expect_ext_buf(mpack_reader_t* reader, int8_t* type, char* buf, size_t bufsize) {
    mpack_assert(buf != NULL, "buf cannot be NULL");

    size_t extsize = mpack_expect_ext(reader, type);
    if (mpack_reader_error(reader))
        return 0;
    if (extsize > bufsize) {
        *type = 0;
        mpack_reader_flag_error(reader, mpack_error_too_big);
        return 0;
    }
    mpack_read_bytes(reader, buf, extsize);
    if (mpack_reader_error(reader)) {
        *type = 0;
        return 0;
    }
    mpack_done_ext(reader);
    return extsize;
}
#endif

void mpack_expect_cstr(mpack_reader_t* reader, char* buf, size_t bufsize) {
    uint32_t length = mpack_expect_str(reader);
    mpack_read_cstr(reader, buf, bufsize, length);
    mpack_done_str(reader);
}

void mpack_expect_utf8_cstr(mpack_reader_t* reader, char* buf, size_t bufsize) {
    uint32_t length = mpack_expect_str(reader);
    mpack_read_utf8_cstr(reader, buf, bufsize, length);
    mpack_done_str(reader);
}

#ifdef MPACK_MALLOC
static char* mpack_expect_cstr_alloc_unchecked(mpack_reader_t* reader, size_t maxsize, size_t* out_length) {
    mpack_assert(out_length != NULL, "out_length cannot be NULL");
    *out_length = 0;

    // make sure argument makes sense
    if (maxsize < 1) {
        mpack_break("maxsize is zero; you must have room for at least a null-terminator");
        mpack_reader_flag_error(reader, mpack_error_bug);
        return NULL;
    }

    if (SIZE_MAX < MPACK_UINT32_MAX) {
        if (maxsize > SIZE_MAX)
            maxsize = SIZE_MAX;
    } else {
        if (maxsize > (size_t)MPACK_UINT32_MAX)
            maxsize = (size_t)MPACK_UINT32_MAX;
    }

    size_t length = mpack_expect_str_max(reader, (uint32_t)maxsize - 1);
    char* str = mpack_read_bytes_alloc_impl(reader, length, true);
    mpack_done_str(reader);

    if (str)
        *out_length = length;
    return str;
}

char* mpack_expect_cstr_alloc(mpack_reader_t* reader, size_t maxsize) {
    size_t length;
    char* str = mpack_expect_cstr_alloc_unchecked(reader, maxsize, &length);

    if (str && !mpack_str_check_no_null(str, length)) {
        MPACK_FREE(str);
        mpack_reader_flag_error(reader, mpack_error_type);
        return NULL;
    }

    return str;
}

char* mpack_expect_utf8_cstr_alloc(mpack_reader_t* reader, size_t maxsize) {
    size_t length;
    char* str = mpack_expect_cstr_alloc_unchecked(reader, maxsize, &length);

    if (str && !mpack_utf8_check_no_null(str, length)) {
        MPACK_FREE(str);
        mpack_reader_flag_error(reader, mpack_error_type);
        return NULL;
    }

    return str;
}
#endif

void mpack_expect_str_match(mpack_reader_t* reader, const char* str, size_t len) {
    mpack_assert(str != NULL, "str cannot be NULL");

    // expect a str the correct length
    if (len > MPACK_UINT32_MAX)
        mpack_reader_flag_error(reader, mpack_error_type);
    mpack_expect_str_length(reader, (uint32_t)len);
    if (mpack_reader_error(reader))
        return;
    mpack_reader_track_bytes(reader, (uint32_t)len);

    // check each byte one by one (matched strings are likely to be very small)
    for (; len > 0; --len) {
        if (mpack_expect_native_u8(reader) != *str++) {
            mpack_reader_flag_error(reader, mpack_error_type);
            return;
        }
    }

    mpack_done_str(reader);
}

void mpack_expect_tag(mpack_reader_t* reader, mpack_tag_t expected) {
    mpack_tag_t actual = mpack_read_tag(reader);
    if (!mpack_tag_equal(actual, expected))
        mpack_reader_flag_error(reader, mpack_error_type);
}

#ifdef MPACK_MALLOC
char* mpack_expect_bin_alloc(mpack_reader_t* reader, size_t maxsize, size_t* size) {
    mpack_assert(size != NULL, "size cannot be NULL");
    *size = 0;

    if (SIZE_MAX < MPACK_UINT32_MAX) {
        if (maxsize > SIZE_MAX)
            maxsize = SIZE_MAX;
    } else {
        if (maxsize > (size_t)MPACK_UINT32_MAX)
            maxsize = (size_t)MPACK_UINT32_MAX;
    }

    size_t length = mpack_expect_bin_max(reader, (uint32_t)maxsize);
    if (mpack_reader_error(reader))
        return NULL;

    char* data = mpack_read_bytes_alloc(reader, length);
    mpack_done_bin(reader);

    if (data)
        *size = length;
    return data;
}
#endif

#if MPACK_EXTENSIONS && defined(MPACK_MALLOC)
char* mpack_expect_ext_alloc(mpack_reader_t* reader, int8_t* type, size_t maxsize, size_t* size) {
    mpack_assert(size != NULL, "size cannot be NULL");
    *size = 0;

    if (SIZE_MAX < MPACK_UINT32_MAX) {
        if (maxsize > SIZE_MAX)
            maxsize = SIZE_MAX;
    } else {
        if (maxsize > (size_t)MPACK_UINT32_MAX)
            maxsize = (size_t)MPACK_UINT32_MAX;
    }

    size_t length = mpack_expect_ext_max(reader, type, (uint32_t)maxsize);
    if (mpack_reader_error(reader))
        return NULL;

    char* data = mpack_read_bytes_alloc(reader, length);
    mpack_done_ext(reader);

    if (data) {
        *size = length;
    } else {
        *type = 0;
    }
    return data;
}
#endif

size_t mpack_expect_enum(mpack_reader_t* reader, const char* strings[], size_t count) {

    // read the string in-place
    size_t keylen = mpack_expect_str(reader);
    const char* key = mpack_read_bytes_inplace(reader, keylen);
    mpack_done_str(reader);
    if (mpack_reader_error(reader) != mpack_ok)
        return count;

    // find what key it matches
    size_t i;
    for (i = 0; i < count; ++i) {
        const char* other = strings[i];
        size_t otherlen = mpack_strlen(other);
        if (keylen == otherlen && mpack_memcmp(key, other, keylen) == 0)
            return i;
    }

    // no matches
    mpack_reader_flag_error(reader, mpack_error_type);
    return count;
}

size_t mpack_expect_enum_optional(mpack_reader_t* reader, const char* strings[], size_t count) {
    if (mpack_reader_error(reader) != mpack_ok)
        return count;

    mpack_assert(count != 0, "count cannot be zero; no strings are valid!");
    mpack_assert(strings != NULL, "strings cannot be NULL");

    // the key is only recognized if it is a string
    if (mpack_peek_tag(reader).type != mpack_type_str) {
        mpack_discard(reader);
        return count;
    }

    // read the string in-place
    size_t keylen = mpack_expect_str(reader);
    const char* key = mpack_read_bytes_inplace(reader, keylen);
    mpack_done_str(reader);
    if (mpack_reader_error(reader) != mpack_ok)
        return count;

    // find what key it matches
    size_t i;
    for (i = 0; i < count; ++i) {
        const char* other = strings[i];
        size_t otherlen = mpack_strlen(other);
        if (keylen == otherlen && mpack_memcmp(key, other, keylen) == 0)
            return i;
    }

    // no matches
    return count;
}

size_t mpack_expect_key_uint(mpack_reader_t* reader, bool found[], size_t count) {
    if (mpack_reader_error(reader) != mpack_ok)
        return count;

    if (count == 0) {
        mpack_break("count cannot be zero; no keys are valid!");
        mpack_reader_flag_error(reader, mpack_error_bug);
        return count;
    }
    mpack_assert(found != NULL, "found cannot be NULL");

    // the key is only recognized if it is an unsigned int
    if (mpack_peek_tag(reader).type != mpack_type_uint) {
        mpack_discard(reader);
        return count;
    }

    // read the key
    uint64_t value = mpack_expect_u64(reader);
    if (mpack_reader_error(reader) != mpack_ok)
        return count;

    // unrecognized keys are fine, we just return count
    if (value >= count)
        return count;

    // check if this key is a duplicate
    if (found[value]) {
        mpack_reader_flag_error(reader, mpack_error_invalid);
        return count;
    }

    found[value] = true;
    return (size_t)value;
}

size_t mpack_expect_key_cstr(mpack_reader_t* reader, const char* keys[], bool found[], size_t count) {
    size_t i = mpack_expect_enum_optional(reader, keys, count);

    // unrecognized keys are fine, we just return count
    if (i == count)
        return count;

    // check if this key is a duplicate
    mpack_assert(found != NULL, "found cannot be NULL");
    if (found[i]) {
        mpack_reader_flag_error(reader, mpack_error_invalid);
        return count;
    }

    found[i] = true;
    return i;
}

#endif

MPACK_SILENCE_WARNINGS_END

/* mpack/mpack-node.c.c */

#define MPACK_INTERNAL 1

/* #include "mpack-node.h" */

MPACK_SILENCE_WARNINGS_BEGIN

#if MPACK_NODE

MPACK_STATIC_INLINE const char* mpack_node_data_unchecked(mpack_node_t node) {
    mpack_assert(mpack_node_error(node) == mpack_ok, "tree is in an error state!");

    mpack_type_t type = node.data->type;
    MPACK_UNUSED(type);
    #if MPACK_EXTENSIONS
    mpack_assert(type == mpack_type_str || type == mpack_type_bin || type == mpack_type_ext,
            "node of type %i (%s) is not a data type!", type, mpack_type_to_string(type));
    #else
    mpack_assert(type == mpack_type_str || type == mpack_type_bin,
            "node of type %i (%s) is not a data type!", type, mpack_type_to_string(type));
    #endif

    return node.tree->data + node.data->value.offset;
}

#if MPACK_EXTENSIONS
MPACK_STATIC_INLINE int8_t mpack_node_exttype_unchecked(mpack_node_t node) {
    mpack_assert(mpack_node_error(node) == mpack_ok, "tree is in an error state!");

    mpack_type_t type = node.data->type;
    MPACK_UNUSED(type);
    mpack_assert(type == mpack_type_ext, "node of type %i (%s) is not an ext type!",
            type, mpack_type_to_string(type));

    // the exttype of an ext node is stored in the byte preceding the data
    return mpack_load_i8(mpack_node_data_unchecked(node) - 1);
}
#endif



/*
 * Tree Parsing
 */

#ifdef MPACK_MALLOC

// fix up the alloc size to make sure it exactly fits the
// maximum number of nodes it can contain (the allocator will
// waste it back anyway, but we round it down just in case)

#define MPACK_NODES_PER_PAGE \
    ((MPACK_NODE_PAGE_SIZE - sizeof(mpack_tree_page_t)) / sizeof(mpack_node_data_t) + 1)

#define MPACK_PAGE_ALLOC_SIZE \
    (sizeof(mpack_tree_page_t) + sizeof(mpack_node_data_t) * (MPACK_NODES_PER_PAGE - 1))

#endif

#ifdef MPACK_MALLOC
/*
 * Fills the tree until we have at least enough bytes for the current node.
 */
static bool mpack_tree_reserve_fill(mpack_tree_t* tree) {
    mpack_assert(tree->parser.state == mpack_tree_parse_state_in_progress);

    size_t bytes = tree->parser.current_node_reserved;
    mpack_assert(bytes > tree->parser.possible_nodes_left,
            "there are already enough bytes! call mpack_tree_ensure() instead.");
    mpack_log("filling to reserve %i bytes\n", (int)bytes);

    // if the necessary bytes would put us over the maximum tree
    // size, fail right away.
    // TODO: check for overflow?
    if (tree->data_length + bytes > tree->max_size) {
        mpack_tree_flag_error(tree, mpack_error_too_big);
        return false;
    }

    // we'll need a read function to fetch more data. if there's
    // no read function, the data should contain an entire message
    // (or messages), so we flag it as invalid.
    if (tree->read_fn == NULL) {
        mpack_log("tree has no read function!\n");
        mpack_tree_flag_error(tree, mpack_error_invalid);
        return false;
    }

    // expand the buffer if needed
    if (tree->data_length + bytes > tree->buffer_capacity) {

        // TODO: check for overflow?
        size_t new_capacity = (tree->buffer_capacity == 0) ? MPACK_BUFFER_SIZE : tree->buffer_capacity;
        while (new_capacity < tree->data_length + bytes)
            new_capacity *= 2;
        if (new_capacity > tree->max_size)
            new_capacity = tree->max_size;

        mpack_log("expanding buffer from %i to %i\n", (int)tree->buffer_capacity, (int)new_capacity);

        char* new_buffer;
        if (tree->buffer == NULL)
            new_buffer = (char*)MPACK_MALLOC(new_capacity);
        else
            new_buffer = (char*)mpack_realloc(tree->buffer, tree->data_length, new_capacity);

        if (new_buffer == NULL) {
            mpack_tree_flag_error(tree, mpack_error_memory);
            return false;
        }

        tree->data = new_buffer;
        tree->buffer = new_buffer;
        tree->buffer_capacity = new_capacity;
    }

    // request as much data as possible, looping until we have
    // all the data we need
    do {
        size_t read = tree->read_fn(tree, tree->buffer + tree->data_length, tree->buffer_capacity - tree->data_length);

        // If the fill function encounters an error, it should flag an error on
        // the tree.
        if (mpack_tree_error(tree) != mpack_ok)
            return false;

        // We guard against fill functions that return -1 just in case.
        if (read == (size_t)(-1)) {
            mpack_tree_flag_error(tree, mpack_error_io);
            return false;
        }

        // If the fill function returns 0, the data is not available yet. We
        // return false to stop parsing the current node.
        if (read == 0) {
            mpack_log("not enough data.\n");
            return false;
        }

        mpack_log("read %" PRIu32 " more bytes\n", (uint32_t)read);
        tree->data_length += read;
        tree->parser.possible_nodes_left += read;
    } while (tree->parser.possible_nodes_left < bytes);

    return true;
}
#endif

/*
 * Ensures there are enough additional bytes in the tree for the current node
 * (including reserved bytes for the children of this node, and in addition to
 * the reserved bytes for children of previous compound nodes), reading more
 * data if needed.
 *
 * extra_bytes is the number of additional bytes to reserve for the current
 * node beyond the type byte (since one byte is already reserved for each node
 * by its parent array or map.)
 *
 * This may reallocate the tree, which means the tree->data pointer may change!
 *
 * Returns false if not enough bytes could be read.
 */
MPACK_STATIC_INLINE bool mpack_tree_reserve_bytes(mpack_tree_t* tree, size_t extra_bytes) {
    mpack_assert(tree->parser.state == mpack_tree_parse_state_in_progress);

    // We guard against overflow here. A compound type could declare more than
    // MPACK_UINT32_MAX contents which overflows SIZE_MAX on 32-bit platforms. We
    // flag mpack_error_invalid instead of mpack_error_too_big since it's far
    // more likely that the message is corrupt than that the data is valid but
    // not parseable on this architecture (see test_read_node_possible() in
    // test-node.c .)
    if ((uint64_t)tree->parser.current_node_reserved + (uint64_t)extra_bytes > SIZE_MAX) {
        mpack_tree_flag_error(tree, mpack_error_invalid);
        return false;
    }

    tree->parser.current_node_reserved += extra_bytes;

    // Note that possible_nodes_left already accounts for reserved bytes for
    // children of previous compound nodes. So even if there are hundreds of
    // bytes left in the buffer, we might need to read anyway.
    if (tree->parser.current_node_reserved <= tree->parser.possible_nodes_left)
        return true;

    #ifdef MPACK_MALLOC
    return mpack_tree_reserve_fill(tree);
    #else
    return false;
    #endif
}

MPACK_STATIC_INLINE size_t mpack_tree_parser_stack_capacity(mpack_tree_t* tree) {
    #ifdef MPACK_MALLOC
    return tree->parser.stack_capacity;
    #else
    return sizeof(tree->parser.stack) / sizeof(tree->parser.stack[0]);
    #endif
}

static bool mpack_tree_push_stack(mpack_tree_t* tree, mpack_node_data_t* first_child, size_t total) {
    mpack_tree_parser_t* parser = &tree->parser;
    mpack_assert(parser->state == mpack_tree_parse_state_in_progress);

    // No need to push empty containers
    if (total == 0)
        return true;

    // Make sure we have enough room in the stack
    if (parser->level + 1 == mpack_tree_parser_stack_capacity(tree)) {
        #ifdef MPACK_MALLOC
        size_t new_capacity = parser->stack_capacity * 2;
        mpack_log("growing parse stack to capacity %i\n", (int)new_capacity);

        // Replace the stack-allocated parsing stack
        if (!parser->stack_owned) {
            mpack_level_t* new_stack = (mpack_level_t*)MPACK_MALLOC(sizeof(mpack_level_t) * new_capacity);
            if (!new_stack) {
                mpack_tree_flag_error(tree, mpack_error_memory);
                return false;
            }
            mpack_memcpy(new_stack, parser->stack, sizeof(mpack_level_t) * parser->stack_capacity);
            parser->stack = new_stack;
            parser->stack_owned = true;

        // Realloc the allocated parsing stack
        } else {
            mpack_level_t* new_stack = (mpack_level_t*)mpack_realloc(parser->stack,
                    sizeof(mpack_level_t) * parser->stack_capacity, sizeof(mpack_level_t) * new_capacity);
            if (!new_stack) {
                mpack_tree_flag_error(tree, mpack_error_memory);
                return false;
            }
            parser->stack = new_stack;
        }
        parser->stack_capacity = new_capacity;
        #else
        mpack_tree_flag_error(tree, mpack_error_too_big);
        return false;
        #endif
    }

    // Push the contents of this node onto the parsing stack
    ++parser->level;
    parser->stack[parser->level].child = first_child;
    parser->stack[parser->level].left = total;
    return true;
}

static bool mpack_tree_parse_children(mpack_tree_t* tree, mpack_node_data_t* node) {
    mpack_tree_parser_t* parser = &tree->parser;
    mpack_assert(parser->state == mpack_tree_parse_state_in_progress);

    mpack_type_t type = node->type;
    size_t total = node->len;

    // Calculate total elements to read
    if (type == mpack_type_map) {
        if ((uint64_t)total * 2 > SIZE_MAX) {
            mpack_tree_flag_error(tree, mpack_error_too_big);
            return false;
        }
        total *= 2;
    }

    // Make sure we are under our total node limit (TODO can this overflow?)
    tree->node_count += total;
    if (tree->node_count > tree->max_nodes) {
        mpack_tree_flag_error(tree, mpack_error_too_big);
        return false;
    }

    // Each node is at least one byte. Count these bytes now to make
    // sure there is enough data left.
    if (!mpack_tree_reserve_bytes(tree, total))
        return false;

    // If there are enough nodes left in the current page, no need to grow
    if (total <= parser->nodes_left) {
        node->value.children = parser->nodes;
        parser->nodes += total;
        parser->nodes_left -= total;

    } else {

        #ifdef MPACK_MALLOC

        // We can't grow if we're using a fixed pool (i.e. we didn't start with a page)
        if (!tree->next) {
            mpack_tree_flag_error(tree, mpack_error_too_big);
            return false;
        }

        // Otherwise we need to grow, and the node's children need to be contiguous.
        // This is a heuristic to decide whether we should waste the remaining space
        // in the current page and start a new one, or give the children their
        // own page. With a fraction of 1/8, this causes at most 12% additional
        // waste. Note that reducing this too much causes less cache coherence and
        // more malloc() overhead due to smaller allocations, so there's a tradeoff
        // here. This heuristic could use some improvement, especially with custom
        // page sizes.

        mpack_tree_page_t* page;

        if (total > MPACK_NODES_PER_PAGE || parser->nodes_left > MPACK_NODES_PER_PAGE / 8) {
            // TODO: this should check for overflow
            page = (mpack_tree_page_t*)MPACK_MALLOC(
                    sizeof(mpack_tree_page_t) + sizeof(mpack_node_data_t) * (total - 1));
            if (page == NULL) {
                mpack_tree_flag_error(tree, mpack_error_memory);
                return false;
            }
            mpack_log("allocated seperate page %p for %i children, %i left in page of %i total\n",
                    (void*)page, (int)total, (int)parser->nodes_left, (int)MPACK_NODES_PER_PAGE);

            node->value.children = page->nodes;

        } else {
            page = (mpack_tree_page_t*)MPACK_MALLOC(MPACK_PAGE_ALLOC_SIZE);
            if (page == NULL) {
                mpack_tree_flag_error(tree, mpack_error_memory);
                return false;
            }
            mpack_log("allocated new page %p for %i children, wasting %i in page of %i total\n",
                    (void*)page, (int)total, (int)parser->nodes_left, (int)MPACK_NODES_PER_PAGE);

            node->value.children = page->nodes;
            parser->nodes = page->nodes + total;
            parser->nodes_left = MPACK_NODES_PER_PAGE - total;
        }

        page->next = tree->next;
        tree->next = page;

        #else
        // We can't grow if we don't have an allocator
        mpack_tree_flag_error(tree, mpack_error_too_big);
        return false;
        #endif
    }

    return mpack_tree_push_stack(tree, node->value.children, total);
}

static bool mpack_tree_parse_bytes(mpack_tree_t* tree, mpack_node_data_t* node) {
    node->value.offset = tree->size + tree->parser.current_node_reserved + 1;
    return mpack_tree_reserve_bytes(tree, node->len);
}

#if MPACK_EXTENSIONS
static bool mpack_tree_parse_ext(mpack_tree_t* tree, mpack_node_data_t* node) {
    // reserve space for exttype
    tree->parser.current_node_reserved += sizeof(int8_t);
    node->type = mpack_type_ext;
    return mpack_tree_parse_bytes(tree, node);
}
#endif

static bool mpack_tree_parse_node_contents(mpack_tree_t* tree, mpack_node_data_t* node) {
    mpack_assert(tree->parser.state == mpack_tree_parse_state_in_progress);
    mpack_assert(node != NULL, "null node?");

    // read the type. we've already accounted for this byte in
    // possible_nodes_left, so we already know it is in bounds, and we don't
    // need to reserve it for this node.
    mpack_assert(tree->data_length > tree->size);
    uint8_t type = mpack_load_u8(tree->data + tree->size);
    mpack_log("node type %x\n", type);
    tree->parser.current_node_reserved = 0;

    // as with mpack_read_tag(), the fastest way to parse a node is to switch
    // on the first byte, and to explicitly list every possible byte. we switch
    // on the first four bits in size-optimized builds.

    #if MPACK_OPTIMIZE_FOR_SIZE
    switch (type >> 4) {

        // positive fixnum
        case 0x0: case 0x1: case 0x2: case 0x3:
        case 0x4: case 0x5: case 0x6: case 0x7:
            node->type = mpack_type_uint;
            node->value.u = type;
            return true;

        // negative fixnum
        case 0xe: case 0xf:
            node->type = mpack_type_int;
            node->value.i = (int8_t)type;
            return true;

        // fixmap
        case 0x8:
            node->type = mpack_type_map;
            node->len = (uint32_t)(type & ~0xf0);
            return mpack_tree_parse_children(tree, node);

        // fixarray
        case 0x9:
            node->type = mpack_type_array;
            node->len = (uint32_t)(type & ~0xf0);
            return mpack_tree_parse_children(tree, node);

        // fixstr
        case 0xa: case 0xb:
            node->type = mpack_type_str;
            node->len = (uint32_t)(type & ~0xe0);
            return mpack_tree_parse_bytes(tree, node);

        // not one of the common infix types
        default:
            break;
    }
    #endif

    switch (type) {

        #if !MPACK_OPTIMIZE_FOR_SIZE
        // positive fixnum
        case 0x00: case 0x01: case 0x02: case 0x03: case 0x04: case 0x05: case 0x06: case 0x07:
        case 0x08: case 0x09: case 0x0a: case 0x0b: case 0x0c: case 0x0d: case 0x0e: case 0x0f:
        case 0x10: case 0x11: case 0x12: case 0x13: case 0x14: case 0x15: case 0x16: case 0x17:
        case 0x18: case 0x19: case 0x1a: case 0x1b: case 0x1c: case 0x1d: case 0x1e: case 0x1f:
        case 0x20: case 0x21: case 0x22: case 0x23: case 0x24: case 0x25: case 0x26: case 0x27:
        case 0x28: case 0x29: case 0x2a: case 0x2b: case 0x2c: case 0x2d: case 0x2e: case 0x2f:
        case 0x30: case 0x31: case 0x32: case 0x33: case 0x34: case 0x35: case 0x36: case 0x37:
        case 0x38: case 0x39: case 0x3a: case 0x3b: case 0x3c: case 0x3d: case 0x3e: case 0x3f:
        case 0x40: case 0x41: case 0x42: case 0x43: case 0x44: case 0x45: case 0x46: case 0x47:
        case 0x48: case 0x49: case 0x4a: case 0x4b: case 0x4c: case 0x4d: case 0x4e: case 0x4f:
        case 0x50: case 0x51: case 0x52: case 0x53: case 0x54: case 0x55: case 0x56: case 0x57:
        case 0x58: case 0x59: case 0x5a: case 0x5b: case 0x5c: case 0x5d: case 0x5e: case 0x5f:
        case 0x60: case 0x61: case 0x62: case 0x63: case 0x64: case 0x65: case 0x66: case 0x67:
        case 0x68: case 0x69: case 0x6a: case 0x6b: case 0x6c: case 0x6d: case 0x6e: case 0x6f:
        case 0x70: case 0x71: case 0x72: case 0x73: case 0x74: case 0x75: case 0x76: case 0x77:
        case 0x78: case 0x79: case 0x7a: case 0x7b: case 0x7c: case 0x7d: case 0x7e: case 0x7f:
            node->type = mpack_type_uint;
            node->value.u = type;
            return true;

        // negative fixnum
        case 0xe0: case 0xe1: case 0xe2: case 0xe3: case 0xe4: case 0xe5: case 0xe6: case 0xe7:
        case 0xe8: case 0xe9: case 0xea: case 0xeb: case 0xec: case 0xed: case 0xee: case 0xef:
        case 0xf0: case 0xf1: case 0xf2: case 0xf3: case 0xf4: case 0xf5: case 0xf6: case 0xf7:
        case 0xf8: case 0xf9: case 0xfa: case 0xfb: case 0xfc: case 0xfd: case 0xfe: case 0xff:
            node->type = mpack_type_int;
            node->value.i = (int8_t)type;
            return true;

        // fixmap
        case 0x80: case 0x81: case 0x82: case 0x83: case 0x84: case 0x85: case 0x86: case 0x87:
        case 0x88: case 0x89: case 0x8a: case 0x8b: case 0x8c: case 0x8d: case 0x8e: case 0x8f:
            node->type = mpack_type_map;
            node->len = (uint32_t)(type & ~0xf0);
            return mpack_tree_parse_children(tree, node);

        // fixarray
        case 0x90: case 0x91: case 0x92: case 0x93: case 0x94: case 0x95: case 0x96: case 0x97:
        case 0x98: case 0x99: case 0x9a: case 0x9b: case 0x9c: case 0x9d: case 0x9e: case 0x9f:
            node->type = mpack_type_array;
            node->len = (uint32_t)(type & ~0xf0);
            return mpack_tree_parse_children(tree, node);

        // fixstr
        case 0xa0: case 0xa1: case 0xa2: case 0xa3: case 0xa4: case 0xa5: case 0xa6: case 0xa7:
        case 0xa8: case 0xa9: case 0xaa: case 0xab: case 0xac: case 0xad: case 0xae: case 0xaf:
        case 0xb0: case 0xb1: case 0xb2: case 0xb3: case 0xb4: case 0xb5: case 0xb6: case 0xb7:
        case 0xb8: case 0xb9: case 0xba: case 0xbb: case 0xbc: case 0xbd: case 0xbe: case 0xbf:
            node->type = mpack_type_str;
            node->len = (uint32_t)(type & ~0xe0);
            return mpack_tree_parse_bytes(tree, node);
        #endif

        // nil
        case 0xc0:
            node->type = mpack_type_nil;
            return true;

        // bool
        case 0xc2: case 0xc3:
            node->type = mpack_type_bool;
            node->value.b = type & 1;
            return true;

        // bin8
        case 0xc4:
            node->type = mpack_type_bin;
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint8_t)))
                return false;
            node->len = mpack_load_u8(tree->data + tree->size + 1);
            return mpack_tree_parse_bytes(tree, node);

        // bin16
        case 0xc5:
            node->type = mpack_type_bin;
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint16_t)))
                return false;
            node->len = mpack_load_u16(tree->data + tree->size + 1);
            return mpack_tree_parse_bytes(tree, node);

        // bin32
        case 0xc6:
            node->type = mpack_type_bin;
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint32_t)))
                return false;
            node->len = mpack_load_u32(tree->data + tree->size + 1);
            return mpack_tree_parse_bytes(tree, node);

        #if MPACK_EXTENSIONS
        // ext8
        case 0xc7:
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint8_t)))
                return false;
            node->len = mpack_load_u8(tree->data + tree->size + 1);
            return mpack_tree_parse_ext(tree, node);

        // ext16
        case 0xc8:
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint16_t)))
                return false;
            node->len = mpack_load_u16(tree->data + tree->size + 1);
            return mpack_tree_parse_ext(tree, node);

        // ext32
        case 0xc9:
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint32_t)))
                return false;
            node->len = mpack_load_u32(tree->data + tree->size + 1);
            return mpack_tree_parse_ext(tree, node);
        #endif

        // float
        case 0xca:
            #if MPACK_FLOAT
            if (!mpack_tree_reserve_bytes(tree, sizeof(float)))
                return false;
            node->value.f = mpack_load_float(tree->data + tree->size + 1);
            #else
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint32_t)))
                return false;
            node->value.f = mpack_load_u32(tree->data + tree->size + 1);
            #endif
            node->type = mpack_type_float;
            return true;

        // double
        case 0xcb:
            #if MPACK_DOUBLE
            if (!mpack_tree_reserve_bytes(tree, sizeof(double)))
                return false;
            node->value.d = mpack_load_double(tree->data + tree->size + 1);
            #else
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint64_t)))
                return false;
            node->value.d = mpack_load_u64(tree->data + tree->size + 1);
            #endif
            node->type = mpack_type_double;
            return true;

        // uint8
        case 0xcc:
            node->type = mpack_type_uint;
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint8_t)))
                return false;
            node->value.u = mpack_load_u8(tree->data + tree->size + 1);
            return true;

        // uint16
        case 0xcd:
            node->type = mpack_type_uint;
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint16_t)))
                return false;
            node->value.u = mpack_load_u16(tree->data + tree->size + 1);
            return true;

        // uint32
        case 0xce:
            node->type = mpack_type_uint;
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint32_t)))
                return false;
            node->value.u = mpack_load_u32(tree->data + tree->size + 1);
            return true;

        // uint64
        case 0xcf:
            node->type = mpack_type_uint;
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint64_t)))
                return false;
            node->value.u = mpack_load_u64(tree->data + tree->size + 1);
            return true;

        // int8
        case 0xd0:
            node->type = mpack_type_int;
            if (!mpack_tree_reserve_bytes(tree, sizeof(int8_t)))
                return false;
            node->value.i = mpack_load_i8(tree->data + tree->size + 1);
            return true;

        // int16
        case 0xd1:
            node->type = mpack_type_int;
            if (!mpack_tree_reserve_bytes(tree, sizeof(int16_t)))
                return false;
            node->value.i = mpack_load_i16(tree->data + tree->size + 1);
            return true;

        // int32
        case 0xd2:
            node->type = mpack_type_int;
            if (!mpack_tree_reserve_bytes(tree, sizeof(int32_t)))
                return false;
            node->value.i = mpack_load_i32(tree->data + tree->size + 1);
            return true;

        // int64
        case 0xd3:
            node->type = mpack_type_int;
            if (!mpack_tree_reserve_bytes(tree, sizeof(int64_t)))
                return false;
            node->value.i = mpack_load_i64(tree->data + tree->size + 1);
            return true;

        #if MPACK_EXTENSIONS
        // fixext1
        case 0xd4:
            node->len = 1;
            return mpack_tree_parse_ext(tree, node);

        // fixext2
        case 0xd5:
            node->len = 2;
            return mpack_tree_parse_ext(tree, node);

        // fixext4
        case 0xd6:
            node->len = 4;
            return mpack_tree_parse_ext(tree, node);

        // fixext8
        case 0xd7:
            node->len = 8;
            return mpack_tree_parse_ext(tree, node);

        // fixext16
        case 0xd8:
            node->len = 16;
            return mpack_tree_parse_ext(tree, node);
        #endif

        // str8
        case 0xd9:
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint8_t)))
                return false;
            node->len = mpack_load_u8(tree->data + tree->size + 1);
            node->type = mpack_type_str;
            return mpack_tree_parse_bytes(tree, node);

        // str16
        case 0xda:
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint16_t)))
                return false;
            node->len = mpack_load_u16(tree->data + tree->size + 1);
            node->type = mpack_type_str;
            return mpack_tree_parse_bytes(tree, node);

        // str32
        case 0xdb:
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint32_t)))
                return false;
            node->len = mpack_load_u32(tree->data + tree->size + 1);
            node->type = mpack_type_str;
            return mpack_tree_parse_bytes(tree, node);

        // array16
        case 0xdc:
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint16_t)))
                return false;
            node->len = mpack_load_u16(tree->data + tree->size + 1);
            node->type = mpack_type_array;
            return mpack_tree_parse_children(tree, node);

        // array32
        case 0xdd:
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint32_t)))
                return false;
            node->len = mpack_load_u32(tree->data + tree->size + 1);
            node->type = mpack_type_array;
            return mpack_tree_parse_children(tree, node);

        // map16
        case 0xde:
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint16_t)))
                return false;
            node->len = mpack_load_u16(tree->data + tree->size + 1);
            node->type = mpack_type_map;
            return mpack_tree_parse_children(tree, node);

        // map32
        case 0xdf:
            if (!mpack_tree_reserve_bytes(tree, sizeof(uint32_t)))
                return false;
            node->len = mpack_load_u32(tree->data + tree->size + 1);
            node->type = mpack_type_map;
            return mpack_tree_parse_children(tree, node);

        // reserved
        case 0xc1:
            mpack_tree_flag_error(tree, mpack_error_invalid);
            return false;

        #if !MPACK_EXTENSIONS
        // ext
        case 0xc7: // fallthrough
        case 0xc8: // fallthrough
        case 0xc9: // fallthrough
        // fixext
        case 0xd4: // fallthrough
        case 0xd5: // fallthrough
        case 0xd6: // fallthrough
        case 0xd7: // fallthrough
        case 0xd8:
            mpack_tree_flag_error(tree, mpack_error_unsupported);
            return false;
        #endif

        #if MPACK_OPTIMIZE_FOR_SIZE
        // any other bytes should have been handled by the infix switch
        default:
            break;
        #endif
    }

    mpack_assert(0, "unreachable");
    return false;
}

static bool mpack_tree_parse_node(mpack_tree_t* tree, mpack_node_data_t* node) {
    mpack_log("parsing a node at position %i in level %i\n",
            (int)tree->size, (int)tree->parser.level);

    if (!mpack_tree_parse_node_contents(tree, node)) {
        mpack_log("node parsing returned false\n");
        return false;
    }

    tree->parser.possible_nodes_left -= tree->parser.current_node_reserved;

    // The reserve for the current node does not include the initial byte
    // previously reserved as part of its parent.
    size_t node_size = tree->parser.current_node_reserved + 1;

    // If the parsed type is a map or array, the reserve includes one byte for
    // each child. We want to subtract these out of possible_nodes_left, but
    // not out of the current size of the tree.
    if (node->type == mpack_type_array)
        node_size -= node->len;
    else if (node->type == mpack_type_map)
        node_size -= node->len * 2;
    tree->size += node_size;

    mpack_log("parsed a node of type %s of %i bytes and "
            "%i additional bytes reserved for children.\n",
            mpack_type_to_string(node->type), (int)node_size,
            (int)tree->parser.current_node_reserved + 1 - (int)node_size);

    return true;
}

/*
 * We read nodes in a loop instead of recursively for maximum performance. The
 * stack holds the amount of children left to read in each level of the tree.
 * Parsing can pause and resume when more data becomes available.
 */
static bool mpack_tree_continue_parsing(mpack_tree_t* tree) {
    if (mpack_tree_error(tree) != mpack_ok)
        return false;

    mpack_tree_parser_t* parser = &tree->parser;
    mpack_assert(parser->state == mpack_tree_parse_state_in_progress);
    mpack_log("parsing tree elements, %i bytes in buffer\n", (int)tree->data_length);

    // we loop parsing nodes until the parse stack is empty. we break
    // by returning out of the function.
    while (true) {
        mpack_node_data_t* node = parser->stack[parser->level].child;
        size_t level = parser->level;
        if (!mpack_tree_parse_node(tree, node))
            return false;
        --parser->stack[level].left;
        ++parser->stack[level].child;

        mpack_assert(mpack_tree_error(tree) == mpack_ok,
                "mpack_tree_parse_node() should have returned false due to error!");

        // pop empty stack levels, exiting the outer loop when the stack is empty.
        // (we could tail-optimize containers by pre-emptively popping empty
        // stack levels before reading the new element, this way we wouldn't
        // have to loop. but we eventually want to use the parse stack to give
        // better error messages that contain the location of the error, so
        // it needs to be complete.)
        while (parser->stack[parser->level].left == 0) {
            if (parser->level == 0)
                return true;
            --parser->level;
        }
    }
}

static void mpack_tree_cleanup(mpack_tree_t* tree) {
    MPACK_UNUSED(tree);

    #ifdef MPACK_MALLOC
    if (tree->parser.stack_owned) {
        MPACK_FREE(tree->parser.stack);
        tree->parser.stack = NULL;
        tree->parser.stack_owned = false;
    }

    mpack_tree_page_t* page = tree->next;
    while (page != NULL) {
        mpack_tree_page_t* next = page->next;
        mpack_log("freeing page %p\n", (void*)page);
        MPACK_FREE(page);
        page = next;
    }
    tree->next = NULL;
    #endif
}

static bool mpack_tree_parse_start(mpack_tree_t* tree) {
    if (mpack_tree_error(tree) != mpack_ok)
        return false;

    mpack_tree_parser_t* parser = &tree->parser;
    mpack_assert(parser->state != mpack_tree_parse_state_in_progress,
            "previous parsing was not finished!");

    if (parser->state == mpack_tree_parse_state_parsed)
        mpack_tree_cleanup(tree);

    mpack_log("starting parse\n");
    tree->parser.state = mpack_tree_parse_state_in_progress;
    tree->parser.current_node_reserved = 0;

    // check if we previously parsed a tree
    if (tree->size > 0) {
        #ifdef MPACK_MALLOC
        // if we're buffered, move the remaining data back to the
        // start of the buffer
        // TODO: This is not ideal performance-wise. We should only move data
        // when we need to call the fill function.
        // TODO: We could consider shrinking the buffer here, especially if we
        // determine that the fill function is providing less than a quarter of
        // the buffer size or if messages take up less than a quarter of the
        // buffer size. Maybe this should be configurable.
        if (tree->buffer != NULL) {
            mpack_memmove(tree->buffer, tree->buffer + tree->size, tree->data_length - tree->size);
        }
        else
        #endif
        // otherwise advance past the parsed data
        {
            tree->data += tree->size;
        }
        tree->data_length -= tree->size;
        tree->size = 0;
        tree->node_count = 0;
    }

    // make sure we have at least one byte available before allocating anything
    parser->possible_nodes_left = tree->data_length;
    if (!mpack_tree_reserve_bytes(tree, sizeof(uint8_t))) {
        tree->parser.state = mpack_tree_parse_state_not_started;
        return false;
    }
    mpack_log("parsing tree at %p starting with byte %x\n", tree->data, (uint8_t)tree->data[0]);
    parser->possible_nodes_left -= 1;
    tree->node_count = 1;

    #ifdef MPACK_MALLOC
    parser->stack = parser->stack_local;
    parser->stack_owned = false;
    parser->stack_capacity = sizeof(parser->stack_local) / sizeof(*parser->stack_local);

    if (tree->pool == NULL) {

        // allocate first page
        mpack_tree_page_t* page = (mpack_tree_page_t*)MPACK_MALLOC(MPACK_PAGE_ALLOC_SIZE);
        mpack_log("allocated initial page %p of size %i count %i\n",
                (void*)page, (int)MPACK_PAGE_ALLOC_SIZE, (int)MPACK_NODES_PER_PAGE);
        if (page == NULL) {
            tree->error = mpack_error_memory;
            return false;
        }
        page->next = NULL;
        tree->next = page;

        parser->nodes = page->nodes;
        parser->nodes_left = MPACK_NODES_PER_PAGE;
    }
    else
    #endif
    {
        // otherwise use the provided pool
        mpack_assert(tree->pool != NULL, "no pool provided?");
        parser->nodes = tree->pool;
        parser->nodes_left = tree->pool_count;
    }

    tree->root = parser->nodes;
    ++parser->nodes;
    --parser->nodes_left;

    parser->level = 0;
    parser->stack[0].child = tree->root;
    parser->stack[0].left = 1;

    return true;
}

void mpack_tree_parse(mpack_tree_t* tree) {
    if (mpack_tree_error(tree) != mpack_ok)
        return;

    if (tree->parser.state != mpack_tree_parse_state_in_progress) {
        if (!mpack_tree_parse_start(tree)) {
            mpack_tree_flag_error(tree, (tree->read_fn == NULL) ?
                    mpack_error_invalid : mpack_error_io);
            return;
        }
    }

    if (!mpack_tree_continue_parsing(tree)) {
        if (mpack_tree_error(tree) != mpack_ok)
            return;

        // We're parsing synchronously on a blocking fill function. If we
        // didn't completely finish parsing the tree, it's an error.
        mpack_log("tree parsing incomplete. flagging error.\n");
        mpack_tree_flag_error(tree, (tree->read_fn == NULL) ?
                mpack_error_invalid : mpack_error_io);
        return;
    }

    mpack_assert(mpack_tree_error(tree) == mpack_ok);
    mpack_assert(tree->parser.level == 0);
    tree->parser.state = mpack_tree_parse_state_parsed;
    mpack_log("parsed tree of %i bytes, %i bytes left\n", (int)tree->size, (int)tree->parser.possible_nodes_left);
    mpack_log("%i nodes in final page\n", (int)tree->parser.nodes_left);
}

bool mpack_tree_try_parse(mpack_tree_t* tree) {
    if (mpack_tree_error(tree) != mpack_ok)
        return false;

    if (tree->parser.state != mpack_tree_parse_state_in_progress)
        if (!mpack_tree_parse_start(tree))
            return false;

    if (!mpack_tree_continue_parsing(tree))
        return false;

    mpack_assert(mpack_tree_error(tree) == mpack_ok);
    mpack_assert(tree->parser.level == 0);
    tree->parser.state = mpack_tree_parse_state_parsed;
    return true;
}



/*
 * Tree functions
 */

mpack_node_t mpack_tree_root(mpack_tree_t* tree) {
    if (mpack_tree_error(tree) != mpack_ok)
        return mpack_tree_nil_node(tree);

    // We check that a tree was parsed successfully and assert if not. You must
    // call mpack_tree_parse() (or mpack_tree_try_parse() with a success
    // result) in order to access the root node.
    if (tree->parser.state != mpack_tree_parse_state_parsed) {
        mpack_break("Tree has not been parsed! "
                "Did you call mpack_tree_parse() or mpack_tree_try_parse()?");
        mpack_tree_flag_error(tree, mpack_error_bug);
        return mpack_tree_nil_node(tree);
    }

    return mpack_node(tree, tree->root);
}

static void mpack_tree_init_clear(mpack_tree_t* tree) {
    mpack_memset(tree, 0, sizeof(*tree));
    tree->nil_node.type = mpack_type_nil;
    tree->missing_node.type = mpack_type_missing;
    tree->max_size = SIZE_MAX;
    tree->max_nodes = SIZE_MAX;
}

#ifdef MPACK_MALLOC
void mpack_tree_init_data(mpack_tree_t* tree, const char* data, size_t length) {
    mpack_tree_init_clear(tree);

    MPACK_STATIC_ASSERT(MPACK_NODE_PAGE_SIZE >= sizeof(mpack_tree_page_t),
            "MPACK_NODE_PAGE_SIZE is too small");

    MPACK_STATIC_ASSERT(MPACK_PAGE_ALLOC_SIZE <= MPACK_NODE_PAGE_SIZE,
            "incorrect page rounding?");

    tree->data = data;
    tree->data_length = length;
    tree->pool = NULL;
    tree->pool_count = 0;
    tree->next = NULL;

    mpack_log("===========================\n");
    mpack_log("initializing tree with data of size %i\n", (int)length);
}
#endif

void mpack_tree_init_pool(mpack_tree_t* tree, const char* data, size_t length,
        mpack_node_data_t* node_pool, size_t node_pool_count)
{
    mpack_tree_init_clear(tree);
    #ifdef MPACK_MALLOC
    tree->next = NULL;
    #endif

    if (node_pool_count == 0) {
        mpack_break("initial page has no nodes!");
        mpack_tree_flag_error(tree, mpack_error_bug);
        return;
    }

    tree->data = data;
    tree->data_length = length;
    tree->pool = node_pool;
    tree->pool_count = node_pool_count;

    mpack_log("===========================\n");
    mpack_log("initializing tree with data of size %i and pool of count %i\n",
            (int)length, (int)node_pool_count);
}

void mpack_tree_init_error(mpack_tree_t* tree, mpack_error_t error) {
    mpack_tree_init_clear(tree);
    tree->error = error;

    mpack_log("===========================\n");
    mpack_log("initializing tree error state %i\n", (int)error);
}

#ifdef MPACK_MALLOC
void mpack_tree_init_stream(mpack_tree_t* tree, mpack_tree_read_t read_fn, void* context,
        size_t max_message_size, size_t max_message_nodes) {
    mpack_tree_init_clear(tree);

    tree->read_fn = read_fn;
    tree->context = context;

    mpack_tree_set_limits(tree, max_message_size, max_message_nodes);
    tree->max_size = max_message_size;
    tree->max_nodes = max_message_nodes;

    mpack_log("===========================\n");
    mpack_log("initializing tree with stream, max size %i max nodes %i\n",
            (int)max_message_size, (int)max_message_nodes);
}
#endif

void mpack_tree_set_limits(mpack_tree_t* tree, size_t max_message_size, size_t max_message_nodes) {
    mpack_assert(max_message_size > 0);
    mpack_assert(max_message_nodes > 0);
    tree->max_size = max_message_size;
    tree->max_nodes = max_message_nodes;
}

#if MPACK_STDIO
typedef struct mpack_file_tree_t {
    char* data;
    size_t size;
    char buffer[MPACK_BUFFER_SIZE];
} mpack_file_tree_t;

static void mpack_file_tree_teardown(mpack_tree_t* tree) {
    mpack_file_tree_t* file_tree = (mpack_file_tree_t*)tree->context;
    MPACK_FREE(file_tree->data);
    MPACK_FREE(file_tree);
}

static bool mpack_file_tree_read(mpack_tree_t* tree, mpack_file_tree_t* file_tree, FILE* file, size_t max_bytes) {

    // get the file size
    errno = 0;
    int error = 0;
    fseek(file, 0, SEEK_END);
    error |= errno;
    long size = ftell(file);
    error |= errno;
    fseek(file, 0, SEEK_SET);
    error |= errno;

    // check for errors
    if (error != 0 || size < 0) {
        mpack_tree_init_error(tree, mpack_error_io);
        return false;
    }
    if (size == 0) {
        mpack_tree_init_error(tree, mpack_error_invalid);
        return false;
    }

    // make sure the size is less than max_bytes
    // (this mess exists to safely convert between long and size_t regardless of their widths)
    if (max_bytes != 0 && (((uint64_t)LONG_MAX > (uint64_t)SIZE_MAX && size > (long)SIZE_MAX) || (size_t)size > max_bytes)) {
        mpack_tree_init_error(tree, mpack_error_too_big);
        return false;
    }

    // allocate data
    file_tree->data = (char*)MPACK_MALLOC((size_t)size);
    if (file_tree->data == NULL) {
        mpack_tree_init_error(tree, mpack_error_memory);
        return false;
    }

    // read the file
    long total = 0;
    while (total < size) {
        size_t read = fread(file_tree->data + total, 1, (size_t)(size - total), file);
        if (read <= 0) {
            mpack_tree_init_error(tree, mpack_error_io);
            MPACK_FREE(file_tree->data);
            return false;
        }
        total += (long)read;
    }

    file_tree->size = (size_t)size;
    return true;
}

static bool mpack_tree_file_check_max_bytes(mpack_tree_t* tree, size_t max_bytes) {

    // the C STDIO family of file functions use long (e.g. ftell)
    if (max_bytes > LONG_MAX) {
        mpack_break("max_bytes of %" PRIu64 " is invalid, maximum is LONG_MAX", (uint64_t)max_bytes);
        mpack_tree_init_error(tree, mpack_error_bug);
        return false;
    }

    return true;
}

static void mpack_tree_init_stdfile_noclose(mpack_tree_t* tree, FILE* stdfile, size_t max_bytes) {

    // allocate file tree
    mpack_file_tree_t* file_tree = (mpack_file_tree_t*) MPACK_MALLOC(sizeof(mpack_file_tree_t));
    if (file_tree == NULL) {
        mpack_tree_init_error(tree, mpack_error_memory);
        return;
    }

    // read all data
    if (!mpack_file_tree_read(tree, file_tree, stdfile, max_bytes)) {
        MPACK_FREE(file_tree);
        return;
    }

    mpack_tree_init_data(tree, file_tree->data, file_tree->size);
    mpack_tree_set_context(tree, file_tree);
    mpack_tree_set_teardown(tree, mpack_file_tree_teardown);
}

void mpack_tree_init_stdfile(mpack_tree_t* tree, FILE* stdfile, size_t max_bytes, bool close_when_done) {
    if (!mpack_tree_file_check_max_bytes(tree, max_bytes))
        return;

    mpack_tree_init_stdfile_noclose(tree, stdfile, max_bytes);

    if (close_when_done)
        fclose(stdfile);
}

void mpack_tree_init_filename(mpack_tree_t* tree, const char* filename, size_t max_bytes) {
    if (!mpack_tree_file_check_max_bytes(tree, max_bytes))
        return;

    // open the file
    FILE* file = fopen(filename, "rb");
    if (!file) {
        mpack_tree_init_error(tree, mpack_error_io);
        return;
    }

    mpack_tree_init_stdfile(tree, file, max_bytes, true);
}
#endif

mpack_error_t mpack_tree_destroy(mpack_tree_t* tree) {
    mpack_tree_cleanup(tree);

    #ifdef MPACK_MALLOC
    if (tree->buffer)
        MPACK_FREE(tree->buffer);
    #endif

    if (tree->teardown)
        tree->teardown(tree);
    tree->teardown = NULL;

    return tree->error;
}

void mpack_tree_flag_error(mpack_tree_t* tree, mpack_error_t error) {
    if (tree->error == mpack_ok) {
        mpack_log("tree %p setting error %i: %s\n", (void*)tree, (int)error, mpack_error_to_string(error));
        tree->error = error;
        if (tree->error_fn)
            tree->error_fn(tree, error);
    }

}



/*
 * Node misc functions
 */

void mpack_node_flag_error(mpack_node_t node, mpack_error_t error) {
    mpack_tree_flag_error(node.tree, error);
}

mpack_tag_t mpack_node_tag(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return mpack_tag_nil();

    mpack_tag_t tag = MPACK_TAG_ZERO;

    tag.type = node.data->type;
    switch (node.data->type) {
        case mpack_type_missing:
            // If a node is missing, I don't know if it makes sense to ask for
            // a tag for it. We'll return a missing tag to match the missing
            // node I guess, but attempting to use the tag for anything (like
            // writing it for example) will flag mpack_error_bug.
            break;
        case mpack_type_nil:                                            break;
        case mpack_type_bool:    tag.v.b = node.data->value.b;          break;
        case mpack_type_float:   tag.v.f = node.data->value.f;          break;
        case mpack_type_double:  tag.v.d = node.data->value.d;          break;
        case mpack_type_int:     tag.v.i = node.data->value.i;          break;
        case mpack_type_uint:    tag.v.u = node.data->value.u;          break;

        case mpack_type_str:     tag.v.l = node.data->len;     break;
        case mpack_type_bin:     tag.v.l = node.data->len;     break;

        #if MPACK_EXTENSIONS
        case mpack_type_ext:
            tag.v.l = node.data->len;
            tag.exttype = mpack_node_exttype_unchecked(node);
            break;
        #endif

        case mpack_type_array:   tag.v.n = node.data->len;  break;
        case mpack_type_map:     tag.v.n = node.data->len;  break;

        default:
            mpack_assert(0, "unrecognized type %i", (int)node.data->type);
            break;
    }
    return tag;
}

#if MPACK_DEBUG && MPACK_STDIO
static void mpack_node_print_element(mpack_node_t node, mpack_print_t* print, size_t depth) {
    mpack_node_data_t* data = node.data;
    size_t i,j;
    switch (data->type) {
        case mpack_type_str:
            {
                mpack_print_append_cstr(print, "\"");
                const char* bytes = mpack_node_data_unchecked(node);
                for (i = 0; i < data->len; ++i) {
                    char c = bytes[i];
                    switch (c) {
                        case '\n': mpack_print_append_cstr(print, "\\n"); break;
                        case '\\': mpack_print_append_cstr(print, "\\\\"); break;
                        case '"': mpack_print_append_cstr(print, "\\\""); break;
                        default: mpack_print_append(print, &c, 1); break;
                    }
                }
                mpack_print_append_cstr(print, "\"");
            }
            break;

        case mpack_type_array:
            mpack_print_append_cstr(print, "[\n");
            for (i = 0; i < data->len; ++i) {
                for (j = 0; j < depth + 1; ++j)
                    mpack_print_append_cstr(print, "    ");
                mpack_node_print_element(mpack_node_array_at(node, i), print, depth + 1);
                if (i != data->len - 1)
                    mpack_print_append_cstr(print, ",");
                mpack_print_append_cstr(print, "\n");
            }
            for (i = 0; i < depth; ++i)
                mpack_print_append_cstr(print, "    ");
            mpack_print_append_cstr(print, "]");
            break;

        case mpack_type_map:
            mpack_print_append_cstr(print, "{\n");
            for (i = 0; i < data->len; ++i) {
                for (j = 0; j < depth + 1; ++j)
                    mpack_print_append_cstr(print, "    ");
                mpack_node_print_element(mpack_node_map_key_at(node, i), print, depth + 1);
                mpack_print_append_cstr(print, ": ");
                mpack_node_print_element(mpack_node_map_value_at(node, i), print, depth + 1);
                if (i != data->len - 1)
                    mpack_print_append_cstr(print, ",");
                mpack_print_append_cstr(print, "\n");
            }
            for (i = 0; i < depth; ++i)
                mpack_print_append_cstr(print, "    ");
            mpack_print_append_cstr(print, "}");
            break;

        default:
            {
                const char* prefix = NULL;
                size_t prefix_length = 0;
                if (mpack_node_type(node) == mpack_type_bin
                        #if MPACK_EXTENSIONS
                        || mpack_node_type(node) == mpack_type_ext
                        #endif
                ) {
                    prefix = mpack_node_data(node);
                    prefix_length = mpack_node_data_len(node);
                }

                char buf[256];
                mpack_tag_t tag = mpack_node_tag(node);
                mpack_tag_debug_pseudo_json(tag, buf, sizeof(buf), prefix, prefix_length);
                mpack_print_append_cstr(print, buf);
            }
            break;
    }
}

void mpack_node_print_to_buffer(mpack_node_t node, char* buffer, size_t buffer_size) {
    if (buffer_size == 0) {
        mpack_assert(false, "buffer size is zero!");
        return;
    }

    mpack_print_t print;
    mpack_memset(&print, 0, sizeof(print));
    print.buffer = buffer;
    print.size = buffer_size;
    mpack_node_print_element(node, &print, 0);
    mpack_print_append(&print, "",  1); // null-terminator
    mpack_print_flush(&print);

    // we always make sure there's a null-terminator at the end of the buffer
    // in case we ran out of space.
    print.buffer[print.size - 1] = '\0';
}

void mpack_node_print_to_callback(mpack_node_t node, mpack_print_callback_t callback, void* context) {
    char buffer[1024];
    mpack_print_t print;
    mpack_memset(&print, 0, sizeof(print));
    print.buffer = buffer;
    print.size = sizeof(buffer);
    print.callback = callback;
    print.context = context;
    mpack_node_print_element(node, &print, 0);
    mpack_print_flush(&print);
}

void mpack_node_print_to_file(mpack_node_t node, FILE* file) {
    mpack_assert(file != NULL, "file is NULL");

    char buffer[1024];
    mpack_print_t print;
    mpack_memset(&print, 0, sizeof(print));
    print.buffer = buffer;
    print.size = sizeof(buffer);
    print.callback = &mpack_print_file_callback;
    print.context = file;

    size_t depth = 2;
    size_t i;
    for (i = 0; i < depth; ++i)
        mpack_print_append_cstr(&print, "    ");
    mpack_node_print_element(node, &print, depth);
    mpack_print_append_cstr(&print, "\n");
    mpack_print_flush(&print);
}
#endif



/*
 * Node Value Functions
 */

#if MPACK_EXTENSIONS
mpack_timestamp_t mpack_node_timestamp(mpack_node_t node) {
    mpack_timestamp_t timestamp = {0, 0};

    // we'll let mpack_node_exttype() do most checks
    if (mpack_node_exttype(node) != MPACK_EXTTYPE_TIMESTAMP) {
        mpack_log("exttype %i\n", mpack_node_exttype(node));
        mpack_node_flag_error(node, mpack_error_type);
        return timestamp;
    }

    const char* p = mpack_node_data_unchecked(node);

    switch (node.data->len) {
        case 4:
            timestamp.nanoseconds = 0;
            timestamp.seconds = mpack_load_u32(p);
            break;

        case 8: {
            uint64_t value = mpack_load_u64(p);
            timestamp.nanoseconds = (uint32_t)(value >> 34);
            timestamp.seconds = value & ((MPACK_UINT64_C(1) << 34) - 1);
            break;
        }

        case 12:
            timestamp.nanoseconds = mpack_load_u32(p);
            timestamp.seconds = mpack_load_i64(p + 4);
            break;

        default:
            mpack_tree_flag_error(node.tree, mpack_error_invalid);
            return timestamp;
    }

    if (timestamp.nanoseconds > MPACK_TIMESTAMP_NANOSECONDS_MAX) {
        mpack_tree_flag_error(node.tree, mpack_error_invalid);
        mpack_timestamp_t zero = {0, 0};
        return zero;
    }

    return timestamp;
}

int64_t mpack_node_timestamp_seconds(mpack_node_t node) {
    return mpack_node_timestamp(node).seconds;
}

uint32_t mpack_node_timestamp_nanoseconds(mpack_node_t node) {
    return mpack_node_timestamp(node).nanoseconds;
}
#endif



/*
 * Node Data Functions
 */

void mpack_node_check_utf8(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return;
    mpack_node_data_t* data = node.data;
    if (data->type != mpack_type_str || !mpack_utf8_check(mpack_node_data_unchecked(node), data->len))
        mpack_node_flag_error(node, mpack_error_type);
}

void mpack_node_check_utf8_cstr(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return;
    mpack_node_data_t* data = node.data;
    if (data->type != mpack_type_str || !mpack_utf8_check_no_null(mpack_node_data_unchecked(node), data->len))
        mpack_node_flag_error(node, mpack_error_type);
}

size_t mpack_node_copy_data(mpack_node_t node, char* buffer, size_t bufsize) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    mpack_assert(bufsize == 0 || buffer != NULL, "buffer is NULL for maximum of %i bytes", (int)bufsize);

    mpack_type_t type = node.data->type;
    if (type != mpack_type_str && type != mpack_type_bin
            #if MPACK_EXTENSIONS
            && type != mpack_type_ext
            #endif
    ) {
        mpack_node_flag_error(node, mpack_error_type);
        return 0;
    }

    if (node.data->len > bufsize) {
        mpack_node_flag_error(node, mpack_error_too_big);
        return 0;
    }

    mpack_memcpy(buffer, mpack_node_data_unchecked(node), node.data->len);
    return (size_t)node.data->len;
}

size_t mpack_node_copy_utf8(mpack_node_t node, char* buffer, size_t bufsize) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    mpack_assert(bufsize == 0 || buffer != NULL, "buffer is NULL for maximum of %i bytes", (int)bufsize);

    mpack_type_t type = node.data->type;
    if (type != mpack_type_str) {
        mpack_node_flag_error(node, mpack_error_type);
        return 0;
    }

    if (node.data->len > bufsize) {
        mpack_node_flag_error(node, mpack_error_too_big);
        return 0;
    }

    if (!mpack_utf8_check(mpack_node_data_unchecked(node), node.data->len)) {
        mpack_node_flag_error(node, mpack_error_type);
        return 0;
    }

    mpack_memcpy(buffer, mpack_node_data_unchecked(node), node.data->len);
    return (size_t)node.data->len;
}

void mpack_node_copy_cstr(mpack_node_t node, char* buffer, size_t bufsize) {

    // we can't break here because the error isn't recoverable; we
    // have to add a null-terminator.
    mpack_assert(buffer != NULL, "buffer is NULL");
    mpack_assert(bufsize >= 1, "buffer size is zero; you must have room for at least a null-terminator");

    if (mpack_node_error(node) != mpack_ok) {
        buffer[0] = '\0';
        return;
    }

    if (node.data->type != mpack_type_str) {
        buffer[0] = '\0';
        mpack_node_flag_error(node, mpack_error_type);
        return;
    }

    if (node.data->len > bufsize - 1) {
        buffer[0] = '\0';
        mpack_node_flag_error(node, mpack_error_too_big);
        return;
    }

    if (!mpack_str_check_no_null(mpack_node_data_unchecked(node), node.data->len)) {
        buffer[0] = '\0';
        mpack_node_flag_error(node, mpack_error_type);
        return;
    }

    mpack_memcpy(buffer, mpack_node_data_unchecked(node), node.data->len);
    buffer[node.data->len] = '\0';
}

void mpack_node_copy_utf8_cstr(mpack_node_t node, char* buffer, size_t bufsize) {

    // we can't break here because the error isn't recoverable; we
    // have to add a null-terminator.
    mpack_assert(buffer != NULL, "buffer is NULL");
    mpack_assert(bufsize >= 1, "buffer size is zero; you must have room for at least a null-terminator");

    if (mpack_node_error(node) != mpack_ok) {
        buffer[0] = '\0';
        return;
    }

    if (node.data->type != mpack_type_str) {
        buffer[0] = '\0';
        mpack_node_flag_error(node, mpack_error_type);
        return;
    }

    if (node.data->len > bufsize - 1) {
        buffer[0] = '\0';
        mpack_node_flag_error(node, mpack_error_too_big);
        return;
    }

    if (!mpack_utf8_check_no_null(mpack_node_data_unchecked(node), node.data->len)) {
        buffer[0] = '\0';
        mpack_node_flag_error(node, mpack_error_type);
        return;
    }

    mpack_memcpy(buffer, mpack_node_data_unchecked(node), node.data->len);
    buffer[node.data->len] = '\0';
}

#ifdef MPACK_MALLOC
char* mpack_node_data_alloc(mpack_node_t node, size_t maxlen) {
    if (mpack_node_error(node) != mpack_ok)
        return NULL;

    // make sure this is a valid data type
    mpack_type_t type = node.data->type;
    if (type != mpack_type_str && type != mpack_type_bin
            #if MPACK_EXTENSIONS
            && type != mpack_type_ext
            #endif
    ) {
        mpack_node_flag_error(node, mpack_error_type);
        return NULL;
    }

    if (node.data->len > maxlen) {
        mpack_node_flag_error(node, mpack_error_too_big);
        return NULL;
    }

    char* ret = (char*) MPACK_MALLOC((size_t)node.data->len);
    if (ret == NULL) {
        mpack_node_flag_error(node, mpack_error_memory);
        return NULL;
    }

    mpack_memcpy(ret, mpack_node_data_unchecked(node), node.data->len);
    return ret;
}

char* mpack_node_cstr_alloc(mpack_node_t node, size_t maxlen) {
    if (mpack_node_error(node) != mpack_ok)
        return NULL;

    // make sure maxlen makes sense
    if (maxlen < 1) {
        mpack_break("maxlen is zero; you must have room for at least a null-terminator");
        mpack_node_flag_error(node, mpack_error_bug);
        return NULL;
    }

    if (node.data->type != mpack_type_str) {
        mpack_node_flag_error(node, mpack_error_type);
        return NULL;
    }

    if (node.data->len > maxlen - 1) {
        mpack_node_flag_error(node, mpack_error_too_big);
        return NULL;
    }

    if (!mpack_str_check_no_null(mpack_node_data_unchecked(node), node.data->len)) {
        mpack_node_flag_error(node, mpack_error_type);
        return NULL;
    }

    char* ret = (char*) MPACK_MALLOC((size_t)(node.data->len + 1));
    if (ret == NULL) {
        mpack_node_flag_error(node, mpack_error_memory);
        return NULL;
    }

    mpack_memcpy(ret, mpack_node_data_unchecked(node), node.data->len);
    ret[node.data->len] = '\0';
    return ret;
}

char* mpack_node_utf8_cstr_alloc(mpack_node_t node, size_t maxlen) {
    if (mpack_node_error(node) != mpack_ok)
        return NULL;

    // make sure maxlen makes sense
    if (maxlen < 1) {
        mpack_break("maxlen is zero; you must have room for at least a null-terminator");
        mpack_node_flag_error(node, mpack_error_bug);
        return NULL;
    }

    if (node.data->type != mpack_type_str) {
        mpack_node_flag_error(node, mpack_error_type);
        return NULL;
    }

    if (node.data->len > maxlen - 1) {
        mpack_node_flag_error(node, mpack_error_too_big);
        return NULL;
    }

    if (!mpack_utf8_check_no_null(mpack_node_data_unchecked(node), node.data->len)) {
        mpack_node_flag_error(node, mpack_error_type);
        return NULL;
    }

    char* ret = (char*) MPACK_MALLOC((size_t)(node.data->len + 1));
    if (ret == NULL) {
        mpack_node_flag_error(node, mpack_error_memory);
        return NULL;
    }

    mpack_memcpy(ret, mpack_node_data_unchecked(node), node.data->len);
    ret[node.data->len] = '\0';
    return ret;
}
#endif


/*
 * Compound Node Functions
 */

static mpack_node_data_t* mpack_node_map_int_impl(mpack_node_t node, int64_t num) {
    if (mpack_node_error(node) != mpack_ok)
        return NULL;

    if (node.data->type != mpack_type_map) {
        mpack_node_flag_error(node, mpack_error_type);
        return NULL;
    }

    mpack_node_data_t* found = NULL;

    size_t i;
    for (i = 0; i < node.data->len; ++i) {
        mpack_node_data_t* key = mpack_node_child(node, i * 2);

        if ((key->type == mpack_type_int && key->value.i == num) ||
            (key->type == mpack_type_uint && num >= 0 && key->value.u == (uint64_t)num))
        {
            if (found) {
                mpack_node_flag_error(node, mpack_error_data);
                return NULL;
            }
            found = mpack_node_child(node, i * 2 + 1);
        }
    }

    if (found)
        return found;

    return NULL;
}

static mpack_node_data_t* mpack_node_map_uint_impl(mpack_node_t node, uint64_t num) {
    if (mpack_node_error(node) != mpack_ok)
        return NULL;

    if (node.data->type != mpack_type_map) {
        mpack_node_flag_error(node, mpack_error_type);
        return NULL;
    }

    mpack_node_data_t* found = NULL;

    size_t i;
    for (i = 0; i < node.data->len; ++i) {
        mpack_node_data_t* key = mpack_node_child(node, i * 2);

        if ((key->type == mpack_type_uint && key->value.u == num) ||
            (key->type == mpack_type_int && key->value.i >= 0 && (uint64_t)key->value.i == num))
        {
            if (found) {
                mpack_node_flag_error(node, mpack_error_data);
                return NULL;
            }
            found = mpack_node_child(node, i * 2 + 1);
        }
    }

    if (found)
        return found;

    return NULL;
}

static mpack_node_data_t* mpack_node_map_str_impl(mpack_node_t node, const char* str, size_t length) {
    if (mpack_node_error(node) != mpack_ok)
        return NULL;

    mpack_assert(length == 0 || str != NULL, "str of length %i is NULL", (int)length);

    if (node.data->type != mpack_type_map) {
        mpack_node_flag_error(node, mpack_error_type);
        return NULL;
    }

    mpack_tree_t* tree = node.tree;
    mpack_node_data_t* found = NULL;

    size_t i;
    for (i = 0; i < node.data->len; ++i) {
        mpack_node_data_t* key = mpack_node_child(node, i * 2);

        if (key->type == mpack_type_str && key->len == length &&
                mpack_memcmp(str, mpack_node_data_unchecked(mpack_node(tree, key)), length) == 0) {
            if (found) {
                mpack_node_flag_error(node, mpack_error_data);
                return NULL;
            }
            found = mpack_node_child(node, i * 2 + 1);
        }
    }

    if (found)
        return found;

    return NULL;
}

static mpack_node_t mpack_node_wrap_lookup(mpack_tree_t* tree, mpack_node_data_t* data) {
    if (!data) {
        if (tree->error == mpack_ok)
            mpack_tree_flag_error(tree, mpack_error_data);
        return mpack_tree_nil_node(tree);
    }
    return mpack_node(tree, data);
}

static mpack_node_t mpack_node_wrap_lookup_optional(mpack_tree_t* tree, mpack_node_data_t* data) {
    if (!data) {
        if (tree->error == mpack_ok)
            return mpack_tree_missing_node(tree);
        return mpack_tree_nil_node(tree);
    }
    return mpack_node(tree, data);
}

mpack_node_t mpack_node_map_int(mpack_node_t node, int64_t num) {
    return mpack_node_wrap_lookup(node.tree, mpack_node_map_int_impl(node, num));
}

mpack_node_t mpack_node_map_int_optional(mpack_node_t node, int64_t num) {
    return mpack_node_wrap_lookup_optional(node.tree, mpack_node_map_int_impl(node, num));
}

mpack_node_t mpack_node_map_uint(mpack_node_t node, uint64_t num) {
    return mpack_node_wrap_lookup(node.tree, mpack_node_map_uint_impl(node, num));
}

mpack_node_t mpack_node_map_uint_optional(mpack_node_t node, uint64_t num) {
    return mpack_node_wrap_lookup_optional(node.tree, mpack_node_map_uint_impl(node, num));
}

mpack_node_t mpack_node_map_str(mpack_node_t node, const char* str, size_t length) {
    return mpack_node_wrap_lookup(node.tree, mpack_node_map_str_impl(node, str, length));
}

mpack_node_t mpack_node_map_str_optional(mpack_node_t node, const char* str, size_t length) {
    return mpack_node_wrap_lookup_optional(node.tree, mpack_node_map_str_impl(node, str, length));
}

mpack_node_t mpack_node_map_cstr(mpack_node_t node, const char* cstr) {
    mpack_assert(cstr != NULL, "cstr is NULL");
    return mpack_node_map_str(node, cstr, mpack_strlen(cstr));
}

mpack_node_t mpack_node_map_cstr_optional(mpack_node_t node, const char* cstr) {
    mpack_assert(cstr != NULL, "cstr is NULL");
    return mpack_node_map_str_optional(node, cstr, mpack_strlen(cstr));
}

bool mpack_node_map_contains_int(mpack_node_t node, int64_t num) {
    return mpack_node_map_int_impl(node, num) != NULL;
}

bool mpack_node_map_contains_uint(mpack_node_t node, uint64_t num) {
    return mpack_node_map_uint_impl(node, num) != NULL;
}

bool mpack_node_map_contains_str(mpack_node_t node, const char* str, size_t length) {
    return mpack_node_map_str_impl(node, str, length) != NULL;
}

bool mpack_node_map_contains_cstr(mpack_node_t node, const char* cstr) {
    mpack_assert(cstr != NULL, "cstr is NULL");
    return mpack_node_map_contains_str(node, cstr, mpack_strlen(cstr));
}

size_t mpack_node_enum_optional(mpack_node_t node, const char* strings[], size_t count) {
    if (mpack_node_error(node) != mpack_ok)
        return count;

    // the value is only recognized if it is a string
    if (mpack_node_type(node) != mpack_type_str)
        return count;

    // fetch the string
    const char* key = mpack_node_str(node);
    size_t keylen = mpack_node_strlen(node);
    mpack_assert(mpack_node_error(node) == mpack_ok, "these should not fail");

    // find what key it matches
    size_t i;
    for (i = 0; i < count; ++i) {
        const char* other = strings[i];
        size_t otherlen = mpack_strlen(other);
        if (keylen == otherlen && mpack_memcmp(key, other, keylen) == 0)
            return i;
    }

    // no matches
    return count;
}

size_t mpack_node_enum(mpack_node_t node, const char* strings[], size_t count) {
    size_t value = mpack_node_enum_optional(node, strings, count);
    if (value == count)
        mpack_node_flag_error(node, mpack_error_type);
    return value;
}

mpack_type_t mpack_node_type(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return mpack_type_nil;
    return node.data->type;
}

bool mpack_node_is_nil(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok) {
        // All nodes are treated as nil nodes when we are in error.
        return true;
    }
    return node.data->type == mpack_type_nil;
}

bool mpack_node_is_missing(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok) {
        // errors still return nil nodes, not missing nodes.
        return false;
    }
    return node.data->type == mpack_type_missing;
}

void mpack_node_nil(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return;
    if (node.data->type != mpack_type_nil)
        mpack_node_flag_error(node, mpack_error_type);
}

void mpack_node_missing(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return;
    if (node.data->type != mpack_type_missing)
        mpack_node_flag_error(node, mpack_error_type);
}

bool mpack_node_bool(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return false;

    if (node.data->type == mpack_type_bool)
        return node.data->value.b;

    mpack_node_flag_error(node, mpack_error_type);
    return false;
}

void mpack_node_true(mpack_node_t node) {
    if (mpack_node_bool(node) != true)
        mpack_node_flag_error(node, mpack_error_type);
}

void mpack_node_false(mpack_node_t node) {
    if (mpack_node_bool(node) != false)
        mpack_node_flag_error(node, mpack_error_type);
}

uint8_t mpack_node_u8(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    if (node.data->type == mpack_type_uint) {
        if (node.data->value.u <= MPACK_UINT8_MAX)
            return (uint8_t)node.data->value.u;
    } else if (node.data->type == mpack_type_int) {
        if (node.data->value.i >= 0 && node.data->value.i <= MPACK_UINT8_MAX)
            return (uint8_t)node.data->value.i;
    }

    mpack_node_flag_error(node, mpack_error_type);
    return 0;
}

int8_t mpack_node_i8(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    if (node.data->type == mpack_type_uint) {
        if (node.data->value.u <= MPACK_INT8_MAX)
            return (int8_t)node.data->value.u;
    } else if (node.data->type == mpack_type_int) {
        if (node.data->value.i >= MPACK_INT8_MIN && node.data->value.i <= MPACK_INT8_MAX)
            return (int8_t)node.data->value.i;
    }

    mpack_node_flag_error(node, mpack_error_type);
    return 0;
}

uint16_t mpack_node_u16(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    if (node.data->type == mpack_type_uint) {
        if (node.data->value.u <= MPACK_UINT16_MAX)
            return (uint16_t)node.data->value.u;
    } else if (node.data->type == mpack_type_int) {
        if (node.data->value.i >= 0 && node.data->value.i <= MPACK_UINT16_MAX)
            return (uint16_t)node.data->value.i;
    }

    mpack_node_flag_error(node, mpack_error_type);
    return 0;
}

int16_t mpack_node_i16(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    if (node.data->type == mpack_type_uint) {
        if (node.data->value.u <= MPACK_INT16_MAX)
            return (int16_t)node.data->value.u;
    } else if (node.data->type == mpack_type_int) {
        if (node.data->value.i >= MPACK_INT16_MIN && node.data->value.i <= MPACK_INT16_MAX)
            return (int16_t)node.data->value.i;
    }

    mpack_node_flag_error(node, mpack_error_type);
    return 0;
}

uint32_t mpack_node_u32(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    if (node.data->type == mpack_type_uint) {
        if (node.data->value.u <= MPACK_UINT32_MAX)
            return (uint32_t)node.data->value.u;
    } else if (node.data->type == mpack_type_int) {
        if (node.data->value.i >= 0 && node.data->value.i <= MPACK_UINT32_MAX)
            return (uint32_t)node.data->value.i;
    }

    mpack_node_flag_error(node, mpack_error_type);
    return 0;
}

int32_t mpack_node_i32(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    if (node.data->type == mpack_type_uint) {
        if (node.data->value.u <= MPACK_INT32_MAX)
            return (int32_t)node.data->value.u;
    } else if (node.data->type == mpack_type_int) {
        if (node.data->value.i >= MPACK_INT32_MIN && node.data->value.i <= MPACK_INT32_MAX)
            return (int32_t)node.data->value.i;
    }

    mpack_node_flag_error(node, mpack_error_type);
    return 0;
}

uint64_t mpack_node_u64(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    if (node.data->type == mpack_type_uint) {
        return node.data->value.u;
    } else if (node.data->type == mpack_type_int) {
        if (node.data->value.i >= 0)
            return (uint64_t)node.data->value.i;
    }

    mpack_node_flag_error(node, mpack_error_type);
    return 0;
}

int64_t mpack_node_i64(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    if (node.data->type == mpack_type_uint) {
        if (node.data->value.u <= (uint64_t)MPACK_INT64_MAX)
            return (int64_t)node.data->value.u;
    } else if (node.data->type == mpack_type_int) {
        return node.data->value.i;
    }

    mpack_node_flag_error(node, mpack_error_type);
    return 0;
}

unsigned int mpack_node_uint(mpack_node_t node) {

    // This should be true at compile-time, so this just wraps the 32-bit function.
    if (sizeof(unsigned int) == 4)
        return (unsigned int)mpack_node_u32(node);

    // Otherwise we use u64 and check the range.
    uint64_t val = mpack_node_u64(node);
    if (val <= MPACK_UINT_MAX)
        return (unsigned int)val;

    mpack_node_flag_error(node, mpack_error_type);
    return 0;
}

int mpack_node_int(mpack_node_t node) {

    // This should be true at compile-time, so this just wraps the 32-bit function.
    if (sizeof(int) == 4)
        return (int)mpack_node_i32(node);

    // Otherwise we use i64 and check the range.
    int64_t val = mpack_node_i64(node);
    if (val >= MPACK_INT_MIN && val <= MPACK_INT_MAX)
        return (int)val;

    mpack_node_flag_error(node, mpack_error_type);
    return 0;
}

#if MPACK_FLOAT
float mpack_node_float(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0.0f;

    if (node.data->type == mpack_type_uint)
        return (float)node.data->value.u;
    if (node.data->type == mpack_type_int)
        return (float)node.data->value.i;
    if (node.data->type == mpack_type_float)
        return node.data->value.f;

    if (node.data->type == mpack_type_double) {
        #if MPACK_DOUBLE
        return (float)node.data->value.d;
        #else
        return mpack_shorten_raw_double_to_float(node.data->value.d);
        #endif
    }

    mpack_node_flag_error(node, mpack_error_type);
    return 0.0f;
}
#endif

#if MPACK_DOUBLE
double mpack_node_double(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0.0;

    if (node.data->type == mpack_type_uint)
        return (double)node.data->value.u;
    else if (node.data->type == mpack_type_int)
        return (double)node.data->value.i;
    else if (node.data->type == mpack_type_float)
        return (double)node.data->value.f;
    else if (node.data->type == mpack_type_double)
        return node.data->value.d;

    mpack_node_flag_error(node, mpack_error_type);
    return 0.0;
}
#endif

#if MPACK_FLOAT
float mpack_node_float_strict(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0.0f;

    if (node.data->type == mpack_type_float)
        return node.data->value.f;

    mpack_node_flag_error(node, mpack_error_type);
    return 0.0f;
}
#endif

#if MPACK_DOUBLE
double mpack_node_double_strict(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0.0;

    if (node.data->type == mpack_type_float)
        return (double)node.data->value.f;
    else if (node.data->type == mpack_type_double)
        return node.data->value.d;

    mpack_node_flag_error(node, mpack_error_type);
    return 0.0;
}
#endif

#if !MPACK_FLOAT
uint32_t mpack_node_raw_float(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    if (node.data->type == mpack_type_float)
        return node.data->value.f;

    mpack_node_flag_error(node, mpack_error_type);
    return 0;
}
#endif

#if !MPACK_DOUBLE
uint64_t mpack_node_raw_double(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    if (node.data->type == mpack_type_double)
        return node.data->value.d;

    mpack_node_flag_error(node, mpack_error_type);
    return 0;
}
#endif

#if MPACK_EXTENSIONS
int8_t mpack_node_exttype(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    if (node.data->type == mpack_type_ext)
        return mpack_node_exttype_unchecked(node);

    mpack_node_flag_error(node, mpack_error_type);
    return 0;
}
#endif

uint32_t mpack_node_data_len(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    mpack_type_t type = node.data->type;
    if (type == mpack_type_str || type == mpack_type_bin
            #if MPACK_EXTENSIONS
            || type == mpack_type_ext
            #endif
            )
        return (uint32_t)node.data->len;

    mpack_node_flag_error(node, mpack_error_type);
    return 0;
}

size_t mpack_node_strlen(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    if (node.data->type == mpack_type_str)
        return (size_t)node.data->len;

    mpack_node_flag_error(node, mpack_error_type);
    return 0;
}

const char* mpack_node_str(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return NULL;

    mpack_type_t type = node.data->type;
    if (type == mpack_type_str)
        return mpack_node_data_unchecked(node);

    mpack_node_flag_error(node, mpack_error_type);
    return NULL;
}

const char* mpack_node_data(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return NULL;

    mpack_type_t type = node.data->type;
    if (type == mpack_type_str || type == mpack_type_bin
            #if MPACK_EXTENSIONS
            || type == mpack_type_ext
            #endif
            )
        return mpack_node_data_unchecked(node);

    mpack_node_flag_error(node, mpack_error_type);
    return NULL;
}

const char* mpack_node_bin_data(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return NULL;

    if (node.data->type == mpack_type_bin)
        return mpack_node_data_unchecked(node);

    mpack_node_flag_error(node, mpack_error_type);
    return NULL;
}

size_t mpack_node_bin_size(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    if (node.data->type == mpack_type_bin)
        return (size_t)node.data->len;

    mpack_node_flag_error(node, mpack_error_type);
    return 0;
}

size_t mpack_node_array_length(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    if (node.data->type != mpack_type_array) {
        mpack_node_flag_error(node, mpack_error_type);
        return 0;
    }

    return (size_t)node.data->len;
}

mpack_node_t mpack_node_array_at(mpack_node_t node, size_t index) {
    if (mpack_node_error(node) != mpack_ok)
        return mpack_tree_nil_node(node.tree);

    if (node.data->type != mpack_type_array) {
        mpack_node_flag_error(node, mpack_error_type);
        return mpack_tree_nil_node(node.tree);
    }

    if (index >= node.data->len) {
        mpack_node_flag_error(node, mpack_error_data);
        return mpack_tree_nil_node(node.tree);
    }

    return mpack_node(node.tree, mpack_node_child(node, index));
}

size_t mpack_node_map_count(mpack_node_t node) {
    if (mpack_node_error(node) != mpack_ok)
        return 0;

    if (node.data->type != mpack_type_map) {
        mpack_node_flag_error(node, mpack_error_type);
        return 0;
    }

    return node.data->len;
}

// internal node map lookup
static mpack_node_t mpack_node_map_at(mpack_node_t node, size_t index, size_t offset) {
    if (mpack_node_error(node) != mpack_ok)
        return mpack_tree_nil_node(node.tree);

    if (node.data->type != mpack_type_map) {
        mpack_node_flag_error(node, mpack_error_type);
        return mpack_tree_nil_node(node.tree);
    }

    if (index >= node.data->len) {
        mpack_node_flag_error(node, mpack_error_data);
        return mpack_tree_nil_node(node.tree);
    }

    return mpack_node(node.tree, mpack_node_child(node, index * 2 + offset));
}

mpack_node_t mpack_node_map_key_at(mpack_node_t node, size_t index) {
    return mpack_node_map_at(node, index, 0);
}

mpack_node_t mpack_node_map_value_at(mpack_node_t node, size_t index) {
    return mpack_node_map_at(node, index, 1);
}

#endif

MPACK_SILENCE_WARNINGS_END
// ended inlining mpack.c 
// start inlining pcg_basic.c 
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
 *       http://www.pcg-random.org
 */

/*
 * This code is derived from the full C implementation, which is in turn
 * derived from the canonical C++ PCG implementation. The C++ version
 * has many additional features and is preferable if you can use C++ in
 * your project.
 */

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

// state for global RNGs

static pcg32_random_t pcg32_global = PCG32_INITIALIZER;

// pcg32_srandom(initstate, initseq)
// pcg32_srandom_r(rng, initstate, initseq):
//     Seed the rng.  Specified in two parts, state initializer and a
//     sequence selection constant (a.k.a. stream id)

void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq)
{
    rng->state = 0U;
    rng->inc = (initseq << 1u) | 1u;
    pcg32_random_r(rng);
    rng->state += initstate;
    pcg32_random_r(rng);
}

void pcg32_srandom(uint64_t seed, uint64_t seq)
{
    pcg32_srandom_r(&pcg32_global, seed, seq);
}

// pcg32_random()
// pcg32_random_r(rng)
//     Generate a uniformly distributed 32-bit random number

uint32_t pcg32_random_r(pcg32_random_t* rng)
{
    uint64_t oldstate = rng->state;
    rng->state = oldstate * 6364136223846793005ULL + rng->inc;
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

uint32_t pcg32_random()
{
    return pcg32_random_r(&pcg32_global);
}


// pcg32_boundedrand(bound):
// pcg32_boundedrand_r(rng, bound):
//     Generate a uniformly distributed number, r, where 0 <= r < bound

uint32_t pcg32_boundedrand_r(pcg32_random_t* rng, uint32_t bound)
{
    // To avoid bias, we need to make the range of the RNG a multiple of
    // bound, which we do by dropping output less than a threshold.
    // A naive scheme to calculate the threshold would be to do
    //
    //     uint32_t threshold = 0x100000000ull % bound;
    //
    // but 64-bit div/mod is slower than 32-bit div/mod (especially on
    // 32-bit platforms).  In essence, we do
    //
    //     uint32_t threshold = (0x100000000ull-bound) % bound;
    //
    // because this version will calculate the same modulus, but the LHS
    // value is less than 2^32.

    uint32_t threshold = -bound % bound;

    // Uniformity guarantees that this loop will terminate.  In practice, it
    // should usually terminate quickly; on average (assuming all bounds are
    // equally likely), 82.25% of the time, we can expect it to require just
    // one iteration.  In the worst case, someone passes a bound of 2^31 + 1
    // (i.e., 2147483649), which invalidates almost 50% of the range.  In 
    // practice, bounds are typically small and only a tiny amount of the range
    // is eliminated.
    for (;;) {
        uint32_t r = pcg32_random_r(rng);
        if (r >= threshold)
            return r % bound;
    }
}


uint32_t pcg32_boundedrand(uint32_t bound)
{
    return pcg32_boundedrand_r(&pcg32_global, bound);
}

// ended inlining pcg_basic.c 
// start inlining thpool.c 
/* ********************************
 * Author:       Johan Hanssen Seferidis
 * License:	     MIT
 * Description:  Library providing a threading pool where you can add
 *               work. For usage, check the thpool.h file or README.md
 *
 *//** @file thpool.h *//*
 *
 ********************************/

#define _POSIX_C_SOURCE 200809L
#include <unistd.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <errno.h>
#include <time.h>
#if defined(__linux__)
#include <sys/prctl.h>
#endif

// start inlining thpool.h 
/**********************************
 * @author      Johan Hanssen Seferidis
 * License:     MIT
 *
 **********************************/

#ifndef _THPOOL_
#define _THPOOL_

#ifdef __cplusplus
extern "C" {
#endif

/* =================================== API ======================================= */


typedef struct thpool_* threadpool;


/**
 * @brief  Initialize threadpool
 *
 * Initializes a threadpool. This function will not return until all
 * threads have initialized successfully.
 *
 * @example
 *
 *    ..
 *    threadpool thpool;                     //First we declare a threadpool
 *    thpool = thpool_init(4);               //then we initialize it to 4 threads
 *    ..
 *
 * @param  num_threads   number of threads to be created in the threadpool
 * @return threadpool    created threadpool on success,
 *                       NULL on error
 */
threadpool thpool_init(int num_threads);


/**
 * @brief Add work to the job queue
 *
 * Takes an action and its argument and adds it to the threadpool's job queue.
 * If you want to add to work a function with more than one arguments then
 * a way to implement this is by passing a pointer to a structure.
 *
 * NOTICE: You have to cast both the function and argument to not get warnings.
 *
 * @example
 *
 *    void print_num(int num){
 *       printf("%d\n", num);
 *    }
 *
 *    int main() {
 *       ..
 *       int a = 10;
 *       thpool_add_work(thpool, (void*)print_num, (void*)a);
 *       ..
 *    }
 *
 * @param  threadpool    threadpool to which the work will be added
 * @param  function_p    pointer to function to add as work
 * @param  arg_p         pointer to an argument
 * @return 0 on success, -1 otherwise.
 */
int thpool_add_work(threadpool, void (*function_p)(void*), void* arg_p);


/**
 * @brief Wait for all queued jobs to finish
 *
 * Will wait for all jobs - both queued and currently running to finish.
 * Once the queue is empty and all work has completed, the calling thread
 * (probably the main program) will continue.
 *
 * Smart polling is used in wait. The polling is initially 0 - meaning that
 * there is virtually no polling at all. If after 1 seconds the threads
 * haven't finished, the polling interval starts growing exponentially
 * until it reaches max_secs seconds. Then it jumps down to a maximum polling
 * interval assuming that heavy processing is being used in the threadpool.
 *
 * @example
 *
 *    ..
 *    threadpool thpool = thpool_init(4);
 *    ..
 *    // Add a bunch of work
 *    ..
 *    thpool_wait(thpool);
 *    puts("All added work has finished");
 *    ..
 *
 * @param threadpool     the threadpool to wait for
 * @return nothing
 */
void thpool_wait(threadpool);


/**
 * @brief Pauses all threads immediately
 *
 * The threads will be paused no matter if they are idle or working.
 * The threads return to their previous states once thpool_resume
 * is called.
 *
 * While the thread is being paused, new work can be added.
 *
 * @example
 *
 *    threadpool thpool = thpool_init(4);
 *    thpool_pause(thpool);
 *    ..
 *    // Add a bunch of work
 *    ..
 *    thpool_resume(thpool); // Let the threads start their magic
 *
 * @param threadpool    the threadpool where the threads should be paused
 * @return nothing
 */
void thpool_pause(threadpool);


/**
 * @brief Unpauses all threads if they are paused
 *
 * @example
 *    ..
 *    thpool_pause(thpool);
 *    sleep(10);              // Delay execution 10 seconds
 *    thpool_resume(thpool);
 *    ..
 *
 * @param threadpool     the threadpool where the threads should be unpaused
 * @return nothing
 */
void thpool_resume(threadpool);


/**
 * @brief Destroy the threadpool
 *
 * This will wait for the currently active threads to finish and then 'kill'
 * the whole threadpool to free up memory.
 *
 * @example
 * int main() {
 *    threadpool thpool1 = thpool_init(2);
 *    threadpool thpool2 = thpool_init(2);
 *    ..
 *    thpool_destroy(thpool1);
 *    ..
 *    return 0;
 * }
 *
 * @param threadpool     the threadpool to destroy
 * @return nothing
 */
void thpool_destroy(threadpool);


/**
 * @brief Show currently working threads
 *
 * Working threads are the threads that are performing work (not idle).
 *
 * @example
 * int main() {
 *    threadpool thpool1 = thpool_init(2);
 *    threadpool thpool2 = thpool_init(2);
 *    ..
 *    printf("Working threads: %d\n", thpool_num_threads_working(thpool1));
 *    ..
 *    return 0;
 * }
 *
 * @param threadpool     the threadpool of interest
 * @return integer       number of threads working
 */
int thpool_num_threads_working(threadpool);

/**
 * @brief Return number of threads alive
 *
 * Total alive threads are maximum number of threads that can perform work.
 *
 * @param threadpool     the threadpool of interest
 * @return integer       number of threads alive
 */
int thpool_num_threads_alive(threadpool);    


#ifdef __cplusplus
}
#endif

#endif
// ended inlining thpool.h 

#ifdef THPOOL_DEBUG
#define THPOOL_DEBUG 1
#else
#define THPOOL_DEBUG 0
#endif

#if !defined(DISABLE_PRINT) || defined(THPOOL_DEBUG)
#define err(str) fprintf(stderr, str)
#else
#define err(str)
#endif

static volatile int threads_keepalive;
static volatile int threads_on_hold;



/* ========================== STRUCTURES ============================ */


/* Binary semaphore */
typedef struct bsem {
	pthread_mutex_t mutex;
	pthread_cond_t   cond;
	int v;
} bsem;


/* Job */
typedef struct job{
	struct job*  prev;                   /* pointer to previous job   */
	void   (*function)(void* arg);       /* function pointer          */
	void*  arg;                          /* function's argument       */
} job;


/* Job queue */
typedef struct jobqueue{
	pthread_mutex_t rwmutex;             /* used for queue r/w access */
	job  *front;                         /* pointer to front of queue */
	job  *rear;                          /* pointer to rear  of queue */
	bsem *has_jobs;                      /* flag as binary semaphore  */
	int   len;                           /* number of jobs in queue   */
} jobqueue;


/* Thread */
typedef struct thread{
	int       id;                        /* friendly id               */
	pthread_t pthread;                   /* pointer to actual thread  */
	struct thpool_* thpool_p;            /* access to thpool          */
} thread;


/* Threadpool */
typedef struct thpool_{
	thread**   threads;                  /* pointer to threads        */
	volatile int num_threads_alive;      /* threads currently alive   */
	volatile int num_threads_working;    /* threads currently working */
	pthread_mutex_t  thcount_lock;       /* used for thread count etc */
	pthread_cond_t  threads_all_idle;    /* signal to thpool_wait     */
	jobqueue  jobqueue;                  /* job queue                 */
} thpool_;





/* ========================== PROTOTYPES ============================ */


static int  thread_init(thpool_* thpool_p, struct thread** thread_p, int id);
static void* thread_do(struct thread* thread_p);
static void  thread_hold(int sig_id);
static void  thread_destroy(struct thread* thread_p);

static int   jobqueue_init(jobqueue* jobqueue_p);
static void  jobqueue_clear(jobqueue* jobqueue_p);
static void  jobqueue_push(jobqueue* jobqueue_p, struct job* newjob_p);
static struct job* jobqueue_pull(jobqueue* jobqueue_p);
static void  jobqueue_destroy(jobqueue* jobqueue_p);

static void  bsem_init(struct bsem *bsem_p, int value);
static void  bsem_reset(struct bsem *bsem_p);
static void  bsem_post(struct bsem *bsem_p);
static void  bsem_post_all(struct bsem *bsem_p);
static void  bsem_wait(struct bsem *bsem_p);

#if defined(__APPLE__) && defined(__MACH__)
// OSX used BSD pthreads and the signature on it is not the same as
// Linux method with the same name.
void pthread_setname_np(const char *name);
#endif



/* ========================== THREADPOOL ============================ */


/* Initialise thread pool */
struct thpool_* thpool_init(int num_threads){

	threads_on_hold   = 0;
	threads_keepalive = 1;

	if (num_threads < 0){
		num_threads = 0;
	}

	/* Make new thread pool */
	thpool_* thpool_p;
	thpool_p = (struct thpool_*)malloc(sizeof(struct thpool_));
	if (thpool_p == NULL){
		err("thpool_init(): Could not allocate memory for thread pool\n");
		return NULL;
	}
	thpool_p->num_threads_alive   = 0;
	thpool_p->num_threads_working = 0;

	/* Initialise the job queue */
	if (jobqueue_init(&thpool_p->jobqueue) == -1){
		err("thpool_init(): Could not allocate memory for job queue\n");
		free(thpool_p);
		return NULL;
	}

	/* Make threads in pool */
	thpool_p->threads = (struct thread**)malloc(num_threads * sizeof(struct thread *));
	if (thpool_p->threads == NULL){
		err("thpool_init(): Could not allocate memory for threads\n");
		jobqueue_destroy(&thpool_p->jobqueue);
		free(thpool_p);
		return NULL;
	}

	pthread_mutex_init(&(thpool_p->thcount_lock), NULL);
	pthread_cond_init(&thpool_p->threads_all_idle, NULL);

	/* Thread init */
	int n;
	for (n=0; n<num_threads; n++){
		thread_init(thpool_p, &thpool_p->threads[n], n);
#if THPOOL_DEBUG
			printf("THPOOL_DEBUG: Created thread %d in pool \n", n);
#endif
	}

	/* Wait for threads to initialize */
	while (thpool_p->num_threads_alive != num_threads) {}

	return thpool_p;
}


/* Add work to the thread pool */
int thpool_add_work(thpool_* thpool_p, void (*function_p)(void*), void* arg_p){
	job* newjob;

	newjob=(struct job*)malloc(sizeof(struct job));
	if (newjob==NULL){
		err("thpool_add_work(): Could not allocate memory for new job\n");
		return -1;
	}

	/* add function and argument */
	newjob->function=function_p;
	newjob->arg=arg_p;

	/* add job to queue */
	jobqueue_push(&thpool_p->jobqueue, newjob);

	return 0;
}


/* Wait until all jobs have finished */
void thpool_wait(thpool_* thpool_p){
	pthread_mutex_lock(&thpool_p->thcount_lock);
	while (thpool_p->jobqueue.len || thpool_p->num_threads_working) {
		pthread_cond_wait(&thpool_p->threads_all_idle, &thpool_p->thcount_lock);
	}
	pthread_mutex_unlock(&thpool_p->thcount_lock);
}


/* Destroy the threadpool */
void thpool_destroy(thpool_* thpool_p){
	/* No need to destory if it's NULL */
	if (thpool_p == NULL) return ;

	volatile int threads_total = thpool_p->num_threads_alive;

	/* End each thread 's infinite loop */
	threads_keepalive = 0;

	/* Give one second to kill idle threads */
	double TIMEOUT = 1.0;
	time_t start, end;
	double tpassed = 0.0;
	time (&start);
	while (tpassed < TIMEOUT && thpool_p->num_threads_alive){
		bsem_post_all(thpool_p->jobqueue.has_jobs);
		time (&end);
		tpassed = difftime(end,start);
	}

	/* Poll remaining threads */
	while (thpool_p->num_threads_alive){
		bsem_post_all(thpool_p->jobqueue.has_jobs);
		sleep(1);
	}

	/* Job queue cleanup */
	jobqueue_destroy(&thpool_p->jobqueue);
	/* Deallocs */
	int n;
	for (n=0; n < threads_total; n++){
		thread_destroy(thpool_p->threads[n]);
	}
	free(thpool_p->threads);
	free(thpool_p);
}


/* Pause all threads in threadpool */
void thpool_pause(thpool_* thpool_p) {
	int n;
	for (n=0; n < thpool_p->num_threads_alive; n++){
		pthread_kill(thpool_p->threads[n]->pthread, SIGUSR1);
	}
}


/* Resume all threads in threadpool */
void thpool_resume(thpool_* thpool_p) {
    // resuming a single threadpool hasn't been
    // implemented yet, meanwhile this supresses
    // the warnings
    (void)thpool_p;

	threads_on_hold = 0;
}


int thpool_num_threads_working(thpool_* thpool_p){
	return thpool_p->num_threads_working;
}

int thpool_num_threads_alive(thpool_* thpool_p){
	return thpool_p->num_threads_alive;
}





/* ============================ THREAD ============================== */


/* Initialize a thread in the thread pool
 *
 * @param thread        address to the pointer of the thread to be created
 * @param id            id to be given to the thread
 * @return 0 on success, -1 otherwise.
 */
static int thread_init (thpool_* thpool_p, struct thread** thread_p, int id){

	*thread_p = (struct thread*)malloc(sizeof(struct thread));
	if (*thread_p == NULL){
		err("thread_init(): Could not allocate memory for thread\n");
		return -1;
	}

	(*thread_p)->thpool_p = thpool_p;
	(*thread_p)->id       = id;

	pthread_create(&(*thread_p)->pthread, NULL, (void * (*)(void *)) thread_do, (*thread_p));
	pthread_detach((*thread_p)->pthread);
	return 0;
}


/* Sets the calling thread on hold */
static void thread_hold(int sig_id) {
    (void)sig_id;
	threads_on_hold = 1;
	while (threads_on_hold){
		sleep(1);
	}
}


/* What each thread is doing
*
* In principle this is an endless loop. The only time this loop gets interuppted is once
* thpool_destroy() is invoked or the program exits.
*
* @param  thread        thread that will run this function
* @return nothing
*/
static void* thread_do(struct thread* thread_p){

	/* Set thread name for profiling and debuging */
	char thread_name[32] = {0};
	snprintf(thread_name, 32, "thread-pool-%d", thread_p->id);

#if defined(__linux__)
	/* Use prctl instead to prevent using _GNU_SOURCE flag and implicit declaration */
	prctl(PR_SET_NAME, thread_name);
#elif defined(__APPLE__) && defined(__MACH__)
	pthread_setname_np(thread_name);
#else
	err("thread_do(): pthread_setname_np is not supported on this system");
#endif

	/* Assure all threads have been created before starting serving */
	thpool_* thpool_p = thread_p->thpool_p;

	/* Register signal handler */
	struct sigaction act;
	sigemptyset(&act.sa_mask);
	act.sa_flags = 0;
	act.sa_handler = thread_hold;
	if (sigaction(SIGUSR1, &act, NULL) == -1) {
		err("thread_do(): cannot handle SIGUSR1");
	}

	/* Mark thread as alive (initialized) */
	pthread_mutex_lock(&thpool_p->thcount_lock);
	thpool_p->num_threads_alive += 1;
	pthread_mutex_unlock(&thpool_p->thcount_lock);

	while(threads_keepalive){

		bsem_wait(thpool_p->jobqueue.has_jobs);

		if (threads_keepalive){

			pthread_mutex_lock(&thpool_p->thcount_lock);
			thpool_p->num_threads_working++;
			pthread_mutex_unlock(&thpool_p->thcount_lock);

			/* Read job from queue and execute it */
			void (*func_buff)(void*);
			void*  arg_buff;
			job* job_p = jobqueue_pull(&thpool_p->jobqueue);
			if (job_p) {
				func_buff = job_p->function;
				arg_buff  = job_p->arg;
				func_buff(arg_buff);
				free(job_p);
			}

			pthread_mutex_lock(&thpool_p->thcount_lock);
			thpool_p->num_threads_working--;
			if (!thpool_p->num_threads_working) {
				pthread_cond_signal(&thpool_p->threads_all_idle);
			}
			pthread_mutex_unlock(&thpool_p->thcount_lock);

		}
	}
	pthread_mutex_lock(&thpool_p->thcount_lock);
	thpool_p->num_threads_alive --;
	pthread_mutex_unlock(&thpool_p->thcount_lock);

	return NULL;
}


/* Frees a thread  */
static void thread_destroy (thread* thread_p){
	free(thread_p);
}





/* ============================ JOB QUEUE =========================== */


/* Initialize queue */
static int jobqueue_init(jobqueue* jobqueue_p){
	jobqueue_p->len = 0;
	jobqueue_p->front = NULL;
	jobqueue_p->rear  = NULL;

	jobqueue_p->has_jobs = (struct bsem*)malloc(sizeof(struct bsem));
	if (jobqueue_p->has_jobs == NULL){
		return -1;
	}

	pthread_mutex_init(&(jobqueue_p->rwmutex), NULL);
	bsem_init(jobqueue_p->has_jobs, 0);

	return 0;
}


/* Clear the queue */
static void jobqueue_clear(jobqueue* jobqueue_p){

	while(jobqueue_p->len){
		free(jobqueue_pull(jobqueue_p));
	}

	jobqueue_p->front = NULL;
	jobqueue_p->rear  = NULL;
	bsem_reset(jobqueue_p->has_jobs);
	jobqueue_p->len = 0;

}


/* Add (allocated) job to queue
 */
static void jobqueue_push(jobqueue* jobqueue_p, struct job* newjob){

	pthread_mutex_lock(&jobqueue_p->rwmutex);
	newjob->prev = NULL;

	switch(jobqueue_p->len){

		case 0:  /* if no jobs in queue */
					jobqueue_p->front = newjob;
					jobqueue_p->rear  = newjob;
					break;

		default: /* if jobs in queue */
					jobqueue_p->rear->prev = newjob;
					jobqueue_p->rear = newjob;

	}
	jobqueue_p->len++;

	bsem_post(jobqueue_p->has_jobs);
	pthread_mutex_unlock(&jobqueue_p->rwmutex);
}


/* Get first job from queue(removes it from queue)
 * Notice: Caller MUST hold a mutex
 */
static struct job* jobqueue_pull(jobqueue* jobqueue_p){

	pthread_mutex_lock(&jobqueue_p->rwmutex);
	job* job_p = jobqueue_p->front;

	switch(jobqueue_p->len){

		case 0:  /* if no jobs in queue */
		  			break;

		case 1:  /* if one job in queue */
					jobqueue_p->front = NULL;
					jobqueue_p->rear  = NULL;
					jobqueue_p->len = 0;
					break;

		default: /* if >1 jobs in queue */
					jobqueue_p->front = job_p->prev;
					jobqueue_p->len--;
					/* more than one job in queue -> post it */
					bsem_post(jobqueue_p->has_jobs);

	}

	pthread_mutex_unlock(&jobqueue_p->rwmutex);
	return job_p;
}


/* Free all queue resources back to the system */
static void jobqueue_destroy(jobqueue* jobqueue_p){
	jobqueue_clear(jobqueue_p);
	free(jobqueue_p->has_jobs);
}





/* ======================== SYNCHRONISATION ========================= */


/* Init semaphore to 1 or 0 */
static void bsem_init(bsem *bsem_p, int value) {
	if (value < 0 || value > 1) {
		err("bsem_init(): Binary semaphore can take only values 1 or 0");
		exit(1);
	}
	pthread_mutex_init(&(bsem_p->mutex), NULL);
	pthread_cond_init(&(bsem_p->cond), NULL);
	bsem_p->v = value;
}


/* Reset semaphore to 0 */
static void bsem_reset(bsem *bsem_p) {
	bsem_init(bsem_p, 0);
}


/* Post to at least one thread */
static void bsem_post(bsem *bsem_p) {
	pthread_mutex_lock(&bsem_p->mutex);
	bsem_p->v = 1;
	pthread_cond_signal(&bsem_p->cond);
	pthread_mutex_unlock(&bsem_p->mutex);
}


/* Post to all threads */
static void bsem_post_all(bsem *bsem_p) {
	pthread_mutex_lock(&bsem_p->mutex);
	bsem_p->v = 1;
	pthread_cond_broadcast(&bsem_p->cond);
	pthread_mutex_unlock(&bsem_p->mutex);
}


/* Wait on semaphore until semaphore has value 0 */
static void bsem_wait(bsem* bsem_p) {
	pthread_mutex_lock(&bsem_p->mutex);
	while (bsem_p->v != 1) {
		pthread_cond_wait(&bsem_p->cond, &bsem_p->mutex);
	}
	bsem_p->v = 0;
	pthread_mutex_unlock(&bsem_p->mutex);
}
// ended inlining thpool.c 

// start inlining alloc.c 
#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// start inlining gkyl_alloc.h 

// start inlining gkyl_util.h 

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// random number generator
// skipping file: pcg_basic.h 

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
// skipping file: gkyl_util.h 

// Output command for use in debugging memory usage
#define GKYL_MEMMSG(fmt, ...) do {              \
      if (gkyl_mem_debug)                       \
        fprintf(stderr, fmt, __VA_ARGS__);      \
  } while (0);

// Output command for use in debugging CUDA memory usage
#define GKYL_CU_MEMMSG(fmt, ...) do {           \
      if (gkyl_cu_dev_mem_debug)                \
        fprintf(stderr, fmt, __VA_ARGS__);      \
    } while (0);

// by default, do not print memory allocation traces
static bool gkyl_mem_debug = false;
static bool gkyl_cu_dev_mem_debug = false;

void gkyl_mem_debug_set(bool flag)
{
  gkyl_mem_debug = flag;
}

void gkyl_cu_dev_mem_debug_set(bool flag)
{
  gkyl_cu_dev_mem_debug = flag;
}

// Compute first 'align' boundary after 'num'
#define align_up(num, align)                    \
    (((num) + ((align) - 1)) & ~((align) - 1))

static const size_t PTR_OFFSET_SZ = sizeof(uint16_t);

void*
gkyl_malloc_(const char *file, int line, const char *func, size_t size)
{
  void *mem = malloc(size);
  GKYL_MEMMSG("%p [%zu] 0.malloc: %s %s:%d\n", mem, size, file, func, line);  
  if (0 == mem) gkyl_exit("malloc failed!");
  return mem;
}

void*
gkyl_calloc_(const char *file, int line, const char *func, size_t num, size_t size)
{
  void *mem = calloc(num, size);
  GKYL_MEMMSG("%p [%zu] 0.calloc: %s %s:%d\n", mem, size, file, func, line);
  if (0 == mem) gkyl_exit("calloc failed!");
  return mem;
}

void*
gkyl_realloc_(const char *file, int line, const char *func, void *ptr, size_t new_size)
{
  void *mem = realloc(ptr, new_size);
  GKYL_MEMMSG("%p [%zu] 0.realloc: %s %s:%d\n", mem, new_size, file, func, line);  
  if (0 == mem) gkyl_exit("realloc failed!");
  return mem;
}

void
gkyl_free_(const char *file, int line, const char *func, void *ptr)
{
  GKYL_MEMMSG("%p 1.free: %s %s:%d\n", ptr, file, func, line);
  free(ptr);
}

void*
gkyl_aligned_alloc_(const char *file, int line, const char *func,
  size_t align, size_t size)
{
  void *ptr = 0;
  assert((align & (align-1)) == 0); // power of 2?

  if (align && size) {
    uint32_t hdr_size = PTR_OFFSET_SZ + (align-1);
    void *p = gkyl_calloc(size+hdr_size, 1);
    if (p) {
      ptr = (void *) align_up(((uintptr_t)p + PTR_OFFSET_SZ), align);
      *((uint16_t *)ptr - 1) = (uint16_t) ((uintptr_t) ptr - (uintptr_t) p);
    }
  }

  GKYL_MEMMSG("%p 0.aligned_alloc: %s %s:%d\n", ptr, file, func, line);
  return ptr;
}

void*
gkyl_aligned_realloc_(const char *file, int line, const char *func,
  void *ptr, size_t align, size_t old_sz, size_t new_sz)
{
  void *nptr = gkyl_aligned_alloc(align, new_sz);
  if (0 == nptr) {
    gkyl_exit("aligned_realloc failed!");
    return 0;
  }
  memcpy(nptr, ptr, old_sz < new_sz ? old_sz : new_sz);
  gkyl_aligned_free(ptr);

  GKYL_MEMMSG("%p 0.aligned_realloc: %s %s:%d\n", ptr, file, func, line);
  return nptr;
}

void
gkyl_aligned_free_(const char *file, int line, const char *func,
  void* ptr)
{
  assert(ptr);
  GKYL_MEMMSG("%p 1.aligned_free: %s %s:%d\n", ptr, file, func, line);
    
  uint16_t offset = *((uint16_t *)ptr - 1);
  gkyl_free((uint8_t *)ptr - offset);
}

struct gkyl_mem_buff_tag {
  bool on_gpu; // is this on GPU?
  size_t count; // size of memory in bytes
  char *data; // Allocated memory
};

gkyl_mem_buff
gkyl_mem_buff_new(size_t count)
{
  struct gkyl_mem_buff_tag *mem = gkyl_malloc(sizeof(*mem));
  mem->on_gpu = false;
  mem->count = count;
  mem->data = gkyl_malloc(count);
  return mem;
}

gkyl_mem_buff
gkyl_mem_buff_cu_new(size_t count)
{
  struct gkyl_mem_buff_tag *mem = gkyl_malloc(sizeof(*mem));
  mem->on_gpu = true;
  mem->count = count;
  mem->data = gkyl_cu_malloc(count);
  return mem;
}

gkyl_mem_buff
gkyl_mem_buff_resize(gkyl_mem_buff mem, size_t count)
{
  if (count > mem->count) {
    if (mem->on_gpu) {
      char *data_new = gkyl_cu_malloc(count);
      gkyl_cu_memcpy(data_new, mem->data, mem->count, GKYL_CU_MEMCPY_D2D);
      gkyl_cu_free(mem->data);
      mem->data = data_new;
    }
    else {
      mem->data = gkyl_realloc(mem->data, count);
    }
    mem->count = count;
  }
  return mem;
}

size_t
gkyl_mem_buff_size(gkyl_mem_buff mem)
{
  return mem->count;
}

char*
gkyl_mem_buff_data(gkyl_mem_buff mem)
{
  return mem->data;
}

void
gkyl_mem_buff_release(gkyl_mem_buff mem)
{
  if (mem->on_gpu)
    gkyl_cu_free(mem->data);
  else
    gkyl_free(mem->data);
  
  gkyl_free(mem);
}

// CUDA specific code

#ifdef GKYL_HAVE_CUDA

#include <cuda_runtime.h>

void*
gkyl_cu_malloc_(const char *file, int line, const char *func, size_t size)
{
  void *ptr;
  cudaError_t err = cudaMalloc(&ptr, size);
  if (err != cudaSuccess)
    gkyl_exit("cudaMalloc failed!");
  
  GKYL_CU_MEMMSG("%p 0.cudaMalloc: %s %s:%d\n", ptr, file, func, line);

  return ptr;
}

void*
gkyl_cu_malloc_host_(const char *file, int line, const char *func, size_t size)
{
  // Allocate pinned host memory.
  void *ptr;
  cudaError_t err = cudaMallocHost(&ptr, size);
  if (err != cudaSuccess)
    gkyl_exit("cudaMallocHost failed!");

  GKYL_CU_MEMMSG("%p 0.cudaMallocHost: %s %s:%d\n", ptr, file, func, line);

  return ptr;
}

void
gkyl_cu_free_(const char *file, int line, const char *func, void *ptr)
{
  GKYL_CU_MEMMSG("%p 1.cudaFree: %s %s:%d\n", ptr, file, func, line);
  cudaFree(ptr);
}

void
gkyl_cu_free_host_(const char *file, int line, const char *func, void *ptr)
{
  GKYL_CU_MEMMSG("%p 1.cudaFreeHost: %s %s:%d\n", ptr, file, func, line);
  cudaFreeHost(ptr);
}

void
gkyl_cu_memcpy(void *dst, const void *src, size_t count, enum gkyl_cu_memcpy_kind kind)
{
  cudaError_t err = cudaMemcpy(dst, src, count, kind);
  if (err != cudaSuccess) {
    char str[1024];
    sprintf(str, "\nCUDA error: %s\n", cudaGetErrorString(err));
    gkyl_exit(str);
  }
}

void
gkyl_cu_memcpy_async(void *dst, const void *src, size_t count, enum gkyl_cu_memcpy_kind kind, cudaStream_t stream)
{
  cudaError_t err = cudaMemcpyAsync(dst, src, count, kind, stream);
  if (err != cudaSuccess) {
    char str[1024];
    sprintf(str, "\nCUDA error: %s\n", cudaGetErrorString(err));
    gkyl_exit(str);
  }
}

void
gkyl_cu_memset(void *data, int val, size_t count)
{
  cudaError_t err = cudaMemset(data, val, count);
  if (err != cudaSuccess)
    gkyl_exit("gkyl_cu_memset failed!");
}

#else

// These non-CUDA functions will simply abort. When not using CUDA
// none of these methods should be called at all.

void*
gkyl_cu_malloc_(const char *file, int line, const char *func, size_t size)
{
  assert(false);
  return 0;
}

void*
gkyl_cu_malloc_host_(const char *file, int line, const char *func, size_t size)
{
  assert(false);
  return 0;
}

void
gkyl_cu_free_(const char *file, int line, const char *func, void *ptr)
{
  assert(false);
}

void
gkyl_cu_free_host_(const char *file, int line, const char *func, void *ptr)
{
  assert(false);
}

void
gkyl_cu_memcpy(void *dst, const void *src, size_t count, enum gkyl_cu_memcpy_kind kind)
{
  assert(false);  
}

void
gkyl_cu_memcpy_async(void *dst, const void *src, size_t count, enum gkyl_cu_memcpy_kind kind, int stream)
{
  assert(false);  
}

void
gkyl_cu_memset(void *data, int val, size_t count)
{
  assert(false);  
}

#endif // CUDA specific code
// ended inlining alloc.c 
// start inlining array.c 
#include <assert.h>
#include <stdlib.h>
#include <string.h>

// skipping file: gkyl_alloc.h 
// start inlining gkyl_alloc_flags_priv.h 
// private header to provide some commonly used flags in allocations

#include <stdint.h>

// flags and corresponding bit-masks
enum gkyl_alloc_flags { GKYL_IS_CU_ALLOC, GKYL_IS_ALLOC_ALIGNED, GKYL_IS_ALLOC_EXTERN };
static const uint32_t gkyl_alloc_flags_masks[] =
{ 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

// NV-GPU flags
#define GKYL_SET_CU_ALLOC(flags) ((flags) |= gkyl_alloc_flags_masks[GKYL_IS_CU_ALLOC])
#define GKYL_CLEAR_CU_ALLOC(flags) ((flags) &= ~gkyl_alloc_flags_masks[GKYL_IS_CU_ALLOC])
#define GKYL_IS_CU_ALLOC(flags) (((flags) & gkyl_alloc_flags_masks[GKYL_IS_CU_ALLOC]) != 0)

// Alignment flags
#define GKYL_SET_ALLOC_ALIGNED(flags) ((flags) |= gkyl_alloc_flags_masks[GKYL_IS_ALLOC_ALIGNED])
#define GKYL_CLEAR_ALLOC_ALIGNED(flags) ((flags) &= ~gkyl_alloc_flags_masks[GKYL_IS_ALLOC_ALIGNED])
#define GKYL_IS_ALLOC_ALIGNED(flags) (((flags) & gkyl_alloc_flags_masks[GKYL_IS_ALLOC_ALIGNED]) != 0)

// Flag to indicate if there was no actual allocation but only external allocated memory was used
#define GKYL_SET_ALLOC_EXTERN(flags) ((flags) |= gkyl_alloc_flags_masks[GKYL_IS_ALLOC_EXTERN])
#define GKYL_CLEAR_ALLOC_EXTERN(flags) ((flags) &= ~gkyl_alloc_flags_masks[GKYL_IS_ALLOC_EXTERN])
#define GKYL_IS_ALLOC_EXTERN(flags) (((flags) & gkyl_alloc_flags_masks[GKYL_IS_ALLOC_EXTERN]) != 0)
// ended inlining gkyl_alloc_flags_priv.h 
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
// skipping file: gkyl_util.h 

// undefine this to use non-aligned memory allocations
#define USE_ALIGNED_ALLOC
// alignment boundary is 32 bytes to be compatible with AVX
static const size_t ARRAY_ALIGN_BND = 32;

#define set_arr_dat_zero_ho(arr, data) \
  for (size_t i=0; i<arr->size*arr->ncomp; ++i) data[i] = 0

#define set_arr_dat_zero_dev(arr, data_ho) \
  for (size_t i=0; i<arr->size*arr->ncomp; ++i) data_ho[i] = 0; \
  gkyl_cu_memcpy(arr->data, data_ho, arr->size*arr->esznc, GKYL_CU_MEMCPY_H2D);

static void*
g_array_alloc(size_t num, size_t sz)
{
#ifdef USE_ALIGNED_ALLOC
  return gkyl_aligned_alloc(ARRAY_ALIGN_BND, num*sz);
#else
  return gkyl_calloc(num, sz);
#endif
}
static void
g_array_free(void* ptr)
{
#ifdef USE_ALIGNED_ALLOC  
  gkyl_aligned_free(ptr);
#else
  gkyl_free(ptr);
#endif  
}

// size in bytes for various data-types
static const size_t array_elem_size[] = {
  [GKYL_INT] = sizeof(int),
  [GKYL_FLOAT] = sizeof(float),
  [GKYL_DOUBLE] = sizeof(double),
  [GKYL_USER] = 1,
};

static void
array_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_array *arr = container_of(ref, struct gkyl_array, ref_count);

  if (false == GKYL_IS_ALLOC_EXTERN(arr->flags)) {
    // only free if we allocated memory ourselves
  
    if (GKYL_IS_CU_ALLOC(arr->flags)) {
#ifdef GKYL_HAVE_CUDA 
      cudaStreamDestroy(arr->iostream);
#endif
      gkyl_cu_free(arr->data);
      gkyl_cu_free(arr->on_dev);
    }
    else {
      g_array_free(arr->data);
    }
    
  }
  gkyl_free(arr);
}

// internal method to allocate array
static struct gkyl_array*
array_new(enum gkyl_elem_type type, size_t ncomp, size_t size, bool is_alloc_extern, void *buff)
{
  struct gkyl_array* arr = gkyl_malloc(sizeof(struct gkyl_array));

  arr->type = type;
  arr->elemsz = array_elem_size[type];
  arr->ncomp = ncomp;
  arr->size = size;
  arr->flags = 0;

  if (is_alloc_extern)
    GKYL_SET_ALLOC_EXTERN(arr->flags);
  else
    GKYL_CLEAR_ALLOC_EXTERN(arr->flags);
  
  GKYL_CLEAR_CU_ALLOC(arr->flags);
#ifdef USE_ALIGNED_ALLOC  
  GKYL_SET_ALLOC_ALIGNED(arr->flags);
#else
  GKYL_CLEAR_ALLOC_ALIGNED(arr->flags);
#endif
  
  arr->esznc = arr->elemsz*arr->ncomp;
  arr->data = buff;  
  if (!is_alloc_extern)
    arr->data = g_array_alloc(arr->size, arr->esznc);
  
  arr->ref_count = gkyl_ref_count_init(array_free);

  arr->nthreads = 1;
  arr->nblocks = 1;

  arr->on_dev = arr; // on_dev reference

  if (!is_alloc_extern) {
    // Zero out array elements (not for user-defined type).
    if (type == GKYL_INT) {
      int *dat_p = arr->data;
      set_arr_dat_zero_ho(arr, dat_p);
    }
    else if (type == GKYL_FLOAT) {
      float *dat_p = arr->data;
      set_arr_dat_zero_ho(arr, dat_p);
    }
    else if (type == GKYL_DOUBLE) {
      double *dat_p = arr->data;
      set_arr_dat_zero_ho(arr, dat_p);
    }
  }

  return arr;
}

struct gkyl_array*
gkyl_array_new(enum gkyl_elem_type type, size_t ncomp, size_t size)
{
  return array_new(type, ncomp, size, false, 0);
}

struct gkyl_array*
gkyl_array_new_from_buff(enum gkyl_elem_type type, size_t ncomp, size_t size, void *buff)
{
  return array_new(type, ncomp, size, true, buff);
}

bool
gkyl_array_is_cu_dev(const struct gkyl_array *arr)
{
  return GKYL_IS_CU_ALLOC(arr->flags);  
}

bool
gkyl_array_is_using_buffer(const struct gkyl_array *arr)
{
  return GKYL_IS_ALLOC_EXTERN(arr->flags);
}

struct gkyl_array*
gkyl_array_copy(struct gkyl_array* dest, const struct gkyl_array* src)
{
  assert(dest->esznc == src->esznc);
  
  long ncopy = src->size < dest->size ? src->size : dest->size;

  bool dest_is_cu_dev = gkyl_array_is_cu_dev(dest);
  bool src_is_cu_dev = gkyl_array_is_cu_dev(src);

  if (src_is_cu_dev) {
    // source is on device
    if (dest_is_cu_dev)
      gkyl_cu_memcpy(dest->data, src->data, ncopy*src->esznc, GKYL_CU_MEMCPY_D2D);
    else
      gkyl_cu_memcpy(dest->data, src->data, ncopy*src->esznc, GKYL_CU_MEMCPY_D2H);
  }
  else {
    // source is on host
    if (dest_is_cu_dev)
      gkyl_cu_memcpy(dest->data, src->data, ncopy*src->esznc, GKYL_CU_MEMCPY_H2D);
    else
      memcpy(dest->data, src->data, ncopy*src->esznc);
  }
  
  return dest;
}

struct gkyl_array*
gkyl_array_copy_async(struct gkyl_array* dest, const struct gkyl_array* src)
{
  assert(dest->esznc == src->esznc);
  
  long ncopy = src->size < dest->size ? src->size : dest->size;

  bool dest_is_cu_dev = gkyl_array_is_cu_dev(dest);
  bool src_is_cu_dev = gkyl_array_is_cu_dev(src);

  if (src_is_cu_dev) {
    // source is on device
    if (dest_is_cu_dev)
      gkyl_cu_memcpy_async(dest->data, src->data, ncopy*src->esznc, GKYL_CU_MEMCPY_D2D, src->iostream);
    else
      gkyl_cu_memcpy_async(dest->data, src->data, ncopy*src->esznc, GKYL_CU_MEMCPY_D2H, src->iostream);
  }
  else {
    // source is on host
    if (dest_is_cu_dev)
      gkyl_cu_memcpy_async(dest->data, src->data, ncopy*src->esznc, GKYL_CU_MEMCPY_H2D, dest->iostream);
    else
      memcpy(dest->data, src->data, ncopy*src->esznc);
  }
  
  return dest;
}

struct gkyl_array*
gkyl_array_clone(const struct gkyl_array* src)
{
  struct gkyl_array* arr = gkyl_malloc(sizeof(struct gkyl_array));

  arr->type = src->type;
  arr->elemsz = src->elemsz;
  arr->ncomp = src->ncomp;
  arr->esznc = src->esznc;
  arr->size = src->size;
  arr->flags = src->flags;

  GKYL_CLEAR_ALLOC_EXTERN(arr->flags);

  if (GKYL_IS_CU_ALLOC(src->flags)) {
    arr->nthreads = src->nthreads;
    arr->nblocks = src->nblocks;
    arr->data = gkyl_cu_malloc(arr->size*arr->esznc);
    arr->on_dev = gkyl_cu_malloc(sizeof(struct gkyl_array));
    gkyl_cu_memcpy(arr->data, src->data, arr->size*arr->esznc, GKYL_CU_MEMCPY_D2D);
    gkyl_cu_memcpy(arr->on_dev, src->on_dev, sizeof(struct gkyl_array), GKYL_CU_MEMCPY_D2D);
    gkyl_cu_memcpy(&((arr->on_dev)->data), &arr->data, sizeof(void*), GKYL_CU_MEMCPY_H2D);
  }
  else {
    arr->data = g_array_alloc(arr->size, arr->esznc);
    memcpy(arr->data, src->data, arr->size*arr->esznc);
  }
  
  arr->ref_count = gkyl_ref_count_init(array_free);
  
  return arr;
}

struct gkyl_array*
gkyl_array_acquire(const struct gkyl_array* arr)
{
  gkyl_ref_count_inc(&arr->ref_count);
  return (struct gkyl_array*) arr;
}

void
gkyl_array_release(const struct gkyl_array* arr)
{
  if (arr)
    gkyl_ref_count_dec(&arr->ref_count);
}

// CUDA specific code

#ifdef GKYL_HAVE_CUDA

struct gkyl_array*
gkyl_array_cu_dev_new(enum gkyl_elem_type type, size_t ncomp, size_t size)
{
  struct gkyl_array* arr = gkyl_malloc(sizeof(struct gkyl_array));

  arr->type = type;
  arr->elemsz = array_elem_size[type];
  arr->ncomp = ncomp;
  arr->size = size;
  arr->flags = 0;

  GKYL_CLEAR_ALLOC_EXTERN(arr->flags);
  GKYL_SET_CU_ALLOC(arr->flags);
  GKYL_CLEAR_ALLOC_ALIGNED(arr->flags);
  
  arr->esznc = arr->elemsz*arr->ncomp;
  arr->ref_count = gkyl_ref_count_init(array_free);
  arr->data = gkyl_cu_malloc(arr->size*arr->esznc);
  arr->nthreads = GKYL_DEFAULT_NUM_THREADS;
  arr->nblocks = gkyl_int_div_up(arr->size*arr->ncomp, arr->nthreads);

  cudaStreamCreate(&arr->iostream);

  // create a clone of the struct arr->on_dev that lives on the device,
  // so that the whole arr->on_dev struct can be passed to a device kernel
  arr->on_dev = gkyl_cu_malloc(sizeof(struct gkyl_array));
  gkyl_cu_memcpy(arr->on_dev, arr, sizeof(struct gkyl_array), GKYL_CU_MEMCPY_H2D);
  // set device-side data pointer in arr->on_dev to arr->data 
  // (which is the host-side pointer to the device data)
  gkyl_cu_memcpy(&((arr->on_dev)->data), &arr->data, sizeof(void*), GKYL_CU_MEMCPY_H2D);

  // Zero out array elements (not for user-defined type).
  if (type == GKYL_INT) {
    int *data_ho = gkyl_malloc(arr->size*arr->esznc);
    set_arr_dat_zero_dev(arr, data_ho);
    gkyl_free(data_ho);
  }
  else if (type == GKYL_FLOAT) {
    float *data_ho = gkyl_malloc(arr->size*arr->esznc);
    set_arr_dat_zero_dev(arr, data_ho);
    gkyl_free(data_ho);
  }
  else if (type == GKYL_DOUBLE) {
    double *data_ho = gkyl_malloc(arr->size*arr->esznc);
    set_arr_dat_zero_dev(arr, data_ho);
    gkyl_free(data_ho);
  }

  return arr;
}

struct gkyl_array*
gkyl_array_cu_host_new(enum gkyl_elem_type type, size_t ncomp, size_t size)
{
  struct gkyl_array* arr = gkyl_cu_malloc_host(sizeof(struct gkyl_array));

  arr->type = type;
  arr->elemsz = array_elem_size[type];
  arr->ncomp = ncomp;
  arr->size = size;
  arr->flags = 0;

  GKYL_CLEAR_ALLOC_EXTERN(arr->flags);
  GKYL_CLEAR_CU_ALLOC(arr->flags);
#ifdef USE_ALIGNED_ALLOC  
  GKYL_SET_ALLOC_ALIGNED(arr->flags);
#else
  GKYL_CLEAR_ALLOC_ALIGNED(arr->flags);
#endif
  
  arr->esznc = arr->elemsz*arr->ncomp;
  arr->data = gkyl_cu_malloc_host(arr->size*arr->esznc);
  arr->ref_count = gkyl_ref_count_init(array_free);

  arr->nthreads = 1;
  arr->nblocks = 1;

  arr->on_dev = arr; // on_dev reference
  
  // Zero out array elements (not for user-defined type).
  if (type == GKYL_INT) {
    int *dat_p = arr->data;
    set_arr_dat_zero_ho(arr, dat_p);
  }
  else if (type == GKYL_FLOAT) {
    float *dat_p = arr->data;
    set_arr_dat_zero_ho(arr, dat_p);
  }
  else if (type == GKYL_DOUBLE) {
    double *dat_p = arr->data;
    set_arr_dat_zero_ho(arr, dat_p);
  }

  return arr;
}


#else

struct gkyl_array*
gkyl_array_cu_dev_new(enum gkyl_elem_type type, size_t ncomp, size_t size)
{
  assert(false);
  return 0;
}

struct gkyl_array*
gkyl_array_cu_host_new(enum gkyl_elem_type type, size_t ncomp, size_t size)
{
  assert(false);
  return 0;
}

#endif // CUDA specific code
// ended inlining array.c 
// start inlining array_ops.c 

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

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
// skipping file: gkyl_util.h 
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
// start inlining gkyl_array_ops_priv.h 

// Private header, not for direct use in user code

// Compute number of elements stored in array 'arr'
#define NELM(arr) (arr->size*arr->ncomp)
// Compute size of 'arr' 
#define NSIZE(arr) (arr->size)
// Compute number of components stored in array 'arr'
#define NCOM(arr) (arr->ncomp)

GKYL_CU_DH
static inline void
array_clear1(long n, double *out, double val)
{
  for (int c=0; c<n; ++c)
    out[c] = val;
}

GKYL_CU_DH
static inline void
array_acc1(long n, double * GKYL_RESTRICT out, double a, const double * GKYL_RESTRICT inp)
{
  for (int c=0; c<n; ++c)
    out[c] += a*inp[c];
}

GKYL_CU_DH
static inline void
array_set1(long n,
  double * GKYL_RESTRICT out, double a, const double * GKYL_RESTRICT inp)
{
  for (int c=0; c<n; ++c)
    out[c] = a*inp[c];
}

GKYL_CU_DH
static inline void
array_set2(long n, long m,
  double * GKYL_RESTRICT out, double a, const double * GKYL_RESTRICT inp)
{
  for (int c=0; c<n; ++c)
    out[c] = a*inp[m+c];
}
// ended inlining gkyl_array_ops_priv.h 
// skipping file: gkyl_alloc_flags_priv.h 

bool
gkyl_array_copy_func_is_cu_dev(const struct gkyl_array_copy_func *bc)
{
  return GKYL_IS_CU_ALLOC(bc->flags);
}

struct gkyl_array*
gkyl_array_clear(struct gkyl_array* out, double val)
{
  assert(out->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {gkyl_array_clear_cu(out, val); return out; }
#endif

  double *out_d = out->data;
  for (size_t i=0; i<NELM(out); ++i)
    out_d[i] = val;
  return out;
}

struct gkyl_array*
gkyl_array_accumulate(struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == inp->size && out->elemsz == inp->elemsz);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out) && gkyl_array_is_cu_dev(inp)) { gkyl_array_accumulate_cu(out, a, inp); return out; }
#endif

  double *out_d = out->data;
  const double *inp_d = inp->data;
  for (size_t i=0; i<NELM(out); ++i)
    out_d[i] += a*inp_d[i];
  return out;
}

struct gkyl_array*
gkyl_array_accumulate_offset(struct gkyl_array* out, double a,
  const struct gkyl_array* inp, int coff)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == inp->size);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out) && gkyl_array_is_cu_dev(inp)) { gkyl_array_accumulate_offset_cu(out, a, inp, coff); return out; }
#endif

  double *out_d = out->data;
  const double *inp_d = inp->data;
  if (NCOM(out) < NCOM(inp)) {
    // Interpret offset as offset in input components.
    for (size_t i=0; i<out->size; ++i)
      for (size_t c=0; c<NCOM(out); ++c)
        out_d[i*NCOM(out)+c] += a*inp_d[i*NCOM(inp)+c+coff];
  } else {
    // Interpret offset as offset in output components.
    for (size_t i=0; i<out->size; ++i)
      for (size_t c=0; c<NCOM(inp); ++c)
        out_d[i*NCOM(out)+c+coff] += a*inp_d[i*NCOM(inp)+c];
  }
  return out;
}

struct gkyl_array*
gkyl_array_set(struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == inp->size && out->elemsz == inp->elemsz);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_set_cu(out, a, inp); return out; }
#endif

  double *out_d = out->data;
  const double *inp_d = inp->data;
  for (size_t i=0; i<NELM(out); ++i)
    out_d[i] = a*inp_d[i];
  return out;
}

struct gkyl_array*
gkyl_array_set_offset(struct gkyl_array* out, double a,
  const struct gkyl_array* inp, int coff)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == inp->size);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_set_offset_cu(out, a, inp, coff); return out; }
#endif

  double *out_d = out->data;
  const double *inp_d = inp->data;
  if (NCOM(out) < NCOM(inp)) {
    // Interpret offset as offset in input components.
    for (size_t i=0; i<out->size; ++i)
      for (size_t c=0; c<NCOM(out); ++c)
        out_d[i*NCOM(out)+c] = a*inp_d[i*NCOM(inp)+c+coff];
  } else {
    // Interpret offset as offset in output components.
    for (size_t i=0; i<out->size; ++i)
      for (size_t c=0; c<NCOM(inp); ++c)
        out_d[i*NCOM(out)+c+coff] = a*inp_d[i*NCOM(inp)+c];
  }
  return out;
}

struct gkyl_array*
gkyl_array_scale(struct gkyl_array* out, double a)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_scale_cu(out, a); return out; }
#endif

  return gkyl_array_set(out, a, out);
}

struct gkyl_array*
gkyl_array_scale_by_cell(struct gkyl_array* out, const struct gkyl_array* a)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == a->size && NCOM(a) == 1);
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_scale_by_cell_cu(out, a); return out; }
#endif

  double *out_d = out->data;
  const double *a_d = a->data;
  for (size_t i=0; i<out->size; ++i)
    for (size_t c=0; c<NCOM(out); ++c)
      out_d[i*NCOM(out)+c] = a_d[i]*out_d[i*NCOM(out)+c];
  return out;
}

struct gkyl_array*
gkyl_array_shiftc(struct gkyl_array* out, double a, unsigned k)
{
  assert(out->type == GKYL_DOUBLE);
  assert(k < NCOM(out));
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_shiftc_cu(out, a, k); return out; }
#endif

  double *out_d = out->data;
  for (size_t i=0; i<out->size; ++i)
    out_d[i*NCOM(out)+k] = a+out_d[i*NCOM(out)+k];
  return out;
}

void 
gkyl_array_reduce(double *out, const struct gkyl_array *arr, enum gkyl_array_op op)
{
  assert(arr->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) {
    switch (op) {
      case GKYL_MAX:
        gkyl_array_reduce_max_cu(out, arr);
        break;
      case GKYL_MIN:
        gkyl_array_reduce_min_cu(out, arr);
        break;
      case GKYL_SUM:
        gkyl_array_reduce_sum_cu(out, arr);
        break;
    }
    return;
  }
#endif

  long nc = NCOM(arr);
  double *arr_d = arr->data;

  switch (op) {
    case GKYL_MIN:
      for (long k=0; k<nc; ++k) out[k] = DBL_MAX;
      for (size_t i=0; i<arr->size; ++i) {
        const double *d = gkyl_array_cfetch(arr, i);
        for (long k=0; k<nc; ++k)
          out[k] = fmin(out[k], d[k]);
      }
      break;

    case GKYL_MAX:
      for (long k=0; k<nc; ++k) out[k] = -DBL_MAX;
      for (size_t i=0; i<arr->size; ++i) {
        const double *d = gkyl_array_cfetch(arr, i);
        for (long k=0; k<nc; ++k)
          out[k] = fmax(out[k], d[k]);
      }
      break;

    case GKYL_SUM:
      for (long k=0; k<nc; ++k) out[k] = 0;
      for (size_t i=0; i<arr->size; ++i) {
        const double *d = gkyl_array_cfetch(arr, i);
        for (long k=0; k<nc; ++k)
          out[k] += d[k];
      }
      break;
  }
}

// range based methods
struct gkyl_array*
gkyl_array_clear_range(struct gkyl_array *out, double val, const struct gkyl_range *range)
{
  assert(out->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_clear_range_cu(out, val, range); return out; }
#endif

  long n = NCOM(out);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    array_clear1(n, gkyl_array_fetch(out, start), val);
  }

  return out;
}

struct gkyl_array*
gkyl_array_accumulate_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, const struct gkyl_range *range)
{
  assert(out->type == GKYL_DOUBLE);

  assert(out->size == inp->size);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_accumulate_range_cu(out, a, inp, range); return out; }
#endif

  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n = outnc<inpnc ? outnc : inpnc;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    array_acc1(n,
      gkyl_array_fetch(out, start), a, gkyl_array_cfetch(inp, start));
  }

  return out;
}

struct gkyl_array*
gkyl_array_accumulate_offset_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, int coff, const struct gkyl_range *range)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == inp->size);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_accumulate_offset_range_cu(out, a, inp, coff, range); return out; }
#endif

  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n;
  int outoff, inoff;
  if (outnc < inpnc) {
    n = outnc;
    outoff = 0;
    inoff = coff;
  } else {
    n = inpnc;
    outoff = coff;
    inoff = 0;
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    double *out_d = gkyl_array_fetch(out, start);
    const double *inp_d = gkyl_array_cfetch(inp, start);
    array_acc1(n, out_d+outoff, a, inp_d+inoff);
  }

  return out;
}

struct gkyl_array*
gkyl_array_set_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, const struct gkyl_range *range)
{
  assert(out->type == GKYL_DOUBLE && inp->type == GKYL_DOUBLE);
  assert(out->size == inp->size);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_set_range_cu(out, a, inp, range); return out; }
#endif

  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n = outnc<inpnc ? outnc : inpnc;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    array_set1(n,
      gkyl_array_fetch(out, start), a, gkyl_array_cfetch(inp, start));
  }

  return out;
}

struct gkyl_array*
gkyl_array_set_range_to_range(struct gkyl_array *out, double a,
  const struct gkyl_array *inp, struct gkyl_range *out_range, struct gkyl_range *inp_range)
{
  assert(out->elemsz == inp->elemsz);
  assert((inp_range->volume < 1) || (out_range->volume == inp_range->volume));

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_set_range_to_range_cu(out, a, inp, out_range, inp_range); return out; }
#endif

  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n = outnc<inpnc ? outnc : inpnc;

  // Setup linear counter offset for output range/array.
  int iloLocal_out[GKYL_MAX_DIM], iloLocal_inp[GKYL_MAX_DIM];
  for (int d=0; d<out_range->ndim; ++d){
    iloLocal_out[d] = out_range->lower[d];
    iloLocal_inp[d] = inp_range->lower[d];
  }

  int idx_out[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, inp_range);
  while (gkyl_range_iter_next(&iter)) {
    for (int d=0; d<out_range->ndim; ++d)
      idx_out[d] = iloLocal_out[d] + (iter.idx[d] - iloLocal_inp[d]);

    long linidx_inp = gkyl_range_idx(inp_range, iter.idx);
    long linidx_out = gkyl_range_idx(out_range, idx_out);
    array_set1(n,
      gkyl_array_fetch(out, linidx_out), a, gkyl_array_cfetch(inp, linidx_inp));
  }

  return out;
}

struct gkyl_array*
gkyl_array_set_offset_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, int coff, const struct gkyl_range *range)
{
  assert(out->type == GKYL_DOUBLE && inp->type == GKYL_DOUBLE);
  assert(out->size == inp->size);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_set_offset_range_cu(out, a, inp, coff, range); return out; }
#endif

  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n;
  int outoff, inoff;
  if (outnc < inpnc) {
    n = outnc;
    outoff = 0;
    inoff = coff;
  } else {
    n = inpnc;
    outoff = coff;
    inoff = 0;
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    double *out_d = gkyl_array_fetch(out, start);
    const double *inp_d = gkyl_array_cfetch(inp, start);
    array_set1(n, out_d+outoff, a, inp_d+inoff);
  }

  return out;
}

struct gkyl_array*
gkyl_array_scale_range(struct gkyl_array *out,
  double a, const struct gkyl_range *range)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_scale_range_cu(out, a, range); return out; }
#endif

  return gkyl_array_set_range(out, a, out, range);
}

struct gkyl_array*
gkyl_array_shiftc_range(struct gkyl_array* out, double a, unsigned k, const struct gkyl_range *range)
{
  assert(out->type == GKYL_DOUBLE);
  assert(k < NCOM(out));
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_shiftc_range_cu(out, a, k, range); return out; }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    double *out_d = gkyl_array_fetch(out, start);
    out_d[k] += a;
  }
  return out;
}

void
gkyl_array_reduce_range(double *res,
  const struct gkyl_array *arr, enum gkyl_array_op op, const struct gkyl_range *range)
{
  assert(arr->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) {
    switch (op) {
      case GKYL_MAX:
        gkyl_array_reduce_range_max_cu(res, arr, range);
        break;
      case GKYL_MIN:
        gkyl_array_reduce_range_min_cu(res, arr, range);
        break;
      case GKYL_SUM:
        gkyl_array_reduce_range_sum_cu(res, arr, range);
        break;
    }
    return;
  }
#endif

  long n = NCOM(arr);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  switch (op) {
    case GKYL_MIN:
      for (long i=0; i<n; ++i) res[i] = DBL_MAX;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *d = gkyl_array_cfetch(arr, start);
        for (long i=0; i<n; ++i)
          res[i] = fmin(res[i], d[i]);
      }
      break;
    case GKYL_MAX:
      for (long i=0; i<n; ++i) res[i] = -DBL_MAX;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *d = gkyl_array_cfetch(arr, start);
        for (long i=0; i<n; ++i)
          res[i] = fmax(res[i], d[i]);
      }
      break;
    case GKYL_SUM:
      for (long i=0; i<n; ++i) res[i] = 0;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(range, iter.idx);
        const double *d = gkyl_array_cfetch(arr, start);
        for (long i=0; i<n; ++i)
          res[i] += d[i];
      }
      break;
  }
}

struct gkyl_array*
gkyl_array_copy_range(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *range)
{
  assert(out->size == inp->size && out->elemsz == inp->elemsz);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_copy_range_cu(out, inp, range); return out; }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    memcpy(gkyl_array_fetch(out, start), gkyl_array_cfetch(inp, start), inp->esznc);
  }
  return out;
}

struct gkyl_array*
gkyl_array_copy_range_to_range(struct gkyl_array *out,
  const struct gkyl_array *inp, const struct gkyl_range *out_range, const struct gkyl_range *inp_range)
{
  assert(out->elemsz == inp->elemsz);
  assert((inp_range->volume < 1) || (out_range->volume == inp_range->volume));

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_copy_range_to_range_cu(out, inp, out_range, inp_range); return out; }
#endif

  // Setup linear counter offset for output range/array.
  int iloLocal_out[GKYL_MAX_DIM], iloLocal_inp[GKYL_MAX_DIM];
  for (int d=0; d<out_range->ndim; ++d){
    iloLocal_out[d] = out_range->lower[d];
    iloLocal_inp[d] = inp_range->lower[d];
  }

  int idx_out[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, inp_range);
  while (gkyl_range_iter_next(&iter)) {
    for (int d=0; d<out_range->ndim; ++d)
      idx_out[d] = iloLocal_out[d] + (iter.idx[d] - iloLocal_inp[d]);

    long linidx_inp = gkyl_range_idx(inp_range, iter.idx);
    long linidx_out = gkyl_range_idx(out_range, idx_out);
    memcpy(gkyl_array_fetch(out, linidx_out), gkyl_array_cfetch(inp, linidx_inp), inp->esznc);
  }
  return out;
}

void
gkyl_array_copy_to_buffer(void *data, const struct gkyl_array *arr,
  const struct gkyl_range *range)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) { gkyl_array_copy_to_buffer_cu(data, arr, range); return; }
#endif

#define _F(loc) gkyl_array_cfetch(arr, loc)

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    memcpy(((char*) data) + arr->esznc*count++, _F(start), arr->esznc);
  }

#undef _F
}

void
gkyl_array_copy_from_buffer(struct gkyl_array *arr,
  const void *data, const struct gkyl_range *range)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) { gkyl_array_copy_from_buffer_cu(arr, data, range); return; }
#endif

#define _F(loc) gkyl_array_fetch(arr, loc)

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    memcpy(_F(start), ((char*) data) + arr->esznc*count++, arr->esznc);
  }

#undef _F
}

void
gkyl_array_copy_to_buffer_fn(void *data, const struct gkyl_array *arr,
  const struct gkyl_range *range, struct gkyl_array_copy_func *cf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) { gkyl_array_copy_to_buffer_fn_cu(data, arr, range, cf); return; }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *inp = gkyl_array_cfetch(arr, loc);
    double *out = gkyl_flat_fetch(data, arr->esznc*count);
    cf->func(NCOM(arr), out, inp, cf->ctx);
    count += 1;
  }
}

void
gkyl_array_flip_copy_to_buffer_fn(void *data, const struct gkyl_array *arr,
  int dir, const struct gkyl_range *range, struct gkyl_array_copy_func *cf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) {
    if (gkyl_array_is_cu_dev(arr)) { gkyl_array_flip_copy_to_buffer_fn_cu(data, arr, dir, range, cf); return; }
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  int fidx[GKYL_MAX_DIM]; // flipped index
  struct gkyl_range buff_range;
  gkyl_range_init(&buff_range, range->ndim, range->lower, range->upper);

  int uplo = range->upper[dir]+range->lower[dir];

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    gkyl_copy_int_arr(range->ndim, iter.idx, fidx);
    fidx[dir] = uplo - iter.idx[dir];
    
    long count = gkyl_range_idx(&buff_range, fidx);

    const double *inp = gkyl_array_cfetch(arr, loc);
    double *out = gkyl_flat_fetch(data, arr->esznc*count);
    cf->func(NCOM(arr), out, inp, cf->ctx);
  }
}

static double
calc_rel_diff(double a, double b)
{
  if (isnan(a) || isnan(b)) return DBL_MAX;
  
  double absa = fabs(a), absb = fabs(b), diff = fabs(a-b);
  if (a == b) return 0;
  if (a == 0 || b == 0 || (absa+absb < DBL_MIN)) return diff;
  return diff/fmin(absa+absb, DBL_MAX);
}

struct gkyl_array_diff
gkyl_array_diff(const struct gkyl_array *arr1, const struct gkyl_array *arr2, const struct gkyl_range *range)
{
  struct gkyl_array_diff incompat = {
    .is_compatible = false,
    .max_abs_diff = DBL_MAX,
    .min_abs_diff = DBL_MAX,
    .max_rel_diff = DBL_MAX,
    .min_rel_diff = DBL_MAX
  };

  if ((arr1->type != GKYL_DOUBLE) && (arr2->type != GKYL_DOUBLE))
    return incompat;

  if (gkyl_array_is_cu_dev(arr1) || gkyl_array_is_cu_dev(arr2))
    return incompat;    

  if (arr1->elemsz != arr2->elemsz)
    return incompat;    

  if (arr1->ncomp != arr2->ncomp)
    return incompat;

  if (arr1->size != arr2->size)
    return incompat;

  double max_abs_diff = -DBL_MAX, max_rel_diff = -DBL_MAX;
  double min_abs_diff = DBL_MAX, min_rel_diff = DBL_MAX;
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    
    long loc = gkyl_range_idx(range, iter.idx);
    const double *a1 = gkyl_array_cfetch(arr1, loc);
    const double *a2 = gkyl_array_cfetch(arr2, loc);

    for (int c=0; c<arr1->ncomp; ++c) {
      max_abs_diff = fmax(max_abs_diff, a1[c]-a2[c]);
      min_abs_diff = fmin(min_abs_diff, a1[c]-a2[c]);
      double rel_diff = calc_rel_diff(a1[c], a2[c]);
      max_rel_diff = fmax(max_rel_diff, rel_diff);
      min_rel_diff = fmin(min_rel_diff, rel_diff);
    }
  }

  return (struct gkyl_array_diff) {
    .is_compatible = true,
    .max_abs_diff = max_abs_diff,
    .min_abs_diff = min_abs_diff,
    .max_rel_diff = max_rel_diff,
    .min_rel_diff = min_rel_diff
  };
}
// ended inlining array_ops.c 
// start inlining array_rio.c 
#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

// skipping file: gkyl_alloc.h 
// start inlining gkyl_array_rio.h 

#include <stdio.h>

// skipping file: gkyl_array.h 
// skipping file: gkyl_range.h 
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
// start inlining gkyl_array_rio_priv.h 

// skipping file: gkyl_array.h 
// skipping file: gkyl_range.h 
// skipping file: gkyl_rect_grid.h 

struct gkyl_array_header_info;

/**
 * Write out header and meta data to file.
 *
 * @param hrd Header data.
 * @return Status flag
 */
int gkyl_header_meta_write_fp(const struct gkyl_array_header_info *hdr, FILE *fp);

/**
 * Read header and meta data from file.
 *
 * @param hrd On output, header data
 * @return Status flag
 */
int gkyl_header_meta_read_fp(struct gkyl_array_header_info *hdr, FILE *fp);

/**
 * Write out grid and array data header data to file. Note that only
 * HEADER is written and NOT the array data itself.
 *
 * @param grid Grid object to write
 * @param hrd Header data.
 * @return Status flag: 0 if write succeeded, 'errno' otherwise
 */
int gkyl_grid_sub_array_header_write_fp(const struct gkyl_rect_grid *grid,
  const struct gkyl_array_header_info *hdr, FILE *fp);

/**
 * Read grid and array data header data from file. Note that only
 * HEADER is read and NOT the array data itself.
 * 
 * If the array header contains meta-data then the 'meta' char * is
 * allocated in the @a hdr struct.
 *
 * YOU MUST FREE 'hdr->meta' using gkyl_free or else this will result
 * in a memory leak!
 *
 * @param grid Grid object to read
 * @param hrd On output, Header data.
 * @return Status flag: 0 if read succeeded, 'errno' otherwise
 */
int gkyl_grid_sub_array_header_read_fp(struct gkyl_rect_grid *grid,
  struct gkyl_array_header_info *hdr, FILE *fp);

/**
 * Release memory for header.
 *
 * @param hdr Header memory to release
 */
void gkyl_grid_sub_array_header_release(struct gkyl_array_header_info *hdr);
// ended inlining gkyl_array_rio_priv.h 
// start inlining gkyl_elem_type_priv.h 

// skipping file: gkyl_elem_type.h 


// code for array datatype for use in IO
static const uint64_t gkyl_array_data_type[] = {
  [GKYL_INT] = 0,
  [GKYL_FLOAT] = 1,
  [GKYL_DOUBLE] = 2,
  [GKYL_INT_64] = 3,
  [GKYL_USER] = 32,
};

// mapping of code to datatype: MUST be consistent with the
// gkyl_array_data_type array above
static const int gkyl_array_code_to_data_type[] = {
  [0] = GKYL_INT,
  [1] = GKYL_FLOAT,
  [2] = GKYL_DOUBLE,
  [3] = GKYL_INT_64,
  [32] = GKYL_USER
};    

// size in bytes for various data-types
static const size_t gkyl_elem_type_size[] = {
  [GKYL_INT] = sizeof(int),
  [GKYL_FLOAT] = sizeof(float),
  [GKYL_DOUBLE] = sizeof(double),
  [GKYL_INT_64] = sizeof(int64_t),
  [GKYL_USER] = 1,
};


static const uint64_t gkyl_file_type_int[] = {
  [GKYL_FIELD_DATA_FILE] = 1,
  [GKYL_DYNVEC_DATA_FILE] = 2,
  [GKYL_MULTI_RANGE_DATA_FILE] = 3,
  [GKYL_BLOCK_TOPO_DATA_FILE] = 4,
  [GKYL_MULTI_BLOCK_DATA_FILE] = 5,
};

// ended inlining gkyl_elem_type_priv.h 

// Error message strings
static const char *array_rio_status_msg[] = {
  [GKYL_ARRAY_RIO_SUCCESS] = "Success",
  [GKYL_ARRAY_RIO_BAD_VERSION] = "Incorrect header version",
  [GKYL_ARRAY_RIO_FOPEN_FAILED] = "File open failed",
  [GKYL_ARRAY_RIO_FREAD_FAILED] = "Data read failed",
  [GKYL_ARRAY_RIO_DATA_MISMATCH] = "Data mismatch",
  [GKYL_ARRAY_RIO_META_FAILED] = "Metadata output failed"
};

const char*
gkyl_array_rio_status_msg(enum gkyl_array_rio_status status)
{
  return array_rio_status_msg[status];
}

static void
sub_array_write_priv(const struct gkyl_range *range,
  const struct gkyl_array *arr, FILE *fp)
{
#define _F(loc) gkyl_array_cfetch(arr, loc)

  // construct skip iterator to allow writing (potentially) in chunks
  // rather than element by element or requiring a copy of data
  struct gkyl_range_skip_iter skip;
  gkyl_range_skip_iter_init(&skip, range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &skip.range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&skip.range, iter.idx);
    fwrite(_F(start), arr->esznc*skip.delta, 1, fp);
  }
#undef _F
}

int
gkyl_header_meta_write_fp(const struct gkyl_array_header_info *hdr, FILE *fp)
{
  const char g0[5] = "gkyl0";

  // Version 1 header
  fwrite(g0, sizeof(char[5]), 1, fp);
  uint64_t version = 1;
  fwrite(&version, sizeof(uint64_t), 1, fp);
  fwrite(&hdr->file_type, sizeof(uint64_t), 1, fp);
  uint64_t meta_size = hdr->meta_size;
  fwrite(&meta_size, sizeof(uint64_t), 1, fp);
  if (meta_size > 0)
    fwrite(hdr->meta, meta_size, 1, fp);

  return GKYL_ARRAY_RIO_SUCCESS;
}

int
gkyl_grid_sub_array_header_write_fp(const struct gkyl_rect_grid *grid,
  const struct gkyl_array_header_info *hdr, FILE *fp)
{
  gkyl_header_meta_write_fp(hdr, fp);  
  
  // Version 0 format is used for rest of the header
  uint64_t real_type = gkyl_array_data_type[hdr->etype];
  fwrite(&real_type, sizeof(uint64_t), 1, fp);
  gkyl_rect_grid_write(grid, fp);

  fwrite(&hdr->esznc, sizeof(uint64_t), 1, fp);
  fwrite(&hdr->tot_cells, sizeof(uint64_t), 1, fp);

  return GKYL_ARRAY_RIO_SUCCESS;
}

int
gkyl_header_meta_read_fp(struct gkyl_array_header_info *hdr, FILE *fp)
{
  size_t frr;
  hdr->meta_size = 0;

  char g0[6];
  frr = fread(g0, sizeof(char[5]), 1, fp); // no trailing '\0'
  g0[5] = '\0';                            // add the NULL
  if (strcmp(g0, "gkyl0") != 0)
    return GKYL_ARRAY_RIO_BAD_VERSION;
  
  uint64_t version;
  frr = fread(&version, sizeof(uint64_t), 1, fp);
  if (version != 1)
    return GKYL_ARRAY_RIO_BAD_VERSION;

  uint64_t file_type;
  frr = fread(&file_type, sizeof(uint64_t), 1, fp);

  uint64_t meta_size;
  frr = fread(&meta_size, sizeof(uint64_t), 1, fp);
  if (1 != frr)
    return GKYL_ARRAY_RIO_FREAD_FAILED;

  hdr->meta = 0;
  if (meta_size > 0) {
    hdr->meta = gkyl_malloc(meta_size);
    if (1 != fread(hdr->meta, meta_size, 1, fp)) {
      gkyl_free(hdr->meta);
      return GKYL_ARRAY_RIO_FREAD_FAILED;
    }
  }

  hdr->file_type = file_type;
  hdr->esznc = 0;
  hdr->tot_cells = 0;
  hdr->meta_size = meta_size;

  return GKYL_ARRAY_RIO_SUCCESS;  
}

static int
grid_sub_array_header_read_fp(struct gkyl_rect_grid *grid,
  struct gkyl_array_header_info *hdr, bool read_meta, FILE *fp)
{
  size_t frr;
  hdr->meta_size = 0;

  char g0[6];
  frr = fread(g0, sizeof(char[5]), 1, fp); // no trailing '\0'
  g0[5] = '\0';                            // add the NULL
  if (strcmp(g0, "gkyl0") != 0)
    return GKYL_ARRAY_RIO_BAD_VERSION;
  
  uint64_t version;
  frr = fread(&version, sizeof(uint64_t), 1, fp);
  if (version != 1)
    return GKYL_ARRAY_RIO_BAD_VERSION;

  uint64_t file_type;
  frr = fread(&file_type, sizeof(uint64_t), 1, fp);
  if (1 != frr)
    return GKYL_ARRAY_RIO_FREAD_FAILED;

  uint64_t meta_size;
  frr = fread(&meta_size, sizeof(uint64_t), 1, fp);
  if (1 != frr)
    return GKYL_ARRAY_RIO_FREAD_FAILED;

  if (meta_size > 0) {
    if (read_meta) {
      hdr->meta = gkyl_malloc(meta_size);
      if (1 != fread(hdr->meta, meta_size, 1, fp)) {
        gkyl_free(hdr->meta);
        return GKYL_ARRAY_RIO_FREAD_FAILED;
      }
    }
    else {
      fseek(fp, meta_size, SEEK_CUR);
    }
  }

  uint64_t real_type = 0;
  if (1 != fread(&real_type, sizeof(uint64_t), 1, fp))
    return GKYL_ARRAY_RIO_FREAD_FAILED;

  gkyl_rect_grid_read(grid, fp);

  uint64_t esznc = 0;
  if (1 != fread(&esznc, sizeof(uint64_t), 1, fp))
    return GKYL_ARRAY_RIO_FREAD_FAILED;;

  uint64_t tot_cells = 0;
  if (1 != fread(&tot_cells, sizeof(uint64_t), 1, fp))
    return GKYL_ARRAY_RIO_FREAD_FAILED;;

  uint64_t nrange = 1;
  if (file_type == gkyl_file_type_int[GKYL_MULTI_RANGE_DATA_FILE])
    if (1 != fread(&nrange, sizeof(uint64_t), 1, fp))
      return GKYL_ARRAY_RIO_FREAD_FAILED;

  hdr->file_type = file_type;
  hdr->etype = gkyl_array_code_to_data_type[real_type];
  hdr->esznc = esznc;
  hdr->tot_cells = tot_cells;
  hdr->meta_size = meta_size;
  hdr->nrange = nrange;

  return GKYL_ARRAY_RIO_SUCCESS;
}

int
gkyl_grid_sub_array_header_read_fp(struct gkyl_rect_grid *grid,
  struct gkyl_array_header_info *hdr, FILE *fp)
{
  return grid_sub_array_header_read_fp(grid, hdr, true, fp);
}

void
gkyl_grid_sub_array_header_release(struct gkyl_array_header_info *hdr)
{
  if (hdr->meta_size>0) {
    gkyl_free(hdr->meta);
    hdr->meta_size = 0;
  }
}

enum gkyl_array_rio_status
gkyl_grid_sub_array_header_read(struct gkyl_rect_grid *grid,
  struct gkyl_array_header_info *hdr, const char *fname)
{
  enum gkyl_array_rio_status status = GKYL_ARRAY_RIO_FOPEN_FAILED;
  FILE *fp = 0;
  with_file(fp, fname, "r") {
    status = gkyl_grid_sub_array_header_read_fp(grid, hdr, fp);
  }
  return status;
}

void
gkyl_array_header_info_release(struct gkyl_array_header_info *info)
{
  if (info->meta_size > 0)
    gkyl_free(info->meta);
}

enum gkyl_array_rio_status
gkyl_grid_sub_array_write(const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  const struct gkyl_msgpack_data *meta,
  const struct gkyl_array *arr, const char *fname)
{
  enum gkyl_array_rio_status status = GKYL_ARRAY_RIO_FOPEN_FAILED;
  FILE *fp = 0;
  int err;
  with_file (fp, fname, "w") {
    
    status = gkyl_grid_sub_array_header_write_fp(grid,
      &(struct gkyl_array_header_info) {
        .file_type = gkyl_file_type_int[GKYL_FIELD_DATA_FILE],
        .etype = arr->type,
        .esznc = arr->esznc,
        .tot_cells = range->volume,
        .meta_size = meta ? meta->meta_sz : 0,
        .meta = meta ? meta->meta : 0 
      },
      fp
    );

    if (status == GKYL_ARRAY_RIO_SUCCESS)
      sub_array_write_priv(range, arr, fp);
  }
  return status;
}

static enum gkyl_array_rio_status
grid_sub_array_read_ft_1(const struct gkyl_rect_grid *grid,
  struct gkyl_array_header_info *hdr, const struct gkyl_range *range,
  struct gkyl_array *arr, FILE *fp)
{
  size_t loc = gkyl_base_hdr_size(hdr->meta_size)
    + gkyl_file_type_1_hrd_size(grid->ndim);
  fseek(fp, loc, SEEK_SET);

  struct gkyl_range blk_rng;
  gkyl_range_init_from_shape1(&blk_rng, grid->ndim, grid->cells);

  struct gkyl_range inter;
  int not_empty = gkyl_range_intersect(&inter, &blk_rng, range);

  if (not_empty) {
    uint64_t sz = hdr->tot_cells;
    gkyl_mem_buff buff = gkyl_mem_buff_new(sz*hdr->esznc);

    if (1 != fread(gkyl_mem_buff_data(buff), sz*hdr->esznc, 1, fp)) {
      gkyl_mem_buff_release(buff);
      return GKYL_ARRAY_RIO_FREAD_FAILED;
    }
    
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &inter);
    while (gkyl_range_iter_next(&iter)) {
      
      char *out = gkyl_array_fetch(arr, gkyl_range_idx(range, iter.idx));
      const char *inp = gkyl_mem_buff_data(buff) + hdr->esznc*gkyl_range_idx(&blk_rng, iter.idx);
      memcpy(out, inp, hdr->esznc);
    }    

    gkyl_mem_buff_release(buff);
  }
    
  return GKYL_ARRAY_RIO_SUCCESS;
}

static enum gkyl_array_rio_status
grid_sub_array_read_ft_3(const struct gkyl_rect_grid *grid,
  struct gkyl_array_header_info *hdr, const struct gkyl_range *range,
  struct gkyl_array *arr, FILE *fp)
{
  size_t rng_sz = gkyl_file_type_3_range_hrd_size(grid->ndim);
  size_t loc = gkyl_base_hdr_size(hdr->meta_size) + gkyl_file_type_3_hrd_size(grid->ndim);

  gkyl_mem_buff buff = gkyl_mem_buff_new(10); // will be reallocated

  for (int r=0; r<hdr->nrange; ++r) {
    uint64_t sz, loidx[GKYL_MAX_DIM], upidx[GKYL_MAX_DIM];
    fseek(fp, loc, SEEK_SET);

    // read lower, upper indices and number of elements stored
    if (1 != fread(loidx, sizeof(uint64_t[grid->ndim]), 1, fp))
      return GKYL_ARRAY_RIO_FREAD_FAILED;
    if (1 != fread(upidx, sizeof(uint64_t[grid->ndim]), 1, fp))
      return GKYL_ARRAY_RIO_FREAD_FAILED;
    if (1 != fread(&sz, sizeof(uint64_t), 1, fp))
      return GKYL_ARRAY_RIO_FREAD_FAILED;

    // construct range of indices corresponding to data in block
    int loidx_i[GKYL_MAX_DIM]= { 0 } , upidx_i[GKYL_MAX_DIM] = { 0 };
    for (int d=0; d<grid->ndim; ++d) {
      loidx_i[d] = loidx[d];
      upidx_i[d] = upidx[d];
    }
    struct gkyl_range blk_rng; // block range
    gkyl_range_init(&blk_rng, grid->ndim, loidx_i, upidx_i);

    struct gkyl_range inter; // intersection
    int not_empty = gkyl_range_intersect(&inter, &blk_rng, range);
    if (not_empty) {

      buff = gkyl_mem_buff_resize(buff, sz*hdr->esznc);
      if (1 != fread(gkyl_mem_buff_data(buff), sz*hdr->esznc, 1, fp)) {
        gkyl_mem_buff_release(buff);
        return GKYL_ARRAY_RIO_FREAD_FAILED;
      }

      struct gkyl_range_iter iter;
      gkyl_range_iter_init(&iter, &inter);
      while (gkyl_range_iter_next(&iter)) {
        char *out = gkyl_array_fetch(arr, gkyl_range_idx(range, iter.idx));
        const char *inp = gkyl_mem_buff_data(buff) + hdr->esznc*gkyl_range_idx(&blk_rng, iter.idx);
        memcpy(out, inp, hdr->esznc);
      }
    }

    loc += rng_sz + sz*hdr->esznc;
  }

  gkyl_mem_buff_release(buff);
  
  return GKYL_ARRAY_RIO_SUCCESS;
}

enum gkyl_array_rio_status
gkyl_grid_sub_array_read(struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  struct gkyl_array *arr, const char *fname)
{
  enum gkyl_array_rio_status status = GKYL_ARRAY_RIO_FOPEN_FAILED;
  struct gkyl_array_header_info hdr;
  FILE *fp = 0;
  with_file (fp, fname, "r") {
    grid_sub_array_header_read_fp(grid, &hdr, false, fp);
    
    if (hdr.file_type == 1)
      status = grid_sub_array_read_ft_1(grid, &hdr, range, arr, fp);
    if (hdr.file_type == 3)
      status = grid_sub_array_read_ft_3(grid, &hdr, range, arr, fp);
  }
  return status;
}

struct gkyl_array*
gkyl_grid_array_new_from_file(struct gkyl_rect_grid *grid, const char* fname)
{
  struct gkyl_array *arr = 0;
  struct gkyl_array_header_info hdr;

  enum gkyl_array_rio_status status = GKYL_ARRAY_RIO_FREAD_FAILED;
  FILE *fp = 0;
  with_file (fp, fname, "r") {
    status = grid_sub_array_header_read_fp(grid, &hdr, false, fp);
  }

  if (status != GKYL_ARRAY_RIO_SUCCESS)
    return 0;

  size_t nc = hdr.esznc/gkyl_elem_type_size[hdr.etype];
  arr = gkyl_array_new(hdr.etype, nc, hdr.tot_cells);
  struct gkyl_range range;
  gkyl_range_init_from_shape1(&range, grid->ndim, grid->cells);

  status = gkyl_grid_sub_array_read(grid, &range, arr, fname);

  if (status != GKYL_ARRAY_RIO_SUCCESS) {
    gkyl_array_release(arr);
    arr = 0;
  }
    
  return arr;
}
// ended inlining array_rio.c 
// start inlining array_rio_format_desc.c 
// skipping file: gkyl_array_rio_format_desc.h 

#include <stdint.h>
#include <stdio.h>
#include <string.h>

// skipping file: gkyl_util.h 

size_t
gkyl_base_hdr_size(size_t meta_sz)
{
  size_t sz = 0;
  // magic string
  sz += 5; // "gkyl0"
  // version
  sz += sizeof(uint64_t);
  // file_type
  sz += sizeof(uint64_t);
  // metadata
  sz += sizeof(uint64_t) + meta_sz;

  return sz;
}

size_t
gkyl_file_type_1_partial_hrd_size(int ndim)
{
  size_t sz = 0;
  // real_type
  sz += sizeof(uint64_t);
  // ndim
  sz += sizeof(uint64_t);
  // cells
  sz += sizeof(uint64_t[ndim]);
  // lower, upper
  sz += sizeof(double[2*ndim]);
  return sz;
}

size_t
gkyl_file_type_1_hrd_size(int ndim)
{
  size_t sz = gkyl_file_type_1_partial_hrd_size(ndim);
  // esznc
  sz += sizeof(uint64_t);
  // size (total number of cells)
  sz += sizeof(uint64_t);
  return sz;
}

size_t
gkyl_file_type_2_hrd_size(void)
{
  size_t sz = 0;
  // real_type
  sz += sizeof(uint64_t);
  // esznc
  sz += sizeof(uint64_t);
  // size (total number of cells)
  sz += sizeof(uint64_t);
  return sz;
}

size_t
gkyl_file_type_3_hrd_size(int ndim)
{
  size_t sz = gkyl_file_type_1_hrd_size(ndim);
  sz += sizeof(uint64_t); // nrange
  return sz;
}

size_t
gkyl_file_type_3_range_hrd_size(int ndim)
{
  size_t sz = 0;
  // loidx and upidx
  sz += sizeof(uint64_t[2*ndim]);
  sz += sizeof(uint64_t);
  return sz;
}

int
gkyl_get_gkyl_file_type(const char *fname)
{
  int file_type = -1;
  FILE *fp = 0;

  with_file(fp, fname, "r") {
    size_t frr;
    char g0[6];
    frr = fread(g0, sizeof(char[5]), 1, fp); // no trailing '\0'
    g0[5] = '\0';                            // add the NULL
    if (strcmp(g0, "gkyl0") != 0) {
      file_type = -1;
      goto finish_with_file;
    }
  
    uint64_t version;
    frr = fread(&version, sizeof(uint64_t), 1, fp);
    if (version != 1) {
      file_type = -1;
      goto finish_with_file;
    }
    
    uint64_t file_type_u64;
    frr = fread(&file_type_u64, sizeof(uint64_t), 1, fp);
    if (1 != frr) {
      file_type = -1;
      goto finish_with_file;      
    }

    file_type = file_type_u64;
    
    finish_with_file:
    ;
  }
  return file_type;
}
// ended inlining array_rio_format_desc.c 
// start inlining comm.c 
// start inlining gkyl_comm_priv.h 

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

// Get local "rank"
typedef int (*get_rank_t)(struct gkyl_comm *comm, int *rank);

// Get number of ranks
typedef int (*get_size_t)(struct gkyl_comm *comm, int *sz);

// "Reduce" all elements of @a type in array @a data and store output in @a out
typedef int (*allreduce_t)(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp, void *out);

// Gather local arrays into global array on each process.
typedef int (*gkyl_array_allgather_t)(struct gkyl_comm *comm,
  const struct gkyl_range *local, const struct gkyl_range *global,
  const struct gkyl_array *array_local, struct gkyl_array *array_global);

// Broadcast array to other processes.
typedef int (*gkyl_array_bcast_t)(struct gkyl_comm *comm,
  const struct gkyl_array *array_send, struct gkyl_array *array_recv, int root);

// "Synchronize" @a array across the regions or blocks.
typedef int (*gkyl_array_sync_t)(struct gkyl_comm *comm,
  const struct gkyl_range *local, const struct gkyl_range *local_ext,
  struct gkyl_array *array);

// "Synchronize" @a array across the periodic directions
typedef int (*gkyl_array_per_sync_t)(struct gkyl_comm *comm,
  const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  int nper_dirs, const int *per_dirs,
  struct gkyl_array *array);

// Write array to specified file
typedef int (*gkyl_array_write_t)(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_msgpack_data *meta,                                
  const struct gkyl_array *arr, const char *fname);

// Read array from specified file
typedef int (*gkyl_array_read_t)(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  struct gkyl_array *arr, const char *fname);

// Create a new communicator that extends the communcator to work on a
// extended domain specified by erange
typedef struct gkyl_comm* (*extend_comm_t)(const struct gkyl_comm *comm,
  const struct gkyl_range *erange);

// Create a new communicator by splitting a comm, and choosing members
// of new communicator according to the color rank. It can be used with
// a new decomp object, or the same one used for the parent comm, depending
// of the use case.
typedef struct gkyl_comm* (*split_comm_t)(const struct gkyl_comm *comm,
  int color, struct gkyl_rect_decomp *new_decomp);

// Create a new communicator from the input comm that takes a list of
// ranks to include in it.
typedef struct gkyl_comm *(*create_comm_from_ranks_t)(
  const struct gkyl_comm *comm, int nranks, const int *ranks,
  struct gkyl_rect_decomp *new_decomp,
  bool *is_valid
);

// Barrier
typedef int (*barrier_t)(struct gkyl_comm *comm);

// Start and end a group call (e.g. in NCCL).
typedef void (*comm_group_call_start_t)();
typedef void (*comm_group_call_end_t)();

// Private structure: not available to public-facing API but to
// specific comm types
struct gkyl_comm_priv {
  struct gkyl_comm pub_comm; // public facing communicator

  // FOLLOWING DO NOT NEED A DECOMP
  
  get_rank_t get_rank; // get local rank function.
  get_size_t get_size; // get number of ranks.
  barrier_t barrier; // barrier.
  allreduce_t allreduce; // all reduce function
  allreduce_t allreduce_host; // all reduce using the host (MPI) communicator
  
  extend_comm_t extend_comm; // extend communcator
  split_comm_t split_comm;   // split communicator.
  create_comm_from_ranks_t  create_comm_from_ranks; // communictor from ranks

  comm_group_call_start_t comm_group_call_start; // start a group call
  comm_group_call_end_t comm_group_call_end; // end a group call

  // FOLLOWING NEED A DECOMP 

  gkyl_array_allgather_t gkyl_array_allgather; // gather local arrays to global array
  gkyl_array_allgather_t gkyl_array_allgather_host; // gather local arrays to global array on host

  gkyl_array_bcast_t gkyl_array_bcast; // broadcast array to other processes
  gkyl_array_bcast_t gkyl_array_bcast_host; // broadcast host side array to other processes

  gkyl_array_sync_t gkyl_array_sync; // sync array
  gkyl_array_per_sync_t gkyl_array_per_sync; // sync array in periodic dirs

  gkyl_array_write_t gkyl_array_write; // array output
  gkyl_array_read_t gkyl_array_read; // array input
};
// ended inlining gkyl_comm_priv.h 

struct gkyl_comm*
gkyl_comm_acquire(const struct gkyl_comm *comm)
{
  gkyl_ref_count_inc(&comm->ref_count);
  return (struct gkyl_comm*) comm;
}

void
gkyl_comm_release(const struct gkyl_comm *comm)
{
  if (comm)
    gkyl_ref_count_dec(&comm->ref_count);
}

int
gkyl_comm_get_rank(struct gkyl_comm *pcomm, int *rank)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);
  return comm->get_rank(pcomm, rank);
}

int
gkyl_comm_get_size(struct gkyl_comm *pcomm, int *sz)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);
  return comm->get_size(pcomm, sz);
}

int
gkyl_comm_allreduce(struct gkyl_comm *pcomm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp, void *out)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  return comm->allreduce(pcomm, type, op, nelem, inp, out);
}

int
gkyl_comm_allreduce_host(struct gkyl_comm *pcomm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp, void *out)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  return comm->allreduce_host(pcomm, type, op, nelem, inp, out);
}

int
gkyl_comm_array_allgather(struct gkyl_comm *pcomm, 
  const struct gkyl_range *local, const struct gkyl_range *global,
  const struct gkyl_array *array_local, struct gkyl_array *array_global)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  return comm->gkyl_array_allgather(pcomm, local, global, array_local, array_global);
}

int
gkyl_comm_array_allgather_host(struct gkyl_comm *pcomm, 
  const struct gkyl_range *local, const struct gkyl_range *global,
  const struct gkyl_array *array_local, struct gkyl_array *array_global)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  return comm->gkyl_array_allgather_host(pcomm, local, global, array_local, array_global);
}

int
gkyl_comm_array_bcast(struct gkyl_comm *pcomm, 
  const struct gkyl_array *array_send, struct gkyl_array *array_recv, int root)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  return comm->gkyl_array_bcast(pcomm, array_send, array_recv, root);
}

int
gkyl_comm_array_bcast_host(struct gkyl_comm *pcomm, 
  const struct gkyl_array *array_send, struct gkyl_array *array_recv, int root)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);
  return comm->gkyl_array_bcast_host(pcomm, array_send, array_recv, root);
}

int
gkyl_comm_array_sync(struct gkyl_comm *pcomm,
  const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  struct gkyl_array *array)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  comm->barrier(pcomm);
  return comm->gkyl_array_sync(pcomm, local, local_ext, array);
}

int
gkyl_comm_array_per_sync(struct gkyl_comm *pcomm,
  const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  int nper_dirs, const int *per_dirs,
  struct gkyl_array *array)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  comm->barrier(pcomm);
  return comm->gkyl_array_per_sync(pcomm, local, local_ext,
    nper_dirs, per_dirs, array);
}

int
gkyl_comm_barrier(struct gkyl_comm *pcomm)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  return comm->barrier(pcomm);
}

void
gkyl_comm_group_call_start(struct gkyl_comm *pcomm)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  comm->comm_group_call_start();
}

void
gkyl_comm_group_call_end(struct gkyl_comm *pcomm)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  comm->comm_group_call_end();
}

int
gkyl_comm_array_write(struct gkyl_comm *pcomm,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_msgpack_data *meta,
  const struct gkyl_array *arr, const char *fname)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  int status = comm->gkyl_array_write(pcomm, grid, range, meta, arr, fname);
  gkyl_comm_barrier(pcomm);
  return status;
}

int
gkyl_comm_array_read(struct gkyl_comm *pcomm,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  struct gkyl_array *arr, const char *fname)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  int status = comm->gkyl_array_read(pcomm, grid, range, arr, fname);
  gkyl_comm_barrier(pcomm);
  return status;
}

struct gkyl_comm*
gkyl_comm_extend_comm(const struct gkyl_comm *pcomm,
  const struct gkyl_range *erange)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);
  return comm->extend_comm(pcomm, erange);
}

struct gkyl_comm*
gkyl_comm_split_comm(const struct gkyl_comm *pcomm, int color,
  struct gkyl_rect_decomp *new_decomp)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  return comm->split_comm(pcomm, color, new_decomp);
}

struct gkyl_comm *
gkyl_comm_create_comm_from_ranks(const struct gkyl_comm *pcomm, int nranks,
  const int *ranks, struct gkyl_rect_decomp *new_decomp,
  bool *is_valid)
{
  struct gkyl_comm_priv *comm = container_of(pcomm, struct gkyl_comm_priv, pub_comm);  
  return comm->create_comm_from_ranks(pcomm, nranks, ranks, new_decomp, is_valid);
}
// ended inlining comm.c 
// start inlining dynvec.c 
// skipping file: gkyl_alloc.h 
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
// skipping file: gkyl_elem_type_priv.h 
// skipping file: gkyl_ref_count.h 
// skipping file: gkyl_util.h 

#include <errno.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

/* Size by which vector grows each time it is reallocated */
static const size_t DYNVEC_ALLOC_SZ = 1024;

struct gkyl_dynvec_tag {
  enum gkyl_elem_type type; // type of data stored in vector
  size_t elemsz, ncomp; // size of elements, number of 'components'

  size_t cloc; // current location to insert data into
  size_t csize; // current number of elements allocated

  size_t esznc; // elemsz*ncomp
  void *data; // pointer to data
  double *tm_mesh; // time stamps
  
  struct gkyl_ref_count ref_count;  
};

static void
dynvec_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dynvec_tag *dv = container_of(ref, struct gkyl_dynvec_tag, ref_count);
  gkyl_free(dv->data);
  gkyl_free(dv->tm_mesh);
  gkyl_free(dv);
}

gkyl_dynvec
gkyl_dynvec_new(enum gkyl_elem_type type, size_t ncomp)
{
  struct gkyl_dynvec_tag *dv = gkyl_malloc(sizeof(struct gkyl_dynvec_tag));
  dv->type = type;
  dv->elemsz = gkyl_elem_type_size[type];
  dv->ncomp = ncomp;
  dv->esznc = dv->elemsz*dv->ncomp;
  dv->csize = DYNVEC_ALLOC_SZ;
  dv->cloc = 0;
  
  dv->data = gkyl_calloc(dv->csize, dv->esznc);
  dv->tm_mesh = gkyl_calloc(dv->csize, sizeof(double));
  
  dv->ref_count = gkyl_ref_count_init(dynvec_free);
  
  return dv;
}

int gkyl_dynvec_elem_type(gkyl_dynvec vec) { return vec->type; }
int gkyl_dynvec_ncomp(gkyl_dynvec vec) { return vec->ncomp; }

void
gkyl_dynvec_reserve_more(gkyl_dynvec dv, size_t rsize)
{
  int n = gkyl_int_div_up(rsize, DYNVEC_ALLOC_SZ);
  dv->csize = n*DYNVEC_ALLOC_SZ + dv->csize;
  dv->data = gkyl_realloc(dv->data, dv->csize*dv->esznc);
  dv->tm_mesh = gkyl_realloc(dv->tm_mesh, dv->csize*sizeof(double));  
}

void
gkyl_dynvec_append(gkyl_dynvec dv, double tm, const void *data)
{
  size_t loc = dv->cloc;
  if (loc >= dv->csize) {
    dv->csize += DYNVEC_ALLOC_SZ;
    dv->data = gkyl_realloc(dv->data, dv->csize*dv->esznc);
    dv->tm_mesh = gkyl_realloc(dv->tm_mesh, dv->csize*sizeof(double));
  }

  // set data
  dv->tm_mesh[loc] = tm;
  memcpy((char*)dv->data + dv->esznc*loc, data, dv->esznc);
  
  dv->cloc += 1;
}

bool
gkyl_dynvec_get(const gkyl_dynvec dv, size_t idx, void *data)
{
  if (idx >= dv->cloc) return false;
  memcpy(data, (char*)dv->data + dv->esznc*idx, dv->esznc);
  return true;  
}

double
gkyl_dynvec_get_tm(const gkyl_dynvec dv, size_t idx)
{
  if (idx >= dv->cloc) return 0.0;
  return dv->tm_mesh[idx];
}

bool
gkyl_dynvec_getlast(const gkyl_dynvec dv, void *data)
{
  size_t loc = dv->cloc;
  if (loc == 0) return false;
  return gkyl_dynvec_get(dv, loc-1, data);
}

double
gkyl_dynvec_getlast_tm(const gkyl_dynvec dv)
{
  size_t loc = dv->cloc;
  return loc == 0 ? 0.0 : dv->tm_mesh[loc-1];
}

size_t
gkyl_dynvec_size(const gkyl_dynvec vec)
{
  return vec->cloc;
}

size_t
gkyl_dynvec_capacity(const gkyl_dynvec vec)
{
  return vec->csize;
}

void
gkyl_dynvec_clear(gkyl_dynvec dv)
{
  dv->cloc = 0;
  dv->csize = DYNVEC_ALLOC_SZ;
  dv->data = gkyl_realloc(dv->data, dv->csize*dv->esznc);
  dv->tm_mesh = gkyl_realloc(dv->tm_mesh, dv->csize*sizeof(double));
}

void
gkyl_dynvec_clear_all_but(gkyl_dynvec dv, size_t num)
{
  if (num>dv->cloc) return;
  
  size_t cloc = dv->cloc;
  dv->csize = DYNVEC_ALLOC_SZ;

  void *data = gkyl_malloc(num*dv->esznc);
  double *tm_mesh = gkyl_malloc(sizeof(double[num]));

  size_t low = num>cloc ? 0 : cloc-num; // lower index to copy from
  size_t ncpy = num>cloc ? cloc : num; // number of elemetns to copy
  dv->cloc = ncpy;

  memcpy(tm_mesh, dv->tm_mesh+low, ncpy*sizeof(double));
  memcpy(data, (char*)dv->data+low*dv->esznc, ncpy*dv->esznc);
  
  dv->data = gkyl_realloc(dv->data, dv->csize*dv->esznc);
  dv->tm_mesh = gkyl_realloc(dv->tm_mesh, dv->csize*sizeof(double));

  memcpy(dv->data, data, ncpy*dv->esznc);
  memcpy(dv->tm_mesh, tm_mesh, ncpy*sizeof(double));

  gkyl_free(data);
  gkyl_free(tm_mesh);
}

gkyl_dynvec
gkyl_dynvec_acquire(const gkyl_dynvec vec)
{
  gkyl_ref_count_inc(&vec->ref_count);
  return (struct gkyl_dynvec_tag*) vec;
}

static int
gkyl_dynvec_write_mode(const gkyl_dynvec vec,
  const char *fname, const char *mode)
{
  const char g0[5] = "gkyl0";

  FILE *fp = 0;
  with_file (fp, fname, mode) {
    fseek(fp, 0, SEEK_END);
    
    // Version 1 header
    fwrite(g0, sizeof(char[5]), 1, fp);
    uint64_t version = 1;
    fwrite(&version, sizeof(uint64_t), 1, fp);
    fwrite(&gkyl_file_type_int[GKYL_DYNVEC_DATA_FILE], sizeof(uint64_t), 1, fp);
    uint64_t meta_size = 0; // THIS WILL CHANGE ONCE METADATA IS EMBEDDED
    fwrite(&meta_size, sizeof(uint64_t), 1, fp);
    
    uint64_t real_type = gkyl_array_data_type[vec->type];
    fwrite(&real_type, sizeof(uint64_t), 1, fp);

    uint64_t esznc = vec->esznc, size = gkyl_dynvec_size(vec);
    fwrite(&esznc, sizeof(uint64_t), 1, fp);
    fwrite(&size, sizeof(uint64_t), 1, fp); 

    fwrite(vec->tm_mesh, sizeof(double)*size, 1, fp);
    fwrite(vec->data, esznc*size, 1, fp);
  }

  return errno;
}

int
gkyl_dynvec_write(const gkyl_dynvec vec, const char *fname)
{
  return gkyl_dynvec_write_mode(vec, fname, "w");
}

int
gkyl_dynvec_awrite(const gkyl_dynvec vec, const char *fname)
{
  return gkyl_dynvec_write_mode(vec, fname, "a");
}

// ncomp returned in 'ncomp'
static bool
gkyl_dynvec_read_ncomp_1(FILE *fp, struct gkyl_dynvec_etype_ncomp *enc)
{
  size_t frr;
  // Version 1 header
  char g0[6];
  if (1 != fread(g0, sizeof(char[5]), 1, fp))
    return false;
  g0[5] = '\0'; // add the NULL
  if (strcmp(g0, "gkyl0") != 0)
    return false;

  uint64_t version;
  frr = fread(&version, sizeof(uint64_t), 1, fp);
  if (version != 1)
    return false;

  uint64_t file_type;
  frr = fread(&file_type, sizeof(uint64_t), 1, fp);
  if (file_type != gkyl_file_type_int[GKYL_DYNVEC_DATA_FILE])
    return false;

  uint64_t meta_size;
  frr = fread(&meta_size, sizeof(uint64_t), 1, fp);

  // read ahead by specified bytes: meta-data is not read in this
  // method
  fseek(fp, meta_size, SEEK_CUR);

  uint64_t real_code = 0;
  if (1 != fread(&real_code, sizeof(uint64_t), 1, fp))
    return false;
  enc->type = gkyl_array_code_to_data_type[real_code];

  uint64_t esznc;
  if (1 != fread(&esznc, sizeof(uint64_t), 1, fp))
    return false;

  int real_type = gkyl_array_code_to_data_type[real_code];
  enc->ncomp = esznc/gkyl_elem_type_size[real_type];
  
  return true;
}


struct gkyl_dynvec_etype_ncomp
gkyl_dynvec_read_ncomp(const char *fname)
{
  struct gkyl_dynvec_etype_ncomp enc = {
    .type = GKYL_DOUBLE,
    .ncomp = 0
  };
  FILE *fp = 0;
  with_file(fp, fname, "r")
    gkyl_dynvec_read_ncomp_1(fp, &enc);
  return enc;
}

static bool
gkyl_dynvec_read_1(gkyl_dynvec vec, FILE *fp) {
  size_t frr;
  // Version 1 header
  char g0[6];
  if (1 != fread(g0, sizeof(char[5]), 1, fp))
    return false;
  g0[5] = '\0'; // add the NULL
  if (strcmp(g0, "gkyl0") != 0)
    return false;

  uint64_t version;
  frr = fread(&version, sizeof(uint64_t), 1, fp);
  if (version != 1)
    return false;

  uint64_t file_type;
  frr = fread(&file_type, sizeof(uint64_t), 1, fp);
  if (file_type != gkyl_file_type_int[GKYL_DYNVEC_DATA_FILE])
    return false;

  uint64_t meta_size;
  frr = fread(&meta_size, sizeof(uint64_t), 1, fp);

  // read ahead by specified bytes: meta-data is not read in this
  // method
  fseek(fp, meta_size, SEEK_CUR);

  uint64_t real_type = 0;
  if (1 != fread(&real_type, sizeof(uint64_t), 1, fp))
    return false;
  if (real_type != gkyl_array_data_type[vec->type])
    return false;

  uint64_t esznc, size;
  if (1 != fread(&esznc, sizeof(uint64_t), 1, fp))
    return false;
  if (vec->esznc != esznc)
    return false;

  if (1 != fread(&size, sizeof(uint64_t), 1, fp))
    return false;

  // resize vector to allow storing new data
  gkyl_dynvec_reserve_more(vec, size);

  // read time-mesh data
  frr = fread(&vec->tm_mesh[vec->cloc], sizeof(double[size]), 1, fp);
  // read dynvec data
  frr = fread((char *)vec->data + vec->esznc * vec->cloc, size * esznc, 1, fp);

  // bump location so further inserts occurs after newly read data
  vec->cloc = vec->cloc + size;

  return true;
}

bool
gkyl_dynvec_read(gkyl_dynvec vec, const char *fname)
{
  bool status = false;
  FILE *fp = fopen(fname, "r");

  // keep reading till we have no more datasets
  while (1) {
    status = gkyl_dynvec_read_1(vec, fp);
    fpos_t curr_pos;
    fgetpos(fp, &curr_pos);
    
    char g0[6];
    if (1 != fread(g0, sizeof(char[5]), 1, fp))
      break;
    fsetpos(fp, &curr_pos);
  }
  fclose(fp);
  
  return status;
}

void
gkyl_dynvec_to_array(const gkyl_dynvec vec, struct gkyl_array *tm_mesh,
  struct gkyl_array *dyndata)
{
  int nv = gkyl_dynvec_size(vec);
  for (int i=0; i<nv; ++i) {
    double *tmm = gkyl_array_fetch(tm_mesh, i);
    tmm[0] = gkyl_dynvec_get_tm(vec, i);

    void *dd = gkyl_array_fetch(dyndata, i);
    gkyl_dynvec_get(vec, i, dd);
  }
}

void
gkyl_dynvec_release(gkyl_dynvec vec)
{
  if (vec)
    gkyl_ref_count_dec(&vec->ref_count);
}
// ended inlining dynvec.c 
// start inlining eval_offset_fd.c 
// skipping file: gkyl_alloc.h 
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

#include <string.h>

struct gkyl_eval_offset_fd {
  
  struct gkyl_rect_grid grid;
  int num_ret_vals; // number of values returned by eval function
  struct gkyl_offset_descr *offsets; // size num_ret_vals
  evalf_t eval; // function to project
  void *ctx; // evaluation context
};

gkyl_eval_offset_fd*
gkyl_eval_offset_fd_new(const struct gkyl_eval_offset_fd_inp *inp)
{
  struct gkyl_eval_offset_fd *up = gkyl_malloc(sizeof(*up));

  up->grid = *inp->grid;
  int num_ret_vals = up->num_ret_vals = inp->num_ret_vals;
  up->eval = inp->eval;
  up->ctx = inp->ctx;

  up->offsets = gkyl_malloc(sizeof(struct gkyl_offset_descr[num_ret_vals]));
  memcpy(up->offsets, inp->offsets, sizeof(struct gkyl_offset_descr[num_ret_vals]));

  return up;
}

static inline void
comp_to_phys(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = dx[d]*eta[d]+xc[d];
}

void
gkyl_eval_offset_fd_advance(const gkyl_eval_offset_fd *up,
  double tm, const struct gkyl_range *update_rng, struct gkyl_array *out)
{
  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  int num_ret_vals = up->num_ret_vals;
  double fvals[num_ret_vals];
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_rng);
  
  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(&up->grid, iter.idx, xc);

    long lidx = gkyl_range_idx(update_rng, iter.idx);
    double *out_p = gkyl_array_fetch(out, lidx);

    double xc[GKYL_MAX_DIM];
    gkyl_rect_grid_cell_center(&up->grid, iter.idx, xc);

    for (int c=0; c<num_ret_vals; ++c) {
      // We need to evaluate the function once for each ret value as
      // each can be on a different location in the cell. This is not
      // too efficient, but likely does not matter.
      double xout[GKYL_MAX_DIM];
      comp_to_phys(up->grid.ndim, up->offsets[c].od_off, up->grid.dx, xc, xout);

      up->eval(tm, xout, fvals, up->ctx);
      out_p[c] = fvals[c];
    }
  }
}

void
gkyl_eval_offset_fd_release(gkyl_eval_offset_fd *up)
{
  gkyl_free(up->offsets);
  gkyl_free(up);
}
// ended inlining eval_offset_fd.c 
// start inlining fv_proj.c 
#include <math.h>
#include <string.h>
#include <assert.h>

// skipping file: gkyl_alloc.h 
// skipping file: gkyl_array_ops.h 
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
// start inlining gkyl_quad_type.h 

// type of quadrature to use
enum gkyl_quad_type {
  GKYL_GAUSS_QUAD = 0, // Gauss-Legendre quadrature
  GKYL_GAUSS_LOBATTO_QUAD, // Gauss-Lobatto quadrature
  GKYL_POSITIVITY_QUAD // Positivity quadrature nodes
};
// ended inlining gkyl_quad_type.h 


GKYL_CU_DH
void static inline
fv_eval_1d_ser_p0(const double *z, double *b )
{
  b[0] = 7.0710678118654757e-01;
}

GKYL_CU_DH
void static inline
fv_eval_2d_ser_p0(const double *z, double *b )
{
  b[0] = 5.0000000000000000e-01;
}

GKYL_CU_DH
void static inline
fv_eval_3d_ser_p0(const double *z, double *b )
{
  b[0] = 3.5355339059327379e-01;
}

GKYL_CU_DH
void static inline
fv_eval_4d_ser_p0(const double *z, double *b )
{
  b[0] = 2.5000000000000000e-01;
}

GKYL_CU_DH
void static inline
fv_eval_5d_ser_p0(const double *z, double *b )
{
  b[0] = 1.7677669529663689e-01;
}

GKYL_CU_DH
void static inline
fv_eval_6d_ser_p0(const double *z, double *b )
{
  b[0] = 1.2500000000000000e-01;
}

// Types for various kernels.
typedef void (*fv_proj_on_basis_c2p_t)(const double *xcomp, double *xphys, void *ctx);

// input packaged as a struct
struct gkyl_fv_proj_inp {
  const struct gkyl_rect_grid *grid; // grid on which to project

  enum gkyl_quad_type qtype; // quadrature to use
  
  int num_quad; // number of quadrature points
  int num_ret_vals; // number of return values in eval function
  evalf_t eval; // function to project
  void *ctx; // function context

  fv_proj_on_basis_c2p_t c2p_func; // Function that transforms a set of ndim
                                // computational coordinates to physical ones.
  void *c2p_func_ctx; // Context for c2p_func.
};

struct gkyl_fv_proj {
  struct gkyl_rect_grid grid;
  int num_quad; // number of quadrature points to use in each direction
  int num_ret_vals; // number of values returned by eval function
  evalf_t eval; // function to project
  void *ctx; // evaluation context

  int num_basis; // number of basis functions
  int tot_quad; // total number of quadrature points
  struct gkyl_array *ordinates; // ordinates for quadrature
  struct gkyl_array *weights; // weights for quadrature
  struct gkyl_array *basis_at_ords; // basis functions at ordinates

  fv_proj_on_basis_c2p_t c2p; // Function transformin comp to phys coords.
  void *c2p_ctx; // Context for the c2p mapping.
};

// Identity comp to phys coord mapping, for when user doesn't provide a map.
static inline void
c2p_identity(const double *xcomp, double *xphys, void *ctx)
{
  struct gkyl_rect_grid *grid = ctx;
  int ndim = grid->ndim;
  for (int d=0; d<ndim; d++) xphys[d] = xcomp[d];
}


static struct gkyl_fv_proj*
gkyl_fv_proj_priv_inew(const struct gkyl_fv_proj_inp *inp)
{
  struct gkyl_fv_proj *up = gkyl_malloc(sizeof(struct gkyl_fv_proj));

  up->grid = *inp->grid;
  int num_quad = up->num_quad = inp->num_quad == 0 ? 1 : inp->num_quad;
  int num_ret_vals = up->num_ret_vals = inp->num_ret_vals;
  up->eval = inp->eval;
  up->ctx = inp->ctx;
  up->num_basis = 1;

  if (inp->c2p_func == 0) {
    up->c2p = c2p_identity;
    up->c2p_ctx = &up->grid; // Use grid as the context since all we need is ndim.
  }
  else {
    up->c2p = inp->c2p_func;
    up->c2p_ctx = inp->c2p_func_ctx;
  }

  double ordinates1[num_quad], weights1[num_quad];

  if (inp->qtype == GKYL_GAUSS_QUAD) {
    if (num_quad <= gkyl_gauss_max) {
      // use pre-computed values if possible (these are more accurate
      // than computing them on the fly)
      memcpy(ordinates1, gkyl_gauss_ordinates[num_quad], sizeof(double[num_quad]));
      memcpy(weights1, gkyl_gauss_weights[num_quad], sizeof(double[num_quad]));
    }
    else {
      gkyl_gauleg(-1, 1, ordinates1, weights1, num_quad);
    }
  }
  else if (inp->qtype == GKYL_GAUSS_LOBATTO_QUAD) {
    assert( (num_quad > 1) && (num_quad <= gkyl_gauss_max) );
    
    // Gauss-Lobatto quadrature
    memcpy(ordinates1, gkyl_gauss_lobatto_ordinates[num_quad], sizeof(double[num_quad]));
    memcpy(weights1, gkyl_gauss_lobatto_weights[num_quad], sizeof(double[num_quad]));
  }
  else {
    fprintf(stderr, "Quadrature type not available. Exiting... \n");
    assert(false);
  }

  // create range to loop over quadrature points
  int qshape[GKYL_MAX_DIM];
  for (int i=0; i<inp->grid->ndim; ++i) qshape[i] = num_quad;
  struct gkyl_range qrange;
  gkyl_range_init_from_shape(&qrange, inp->grid->ndim, qshape);

  int tot_quad = up->tot_quad = qrange.volume;

  // create ordinates and weights for multi-D quadrature  
  up->ordinates = gkyl_array_new(GKYL_DOUBLE, inp->grid->ndim, tot_quad);
  up->weights = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &qrange);

  while (gkyl_range_iter_next(&iter)) {
    long node = gkyl_range_idx(&qrange, iter.idx);
    
    // set ordinates
    double *ord = gkyl_array_fetch(up->ordinates, node);
    for (int i=0; i<inp->grid->ndim; ++i)
      ord[i] = ordinates1[iter.idx[i]-qrange.lower[i]];
    
    // set weights
    double *wgt = gkyl_array_fetch(up->weights, node);
    wgt[0] = 1.0;
    for (int i=0; i<qrange.ndim; ++i)
      wgt[0] *= weights1[iter.idx[i]-qrange.lower[i]];
  }

  void (*eval_func[])(const double *, double *) = {
    0,
    fv_eval_1d_ser_p0,
    fv_eval_2d_ser_p0,
    fv_eval_3d_ser_p0,
    fv_eval_4d_ser_p0,
    fv_eval_5d_ser_p0,
    fv_eval_6d_ser_p0
  };

  // pre-compute basis functions at ordinates
  up->basis_at_ords = gkyl_array_new(GKYL_DOUBLE, 1, tot_quad);
  for (int n=0; n<tot_quad; ++n)
    eval_func[inp->grid->ndim](gkyl_array_fetch(up->ordinates, n),
      gkyl_array_fetch(up->basis_at_ords, n));

  return up;
}

static struct gkyl_fv_proj*
gkyl_fv_proj_priv_new(const struct gkyl_rect_grid *grid,
  int num_quad, int num_ret_vals, evalf_t eval, void *ctx)
{
  return gkyl_fv_proj_priv_inew( &(struct gkyl_fv_proj_inp) {
      .grid = grid,
      .qtype = GKYL_GAUSS_QUAD,
      .num_quad = num_quad,
      .num_ret_vals = num_ret_vals,
      .eval = eval,
      .ctx = ctx,
      .c2p_func = 0,
      .c2p_func_ctx = NULL,
    }
  );
}


static inline void
log_to_comp(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  // Convert logical to computational coordinates.
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
}

static void
gkyl_fv_proj_priv_quad(const struct gkyl_fv_proj *up, const struct gkyl_array *fun_at_ords, double* f)
{
  int num_basis = up->num_basis;
  int tot_quad = up->tot_quad;
  int num_ret_vals = up->num_ret_vals;

  const double* GKYL_RESTRICT weights = up->weights->data;
  const double* GKYL_RESTRICT basis_at_ords = up->basis_at_ords->data;
  const double* GKYL_RESTRICT func_at_ords = fun_at_ords->data;

  // arrangement of f is as:
  // c0[0], c0[1], ... c1[0], c1[1], ....
  // where c0, c1, ... are components of f (num_ret_vals)
  int offset = 0;
  for (int n=0; n<num_ret_vals; ++n) {
    for (int k=0; k<num_basis; ++k) f[offset+k] = 0.0;

    for (int imu=0; imu<tot_quad; ++imu) {
      double tmp = weights[imu]*func_at_ords[n+num_ret_vals*imu];
      for (int k=0; k<num_basis; ++k)
        f[offset+k] += tmp*basis_at_ords[k+num_basis*imu];
    }
    offset += num_basis;
  }
}

static void
gkyl_fv_proj_priv_advance(const struct gkyl_fv_proj *up,
  double tm, const struct gkyl_range *update_range, struct gkyl_array *arr)
{
  double xc[GKYL_MAX_DIM], xmu[GKYL_MAX_DIM];

  int num_ret_vals = up->num_ret_vals;
  int tot_quad = up->tot_quad;
  struct gkyl_array *fun_at_ords = gkyl_array_new(GKYL_DOUBLE, num_ret_vals, tot_quad);
  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  
  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(&up->grid, iter.idx, xc);

    for (int i=0; i<tot_quad; ++i) {
      log_to_comp(up->grid.ndim, gkyl_array_cfetch(up->ordinates, i),
        up->grid.dx, xc, xmu);
      up->c2p(xmu, xmu, up->c2p_ctx);
      up->eval(tm, xmu, gkyl_array_fetch(fun_at_ords, i), up->ctx);
    }

    long lidx = gkyl_range_idx(update_range, iter.idx);
    gkyl_fv_proj_priv_quad(up, fun_at_ords, gkyl_array_fetch(arr, lidx));
  }

  gkyl_array_release(fun_at_ords);
}

static void
gkyl_fv_proj_priv_release(struct gkyl_fv_proj* up)
{
  gkyl_array_release(up->ordinates);
  gkyl_array_release(up->weights);
  gkyl_array_release(up->basis_at_ords);
  gkyl_free(up);
}

// functions below are not static
gkyl_fv_proj*
gkyl_fv_proj_new(const struct gkyl_rect_grid *grid,
  int num_quad, int num_ret_vals, evalf_t eval, void *ctx)
{
  return gkyl_fv_proj_priv_new(grid, num_quad, num_ret_vals, eval, ctx);
}

void
gkyl_fv_proj_advance(const gkyl_fv_proj *pob,
  double tm, const struct gkyl_range *update_rng, struct gkyl_array *out)
{
  gkyl_fv_proj_priv_advance(pob, tm, update_rng, out);

  // from projections, compute cell average
  double denorm = 1.0/sqrt(pow(2, update_rng->ndim));
  gkyl_array_scale_range(out, denorm, update_rng);
}

void
gkyl_fv_proj_release(gkyl_fv_proj* pob)
{
  gkyl_fv_proj_priv_release(pob);
}
// ended inlining fv_proj.c 
// start inlining gauss_quad_data.c 
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
// skipping file: gkyl_gauss_quad_data.h 

#define GUASS_QUAD_EPS 3.0e-15

// This is based on an implementation in Numerical Recipes in C book
static void
priv_gkyl_gauleg( double x1, double x2,  double x[], double w[], int n)
{
  double z1, xm, xl, pp, p3, p2, p1;
  int m = (n+1)/2;
  
  xm = 0.5*(x2+x1);
  xl = 0.5*(x2-x1);
  for (int i = 1; i <= m; i++) {
    double z = cos(GKYL_PI*(i-0.25)/(n+0.5));
    do {
      p1 = 1.0;
      p2 = 0.0;
      for (int j = 1; j <= n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp = n*(z*p1-p2)/(z*z-1.0);
      z1 = z;
      z = z1-p1/pp;
    } while( fabs(z-z1) > GUASS_QUAD_EPS );
    x[i] = xm-xl*z;
    x[n+1-i] = xm+xl*z;
    w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
    w[n+1-i] = w[i];
  }
}

void
gkyl_gauleg(double x1, double x2,  double x[], double w[], int n)
{
  priv_gkyl_gauleg(x1, x2, x-1, w-1, n); // actual routine assumes 1-offset arrays
}
// ended inlining gauss_quad_data.c 
// start inlining job_pool.c 
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

bool
gkyl_job_pool_add_work(const struct gkyl_job_pool *jp, jp_work_func func, void *ctx)
{
  return jp->add_work(jp, func, ctx);
}

void
gkyl_job_pool_wait(const struct gkyl_job_pool *jp)
{
  jp->wait(jp);
}

struct gkyl_job_pool*
gkyl_job_pool_acquire(const struct gkyl_job_pool *jp)
{
  gkyl_ref_count_inc(&jp->ref_count);
  return (struct gkyl_job_pool*) jp;
}

void
gkyl_job_pool_release(const struct gkyl_job_pool* jp)
{
  gkyl_ref_count_dec(&jp->ref_count);
}


// ended inlining job_pool.c 
// start inlining null_comm.c 
// skipping file: gkyl_alloc.h 
// skipping file: gkyl_array_rio.h 
// skipping file: gkyl_comm_priv.h 
// skipping file: gkyl_elem_type_priv.h 
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
// start inlining gkyl_null_comm_priv.h 

// Private header for mpi_comm. Do not include in user-facing header files.

// skipping file: gkyl_null_comm.h 
// skipping file: gkyl_alloc.h 

// ranges for use in BCs
struct skin_ghost_ranges {
  struct gkyl_range lower_skin[GKYL_MAX_DIM];
  struct gkyl_range lower_ghost[GKYL_MAX_DIM];

  struct gkyl_range upper_skin[GKYL_MAX_DIM];
  struct gkyl_range upper_ghost[GKYL_MAX_DIM];

  long max_vol; // maximum vol of send/recv region
};

// define long -> skin_ghost_ranges ...
#define i_key long
#define i_val struct skin_ghost_ranges
#define i_tag l2sgr
// start inlining stc/cmap.h 
/* MIT License
 *
 * Copyright (c) 2022 Tyge Løvset, NORCE, www.norceresearch.no
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// Unordered set/map - implemented as closed hashing with linear probing and no tombstones.
/*
#include <stdio.h>

#define i_tag ichar  // Map int => char
#define i_key int
#define i_val char
// skipping file: stc/cmap.h 

int main(void) {
    c_autovar (cmap_ichar m = cmap_ichar_init(), cmap_ichar_drop(&m))
    {
        cmap_ichar_emplace(&m, 5, 'a');
        cmap_ichar_emplace(&m, 8, 'b');
        cmap_ichar_emplace(&m, 12, 'c');

        cmap_ichar_value* v = cmap_ichar_get(&m, 10); // NULL
        char val = *cmap_ichar_at(&m, 5);               // 'a'
        cmap_ichar_emplace_or_assign(&m, 5, 'd');       // update
        cmap_ichar_erase(&m, 8);

        c_foreach (i, cmap_ichar, m)
            printf("map %d: %c\n", i.ref->first, i.ref->second);
    }
}
*/
// start inlining ccommon.h 
/* MIT License
 *
 * Copyright (c) 2022 Tyge Løvset, NORCE, www.norceresearch.no
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef CCOMMON_H_INCLUDED
#define CCOMMON_H_INCLUDED

#define _CRT_SECURE_NO_WARNINGS
#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>

#if defined(_MSC_VER)
#  pragma warning(disable: 4116 4996) // unnamed type definition in parentheses
#  define STC_FORCE_INLINE static __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#  define STC_FORCE_INLINE static inline __attribute((always_inline))
#else
#  define STC_FORCE_INLINE static inline
#endif
#define STC_INLINE static inline

/* Macro overloading feature support based on: https://rextester.com/ONP80107 */
#define c_MACRO_OVERLOAD(name, ...) \
    c_PASTE(name, c_NUM_ARGS(__VA_ARGS__))(__VA_ARGS__)
#define c_CONCAT(a, b) a ## b
#define c_PASTE(a, b) c_CONCAT(a, b)
#define c_EXPAND(...) __VA_ARGS__
#define c_NUM_ARGS(...) _c_APPLY_ARG_N((__VA_ARGS__, _c_RSEQ_N))

#define _c_APPLY_ARG_N(args) c_EXPAND(_c_ARG_N args)
#define _c_RSEQ_N 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
#define _c_ARG_N(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, \
                 _14, _15, _16, N, ...) N

#define c_static_assert(cond) \
    typedef char c_PASTE(_static_assert_line_, __LINE__)[(cond) ? 1 : -1]
#define c_container_of(ptr, type, member) \
    ((type *)((char *)(ptr) - offsetof(type, member)))

#ifndef __cplusplus
#  define c_alloc(T)            c_malloc(sizeof(T))
#  define c_alloc_n(T, n)       c_malloc(sizeof(T)*(n))
#  define c_make(T)             (T)
#  define c_new(T, ...)         (T*)memcpy(c_alloc(T), (T[]){__VA_ARGS__}, sizeof(T))
#else
#  include <new>
#  define c_alloc(T)            static_cast<T*>(c_malloc(sizeof(T)))
#  define c_alloc_n(T, n)       static_cast<T*>(c_malloc(sizeof(T)*(n)))
#  define c_make(T)             T
#  define c_new(T, ...)         new (c_alloc(T)) T(__VA_ARGS__)
#endif
#ifndef c_malloc
#  define c_malloc(sz)          malloc(sz)
#  define c_calloc(n, sz)       calloc(n, sz)
#  define c_realloc(p, sz)      realloc(p, sz)
#  define c_free(p)             free(p)
#endif

typedef const char              c_strlit[];
#define c_delete(T, ptr)        do { T *_c_p = ptr; T##_drop(_c_p); c_free(_c_p); } while (0)
#define c_swap(T, x, y)         do { T _c_t = x; x = y; y = _c_t; } while (0)
#define c_arraylen(a)           (sizeof (a)/sizeof (a)[0])
#define c_less_cmp(less, x, y)  (less(y, x) - less(x, y))
#define c_default_less(x, y)    (*(x) < *(y))

#define c_default_cmp(x, y)     c_less_cmp(c_default_less, x, y)
#define c_default_eq(x, y)      (*(x) == *(y))
#define c_memcmp_eq(x, y)       (memcmp(x, y, sizeof *(x)) == 0)

#define c_default_from(x)       (x)
#define c_default_toraw(ptr)    (*(ptr))
#define c_default_drop(ptr)     ((void) (ptr))

#define c_option(flag)          ((i_opt) & (flag))
#define c_is_fwd                1
#define c_no_atomic             2
#define c_no_clone              4
#define c_no_cmp                8
#define c_static                16
#define c_header                32
#define c_implement             64

/* Generic algorithms */

typedef const char* crawstr;
#define crawstr_cmp(xp, yp) strcmp(*(xp), *(yp))
#define crawstr_eq(xp, yp) (!strcmp(*(xp), *(yp)))
#define crawstr_hash(p, dummy) c_strhash(*(p))

#define _c_ROTL(x, k) (x << (k) | x >> (8*sizeof(x) - (k)))

STC_INLINE uint64_t c_strhash(const char *s) {
    int c; uint64_t h = *s++;
    if (h) while ((c = *s++)) h = (h << 10) - h + c;
    return _c_ROTL(h, 26) ^ h;
}
STC_INLINE uint64_t c_default_hash(const void* key, size_t len) {
    if (!len) return 1;
    const uint8_t *x = (const uint8_t*) key; 
    uint64_t h = *x++;
    while (--len) h = (h << 10) - h + *x++;
    return _c_ROTL(h, 26) ^ h;
}
STC_INLINE uint64_t c_hash32(const void* key, size_t len) {
    uint32_t x; memcpy(&x, key, 4);
    return x*0xc6a4a7935bd1e99d >> 15;
}
STC_INLINE uint64_t c_hash64(const void* key, size_t len) {
    uint64_t x; memcpy(&x, key, 8);
    return x*0xc6a4a7935bd1e99d;
}

STC_INLINE char* c_strnstrn(const char *s, const char *needle, size_t slen, const size_t nlen) {
    if (!nlen) return (char *)s;
    if (nlen > slen) return NULL;
    slen -= nlen;
    do {
        if (*s == *needle && !memcmp(s, needle, nlen)) 
            return (char *)s;
        ++s;
    } while (slen--);
    return NULL;
}

#define c_foreach(...) c_MACRO_OVERLOAD(c_foreach, __VA_ARGS__)

#define c_foreach3(it, C, cnt) \
    for (C##_iter it = C##_begin(&cnt), it##_end_ = C##_end(&cnt) \
         ; it.ref != it##_end_.ref; C##_next(&it))

#define c_foreach4(it, C, start, finish) \
    for (C##_iter it = start, it##_end_ = finish \
         ; it.ref != it##_end_.ref; C##_next(&it))

#define c_forpair(key, val, C, cnt) /* structured binding */ \
    for (struct {C##_iter _it; C##_value* _endref; C##_key key; C##_mapped val;} \
         _ = {C##_begin(&cnt), C##_end(&cnt).ref} \
         ; _._it.ref != _._endref && (_.key = _._it.ref->first, _.val = _._it.ref->second, true) \
         ; C##_next(&_._it))

#define c_forrange(...) c_MACRO_OVERLOAD(c_forrange, __VA_ARGS__)
#define c_forrange1(stop) for (size_t _c_ii=0, _c_end=stop; _c_ii < _c_end; ++_c_ii)
#define c_forrange2(i, stop) for (size_t i=0, _c_end=stop; i < _c_end; ++i)
#define c_forrange3(i, type, stop) for (type i=0, _c_end=stop; i < _c_end; ++i)
#define c_forrange4(i, type, start, stop) for (type i=start, _c_end=stop; i < _c_end; ++i)
#define c_forrange5(i, type, start, stop, step) \
    for (type i=start, _c_inc=step, _c_end=(stop) - (0 < _c_inc) \
         ; (i <= _c_end) == (0 < _c_inc); i += _c_inc)

#define c_autovar(declvar, ...) for (declvar, **_c_ii = NULL; !_c_ii; ++_c_ii, __VA_ARGS__)
#define c_autoscope(init, ...) for (int _c_ii = (init, 0); !_c_ii; ++_c_ii, __VA_ARGS__)
#define c_autodefer(...) for (int _c_ii = 0; !_c_ii; ++_c_ii, __VA_ARGS__)
#define c_breakauto continue

#define c_auto(...) c_MACRO_OVERLOAD(c_auto, __VA_ARGS__)
#define c_auto2(C, a) \
    c_autovar(C a = C##_init(), C##_drop(&a))
#define c_auto3(C, a, b) \
    c_autovar(c_EXPAND(C a = C##_init(), b = C##_init()), \
              C##_drop(&b), C##_drop(&a))
#define c_auto4(C, a, b, c) \
    c_autovar(c_EXPAND(C a = C##_init(), b = C##_init(), c = C##_init()), \
              C##_drop(&c), C##_drop(&b), C##_drop(&a))
#define c_auto5(C, a, b, c, d) \
    c_autovar(c_EXPAND(C a = C##_init(), b = C##_init(), c = C##_init(), d = C##_init()), \
              C##_drop(&d), C##_drop(&c), C##_drop(&b), C##_drop(&a))

#define c_autobuf(b, type, n) c_autobuf_N(b, type, n, 256)
#define c_autobuf_N(b, type, n, BYTES) \
    for (type _c_b[((BYTES) - 1) / sizeof(type) + 1], \
         *b = (n)*sizeof *b > (BYTES) ? c_alloc_n(type, n) : _c_b \
         ; b; b != _c_b ? c_free(b) : (void)0, b = NULL)

#define c_apply(v, method, T, ...) do { \
    T _c_arr[] = __VA_ARGS__; \
    for (size_t index = 0; index < c_arraylen(_c_arr); ++index) \
        { T v = _c_arr[index]; method; } \
} while (0)
#define c_apply_arr(v, method, T, arr, n) do { \
    T* _c_arr = arr; size_t _n = n; \
    for (size_t index = 0; index < _n; ++index) \
        { T v = _c_arr[index]; method; } \
} while (0)
#define c_apply_cnt(v, method, C, ...) do { \
    size_t index = 0; \
    c_foreach (_it, C, __VA_ARGS__) \
        { const C##_value v = *_it.ref; method; ++index; } \
} while (0)
#define c_pair(v) (v).first, (v).second

#define c_drop(C, ...) do { \
    C* _c_arr[] = {__VA_ARGS__}; \
    for (size_t _c_i = 0; _c_i < c_arraylen(_c_arr); ++_c_i) \
        C##_drop(_c_arr[_c_i]); \
} while (0)

#if defined(__SIZEOF_INT128__)
    #define c_umul128(a, b, lo, hi) \
        do { __uint128_t _z = (__uint128_t)(a)*(b); \
             *(lo) = (uint64_t)_z, *(hi) = _z >> 64; } while(0)
#elif defined(_MSC_VER) && defined(_WIN64)
    #include <intrin.h>
    #define c_umul128(a, b, lo, hi) ((void)(*(lo) = _umul128(a, b, hi)))
#elif defined(__x86_64__)
    #define c_umul128(a, b, lo, hi) \
        asm("mulq %3" : "=a"(*(lo)), "=d"(*(hi)) : "a"(a), "rm"(b))
#endif
#endif // CCOMMON_H_INCLUDED

#undef STC_API
#undef STC_DEF
#undef _i_static
#undef _i_implement

#if !c_option(c_static) && (c_option(c_header) || c_option(c_implement) || \
                            defined(STC_HEADER) || defined(STC_IMPLEMENTATION))
#  define STC_API extern
#  define STC_DEF
#else
#  define _i_static
#  define STC_API static inline
#  define STC_DEF static inline
#endif
#if (c_option(c_implement) || defined(STC_IMPLEMENTATION)) ^ defined(_i_static)
#  define _i_implement
#endif
// ended inlining ccommon.h 

#ifndef CMAP_H_INCLUDED
// start inlining forward.h 
/* MIT License
 *
 * Copyright (c) 2022 Tyge Løvset, NORCE, www.norceresearch.no
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef STC_FORWARD_H_INCLUDED
#define STC_FORWARD_H_INCLUDED

#include <stddef.h>

#define forward_carr2(CX, VAL) _c_carr2_types(CX, VAL)
#define forward_carr3(CX, VAL) _c_carr3_types(CX, VAL)
#define forward_cdeq(CX, VAL) _c_cdeq_types(CX, VAL)
#define forward_clist(CX, VAL) _c_clist_types(CX, VAL)
#define forward_cmap(CX, KEY, VAL) _c_chash_types(CX, KEY, VAL, c_true, c_false)
#define forward_csmap(CX, KEY, VAL) _c_aatree_types(CX, KEY, VAL, c_true, c_false)
#define forward_cset(CX, KEY) _c_chash_types(CX, cset, KEY, KEY, c_false, c_true)
#define forward_csset(CX, KEY) _c_aatree_types(CX, KEY, KEY, c_false, c_true)
#define forward_cbox(CX, VAL) _c_cbox_types(CX, VAL)
#define forward_carc(CX, VAL) _c_carc_types(CX, VAL)
#define forward_cpque(CX, VAL) _c_cpque_types(CX, VAL)
#define forward_cstack(CX, VAL) _c_cstack_types(CX, VAL)
#define forward_cqueue(CX, VAL) _c_cdeq_types(CX, VAL)
#define forward_cvec(CX, VAL) _c_cvec_types(CX, VAL)

typedef struct cstr { char* str; } cstr;
typedef char cstr_value;

typedef struct csview { const char* str; size_t size; } csview;
typedef union csview_iter { const char *ref; csview cp; } csview_iter;
typedef char csview_value;

#ifndef MAP_SIZE_T
#define MAP_SIZE_T uint32_t
#endif
#define c_true(...) __VA_ARGS__
#define c_false(...)

#define _c_carr2_types(SELF, VAL) \
    typedef VAL SELF##_value; \
    typedef struct { SELF##_value *ref; } SELF##_iter; \
    typedef struct { SELF##_value **data; size_t xdim, ydim; } SELF

#define _c_carr3_types(SELF, VAL) \
    typedef VAL SELF##_value; \
    typedef struct { SELF##_value *ref; } SELF##_iter; \
    typedef struct { SELF##_value ***data; size_t xdim, ydim, zdim; } SELF

#define _c_cdeq_types(SELF, VAL) \
    typedef VAL SELF##_value; \
    typedef struct {SELF##_value *ref; } SELF##_iter; \
    typedef struct {SELF##_value *_base, *data;} SELF

#define _c_clist_types(SELF, VAL) \
    typedef VAL SELF##_value; \
    typedef struct SELF##_node SELF##_node; \
\
    typedef struct { \
        SELF##_value *ref; \
        SELF##_node *const *_last, *prev; \
    } SELF##_iter; \
\
    typedef struct { \
        SELF##_node *last; \
    } SELF

#define _c_chash_types(SELF, KEY, VAL, MAP_ONLY, SET_ONLY) \
    typedef KEY SELF##_key; \
    typedef VAL SELF##_mapped; \
    typedef MAP_SIZE_T SELF##_size_t; \
\
    typedef SET_ONLY( SELF##_key ) \
            MAP_ONLY( struct SELF##_value ) \
    SELF##_value; \
\
    typedef struct { \
        SELF##_value *ref; \
        bool inserted; \
    } SELF##_result; \
\
    typedef struct { \
        SELF##_value *ref; \
        uint8_t* _hx; \
    } SELF##_iter; \
\
    typedef struct { \
        SELF##_value* table; \
        uint8_t* _hashx; \
        SELF##_size_t size, bucket_count; \
        float max_load_factor; \
    } SELF

#define _c_aatree_types(SELF, KEY, VAL, MAP_ONLY, SET_ONLY) \
    typedef KEY SELF##_key; \
    typedef VAL SELF##_mapped; \
    typedef MAP_SIZE_T SELF##_size_t; \
    typedef struct SELF##_node SELF##_node; \
\
    typedef SET_ONLY( SELF##_key ) \
            MAP_ONLY( struct SELF##_value ) \
    SELF##_value; \
\
    typedef struct { \
        SELF##_value *ref; \
        bool inserted; \
    } SELF##_result; \
\
    typedef struct { \
        SELF##_value *ref; \
        SELF##_node *_d; \
        int _top; \
        SELF##_size_t _tn, _st[36]; \
    } SELF##_iter; \
\
    typedef struct { \
        SELF##_node *nodes; \
    } SELF

#define _c_cbox_types(SELF, VAL) \
    typedef VAL SELF##_value; \
    typedef struct { \
        SELF##_value* get; \
    } SELF

#define _c_carc_types(SELF, VAL) \
    typedef VAL SELF##_value; \
\
    typedef struct { \
        SELF##_value* get; \
        long* use_count; \
    } SELF

#define _c_cstack_types(SELF, VAL) \
    typedef VAL SELF##_value; \
    typedef struct { SELF##_value *ref; } SELF##_iter; \
    typedef struct SELF { \
        SELF##_value* data; \
        size_t size, capacity; \
    } SELF

#define _c_cpque_types(SELF, VAL) \
    typedef VAL SELF##_value; \
    typedef struct SELF { \
        SELF##_value* data; \
        size_t size, capacity; \
    } SELF

#define _c_cvec_types(SELF, VAL) \
    typedef VAL SELF##_value; \
    typedef struct { SELF##_value *ref; } SELF##_iter; \
    typedef struct { SELF##_value *data; } SELF

#endif // STC_FORWARD_H_INCLUDED
// ended inlining forward.h 
#include <stdlib.h>
#include <string.h>
#define _cmap_inits {NULL, NULL, 0, 0, 0.85f}
typedef struct      { MAP_SIZE_T idx; uint_fast8_t hx; } chash_bucket_t;
#endif // CMAP_H_INCLUDED

#ifndef _i_prefix
#define _i_prefix cmap_
#endif
#ifdef _i_isset
  #define _i_MAP_ONLY c_false
  #define _i_SET_ONLY c_true
  #define _i_keyref(vp) (vp)
#else
  #define _i_MAP_ONLY c_true
  #define _i_SET_ONLY c_false
  #define _i_keyref(vp) (&(vp)->first)
#endif
// start inlining template.h 
/* MIT License
 *
 * Copyright (c) 2022 Tyge Løvset, NORCE, www.norceresearch.no
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef _i_template
#define _i_template

#ifndef STC_TEMPLATE_H_INCLUDED
#define STC_TEMPLATE_H_INCLUDED
  #define _cx_self c_PASTE(_i_prefix, i_tag)
  #define _cx_memb(name) c_PASTE(_cx_self, name)
  #define _cx_deftypes(macro, SELF, ...) c_EXPAND(macro(SELF, __VA_ARGS__))
  #define _cx_value _cx_memb(_value)
  #define _cx_key _cx_memb(_key)
  #define _cx_mapped _cx_memb(_mapped)
  #define _cx_raw _cx_memb(_raw)
  #define _cx_rawkey _cx_memb(_rawkey)
  #define _cx_rawmapped _cx_memb(_rawmapped)
  #define _cx_iter _cx_memb(_iter)
  #define _cx_result _cx_memb(_result)
  #define _cx_node _cx_memb(_node)
  #define _cx_size _cx_memb(_size_t)
#endif

#if defined i_cnt || defined i_equ // [deprecated]
  #define i_type i_cnt
  #error "i_cnt and i_equ no longer supported: use new name i_type / i_eq"
#endif

#ifdef i_type
  #define i_tag i_type
  #undef _i_prefix
  #define _i_prefix
#endif

#if defined i_key_str || defined i_val_str
// start inlining cstr.h 
/* MIT License
 *
 * Copyright (c) 2022 Tyge Løvset, NORCE, www.norceresearch.no
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
#ifndef CSTR_H_INCLUDED
#define CSTR_H_INCLUDED

// skipping file: ccommon.h 
// skipping file: forward.h 
#include <stdlib.h> /* malloc */
#include <string.h>
#include <stdarg.h>
#include <stdio.h> /* vsnprintf */
#include <ctype.h>

#define cstr_npos (SIZE_MAX >> 1)
typedef struct { size_t size, cap; char str[]; } _cstr_rep_t;
#define _cstr_rep(self) c_container_of((self)->str, _cstr_rep_t, str)
#ifdef _i_static 
    static const cstr cstr_null;
#else
    extern const cstr cstr_null;
#endif

/* optimal memory: based on malloc_usable_size() sequence: 24, 40, 56, ... */
#define _cstr_opt_mem(cap)  ((((offsetof(_cstr_rep_t, str) + (cap) + 8)>>4)<<4) + 8)
/* optimal string capacity: 7, 23, 39, ... */
#define _cstr_opt_cap(cap)  (_cstr_opt_mem(cap) - offsetof(_cstr_rep_t, str) - 1)

STC_API cstr            cstr_from_n(const char* str, size_t n);
STC_API cstr            cstr_from_fmt(const char* fmt, ...);
STC_API cstr            cstr_from_replace_all(const char* str, size_t str_len,
                                              const char* find, size_t find_len,
                                              const char* repl, size_t repl_len);
STC_API size_t          cstr_reserve(cstr* self, size_t cap);
STC_API void            cstr_resize(cstr* self, size_t len, char fill);
STC_API cstr*           cstr_assign_n(cstr* self, const char* str, size_t n);
STC_API int             cstr_printf(cstr* self, const char* fmt, ...);
STC_API cstr*           cstr_append_n(cstr* self, const char* str, size_t n);
STC_API void            cstr_replace_n(cstr* self, size_t pos, size_t len, const char* str, size_t n);
STC_API void            cstr_replace_all(cstr* self, const char* find, const char* replace);
STC_API void            cstr_erase_n(cstr* self, size_t pos, size_t n);
STC_API size_t          cstr_find(cstr s, const char* needle);
STC_API size_t          cstr_find_n(cstr s, const char* needle, size_t pos, size_t nmax);
STC_API bool            cstr_getdelim(cstr *self, int delim, FILE *stream);

STC_API char*           c_strnstrn(const char* s, const char* needle, size_t slen, size_t nlen);

STC_INLINE cstr         cstr_init() { return cstr_null; }
#define                 cstr_str(self) (self)->str
#define                 cstr_toraw(self) (self)->str
#define                 cstr_new(literal) \
                            cstr_from_n(literal, sizeof c_make(c_strlit){literal} - 1)
STC_INLINE cstr         cstr_from(const char* str)
                            { return cstr_from_n(str, strlen(str)); }
STC_INLINE char*        cstr_data(cstr* self) { return self->str; }
STC_INLINE size_t       cstr_size(cstr s) { return _cstr_rep(&s)->size; }
STC_INLINE size_t       cstr_length(cstr s) { return _cstr_rep(&s)->size; }
STC_INLINE size_t       cstr_capacity(cstr s) { return _cstr_rep(&s)->cap; }
STC_INLINE bool         cstr_empty(cstr s) { return _cstr_rep(&s)->size == 0; }
STC_INLINE void         cstr_drop(cstr* self)
                            { if (_cstr_rep(self)->cap) c_free(_cstr_rep(self)); }
STC_INLINE cstr         cstr_clone(cstr s)
                            { return cstr_from_n(s.str, _cstr_rep(&s)->size); }
STC_INLINE void         cstr_clear(cstr* self)
                            { self->str[_cstr_rep(self)->size = 0] = '\0'; }
STC_INLINE cstr*        cstr_assign(cstr* self, const char* str)
                            { return cstr_assign_n(self, str, strlen(str)); }
STC_INLINE cstr*        cstr_copy(cstr* self, cstr s)
                            { return cstr_assign_n(self, s.str, _cstr_rep(&s)->size); }
STC_INLINE cstr*        cstr_append(cstr* self, const char* str)
                            { return cstr_append_n(self, str, strlen(str)); }
STC_INLINE cstr*        cstr_append_s(cstr* self, cstr s)
                            { return cstr_append_n(self, s.str, _cstr_rep(&s)->size); }
STC_INLINE void         cstr_push_back(cstr* self, char value)
                            { cstr_append_n(self, &value, 1); }
STC_INLINE void         cstr_pop_back(cstr* self)
                            { self->str[ --_cstr_rep(self)->size ] = '\0'; }
STC_INLINE void         cstr_insert_n(cstr* self, const size_t pos, const char* str, const size_t n)
                            { cstr_replace_n(self, pos, 0, str, n); }
STC_INLINE void         cstr_insert(cstr* self, const size_t pos, const char* str)
                            { cstr_replace_n(self, pos, 0, str, strlen(str)); }
STC_INLINE void         cstr_insert_s(cstr* self, const size_t pos, cstr s)
                            { cstr_replace_n(self, pos, 0, s.str, _cstr_rep(&s)->size); }
STC_INLINE void         cstr_replace(cstr* self, const size_t pos, const size_t len, const char* str)
                            { cstr_replace_n(self, pos, len, str, strlen(str)); }
STC_INLINE void         cstr_replace_s(cstr* self, const size_t pos, const size_t len, cstr s)
                            { cstr_replace_n(self, pos, len, s.str, _cstr_rep(&s)->size); }
STC_INLINE void         cstr_erase(cstr* self, const size_t pos)
                            { cstr_erase_n(self, pos, 1); }
STC_INLINE char*        cstr_front(cstr* self) { return self->str; }
STC_INLINE char*        cstr_back(cstr* self)
                            { return self->str + _cstr_rep(self)->size - 1; }
STC_INLINE bool         cstr_equals(cstr s, const char* str)
                            { return strcmp(s.str, str) == 0; }
STC_INLINE bool         cstr_equals_s(cstr s1, cstr s2)
                            { return strcmp(s1.str, s2.str) == 0; }
STC_INLINE bool         cstr_contains(cstr s, const char* needle)
                            { return strstr(s.str, needle) != NULL; }
STC_INLINE bool         cstr_getline(cstr *self, FILE *stream)
                            { return cstr_getdelim(self, '\n', stream); }

STC_INLINE cstr
cstr_with_capacity(const size_t cap) {
    cstr s = cstr_null;
    cstr_reserve(&s, cap);
    return s;
}

STC_INLINE cstr
cstr_with_size(const size_t len, const char fill) {
    cstr s = cstr_null;
    cstr_resize(&s, len, fill);
    return s;
}

STC_INLINE cstr*
cstr_take(cstr* self, cstr s) {
    if (self->str != s.str && _cstr_rep(self)->cap)
        c_free(_cstr_rep(self));
    self->str = s.str;
    return self;
}

STC_INLINE cstr
cstr_move(cstr* self) {
    cstr tmp = *self;
    *self = cstr_null;
    return tmp;
}

STC_INLINE bool
cstr_starts_with(cstr s, const char* sub) {
    while (*sub && *s.str == *sub) ++s.str, ++sub;
    return *sub == 0;
}

STC_INLINE bool
cstr_ends_with(cstr s, const char* sub) {
    const size_t n = strlen(sub), sz = _cstr_rep(&s)->size;
    return n <= sz && !memcmp(s.str + sz - n, sub, n);
}

STC_INLINE int
c_strncasecmp(const char* s1, const char* s2, size_t nmax) {
    int ret = 0;
    while (nmax-- && (ret = tolower(*s1++) - tolower(*s2)) == 0 && *s2++)
        ;
    return ret;
}

/* container adaptor functions: */
#define  cstr_cmp(xp, yp)     strcmp((xp)->str, (yp)->str)
#define  cstr_eq(xp, yp)      (!cstr_cmp(xp, yp))
#define  cstr_hash(xp, dummy) c_strhash((xp)->str)

/* -------------------------- IMPLEMENTATION ------------------------- */
#if defined(_i_implement)

static struct { size_t size, cap; char str[1]; } _cstr_nullrep = {0, 0, {0}};
#ifdef _i_static
static
#endif
const cstr cstr_null = {_cstr_nullrep.str};

STC_DEF size_t
cstr_reserve(cstr* self, const size_t cap) {
    _cstr_rep_t* rep = _cstr_rep(self);
    const size_t oldcap = rep->cap;
    if (cap > oldcap) {
        rep = (_cstr_rep_t*) c_realloc(oldcap ? rep : NULL, _cstr_opt_mem(cap));
        self->str = rep->str;
        if (oldcap == 0) self->str[rep->size = 0] = '\0';
        return (rep->cap = _cstr_opt_cap(cap));
    }
    return oldcap;
}

STC_DEF void
cstr_resize(cstr* self, const size_t len, const char fill) {
    const size_t n =  _cstr_rep(self)->size;
    cstr_reserve(self, len);
    if (len > n) memset(self->str + n, fill, len - n);
    if (len | n) self->str[_cstr_rep(self)->size = len] = '\0';
}

STC_DEF cstr
cstr_from_n(const char* str, const size_t n) {
    if (n == 0) return cstr_null;
    _cstr_rep_t* rep = (_cstr_rep_t*) c_malloc(_cstr_opt_mem(n));
    rep->str[rep->size = n] = '\0';
    rep->cap = _cstr_opt_cap(n);
    cstr s = {(char *) memcpy(rep->str, str, n)};
    return s;
}

#if defined(__clang__)
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Wdeprecated-declarations"
#elif defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable: 4996)
#endif

STC_DEF int
cstr_vfmt(cstr* self, const char* fmt, va_list args) {
    va_list args2;
    va_copy(args2, args);
    int len = vsnprintf(NULL, (size_t)0, fmt, args);
    cstr_reserve(self, len);
    vsprintf(self->str, fmt, args2);
    va_end(args2);
    return _cstr_rep(self)->size = len;
}

#if defined(__clang__)
#  pragma clang diagnostic pop
#elif defined(_MSC_VER)
#  pragma warning(pop)
#endif

STC_DEF cstr
cstr_from_fmt(const char* fmt, ...) {
    cstr ret = cstr_null;
    va_list args; va_start(args, fmt);
    cstr_vfmt(&ret, fmt, args);
    va_end(args);
    return ret;
}

STC_DEF int
cstr_printf(cstr* self, const char* fmt, ...) {
    cstr ret = cstr_null;
    va_list args;
    va_start(args, fmt);
    int n = cstr_vfmt(&ret, fmt, args);
    va_end(args);
    cstr_drop(self);
    *self = ret;
    return n;
}

STC_DEF cstr*
cstr_assign_n(cstr* self, const char* str, const size_t n) {
    if (n || _cstr_rep(self)->cap) {
        cstr_reserve(self, n);
        memmove(self->str, str, n);
        self->str[_cstr_rep(self)->size = n] = '\0';
    }
    return self;
}

STC_DEF cstr*
cstr_append_n(cstr* self, const char* str, const size_t n) {
    if (n == 0) return self;
    const size_t oldlen = _cstr_rep(self)->size, newlen = oldlen + n;
    if (newlen > _cstr_rep(self)->cap) {
        const size_t off = (size_t) (str - self->str); /* handle self append */
        cstr_reserve(self, (oldlen*3 >> 1) + n);
        if (off <= oldlen) str = self->str + off;
    }
    memcpy(&self->str[oldlen], str, n);
    self->str[_cstr_rep(self)->size = newlen] = '\0';
    return self;
}

STC_INLINE void _cstr_internal_move(cstr* self, const size_t pos1, const size_t pos2) {
    if (pos1 == pos2)
        return;
    const size_t len = _cstr_rep(self)->size, newlen = len + pos2 - pos1;
    if (newlen > _cstr_rep(self)->cap)
        cstr_reserve(self, (len*3 >> 1) + pos2 - pos1);
    memmove(&self->str[pos2], &self->str[pos1], len - pos1);
    self->str[_cstr_rep(self)->size = newlen] = '\0';
}

STC_DEF void
cstr_replace_n(cstr* self, const size_t pos, size_t len, const char* str, const size_t n) {
    const size_t sz = cstr_size(*self);
    if (len > sz - pos) len = sz - pos;
    c_autobuf (xstr, char, n) {
        memcpy(xstr, str, n);
        _cstr_internal_move(self, pos + len, pos + n);
        memcpy(&self->str[pos], xstr, n);
    }
}

STC_DEF cstr
cstr_from_replace_all(const char* str, const size_t str_len,
                      const char* find, const size_t find_len,
                      const char* repl, const size_t repl_len) {
    cstr out = cstr_null;
    size_t from = 0; char* res;
    if (find_len)
        while ((res = c_strnstrn(str + from, find, str_len - from, find_len))) {
            const size_t pos = res - str;
            cstr_append_n(&out, str + from, pos - from);
            cstr_append_n(&out, repl, repl_len);
            from = pos + find_len;
        }
    cstr_append_n(&out, str + from, str_len - from);
    return out;
}

STC_DEF void
cstr_replace_all(cstr* self, const char* find, const char* repl) {
    cstr_take(self, cstr_from_replace_all(self->str, _cstr_rep(self)->size,
                                          find, strlen(find), repl, strlen(repl)));
}

STC_DEF void
cstr_erase_n(cstr* self, const size_t pos, size_t n) {
    const size_t len = _cstr_rep(self)->size;
    if (n > len - pos) n = len - pos;
    if (len) {
        memmove(&self->str[pos], &self->str[pos + n], len - (pos + n));
        self->str[_cstr_rep(self)->size -= n] = '\0';
    }
}

STC_DEF bool
cstr_getdelim(cstr *self, const int delim, FILE *fp) {
    size_t pos = 0, cap = _cstr_rep(self)->cap;
    int c = fgetc(fp);
    if (c == EOF)
        return false;
    for (;;) {
        if (c == delim || c == EOF) {
            if (cap) self->str[_cstr_rep(self)->size = pos] = '\0';
            return true;
        }
        if (pos == cap)
            cap = cstr_reserve(self, (cap*3 >> 1) + 16);
        self->str[pos++] = (char) c;
        c = fgetc(fp);
    }
}

STC_DEF size_t
cstr_find(cstr s, const char* needle) {
    char* res = strstr(s.str, needle);
    return res ? res - s.str : cstr_npos;
}

STC_DEF size_t
cstr_find_n(cstr s, const char* needle, const size_t pos, const size_t nmax) {
    if (pos > _cstr_rep(&s)->size) return cstr_npos;
    const size_t nlen = strlen(needle);
    char* res = c_strnstrn(s.str + pos, needle, _cstr_rep(&s)->size - pos, nmax < nlen ? nmax : nlen);
    return res ? res - s.str : cstr_npos;
}

#endif
#endif
#undef i_opt
// ended inlining cstr.h 
#endif

#ifdef i_key_str
  #define i_key_bind cstr
  #define i_keyraw crawstr
  #ifndef i_tag
    #define i_tag str
  #endif
#elif defined i_key_sptr
  #define i_key_bind i_key_sptr
  #define i_keyraw c_PASTE(i_key_sptr, _value)
#endif

#ifdef i_key_bind
  #define i_key i_key_bind
  #ifndef i_keyraw
    #ifndef i_keyfrom
      #define i_keyfrom c_PASTE(i_key, _clone)
    #endif
  #else
    #ifndef i_keyfrom
      #define i_keyfrom c_PASTE(i_key, _from)
    #endif
    #ifndef i_keyto
      #define i_keyto c_PASTE(i_key, _toraw)
    #endif
  #endif
  #ifndef i_cmp
    #define i_cmp c_PASTE(i_keyraw, _cmp)
  #endif
  #ifndef i_eq
    #define i_eq c_PASTE(i_keyraw, _eq)
  #endif
  #ifndef i_hash
    #define i_hash c_PASTE(i_keyraw, _hash)
  #endif
  #ifndef i_keydrop
    #define i_keydrop c_PASTE(i_key, _drop)
  #endif
#endif

#if defined i_keyraw && !(defined i_keyto && defined i_keyfrom)
  #error "if i_keyraw defined, i_keyfrom and i_keyto must be defined"
#endif

/* Resolve i_drop and i_from here */
#if defined i_drop && defined i_isset
  #define i_keydrop i_drop
#elif defined i_drop && !defined i_key
  #define i_valdrop i_drop
#elif defined i_drop
  #error "i_drop not supported for maps, define i_keydrop / i_valdrop instead."
#endif
#if defined i_from && defined i_isset
  #define i_keyfrom i_from
#elif defined i_from && !defined i_key
  #define i_valfrom i_from
#elif defined i_from
  #error "i_from not supported for maps, define i_keyfrom / i_valfrom instead."
#endif

#ifdef i_val_str
  #define i_val_bind cstr
  #define i_valraw crawstr
  #if !defined i_tag && !defined i_key
    #define i_tag str
  #endif
#elif defined i_val_sptr
  #define i_val_bind i_val_sptr
  #define i_valraw c_PASTE(i_val_sptr, _value)
#endif

#ifdef i_val_bind
  #define i_val i_val_bind
  #ifndef i_valraw
    #ifndef i_valfrom
      #define i_valfrom c_PASTE(i_val, _clone)
    #endif
  #else
    #ifndef i_valfrom
      #define i_valfrom c_PASTE(i_val, _from)
    #endif
    #ifndef i_valto
      #define i_valto c_PASTE(i_val, _toraw)
    #endif
  #endif
  #if !defined i_cmp && !defined i_key
    #define i_cmp c_PASTE(i_valraw, _cmp)
  #endif
  #ifndef i_valdrop
    #define i_valdrop c_PASTE(i_val, _drop)
  #endif
#endif

#if defined i_valraw && !(defined i_valto && defined i_valfrom)
  #error "if i_valraw defined, i_valfrom and i_valto must be defined"
#endif

#ifdef i_key
  #ifdef _i_isset
    #define i_val i_key
  #endif
  #ifndef i_tag
    #define i_tag i_key  
  #endif
  #if !defined _i_has_internal_clone && defined i_keydrop && !defined i_keyfrom && !c_option(c_no_clone)
    #error "i_keydrop defined but not i_keyfrom (e.g. as c_default_from), or no 'i_opt c_no_clone'"
  #endif
  #if !defined i_keyfrom
    #define i_keyfrom c_default_from
  #endif
  #ifndef i_keyraw
    #define _i_no_keyraw
    #define i_keyraw i_key
    #define i_keyto c_default_toraw
  #endif
  #ifndef i_keydrop
    #define i_keydrop c_default_drop
  #endif
#elif defined _i_isset
  #error "i_key define is missing."
#endif

#ifndef i_tag
  #define i_tag i_val
#endif
#if !defined _i_has_internal_clone && defined i_valdrop && !defined i_valfrom && !c_option(c_no_clone)
  #error "i_valdrop/i_drop defined but not i_valfrom (e.g. as c_default_from), or no 'i_opt c_no_clone'"
#endif
#if !defined i_valfrom
  #define i_valfrom c_default_from
#endif
#ifndef i_valraw
  #if !defined i_key || defined _i_no_keyraw
    #define _i_no_raw
  #endif
  #define i_valraw i_val
  #define i_valto c_default_toraw
#endif
#ifndef i_valdrop
  #define i_valdrop c_default_drop
#endif
#if !defined i_eq && defined i_cmp
  #define i_eq !i_cmp
#elif !defined i_eq
  #define i_eq c_default_eq
#endif
#ifndef i_cmp
  #define i_cmp c_default_cmp
#endif
#ifndef i_hash
  #define i_hash c_default_hash
#endif

#else // -------------------------------------------------------

#undef i_type
#undef i_tag
#undef i_imp
#undef i_opt
#undef i_cmp
#undef i_eq
#undef i_hash
#undef i_from
#undef i_drop

#undef i_val
#undef i_val_str
#undef i_val_sptr
#undef i_val_bind
#undef i_valraw
#undef i_valfrom
#undef i_valto
#undef i_valdrop

#undef i_key
#undef i_key_str
#undef i_key_sptr
#undef i_key_bind
#undef i_keyraw
#undef i_keyfrom
#undef i_keyto
#undef i_keydrop

#undef _i_prefix
#undef _i_no_raw
#undef _i_no_keyraw
#undef _i_has_internal_clone
#undef _i_template
#endif
// ended inlining template.h 
#if !c_option(c_is_fwd)
_cx_deftypes(_c_chash_types, _cx_self, i_key, i_val, _i_MAP_ONLY, _i_SET_ONLY);
#endif

_i_MAP_ONLY( struct _cx_value {
    _cx_key first;
    _cx_mapped second;
}; )

typedef i_keyraw _cx_rawkey;
typedef i_valraw _cx_memb(_rawmapped);
typedef _i_SET_ONLY( i_keyraw )
        _i_MAP_ONLY( struct { i_keyraw first;
                              i_valraw second; } )
_cx_raw;

STC_API _cx_self        _cx_memb(_with_capacity)(size_t cap);
#if !c_option(c_no_clone)
STC_API _cx_self        _cx_memb(_clone)(_cx_self map);
#endif
STC_API void            _cx_memb(_drop)(_cx_self* self);
STC_API void            _cx_memb(_clear)(_cx_self* self);
STC_API bool            _cx_memb(_reserve)(_cx_self* self, size_t capacity);
STC_API chash_bucket_t  _cx_memb(_bucket_)(const _cx_self* self, const _cx_rawkey* rkeyptr);
STC_API _cx_result      _cx_memb(_insert_entry_)(_cx_self* self, i_keyraw rkey);
STC_API void            _cx_memb(_erase_entry)(_cx_self* self, _cx_value* val);

STC_INLINE _cx_self     _cx_memb(_init)(void) { return c_make(_cx_self)_cmap_inits; }
STC_INLINE void         _cx_memb(_shrink_to_fit)(_cx_self* self) { _cx_memb(_reserve)(self, self->size); }
STC_INLINE void         _cx_memb(_max_load_factor)(_cx_self* self, float ml) {self->max_load_factor = ml; }
STC_INLINE bool         _cx_memb(_empty)(_cx_self m) { return m.size == 0; }
STC_INLINE size_t       _cx_memb(_size)(_cx_self m) { return m.size; }
STC_INLINE size_t       _cx_memb(_bucket_count)(_cx_self map) { return map.bucket_count; }
STC_INLINE size_t       _cx_memb(_capacity)(_cx_self map)
                            { return (size_t)(map.bucket_count ? (map.bucket_count - 2)*map.max_load_factor : 0.f); }
STC_INLINE void         _cx_memb(_swap)(_cx_self *map1, _cx_self *map2) {c_swap(_cx_self, *map1, *map2); }
STC_INLINE bool         _cx_memb(_contains)(const _cx_self* self, i_keyraw rkey)
                            { return self->size && self->_hashx[_cx_memb(_bucket_)(self, &rkey).idx]; }

#ifndef _i_isset
    #if !c_option(c_no_clone) && !defined _i_no_raw
    STC_API _cx_result  _cx_memb(_emplace_or_assign)(_cx_self* self, i_keyraw rkey, i_valraw rmapped);
    #endif
    STC_API _cx_result  _cx_memb(_insert_or_assign)(_cx_self* self, i_key _key, i_val _mapped);

    STC_INLINE _cx_result  /* short-form, like operator[]: */
    _cx_memb(_put)(_cx_self* self, i_key key, i_val mapped) {
        return _cx_memb(_insert_or_assign)(self, key, mapped);
    }

    STC_INLINE const _cx_mapped*
    _cx_memb(_at)(const _cx_self* self, i_keyraw rkey) {
        chash_bucket_t b = _cx_memb(_bucket_)(self, &rkey);
        assert(self->_hashx[b.idx]);
        return &self->table[b.idx].second;
    }
#endif

#if !c_option(c_no_clone)
STC_INLINE void _cx_memb(_copy)(_cx_self *self, _cx_self other) {
    if (self->table == other.table) return;
    _cx_memb(_drop)(self); *self = _cx_memb(_clone)(other);
}

STC_INLINE _cx_value
_cx_memb(_value_clone)(_cx_value _val) {
    *_i_keyref(&_val) = i_keyfrom(i_keyto(_i_keyref(&_val)));
    _i_MAP_ONLY( _val.second = i_valfrom(i_valto(&_val.second)); )
    return _val;
}

#if !defined _i_no_raw
STC_INLINE _cx_result
_cx_memb(_emplace)(_cx_self* self, i_keyraw rkey _i_MAP_ONLY(, i_valraw rmapped)) {
    _cx_result _res = _cx_memb(_insert_entry_)(self, rkey);
    if (_res.inserted) {
        *_i_keyref(_res.ref) = i_keyfrom(rkey);
        _i_MAP_ONLY( _res.ref->second = i_valfrom(rmapped); )
    }
    return _res;
}
#endif
#endif // !c_no_clone

STC_INLINE _cx_raw
_cx_memb(_value_toraw)(_cx_value* val) {
    return _i_SET_ONLY( i_keyto(val) )
           _i_MAP_ONLY( c_make(_cx_raw){i_keyto(&val->first), i_valto(&val->second)} );
}

STC_INLINE void
_cx_memb(_value_drop)(_cx_value* _val) {
    i_keydrop(_i_keyref(_val));
    _i_MAP_ONLY( i_valdrop(&_val->second); )
}

STC_INLINE _cx_result
_cx_memb(_insert)(_cx_self* self, i_key _key _i_MAP_ONLY(, i_val _mapped)) {
    _cx_result _res = _cx_memb(_insert_entry_)(self, i_keyto(&_key));
    if (_res.inserted) { *_i_keyref(_res.ref) = _key; _i_MAP_ONLY( _res.ref->second = _mapped; )}
    else               { i_keydrop(&_key); _i_MAP_ONLY( i_valdrop(&_mapped); )}
    return _res;
}

STC_INLINE _cx_iter
_cx_memb(_find)(const _cx_self* self, i_keyraw rkey) {
    _cx_size idx;
    if (!(self->size && self->_hashx[idx = _cx_memb(_bucket_)(self, &rkey).idx]))
        idx = self->bucket_count;
    return c_make(_cx_iter){self->table+idx, self->_hashx+idx};
}

STC_INLINE const _cx_value*
_cx_memb(_get)(const _cx_self* self, i_keyraw rkey) {
    _cx_size idx;
    return self->size && self->_hashx[idx = _cx_memb(_bucket_)(self, &rkey).idx] ?
           self->table + idx : NULL;
}

STC_INLINE _cx_value*
_cx_memb(_get_mut)(const _cx_self* self, i_keyraw rkey)
    { return (_cx_value*) _cx_memb(_get)(self, rkey); }

STC_INLINE _cx_iter
_cx_memb(_begin)(const _cx_self* self) {
    _cx_iter it = {self->table, self->_hashx};
    if (it._hx) while (*it._hx == 0) ++it.ref, ++it._hx;
    return it;
}

STC_INLINE _cx_iter
_cx_memb(_end)(const _cx_self* self)
    { return c_make(_cx_iter){self->table + self->bucket_count}; }

STC_INLINE void
_cx_memb(_next)(_cx_iter* it)
    { while ((++it->ref, *++it->_hx == 0)) ; }

STC_INLINE _cx_iter
_cx_memb(_advance)(_cx_iter it, size_t n) {
    // UB if n > elements left
    while (n--) _cx_memb(_next)(&it);
    return it;
}

STC_INLINE size_t
_cx_memb(_erase)(_cx_self* self, i_keyraw rkey) {
    if (self->size == 0) return 0;
    chash_bucket_t b = _cx_memb(_bucket_)(self, &rkey);
    return self->_hashx[b.idx] ? _cx_memb(_erase_entry)(self, self->table + b.idx), 1 : 0;
}

STC_INLINE _cx_iter
_cx_memb(_erase_at)(_cx_self* self, _cx_iter it) {
    _cx_memb(_erase_entry)(self, it.ref);
    if (*it._hx == 0) _cx_memb(_next)(&it);
    return it;
}

/* -------------------------- IMPLEMENTATION ------------------------- */
#if defined(_i_implement)

#ifndef CMAP_H_INCLUDED
//STC_INLINE size_t fastrange_uint64_t(uint64_t x, uint64_t n)
//    { uint64_t lo, hi; c_umul128(x, n, &lo, &hi); return hi; }
#define fastrange_uint32_t(x, n) (uint32_t)((uint32_t)(x)*(uint64_t)(n) >> 32)
#define chash_index_(h, entryPtr) ((entryPtr) - (h).table)
#endif // CMAP_H_INCLUDED

STC_DEF _cx_self
_cx_memb(_with_capacity)(const size_t cap) {
    _cx_self h = _cmap_inits;
    _cx_memb(_reserve)(&h, cap);
    return h;
}

STC_INLINE void _cx_memb(_wipe_)(_cx_self* self) {
    if (self->size == 0) return;
    _cx_value* e = self->table, *end = e + self->bucket_count;
    uint8_t *hx = self->_hashx;
    for (; e != end; ++e) if (*hx++) _cx_memb(_value_drop)(e);
}

STC_DEF void _cx_memb(_drop)(_cx_self* self) {
    _cx_memb(_wipe_)(self);
    c_free(self->_hashx);
    c_free((void *) self->table);
}

STC_DEF void _cx_memb(_clear)(_cx_self* self) {
    _cx_memb(_wipe_)(self);
    self->size = 0;
    memset(self->_hashx, 0, self->bucket_count);
}

#if !defined _i_isset
    STC_DEF _cx_result
    _cx_memb(_insert_or_assign)(_cx_self* self, i_key _key, i_val _mapped) {
        _cx_result _res = _cx_memb(_insert_entry_)(self, i_keyto(&_key));
        if (_res.inserted) _res.ref->first = _key;
        else { i_keydrop(&_key); i_valdrop(&_res.ref->second); }
        _res.ref->second = _mapped; return _res;
    }

    #if !c_option(c_no_clone) && !defined _i_no_raw
    STC_DEF _cx_result
    _cx_memb(_emplace_or_assign)(_cx_self* self, i_keyraw rkey, i_valraw rmapped) {
        _cx_result _res = _cx_memb(_insert_entry_)(self, rkey);
        if (_res.inserted) _res.ref->first = i_keyfrom(rkey);
        else i_valdrop(&_res.ref->second);
        _res.ref->second = i_valfrom(rmapped); return _res;
    }
    #endif
#endif

STC_DEF chash_bucket_t
_cx_memb(_bucket_)(const _cx_self* self, const _cx_rawkey* rkeyptr) {
    const uint64_t _hash = i_hash(rkeyptr, sizeof *rkeyptr);
    uint_fast8_t _hx; _cx_size _cap = self->bucket_count;
    chash_bucket_t b = {c_PASTE(fastrange_,MAP_SIZE_T)(_hash, _cap), (uint_fast8_t)(_hash | 0x80)};
    const uint8_t* _hashx = self->_hashx;
    while ((_hx = _hashx[b.idx])) {
        if (_hx == b.hx) {
            _cx_rawkey _raw = i_keyto(_i_keyref(self->table + b.idx));
            if (i_eq(&_raw, rkeyptr)) break;
        }
        _cx_size _mask = (_cx_size) -(++b.idx != _cap);
        b.idx &= _mask; // b.idx = (b.idx + 1) % _cap
    }
    return b;
}

STC_DEF _cx_result
_cx_memb(_insert_entry_)(_cx_self* self, i_keyraw rkey) {
    if (self->size + 1 >= (_cx_size) (self->bucket_count * self->max_load_factor))
        _cx_memb(_reserve)(self, ((size_t)self->size*3 >> 1) + 4);
    chash_bucket_t b = _cx_memb(_bucket_)(self, &rkey);
    _cx_result res = {&self->table[b.idx], !self->_hashx[b.idx]};
    if (res.inserted) {
        self->_hashx[b.idx] = b.hx;
        ++self->size;
    }
    return res;
}

#if !c_option(c_no_clone)
STC_DEF _cx_self
_cx_memb(_clone)(_cx_self m) {
    _cx_self clone = {
        c_alloc_n(_cx_value, m.bucket_count),
        (uint8_t *) memcpy(c_malloc(m.bucket_count + 1), m._hashx, m.bucket_count + 1),
        m.size, m.bucket_count,
        m.max_load_factor
    };
    _cx_value *e = m.table, *end = e + m.bucket_count, *dst = clone.table;
    for (uint8_t *hx = m._hashx; e != end; ++hx, ++e, ++dst)
        if (*hx) *dst = _cx_memb(_value_clone)(*e);
    return clone;
}
#endif

STC_DEF bool
_cx_memb(_reserve)(_cx_self* self, const size_t _newcap) {
    if (_newcap < self->size) return true;
    const _cx_size _oldbuckets = self->bucket_count;
    const _cx_size _nbuckets = ((_cx_size)(_newcap/self->max_load_factor) + 2) | 1;
    _cx_self _tmp = {
        c_alloc_n(_cx_value, _nbuckets),
        (uint8_t *) c_calloc(_nbuckets + 1, sizeof(uint8_t)),
        self->size, (_cx_size) _nbuckets,
        self->max_load_factor
    };
    bool ret; /* Rehash: */
    if ((ret = _tmp.table && _tmp._hashx)) {
        _tmp._hashx[_nbuckets] = 0xff;
        c_swap(_cx_self, *self, _tmp);
        _cx_value* e = _tmp.table, *_slot = self->table;
        uint8_t* _hashx = self->_hashx;
        for (size_t i = 0; i < _oldbuckets; ++i, ++e) if (_tmp._hashx[i]) {
            _cx_rawkey _raw = i_keyto(_i_keyref(e));
            chash_bucket_t b = _cx_memb(_bucket_)(self, &_raw);
            _slot[b.idx] = *e;
            _hashx[b.idx] = (uint8_t) b.hx;
        }
    }
    c_free(_tmp._hashx);
    c_free((void *) _tmp.table);
    return ret;
}

STC_DEF void
_cx_memb(_erase_entry)(_cx_self* self, _cx_value* _val) {
    _cx_size i = chash_index_(*self, _val), j = i, k;
    const _cx_size _cap = self->bucket_count;
    _cx_value* _slot = self->table;
    uint8_t* _hashx = self->_hashx;
    _cx_memb(_value_drop)(&_slot[i]);
    for (;;) { /* delete without leaving tombstone */
        _cx_size _mask = (_cx_size) -(++j != _cap);
        j &= _mask;
        if (! _hashx[j])
            break;
        _cx_rawkey _raw = i_keyto(_i_keyref(_slot + j));
        k = c_PASTE(fastrange_,MAP_SIZE_T)(i_hash(&_raw, sizeof _raw), _cap);
        if ((j < i) ^ (k <= i) ^ (k > j)) /* is k outside (i, j]? */
            _slot[i] = _slot[j], _hashx[i] = _hashx[j], i = j;
    }
    _hashx[i] = 0;
    --self->size;
}

#endif // _i_implement
#undef _i_isset
#undef _i_keyref
#undef _i_MAP_ONLY
#undef _i_SET_ONLY
#define CMAP_H_INCLUDED
// skipping file: template.h 
// ended inlining stc/cmap.h 
#undef i_key
#undef i_val
#undef i_tag
// ... done with map definition

// Private struct
struct null_comm {
  struct gkyl_comm_priv priv_comm; // base communicator
  struct gkyl_rect_decomp *decomp; // pre-computed decomposition

  bool use_gpu; // flag to use if this communicator is on GPUs
  bool sync_corners; // should we sync corners?
  
  struct gkyl_range grange; // range to "hash" ghost layout

  cmap_l2sgr l2sgr; // map from long -> skin_ghost_ranges
  cmap_l2sgr l2sgr_wc; // map from long -> skin_ghost_ranges with corners
  
  gkyl_mem_buff pbuff; // CUDA buffer for periodic BCs
};

// ended inlining gkyl_null_comm_priv.h 

#include <string.h>
#include <math.h>

// skipping file: gkyl_range.h 

// Create ghost and skin sub-ranges given a parent range
static void
skin_ghost_ranges_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
#define G_MAX(a,b) (a)>(b)?(a):(b)
  
  int ndim = parent->ndim;
  long max_vol = 0;
    
  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);

    max_vol = G_MAX(max_vol, sgr->lower_skin[d].volume);
    max_vol = G_MAX(max_vol, sgr->lower_ghost[d].volume);
    
    gkyl_skin_ghost_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);

    max_vol = G_MAX(max_vol, sgr->upper_skin[d].volume);
    max_vol = G_MAX(max_vol, sgr->upper_ghost[d].volume);
  }

  sgr->max_vol = max_vol;
#undef G_MAX
}

// Create ghost and skin sub-ranges given a parent range: includes
// corners
static void
skin_ghost_ranges_with_corners_init(struct skin_ghost_ranges *sgr,
  const struct gkyl_range *parent, const int *ghost)
{
#define G_MAX(a,b) (a)>(b)?(a):(b)
  
  int ndim = parent->ndim;
  long max_vol = 0;
    
  for (int d=0; d<ndim; ++d) {
    gkyl_skin_ghost_with_corners_ranges(&sgr->lower_skin[d], &sgr->lower_ghost[d],
      d, GKYL_LOWER_EDGE, parent, ghost);

    max_vol = G_MAX(max_vol, sgr->lower_skin[d].volume);
    max_vol = G_MAX(max_vol, sgr->lower_ghost[d].volume);
    
    gkyl_skin_ghost_with_corners_ranges(&sgr->upper_skin[d], &sgr->upper_ghost[d],
      d, GKYL_UPPER_EDGE, parent, ghost);

    max_vol = G_MAX(max_vol, sgr->upper_skin[d].volume);
    max_vol = G_MAX(max_vol, sgr->upper_ghost[d].volume);
  }

  sgr->max_vol = max_vol;
#undef G_MAX
}

static void
comm_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_comm *comm = container_of(ref, struct gkyl_comm, ref_count);  
  struct null_comm *null_comm = container_of(comm, struct null_comm, priv_comm.pub_comm);

  cmap_l2sgr_drop(&null_comm->l2sgr);
  cmap_l2sgr_drop(&null_comm->l2sgr_wc);
  gkyl_rect_decomp_release(null_comm->decomp);
  gkyl_mem_buff_release(null_comm->pbuff);
  gkyl_free(null_comm);
}

static int
get_rank(struct gkyl_comm *comm, int *rank)
{
  *rank = 0;
  return 0;
}

static int
get_size(struct gkyl_comm *comm, int *sz)
{
  *sz = 1;
  return 0;
}

static int
allreduce(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp, void *out)
{
  struct null_comm *null_comm = container_of(comm, struct null_comm, priv_comm.pub_comm);
  if (null_comm->use_gpu)
    gkyl_cu_memcpy(out, inp, gkyl_elem_type_size[type]*nelem, GKYL_CU_MEMCPY_D2D);
  else
    memcpy(out, inp, gkyl_elem_type_size[type]*nelem);
  return 0;
}

static int
allreduce_host(struct gkyl_comm *comm, enum gkyl_elem_type type,
  enum gkyl_array_op op, int nelem, const void *inp, void *out)
{
  struct null_comm *null_comm = container_of(comm, struct null_comm, priv_comm.pub_comm);
  memcpy(out, inp, gkyl_elem_type_size[type]*nelem);
  return 0;
}

static int
array_allgather(struct gkyl_comm *comm,
  const struct gkyl_range *local, const struct gkyl_range *global,
  const struct gkyl_array *array_local, struct gkyl_array *array_global)
{
  gkyl_array_copy(array_global, array_local);
  return 0;
}

static int
array_bcast(struct gkyl_comm *comm, const struct gkyl_array *asend,
  struct gkyl_array *arecv, int root)
{
  gkyl_array_copy(arecv, asend);
  return 0;
}

static int
array_sync(struct gkyl_comm *comm,
  const struct gkyl_range *local, const struct gkyl_range *local_ext,
  struct gkyl_array *array)
{
  return 0;
}

// apply periodic BCs
static void
apply_periodic_bc(const struct skin_ghost_ranges *sgr, char *data,
  int dir, struct gkyl_array *f)
{
  gkyl_array_copy_to_buffer(data, f, &(sgr->lower_skin[dir]));
  gkyl_array_copy_from_buffer(f, data, &(sgr->upper_ghost[dir]));

  gkyl_array_copy_to_buffer(data, f, &(sgr->upper_skin[dir]));
  gkyl_array_copy_from_buffer(f, data, &(sgr->lower_ghost[dir]));
}

static int
array_per_no_corners_sync(struct gkyl_comm *comm, const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  int nper_dirs, const int *per_dirs, struct gkyl_array *array)
{
  struct null_comm *null_comm = container_of(comm, struct null_comm, priv_comm.pub_comm);

  int nghost[GKYL_MAX_DIM];
  for (int d=0; d<null_comm->decomp->ndim; ++d)
    nghost[d] = local_ext->upper[d]-local->upper[d];
  
  long lkey = gkyl_range_idx(&null_comm->grange, nghost);

  if (!cmap_l2sgr_contains(&null_comm->l2sgr, lkey)) {
    struct skin_ghost_ranges sgr;
    skin_ghost_ranges_init(&sgr, local_ext, nghost);
    cmap_l2sgr_insert(&null_comm->l2sgr, lkey, sgr);
  }

  const cmap_l2sgr_value *val = cmap_l2sgr_get(&null_comm->l2sgr, lkey);
  long max_vol_esnz = val->second.max_vol*array->esznc;
  
  if (max_vol_esnz > gkyl_mem_buff_size(null_comm->pbuff))
    gkyl_mem_buff_resize(null_comm->pbuff, max_vol_esnz);

  char *data = gkyl_mem_buff_data(null_comm->pbuff);

  for (int d=0; d<nper_dirs; ++d)
    apply_periodic_bc(&val->second, data, per_dirs[d], array);

  return 0;
}

static int
array_per_with_corners_sync(struct gkyl_comm *comm, const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  int nper_dirs, const int *per_dirs, struct gkyl_array *array)
{
  struct null_comm *null_comm = container_of(comm, struct null_comm, priv_comm.pub_comm);

  int nghost[GKYL_MAX_DIM];
  for (int d=0; d<null_comm->decomp->ndim; ++d)
    nghost[d] = local_ext->upper[d]-local->upper[d];
  
  long lkey = gkyl_range_idx(&null_comm->grange, nghost);

  if (!cmap_l2sgr_contains(&null_comm->l2sgr_wc, lkey)) {
    struct skin_ghost_ranges sgr;
    skin_ghost_ranges_with_corners_init(&sgr, local_ext, nghost);
    cmap_l2sgr_insert(&null_comm->l2sgr_wc, lkey, sgr);
  }

  const cmap_l2sgr_value *val = cmap_l2sgr_get(&null_comm->l2sgr_wc, lkey);
  long max_vol_esnz = val->second.max_vol*array->esznc;
  
  if (max_vol_esnz > gkyl_mem_buff_size(null_comm->pbuff))
    gkyl_mem_buff_resize(null_comm->pbuff, max_vol_esnz);

  char *data = gkyl_mem_buff_data(null_comm->pbuff);

  for (int d=0; d<nper_dirs; ++d)
    apply_periodic_bc(&val->second, data, per_dirs[d], array);

  return 0;
}

static int
array_per_sync(struct gkyl_comm *comm, const struct gkyl_range *local,
  const struct gkyl_range *local_ext,
  int nper_dirs, const int *per_dirs, struct gkyl_array *array)
{
  struct null_comm *null_comm = container_of(comm, struct null_comm, priv_comm.pub_comm);  
  array_per_no_corners_sync(comm, local, local_ext, nper_dirs, per_dirs, array);
  if (null_comm->sync_corners)
    array_per_with_corners_sync(comm, local, local_ext, nper_dirs, per_dirs, array);

  return 0;
}

static int
barrier(struct gkyl_comm *comm)
{
  return 0;
}

static int array_write(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid,
  const struct gkyl_range *range,
  const struct gkyl_msgpack_data *meta,
  const struct gkyl_array *arr, const char *fname)
{
  return gkyl_grid_sub_array_write(grid, range, meta, arr, fname);
}

static int
array_read(struct gkyl_comm *comm,
  const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  struct gkyl_array *arr, const char *fname)
{
  struct gkyl_rect_grid fgrid;
  int status = gkyl_grid_sub_array_read(&fgrid, range, arr, fname);
  if (status == 0) {
    if (!gkyl_rect_grid_cmp(grid, &fgrid))
      status = 1;
  }
  return status;
}

static struct gkyl_comm*
extend_comm(const struct gkyl_comm *comm, const struct gkyl_range *erange)
{
  struct null_comm *null_comm = container_of(comm, struct null_comm, priv_comm.pub_comm);
  // extend internal decomp object and create a new communicator
  struct gkyl_rect_decomp *ext_decomp = gkyl_rect_decomp_extended_new(erange, null_comm->decomp);
  struct gkyl_comm *ext_comm = gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .decomp = ext_decomp,
      .use_gpu = null_comm->use_gpu,
      .sync_corners = null_comm->sync_corners
    }
  );
  gkyl_rect_decomp_release(ext_decomp);
  
  return ext_comm;
}

static struct gkyl_comm*
split_comm(const struct gkyl_comm *comm, int color, struct gkyl_rect_decomp *new_decomp)
{
  struct null_comm *null_comm = container_of(comm, struct null_comm, priv_comm.pub_comm);  

  return gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .use_gpu = null_comm->use_gpu,
      .sync_corners = null_comm->sync_corners,
      .decomp = new_decomp
    }
  );
}

static struct gkyl_comm*
create_comm_from_ranks(const struct gkyl_comm *comm,
  int nranks, const int *ranks, struct gkyl_rect_decomp *new_decomp,
  bool *is_valid)
{
  if (nranks > 1) {
    *is_valid = false;
    return 0;
  }
  
  *is_valid = true;
  
  struct null_comm *null_comm = container_of(comm, struct null_comm, priv_comm.pub_comm);  
  return gkyl_null_comm_inew( &(struct gkyl_null_comm_inp) {
      .use_gpu = null_comm->use_gpu,
      .sync_corners = null_comm->sync_corners,
      .decomp = new_decomp
    }
  );  
}

struct gkyl_comm*
gkyl_null_comm_inew(const struct gkyl_null_comm_inp *inp)
{
  struct null_comm *comm = gkyl_malloc(sizeof *comm);
  strcpy(comm->priv_comm.pub_comm.id, "null_comm");

  comm->priv_comm.pub_comm.has_decomp = true;  
  if (0 == inp->decomp) {
    comm->priv_comm.pub_comm.has_decomp = false;
    
    // construct a dummy decomposition
    comm->decomp =
      gkyl_rect_decomp_new_from_cuts_and_cells(1, (int[]) { 1 }, (int[]) { 1 });
  }
  else {
    comm->decomp = gkyl_rect_decomp_acquire(inp->decomp);
  }

  // construct range to hash ghost layout
  int lower[GKYL_MAX_DIM] = { 0 };
  int upper[GKYL_MAX_DIM];
  for (int d=0; d<comm->decomp->ndim; ++d)
    upper[d] = GKYL_MAX_NGHOST;
  gkyl_range_init(&comm->grange, comm->decomp->ndim, lower, upper);

  comm->use_gpu = inp->use_gpu;
  comm->sync_corners = inp->sync_corners;
  comm->l2sgr = cmap_l2sgr_init();
  comm->l2sgr_wc = cmap_l2sgr_init();

  if (comm->use_gpu)
    comm->pbuff = gkyl_mem_buff_cu_new(1024); // will be reallocated
  else
    comm->pbuff = gkyl_mem_buff_new(1024); // will be reallocated

  comm->priv_comm.get_rank = get_rank;
  comm->priv_comm.get_size = get_size;
  comm->priv_comm.allreduce = allreduce;
  comm->priv_comm.allreduce_host = allreduce_host;
  comm->priv_comm.gkyl_array_allgather = array_allgather;
  comm->priv_comm.gkyl_array_allgather_host = array_allgather;
  comm->priv_comm.gkyl_array_bcast = array_bcast;
  comm->priv_comm.gkyl_array_bcast_host = array_bcast;
  comm->priv_comm.gkyl_array_sync = array_sync;
  comm->priv_comm.gkyl_array_per_sync = array_per_sync;
  comm->priv_comm.barrier = barrier;
  comm->priv_comm.gkyl_array_write = array_write;
  comm->priv_comm.gkyl_array_read = array_read;
  comm->priv_comm.extend_comm = extend_comm;
  comm->priv_comm.split_comm = split_comm;
  comm->priv_comm.create_comm_from_ranks = create_comm_from_ranks;

  comm->priv_comm.pub_comm.ref_count = gkyl_ref_count_init(comm_free);

  return &comm->priv_comm.pub_comm;
}
// ended inlining null_comm.c 
// start inlining range.c 
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// skipping file: gkyl_range.h 
// skipping file: gkyl_util.h 

// flags and corresponding bit-masks
enum range_flags { R_IS_SUB_RANGE };
static const uint32_t masks[] =
{ 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

// sub-range flags
#define SET_SUB_RANGE(flags) (flags) |= masks[R_IS_SUB_RANGE]
#define CLEAR_SUB_RANGE(flags) (flags) &= ~masks[R_IS_SUB_RANGE]
#define IS_SUB_RANGE(flags) (((flags) & masks[R_IS_SUB_RANGE]) != 0)

// Computes coefficients for mapping indices in row-major order
static void
calc_rowmajor_ac(struct gkyl_range* range, long ac[])
{
  int ndim = range->ndim;
  ac[ndim] = 1L;
  for (int i=ndim-1; i>=1; --i)
    ac[i] = ac[i+1]*gkyl_range_shape(range, i);
  long start = 0L;
  for (int i=0; i<ndim; ++i)
    start += ac[i+1]*range->lower[i];
  ac[0] = -start;
}

// Computes stuff needed for "skip iterator"
static long
calc_skip_iter(const struct gkyl_range *rng, int *remDir)
{
  int up[GKYL_MAX_DIM];
  for (int d=0; d<rng->ndim; ++d) {
    remDir[d] = 1;
    up[d] = rng->upper[d];
  }
  int d = 0;
  long vol = rng->volume;
  long loidx = gkyl_range_idx(rng, rng->lower);
  long del = gkyl_range_idx(rng,up)-loidx+1;
  while (del != vol) {
    up[d] = rng->lower[d];
    vol /= gkyl_range_shape(rng, d);
    del = gkyl_range_idx(rng,up)-loidx+1;
    remDir[d] = 0;
    d += 1;
  }
  return del;
}

// compute volume, safely (for malformed ranges
static long
calc_volume_safely(int ndim, const int *lower, const int *upper)
{
  int is_zero_vol = 0;
  long vol = 1L;
  for (int i=0; i<ndim; ++i) {
    vol *= upper[i]-lower[i]+1;
    is_zero_vol = GKYL_MAX2(is_zero_vol, upper[i]<lower[i] ? 1 : 0);
  }
  if (is_zero_vol) vol = 0;
  return vol;
}

void
gkyl_range_init(struct gkyl_range *rng, int ndim,
  const int *lower, const int *upper)
{
//  // MF 2023/07/07: commenting this out because it causes seg faults in g2.
//  *rng = (struct gkyl_range) { };
  
  int is_zero_vol = 0;
  rng->ndim = ndim;
  rng->volume = 1L;
  for (int i=0; i<ndim; ++i) {
    rng->ilo[i] = rng->lower[i] = lower[i];
    rng->upper[i] = upper[i];
    rng->volume *= upper[i]-lower[i]+1;
    // need to handle case when upper[i]<lower[i]
    is_zero_vol = GKYL_MAX2(is_zero_vol, upper[i]<lower[i] ? 1 : 0);
  }
  // reset volume if any lower[d] <= upper[d]
  if (is_zero_vol) rng->volume = 0;
  
  calc_rowmajor_ac(rng, rng->ac);
  gkyl_copy_long_arr(GKYL_MAX_DIM+1, rng->ac, rng->iac);

  int idxZero[GKYL_MAX_DIM];
  for (int i=0; i<ndim; ++i) idxZero[i] = 0;
  rng->linIdxZero = gkyl_range_idx(rng, idxZero);

  rng->nsplit = 1;
  rng->tid = 0;

  rng->flags = 0;

  // for CUDA ops
  rng->nthreads = GKYL_DEFAULT_NUM_THREADS;
  rng->nblocks = rng->volume/rng->nthreads + 1;
}

void
gkyl_range_init_from_shape(struct gkyl_range *rng, int ndim, const int *shape)
{
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  for (int i=0; i<ndim; ++i) {
    lo[i] = 0; // lower-left corner has index (0,0,...)
    up[i] = shape[i]-1;
  }
  gkyl_range_init(rng, ndim, lo, up);
}

void
gkyl_range_init_from_shape1(struct gkyl_range *rng, int ndim, const int *shape)
{
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  for (int i=0; i<ndim; ++i) {
    lo[i] = 1; // lower-left corner has index (1,1,...)
    up[i] = shape[i];
  }
  gkyl_range_init(rng, ndim, lo, up);
}

void
gkyl_range_ten_prod(struct gkyl_range *rng, const struct gkyl_range *a, const struct gkyl_range *b)
{
  int adim = a->ndim, bdim = b->ndim;
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  for (int d=0; d<adim; ++d) {
    lower[d] = a->lower[d];
    upper[d] = a->upper[d];
  }
  for (int d=0; d<bdim; ++d) {
    lower[adim+d] = b->lower[d];
    upper[adim+d] = b->upper[d];
  }
  gkyl_range_init(rng, adim+bdim, lower, upper);
}

void
gkyl_range_shift(struct gkyl_range *rng, const struct gkyl_range *inp,
  const int *delta)
{
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  for (int d=0; d<inp->ndim; ++d) {
    lower[d] = inp->lower[d] + delta[d];
    upper[d] = inp->upper[d] + delta[d];
  }
  gkyl_range_init(rng, inp->ndim, lower, upper);
}

void
gkyl_range_reset_lower(struct gkyl_range *rng, const struct gkyl_range *inp,
  const int *new_lower)
{
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  for (int d=0; d<inp->ndim; ++d) {
    lower[d] = new_lower[d];
    upper[d] = new_lower[d] + gkyl_range_shape(inp, d) - 1;
  }
  gkyl_range_init(rng, inp->ndim, lower, upper);  
}

int
gkyl_range_is_sub_range(const struct gkyl_range *rng)
{
  return IS_SUB_RANGE(rng->flags);
}

int
gkyl_range_contains_idx(const struct gkyl_range *rng, const int *idx)
{
  for (int i=0; i<rng->ndim; ++i) {
    if ( (idx[i] < rng->lower[i]) || (idx[i] > rng->upper[i]) )
      return 0;
  }
  return 1;
}

void
gkyl_sub_range_init(struct gkyl_range *rng,
  const struct gkyl_range *bigrng, const int *sublower, const int *subupper)
{
  rng->ndim = bigrng->ndim;
  rng->volume = 1L;
  for (int i=0; i<rng->ndim; ++i) {
    rng->lower[i] = sublower[i] >= bigrng->lower[i] ? sublower[i] : bigrng->lower[i];
    rng->upper[i] = subupper[i] <= bigrng->upper[i] ? subupper[i] : bigrng->upper[i];
    rng->ilo[i] = bigrng->ilo[i]; // so inv indexer works correctly
    rng->volume *= rng->upper[i]-rng->lower[i]+1;
  }
  for (int i=0; i<rng->ndim+1; ++i)
    rng->ac[i] = bigrng->ac[i];
  rng->linIdxZero = bigrng->linIdxZero;

  rng->nsplit = bigrng->nsplit;
  rng->tid = bigrng->tid;
  
  rng->flags = bigrng->flags;
  SET_SUB_RANGE(rng->flags);

  // we need to construct iac such that sub_range_inv_idx works
  // properly
  struct gkyl_range sub_range;
  gkyl_range_init(&sub_range, rng->ndim, rng->lower, rng->upper);
  gkyl_copy_long_arr(GKYL_MAX_DIM+1, sub_range.ac, rng->iac);

  // for CUDA ops
  rng->nthreads = GKYL_DEFAULT_NUM_THREADS;
  rng->nblocks = rng->volume/rng->nthreads + 1;
}

struct gkyl_range
gkyl_range_split(struct gkyl_range *rng, int nsplit, int tid)
{
  struct gkyl_range r = *rng;
  r.nsplit = nsplit;
  r.tid = tid;
  return r;
}

// Computes split and returns number of elements handled locally and
// the initial index into the range. Number of elements is returned
// and start index set in 'lower'
static long
range_calc_split(const struct gkyl_range *rng, int *lower)
{
  const int nsplit = rng->nsplit, tid = rng->tid;
  long quot = rng->volume/nsplit, rem = rng->volume % nsplit;
  
  long len = gkyl_range_split_len(rng);
  long start = tid < rem ? tid*(quot+1) : rem*(quot+1) + (tid-rem)*quot;

  if (IS_SUB_RANGE(rng->flags)) {
    // as 'start' in sub-range we need to use an additional
    // indirection to compute the 'lower' bounds
    struct gkyl_range subrange;
    gkyl_range_init(&subrange, rng->ndim, rng->lower, rng->upper);
    gkyl_range_inv_idx(&subrange, start, lower);
  }
  else {
    gkyl_range_inv_idx(rng, start, lower);
  }

  return len;
}

long
gkyl_range_split_len(const struct gkyl_range *rng)
{
  const long quot = rng->volume/rng->nsplit, rem = rng->volume % rng->nsplit;
  return rng->tid < rem ? quot+1 : quot;
}

void
gkyl_range_deflate(struct gkyl_range* srng,
  const struct gkyl_range* rng, const int *remDir, const int *locDir)
{
  srng->linIdxZero = rng->linIdxZero;
  srng->ndim = 0;
  srng->volume = 1;  
  for (int i=0, j=0; i<rng->ndim; ++i) {
    if (!remDir[i]) {
      srng->lower[j] = rng->lower[i];
      srng->upper[j] = rng->upper[i];
      srng->ilo[j] = rng->ilo[i];
      srng->ac[j+1] = rng->ac[i+1];      
      srng->ndim += 1;
      srng->volume *= gkyl_range_shape(rng, i);
      j += 1;
    }
  }
  long adel = 0; // need to adjust ac[0]
  for (int i=0; i<rng->ndim; ++i)
    if (remDir[i])
      adel += locDir[i]*rng->ac[i+1];
  srng->ac[0] = rng->ac[0] + adel;

  srng->nsplit = rng->nsplit;
  srng->tid = rng->tid;

  srng->flags = rng->flags;
  SET_SUB_RANGE(srng->flags);

  // for CUDA ops
  srng->nthreads = GKYL_DEFAULT_NUM_THREADS;
  srng->nblocks = srng->volume/srng->nthreads + 1;
}

void
gkyl_range_shorten_from_above(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int len)
{
  int ndim = range->ndim;
  int lo[GKYL_MAX_DIM] = {0}, up[GKYL_MAX_DIM] = {0};
  
  for (int i=0; i<ndim; ++i) {
    lo[i] = range->lower[i];
    up[i] = range->upper[i];
  }
  up[dir] = lo[dir]+len-1;
  gkyl_sub_range_init(rng, range, lo, up);
}

void
gkyl_range_shorten_from_below(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int len)
{
  int ndim = range->ndim;
  int lo[GKYL_MAX_DIM] = {0}, up[GKYL_MAX_DIM] = {0};
  
  for (int i=0; i<ndim; ++i) {
    lo[i] = range->lower[i];
    up[i] = range->upper[i];
  }
  lo[dir] = up[dir]-len+1;
  gkyl_sub_range_init(rng, range, lo, up);
}

void
gkyl_range_extend(struct gkyl_range *erng,
  const struct gkyl_range* range, const int *elo, const int *eup)
{
  int ndim = range->ndim;
  int lo[GKYL_MAX_DIM] = {0}, up[GKYL_MAX_DIM] = {0};

  for (int i=0; i<ndim; ++i) {
    lo[i] = range->lower[i]-elo[i];
    up[i] = range->upper[i]+eup[i];
  }
  gkyl_range_init(erng, ndim, lo, up);
}

void
gkyl_range_perp_extend(struct gkyl_range *erng, int dir,
  const struct gkyl_range* rng, const int *elo, const int *eup)
{
  int ndim = rng->ndim;
  int elo_p[GKYL_MAX_DIM] = {0}, eup_p[GKYL_MAX_DIM] = {0};
  for (int i=0; i<ndim; ++i) {
    elo_p[i] = elo[i];
    eup_p[i] = eup[i];
  }
  elo_p[dir] = 0; eup_p[dir] = 0;
  gkyl_range_extend(erng, rng, elo_p, eup_p);
}

void
gkyl_range_lower_skin(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int nskin)
{
  int ndim = range->ndim;
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  
  for (int i=0; i<ndim; ++i) {
    lo[i] = range->lower[i];
    up[i] = range->upper[i];
  }
  up[dir] = range->lower[dir]+nskin-1;
  gkyl_sub_range_init(rng, range, lo, up);
}

void
gkyl_range_upper_skin(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int nskin)
{
  int ndim = range->ndim;
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  
  for (int i=0; i<ndim; ++i) {
    lo[i] = range->lower[i];
    up[i] = range->upper[i];
  }
  lo[dir] = range->upper[dir]-nskin+1;
  gkyl_sub_range_init(rng, range, lo, up);
}

// Increment an int vector by fact*del[d] in each direction d.
static inline void
incr_int_array(int ndim, int fact, const int * GKYL_RESTRICT del,
  const int * GKYL_RESTRICT inp, int *GKYL_RESTRICT out)
{
  for (int i=0; i<ndim; ++i)
    out[i] = inp[i] + fact*del[i];
}

/**
 * Create ghost and skin sub-ranges given parent (extended
 * range). This code is somewhat convoluted as the skin and ghost
 * ranges need to be sub-ranges of the extended range on the grid and
 * not include corners. I am not sure how to handle corners on
 * physical boundaries. Also, perhaps this code could be simplified.
 */
void
gkyl_skin_ghost_ranges(struct gkyl_range *skin, struct gkyl_range *ghost,
  int dir, enum gkyl_edge_loc edge, const struct gkyl_range *parent, const int *nghost)
{
  int ndim = parent->ndim;
  int lo[GKYL_MAX_DIM] = {0}, up[GKYL_MAX_DIM] = {0};

  if (edge == GKYL_LOWER_EDGE) {

    incr_int_array(ndim, 1, nghost, parent->lower, lo);
    incr_int_array(ndim, -1, nghost, parent->upper, up);
    
    up[dir] = lo[dir]+nghost[dir]-1;
    gkyl_sub_range_init(skin, parent, lo, up);

    incr_int_array(ndim, 1, nghost, parent->lower, lo);
    incr_int_array(ndim, -1, nghost, parent->upper, up);
    
    lo[dir] = lo[dir]-nghost[dir];
    up[dir] = lo[dir]+nghost[dir]-1;
    gkyl_sub_range_init(ghost, parent, lo, up);
  }
  else {

    incr_int_array(ndim, 1, nghost, parent->lower, lo);
    incr_int_array(ndim, -1, nghost, parent->upper, up);
    
    lo[dir] = up[dir]-nghost[dir]+1;
    gkyl_sub_range_init(skin, parent, lo, up);

    incr_int_array(ndim, 1, nghost, parent->lower, lo);
    incr_int_array(ndim, -1, nghost, parent->upper, up);
    
    up[dir] = up[dir]+nghost[dir]+1;
    lo[dir] = up[dir]-nghost[dir];
    gkyl_sub_range_init(ghost, parent, lo, up);
  }
}

void
gkyl_skin_ghost_with_corners_ranges(struct gkyl_range *skin, struct gkyl_range *ghost,
  int dir, enum gkyl_edge_loc edge, const struct gkyl_range *parent, const int *nghost)
{
  int ndim = parent->ndim;
  int lo[GKYL_MAX_DIM] = {0}, up[GKYL_MAX_DIM] = {0};

  for (int i=0; i<ndim; ++i) {
    lo[i] = parent->lower[i];
    up[i] = parent->upper[i];
  }

  if (edge == GKYL_LOWER_EDGE) {

    lo[dir] = parent->lower[dir]+nghost[dir];
    up[dir] = lo[dir]+nghost[dir]-1;
    gkyl_sub_range_init(skin, parent, lo, up);    

    lo[dir] = parent->lower[dir];
    up[dir] = lo[dir]+nghost[dir]-1;
    gkyl_sub_range_init(ghost, parent, lo, up);

  }
  else {

    up[dir] = parent->upper[dir]-nghost[dir];
    lo[dir] = up[dir]-nghost[dir]+1;
    gkyl_sub_range_init(skin, parent, lo, up);

    up[dir] = parent->upper[dir];
    lo[dir] = up[dir]-nghost[dir]+1;
    gkyl_sub_range_init(ghost, parent, lo, up);
  }
}

int
gkyl_range_intersect(struct gkyl_range* irng,
  const struct gkyl_range *r1, const struct gkyl_range *r2)
{
  int ndim = r1->ndim;
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d) {
    lo[d] = r1->lower[d] > r2->lower[d] ? r1->lower[d] : r2->lower[d];
    up[d] = r1->upper[d] < r2->upper[d] ? r1->upper[d] : r2->upper[d];
  }
  gkyl_range_init(irng, ndim, lo, up);
  return irng->volume > 0 ? 1 : 0;
}

int
gkyl_sub_range_intersect(struct gkyl_range* irng,
  const struct gkyl_range *r1, const struct gkyl_range *r2)
{
  int ndim = r1->ndim;
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d) {
    lo[d] = r1->lower[d] > r2->lower[d] ? r1->lower[d] : r2->lower[d];
    up[d] = r1->upper[d] < r2->upper[d] ? r1->upper[d] : r2->upper[d];
  }
  
  long vol = irng->volume = calc_volume_safely(ndim, lo, up);
  if (vol > 0)
    gkyl_sub_range_init(irng, r1, lo, up);
  else
    gkyl_range_init(irng, ndim, lo, up);
  return irng->volume > 0 ? 1 : 0;
}

bool
gkyl_range_is_on_lower_edge(int dir, const struct gkyl_range *range,
  const struct gkyl_range *parent)
{
  if (range->lower[dir] == parent->lower[dir])
    return true;
  return false;
  
}

bool
gkyl_range_is_on_upper_edge(int dir, const struct gkyl_range *range,
  const struct gkyl_range *parent)
{
  if (range->upper[dir] == parent->upper[dir])
    return true;
  return false;  
}

struct gkyl_range_dir_edge
gkyl_range_edge_match(const struct gkyl_range *base,
  const struct gkyl_range *targ)
{
  struct gkyl_range_dir_edge no_dir_ed = {
    .dir = 0,
    .eloc = GKYL_NO_EDGE
  };

  if (base->ndim != targ->ndim)
    return no_dir_ed; // different dimensions do not count

  struct gkyl_range irng;
  if (gkyl_range_intersect(&irng, base, targ))
    return no_dir_ed; // overlapping ranges do not count

  for (int d=0; d<base->ndim; ++d) {

    do {
      int elo[GKYL_MAX_DIM] = { 0 }, eup[GKYL_MAX_DIM] = { 0 };

      // check lower-edge overlap
      elo[d] = 1;
      struct gkyl_range erng;
      gkyl_range_extend(&erng, base, elo, eup);
      if (gkyl_range_intersect(&irng, &erng, targ))
        return (struct gkyl_range_dir_edge) { .dir = d, .eloc = GKYL_LOWER_EDGE };
    } while (0);

    do {
      int elo[GKYL_MAX_DIM] = { 0 }, eup[GKYL_MAX_DIM] = { 0 };    
      // check upper-edge overlap
      eup[d] = 1;
      struct gkyl_range erng;
      gkyl_range_extend(&erng, base, elo, eup);
      if (gkyl_range_intersect(&irng, &erng, targ))
        return (struct gkyl_range_dir_edge) { .dir = d, .eloc = GKYL_UPPER_EDGE };
    } while (0);
  }

  return no_dir_ed;
}

void
gkyl_range_iter_init(struct gkyl_range_iter *iter,
  const struct gkyl_range* range)
{
  iter->is_first = 1;
  iter->ndim = range->ndim;
  iter->bumps_left = range->volume > 0? range_calc_split(range, iter->idx) : 0;
  
  for (int i=0; i<range->ndim; ++i) {
    iter->lower[i] = range->lower[i];
    iter->upper[i] = range->upper[i];
  }
}

void
gkyl_range_iter_no_split_init(struct gkyl_range_iter *iter,
  const struct gkyl_range* range)
{
  iter->is_first = 1;
  iter->ndim = range->ndim;
  iter->bumps_left = range->volume > 0? range_calc_split(range, iter->idx) : 0;
  
  for (int i=0; i<range->ndim; ++i) {
    iter->idx[i] = iter->lower[i] = range->lower[i];
    iter->upper[i] = range->upper[i];
  }  
}

int
gkyl_range_iter_next(struct gkyl_range_iter *iter)
{
  if (iter->bumps_left-- < 1) return 0;
  
  if (iter->is_first) {
    iter->is_first = 0;
    return 1;
  }
  for (int dir=iter->ndim-1; dir>=0; --dir) {
    iter->idx[dir] += 1;
    if (iter->idx[dir] > iter->upper[dir])
      iter->idx[dir] = iter->lower[dir];
    else
      return 1;
  }
  return 0;
}

void
gkyl_range_skip_iter_init(struct gkyl_range_skip_iter *iter,
  const struct gkyl_range* range)
{
  int remDir[GKYL_MAX_DIM];
  iter->delta = calc_skip_iter(range, remDir);
  gkyl_range_deflate(&iter->range, range, remDir, range->lower);
}

void
gkyl_print_range(const struct gkyl_range* range, const char *nm, FILE *fp)
{
  fprintf(fp, "%s = { ndim = %d, ", nm, range->ndim);

  fprintf(fp, " lower = { ");
  for (int d=0; d<range->ndim; ++d)
    fprintf(fp, "%d%c ", range->lower[d], d==range->ndim-1 ? ' ' : ',');
  fprintf(fp, "}, ");

  fprintf(fp, "upper = { ");
  for (int d=0; d<range->ndim; ++d)
    fprintf(fp, "%d%c ", range->upper[d] , d==range->ndim-1 ? ' ' : ',');
  fprintf(fp, "}, ");

  fprintf(fp, " volume = %ld, ", range->volume );
  fprintf(fp, " is_sub_range = %d", gkyl_range_is_sub_range(range) );
  
  fprintf(fp, " }\n");
  fflush(fp);
}

bool
gkyl_range_compare(const struct gkyl_range* r1, const struct gkyl_range* r2)
{
  if (r1->ndim != r2->ndim)
    return false;
  for (int i=0; i<r1->ndim; ++i) {
    if (r1->lower[i] != r2->lower[i])
      return false;
    if (r1->upper[i] != r2->upper[i])
      return false;    
  }
  return true;
}
// ended inlining range.c 
// start inlining rect_decomp.c 
// skipping file: gkyl_alloc.h 
// skipping file: gkyl_array.h 
// skipping file: gkyl_array_ops.h 
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
// skipping file: gkyl_rect_decomp.h 
// skipping file: gkyl_util.h 

#include <string.h>

// Private struct to manage the neighbor data struct
struct rect_decomp_neigh_cont {
  struct gkyl_rect_decomp_neigh neigh;
  int* l_neigh;
  int* l_dir;
  int* l_edge;
};

static void
rect_decomp_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_rect_decomp *decomp = container_of(ref, struct gkyl_rect_decomp, ref_count);
  gkyl_free(decomp->ranges);
  gkyl_free(decomp);
}    

struct gkyl_rect_decomp*
gkyl_rect_decomp_new_from_cuts(int ndim, const int cuts[], const struct gkyl_range *range)
{
  struct gkyl_rect_decomp *decomp = gkyl_malloc(sizeof(*decomp));

  int ndecomp = 1;
  decomp->ndim = ndim;  

  for (int d=0; d<ndim; ++d) ndecomp *= cuts[d];  
  decomp->ndecomp = ndecomp;
  decomp->ranges = gkyl_malloc(sizeof(struct gkyl_range[ndecomp]));

  memcpy(&decomp->parent_range, range, sizeof(struct gkyl_range));

  div_t qr[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d)
    qr[d] = div(gkyl_range_shape(range, d), cuts[d]);

  int *sidx[GKYL_MAX_DIM], *eidx[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d) {
    sidx[d] = gkyl_malloc(sizeof(int[cuts[d]]));
    eidx[d] = gkyl_malloc(sizeof(int[cuts[d]]));

    int *shape = gkyl_malloc(sizeof(int[cuts[d]]));

    // compute shape in direction 'd'
    for (int i=0; i<cuts[d]; ++i)
      shape[i] = i<qr[d].rem ? qr[d].quot+1 : qr[d].quot;

    sidx[d][0] = range->lower[d];
    eidx[d][0] = sidx[d][0]+shape[0]-1;
    for (int i=1; i<cuts[d]; ++i) {
      sidx[d][i] = eidx[d][i-1]+1;
      eidx[d][i] = sidx[d][i]+shape[i]-1;
    }

    gkyl_free(shape);
  }

  struct gkyl_range rcuts;
  gkyl_range_init_from_shape(&rcuts, ndim, cuts);
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, &rcuts);

  int dnum = 0;
  // loop over cuts range, constructing each of the sub-ranges in the
  // decomposition
  while( gkyl_range_iter_next(&citer) ) {
    int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

    for (int d=0; d<ndim; ++d) {
      lower[d] = sidx[d][citer.idx[d]];
      upper[d] = eidx[d][citer.idx[d]];
    }

    gkyl_range_init(&decomp->ranges[dnum++], range->ndim, lower, upper);
  }

  for (int d=0; d<ndim; ++d) {
    gkyl_free(sidx[d]);
    gkyl_free(eidx[d]);
  }

  decomp->ref_count = gkyl_ref_count_init(rect_decomp_free);
  
  return decomp;
}

struct gkyl_rect_decomp*
gkyl_rect_decomp_new_from_cuts_and_cells(int ndim, const int cuts[], const int cells[])
{
  struct gkyl_range range;
  gkyl_create_global_range(ndim, cells, &range);
  return gkyl_rect_decomp_new_from_cuts(ndim, cuts, &range);
}

// ext_range = a X b 
static void
init_extend_range(struct gkyl_range *ext_range,
  const struct gkyl_range *a, const struct gkyl_range *b)
{
  int adim = a->ndim, bdim = b->ndim;
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  for (int d=0; d<adim; ++d) {
    lower[d] = a->lower[d];
    upper[d] = a->upper[d];
  }
  for (int d=0; d<bdim; ++d) {
    lower[adim+d] = b->lower[d];
    upper[adim+d] = b->upper[d];
  }

  gkyl_range_init(ext_range, adim+bdim, lower, upper);
}

struct gkyl_rect_decomp*
gkyl_rect_decomp_extended_new(const struct gkyl_range *arange,
  const struct gkyl_rect_decomp *decomp)
{
  struct gkyl_rect_decomp *extd = gkyl_malloc(sizeof(*extd));

  int ndecomp =  extd->ndecomp = decomp->ndecomp;  
  int ndim = extd->ndim = arange->ndim + decomp->ndim;
  extd->ranges = gkyl_malloc(sizeof(struct gkyl_range[ndecomp]));

  gkyl_range_ten_prod(&extd->parent_range, &decomp->parent_range, arange);
  for (int n=0; n<ndecomp; ++n)
    gkyl_range_ten_prod(&extd->ranges[n], &decomp->ranges[n], arange);

  extd->ref_count = gkyl_ref_count_init(rect_decomp_free);
  
  return extd;
}

struct gkyl_rect_decomp*
gkyl_rect_decomp_acquire(const struct gkyl_rect_decomp *decomp)
{
  gkyl_ref_count_inc(&decomp->ref_count);
  return (struct gkyl_rect_decomp*) decomp;
}

bool
gkyl_rect_decomp_check_covering(const struct gkyl_rect_decomp *decomp)
{
  struct gkyl_array *arr = gkyl_array_new(GKYL_DOUBLE, 1, decomp->parent_range.volume);
  gkyl_array_clear(arr, 0.0);

  // following loops over each sub-range and increments the region it
  // indexes in 'arr'. Each index should be visited exactly once.
  for (int i=0; i<decomp->ndecomp; ++i) {
    // construct a sub-range so indexing into global array works fine
    struct gkyl_range lrange;
    gkyl_sub_range_intersect(&lrange, &decomp->parent_range, &decomp->ranges[i]);
    
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &lrange);

    while (gkyl_range_iter_next(&iter)) {
      double *d = gkyl_array_fetch(arr, gkyl_range_idx(&decomp->parent_range, iter.idx));
      d[0] += 1.0;
    }
  }

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &decomp->parent_range);
  while (gkyl_range_iter_next(&iter)) {
    const double *d = gkyl_array_cfetch(arr, gkyl_range_idx(&decomp->parent_range, iter.idx));
    if (d[0] != 1.0)
      return false;
  }

  gkyl_array_release(arr);
  
  return true;
}

// compute neighbors accounting for corner neighbors
static struct gkyl_rect_decomp_neigh*
calc_neigh_with_corners(const struct gkyl_rect_decomp *decomp, int nidx)
{
  struct rect_decomp_neigh_cont *cont = gkyl_malloc(sizeof(*cont));
  cont->l_neigh = 0;
  cont->l_dir = 0;
  cont->l_edge = 0;
  
  int elo[GKYL_MAX_DIM], eup[GKYL_MAX_DIM];
  for (int i=0; i<decomp->ndim; ++i)
    elo[i] = eup[i] = 1;
  
  struct gkyl_range erng;
  gkyl_range_extend(&erng, &decomp->ranges[nidx], elo, eup);

  for (int i=0; i<decomp->ndecomp; ++i)
    if (i != nidx) {
      struct gkyl_range irng;
      int is_inter = gkyl_range_intersect(&irng, &erng,
        &decomp->ranges[i]);
      if (is_inter) {
        gkyl_fvec_push(cont->l_neigh, i);
        
        struct gkyl_range_dir_edge dir_ed =
          gkyl_range_edge_match(&decomp->ranges[nidx], &decomp->ranges[i]);

        gkyl_fvec_push(cont->l_dir, dir_ed.dir);
        gkyl_fvec_push(cont->l_edge, dir_ed.eloc);
      }
    }

  cont->neigh.num_neigh = gkyl_fvec_size(cont->l_neigh);
  cont->neigh.neigh = &cont->l_neigh[0];
  cont->neigh.dir = &cont->l_dir[0];
  cont->neigh.edge = &cont->l_edge[0];
  
  return &cont->neigh;
}

// compute neighbors leaving out corner neighbors: only face neighbors
// are included
static struct gkyl_rect_decomp_neigh*
calc_neigh_no_corners(const struct gkyl_rect_decomp *decomp, int nidx)
{
  struct rect_decomp_neigh_cont *cont = gkyl_malloc(sizeof(*cont));
  cont->l_neigh = 0;
  cont->l_dir = 0;
  cont->l_edge = 0;  
  
  struct gkyl_range erng;

  for (int n=0; n<decomp->ndim; ++n) {
    
    int elo[GKYL_MAX_DIM] = { 0 }, eup[GKYL_MAX_DIM] = { 0 };
    elo[n] = eup[n] = 1; // only extend in 1 direction
    gkyl_range_extend(&erng, &decomp->ranges[nidx], elo, eup);

    for (int i=0; i<decomp->ndecomp; ++i)
      if (i != nidx) {
        struct gkyl_range irng;
        int is_inter = gkyl_range_intersect(&irng, &erng,
          &decomp->ranges[i]);
        if (is_inter) {
          gkyl_fvec_push(cont->l_neigh, i);

          struct gkyl_range_dir_edge dir_ed =
            gkyl_range_edge_match(&decomp->ranges[nidx], &decomp->ranges[i]);
          
          gkyl_fvec_push(cont->l_dir, dir_ed.dir);
          gkyl_fvec_push(cont->l_edge, dir_ed.eloc);
        }
      }
  }
  
  cont->neigh.num_neigh = gkyl_fvec_size(cont->l_neigh);
  cont->neigh.neigh = &cont->l_neigh[0];
  cont->neigh.dir = &cont->l_dir[0];
  cont->neigh.edge = &cont->l_edge[0];
  
  return &cont->neigh;
}

struct gkyl_rect_decomp_neigh*
gkyl_rect_decomp_calc_neigh(const struct gkyl_rect_decomp *decomp,
  bool inc_corners, int nidx)
{
  if (inc_corners)
    return calc_neigh_with_corners(decomp, nidx);
  return calc_neigh_no_corners(decomp, nidx);
}

struct gkyl_rect_decomp_neigh*
gkyl_rect_decomp_calc_periodic_neigh(const struct gkyl_rect_decomp *decomp,
  int dir, bool inc_corners, int nidx)
{
  struct rect_decomp_neigh_cont *cont = gkyl_malloc(sizeof(*cont));
  cont->l_neigh = 0;
  cont->l_dir = 0;
  cont->l_edge = 0;  

  const struct gkyl_range *curr = &decomp->ranges[nidx];

  int elo[GKYL_MAX_DIM] = { 0 }, eup[GKYL_MAX_DIM] = { 0 };
  if (inc_corners)
    for (int i=0; i<decomp->ndim; ++i)
      elo[i] = eup[i] = 1;
  else
    elo[dir] = eup[dir] = 1;
  
  if (gkyl_range_is_on_lower_edge(dir, curr, &decomp->parent_range)) {
    int delta[GKYL_MAX_DIM] = { 0 };
    delta[dir] = gkyl_range_shape(&decomp->parent_range, dir);
    
    struct gkyl_range curr_shift;
    gkyl_range_shift(&curr_shift, curr, delta);
      
    struct gkyl_range shift_erng;
    gkyl_range_extend(&shift_erng, &curr_shift, elo, eup);

    for (int i=0; i<decomp->ndecomp; ++i)
      if (gkyl_range_is_on_upper_edge(dir, &decomp->ranges[i], &decomp->parent_range)) {
        struct gkyl_range irng;
        int is_inter = gkyl_range_intersect(&irng, &shift_erng,
          &decomp->ranges[i]);
        if (is_inter) {
          gkyl_fvec_push(cont->l_neigh, i);
          gkyl_fvec_push(cont->l_dir, dir);
          // this is not exactly correct: corner neighbors should not
          // be on any edge
          gkyl_fvec_push(cont->l_edge, GKYL_LOWER_EDGE);
        }
      }
  }
  else if (gkyl_range_is_on_upper_edge(dir, curr, &decomp->parent_range)) {
    int delta[GKYL_MAX_DIM] = { 0 };
    delta[dir] = -gkyl_range_shape(&decomp->parent_range, dir);
    
    struct gkyl_range curr_shift;
    gkyl_range_shift(&curr_shift, curr, delta);
      
    struct gkyl_range shift_erng;
    gkyl_range_extend(&shift_erng, &curr_shift, elo, eup);

    for (int i=0; i<decomp->ndecomp; ++i)
      if (gkyl_range_is_on_lower_edge(dir, &decomp->ranges[i], &decomp->parent_range)) {
        struct gkyl_range irng;
        int is_inter = gkyl_range_intersect(&irng, &shift_erng,
          &decomp->ranges[i]);
        if (is_inter) {
          gkyl_fvec_push(cont->l_neigh, i);
          gkyl_fvec_push(cont->l_dir, dir);
          // this is not exactly correct: corner neighbors should not
          // be on any edge          
          gkyl_fvec_push(cont->l_edge, GKYL_UPPER_EDGE);
        }
      }
  }

  cont->neigh.num_neigh = gkyl_fvec_size(cont->l_neigh);  
  cont->neigh.neigh = &cont->l_neigh[0];
  cont->neigh.dir = &cont->l_dir[0];
  cont->neigh.edge = &cont->l_edge[0];
  
  return &cont->neigh;
}

void
gkyl_rect_decomp_neigh_release(struct gkyl_rect_decomp_neigh *ng)
{
  struct rect_decomp_neigh_cont *cont = container_of(ng,
    struct rect_decomp_neigh_cont, neigh);
  gkyl_fvec_free(cont->l_neigh);
  gkyl_fvec_free(cont->l_dir);
  gkyl_fvec_free(cont->l_edge);
  gkyl_free(cont);
}

long
gkyl_rect_decomp_calc_offset(const struct gkyl_rect_decomp *decomp, int nidx)
{
  long offset = 0;
  for (int i=0; i<nidx; ++i)
    offset += decomp->ranges[i].volume;
  return offset;
}

void      
gkyl_rect_decomp_release(struct gkyl_rect_decomp *decomp)
{
  gkyl_ref_count_dec(&decomp->ref_count);
}

// Utility functions

void
gkyl_create_global_range(int ndim, const int *cells, struct gkyl_range *range)
{
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  for (int i=0; i<ndim; ++i) {
    // this needs to be consistent with gkyl_create_grid_ranges below
    lower[i] = 1;
    upper[i] = cells[i];
  }
  gkyl_range_init(range, ndim, lower, upper);
}

void
gkyl_create_grid_ranges(const struct gkyl_rect_grid *grid,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  
  for (int i=0; i<grid->ndim; ++i) {
    lower_ext[i] = 1-nghost[i];
    upper_ext[i] = grid->cells[i]+nghost[i];

    // this needs to be consistent with gkyl_create_global_range above
    lower[i] = 1;
    upper[i] = grid->cells[i];
  }
  gkyl_range_init(ext_range, grid->ndim, lower_ext, upper_ext);
  gkyl_sub_range_init(range, ext_range, lower, upper);
}

void
gkyl_create_ranges(const struct gkyl_range *inrange,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  
  for (int i=0; i<inrange->ndim; ++i) {
    lower_ext[i] = inrange->lower[i]-nghost[i];
    upper_ext[i] = inrange->upper[i]+nghost[i];

    lower[i] = inrange->lower[i];
    upper[i] = inrange->upper[i];
  }
  gkyl_range_init(ext_range, inrange->ndim, lower_ext, upper_ext);
  gkyl_sub_range_init(range, ext_range, lower, upper);  
}

void
gkyl_create_vertex_ranges(const struct gkyl_range *inrange,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range)
{
  int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  for (int i=0; i<inrange->ndim; ++i) {
    lower_ext[i] = inrange->lower[i] - nghost[i];
    upper_ext[i] = inrange->upper[i] + 1 + nghost[i];

    lower[i] = inrange->lower[i];
    upper[i] = inrange->upper[i] + 1;
  }
  gkyl_range_init(ext_range, inrange->ndim, lower_ext, upper_ext);
  gkyl_sub_range_init(range, ext_range, lower, upper);  
}

void
gkyl_rect_decomp_get_cuts(struct gkyl_rect_decomp* decomp, int* cuts)
{
  int ndim = decomp->ndim;

  for (int d=0; d<ndim; d++) {
    int other_dim_lo[GKYL_MAX_DIM] = {0}, other_dim_up[GKYL_MAX_DIM] = {0};
    for (int i=0; i<ndim; i++) {
      if (i != d) {
        other_dim_lo[i] = decomp->ranges[0].lower[i];
        other_dim_up[i] = decomp->ranges[0].upper[i];
      }
    }

    int cuts_curr = 0, range_idx = 0;

    bool not_reached_upper = true;
    while (not_reached_upper) {
      struct gkyl_range range_curr = decomp->ranges[range_idx];
      bool same_other_lims = true;
      for (int i=0; i<ndim; i++) {
        if (i != d) {
          same_other_lims = same_other_lims &&
            ((range_curr.lower[i] == other_dim_lo[i]) && (range_curr.upper[i] == other_dim_up[i]));
        }
      }

      if (same_other_lims) {
        cuts_curr++;
        if (range_curr.upper[d] == decomp->parent_range.upper[d])
          not_reached_upper = false;
      }

      range_idx++;
    }
    cuts[d] = cuts_curr;
  }
}
// ended inlining rect_decomp.c 
// start inlining rect_grid.c 
#include <assert.h>
#include <stdint.h>
#include <math.h>

// skipping file: gkyl_alloc.h 
// skipping file: gkyl_rect_grid.h 
// start inlining gkyl_rect_grid_priv.h 

// skipping file: gkyl_rect_grid.h 

/* Find upper and lower boundaries of given cell
 * @params grid: grid struct which contains the cells
 * @params cell_in: given cell
 * @params dim_trans: The dimensions to check
 * @params known_index: Any already known indices
 * @params lower_boundaries: lower sides of given cell (output)
 * @params upper_boundaries: upper sides of given cell (output)
 */
GKYL_CU_DH
void in_dir(const struct gkyl_rect_grid *grid, int *cell_in, const int *dim_trans,
  const int *known_index, double lower_boundaries[], double upper_boundaries[])
{
  int ndim, check_index;
  const double *dx, *lower;
  ndim = grid -> ndim;
  dx = grid -> dx;
  lower = grid -> lower;
  for (int d=0; d<ndim; d++) {
    check_index = known_index[d] < 0 ? cell_in[dim_trans[d]] : known_index[d];
    lower_boundaries[d] = lower[d]+(check_index-1)*dx[d];
    upper_boundaries[d] = lower[d]+(check_index)*dx[d];
  }
}

/* Checks if given point is in given cell
 * @params grid: grid struct which contains the cells
 * @params point: coordinates of given point
 * @params cell_in: given cell
 * @params dim_trans: The dimensions to check
 * @params known_index: Any already known indices
 * @returns bool: true if point is in given cell
 */
GKYL_CU_DH
bool is_in_cell(const struct gkyl_rect_grid *grid, const double *point,
  int *cell_in, const int *dim_trans, const int *known_index)
{
  int ndim;
  ndim = grid -> ndim;
  double lower_boundaries[ndim], upper_boundaries[ndim];
  for (int d=0; d<ndim; d++) {
    lower_boundaries[d]=0;
    upper_boundaries[d]=0;
  }
  in_dir(grid, cell_in, dim_trans, known_index, lower_boundaries, upper_boundaries);
  bool in_cell = true;
  double eps = 1.0e-14;
  for (int d=0; d<ndim; d++) {
    if (lower_boundaries[d]-eps>point[d] || upper_boundaries[d]+eps<point[d]) {
      in_cell = false;
      break;
    }
  }
  return in_cell;
}

// ended inlining gkyl_rect_grid_priv.h 
// skipping file: gkyl_util.h 
void
gkyl_rect_grid_init(struct gkyl_rect_grid *grid, int ndim,
  const double *lower, const double *upper, const int *cells)
{
//  // MF 2023/07/07: commenting this out because it causes seg faults in g2.
//  *grid = (struct gkyl_rect_grid) { };
  
  grid->ndim = ndim;  
  grid->cellVolume = 1.0;
  for (int i=0; i<ndim; ++i) {
    grid->lower[i] = lower[i];
    grid->upper[i] = upper[i];
    grid->cells[i] = cells[i];
    grid->dx[i] = (upper[i]-lower[i])/cells[i];
    grid->cellVolume *= grid->dx[i];
  }
}

bool
gkyl_rect_grid_cmp(const struct gkyl_rect_grid *grid1, struct gkyl_rect_grid *grid2)
{
  if (grid1->ndim != grid2->ndim)
    return false;

  for (int i=0; i<grid1->ndim; ++i) {
    if (grid1->cells[i] != grid2->cells[i])
      return false;
    if (!gkyl_compare_double(grid1->lower[i], grid2->lower[i], 1e-14))
      return false;
    if (!gkyl_compare_double(grid1->upper[i], grid2->upper[i], 1e-14))
      return false;    
  }
  return true;
}

GKYL_CU_DH
void
gkyl_rect_grid_find_cell(const struct gkyl_rect_grid *grid, const double *point,
  bool pick_lower, const int *known_index, int *cell_index){

  int nDim = grid->ndim;
  int search_num = 0;
  int search_dim[GKYL_MAX_DIM];
  int dim_trans[GKYL_MAX_DIM];
  double low, high;
  
  for (int d=0; d<nDim; d++) {
    if (known_index[d] < 0) {
      search_dim[search_num] = d;
      dim_trans[d] = search_num;
      search_num = search_num + 1;
    } else {
      dim_trans[d] = -1;
      cell_index[d] = known_index[d];
      low = grid->lower[d]+(known_index[d]-1)*grid->dx[d];
      high = grid->lower[d]+(known_index[d])*grid->dx[d];
      assert(low<point[d]&&high>point[d]);
    }
  }

  int start_index[GKYL_MAX_DIM], end_index[GKYL_MAX_DIM], mid_index[GKYL_MAX_DIM], new_index[GKYL_MAX_DIM];
  const int *cells = grid -> cells;
  for (int d=0; d<search_num; d++) {
    start_index[d] = 1;
    end_index[d] = cells[search_dim[d]];
    mid_index[d] = 0;
    new_index[d] = 0;
  }

  int plusminus[2] = {-1,1}, low_high_index[2*GKYL_MAX_DIM];
  bool all_less_eq = true;
  double lower_dir[nDim], upper_dir[nDim];
  /* Below we use a binary search. That is, if the i-th coordinate of the point in
   * question is below(above) the ith-coordinate of the current mid-point, we search
   * the lower(upper) half along that direction in the next iteration.
   */
  while (all_less_eq) {
    for (int d=0; d<search_num; d++)
      mid_index[d] = (start_index[d]+end_index[d])/2;  // Integer division intentional

    if (is_in_cell(grid, point, mid_index, dim_trans, known_index)) {
      // Check if neighboring cells also contain this point.
      for (int k=0; k<search_num; k++){
	new_index[k] = mid_index[k];
	low_high_index[search_dim[k]] = mid_index[k];
	low_high_index[nDim+search_dim[k]] = mid_index[k];
      }
      for (int i=0; i<search_num; i++){	
	for (int j=0; j<2; j++){
	  new_index[i] = GKYL_MAX2(GKYL_MIN2(mid_index[i] + plusminus[j],cells[i]),0);
	  if (is_in_cell(grid, point, new_index, dim_trans, known_index))
	    low_high_index[j*nDim+search_dim[i]] = new_index[i];	    	  
	}
	new_index[i] = mid_index[i];
      }
      if (pick_lower) {
	for (int d=0; d<search_num; d++) {
	  cell_index[search_dim[d]] = low_high_index[search_dim[d]];
	}
      } else {
	for (int d=0; d<search_num; d++) {
	  cell_index[search_dim[d]] = low_high_index[nDim+search_dim[d]];
	}
      }
      break;
    } else {
      in_dir(grid, mid_index, dim_trans, known_index, lower_dir, upper_dir);
      for (int d=0; d<search_num; d++) {
	if (point[search_dim[d]] < lower_dir[search_dim[d]]) {
	  end_index[d] = mid_index[d]-1;
	} else if(point[search_dim[d]] > upper_dir[search_dim[d]]) {
	  start_index[d] = mid_index[d]+1;
	}
      }
    }

    for (int d=0; d<search_num; d++) {
      if (start_index[d]>end_index[d]) {
	all_less_eq = false;
	break;
      }
    }
  }
}

void
gkyl_rect_grid_write(const struct gkyl_rect_grid *grid, FILE *fp)
{
  // dimension and shape are written as 64 bit integers
  uint64_t ndim = grid->ndim;
  uint64_t cells[GKYL_MAX_DIM];
  for (int d=0; d<grid->ndim; ++d)
    cells[d] = grid->cells[d];

  fwrite(&ndim, sizeof(uint64_t), 1, fp);
  fwrite(cells, sizeof(uint64_t), grid->ndim, fp);
  fwrite(grid->lower, sizeof(double), grid->ndim, fp);
  fwrite(grid->upper, sizeof(double), grid->ndim, fp);
}

bool
gkyl_rect_grid_read(struct gkyl_rect_grid *grid, FILE *fp)
{
  uint64_t ndim = grid->ndim;
  uint64_t cells64[GKYL_MAX_DIM];

  if (1 != fread(&ndim, sizeof(uint64_t), 1, fp))
    return false;
  if (ndim != fread(cells64, sizeof(uint64_t), ndim, fp))
    return false;

  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  if (ndim != fread(lower, sizeof(double), ndim, fp))
    return false;
  if (ndim != fread(upper, sizeof(double), ndim, fp))
    return false;

  // copy into regular int array
  int cells[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d) cells[d] = cells64[d];

  gkyl_rect_grid_init(grid, ndim, lower, upper, cells);

  return true;
}
// ended inlining rect_grid.c 
// start inlining thread_pool.c 
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
// skipping file: gkyl_alloc.h 

// skipping file: thpool.h 

struct jp_thread_pool {
  struct gkyl_job_pool jp; // base job-pool object
  threadpool thpool; // thread-pool object
};

static void
thread_pool_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_job_pool *base = container_of(ref, struct gkyl_job_pool, ref_count);
  struct jp_thread_pool *th = container_of(base, struct jp_thread_pool, jp);
  thpool_destroy(th->thpool);
  gkyl_free(th);
}

static bool
thread_pool_add_work(const struct gkyl_job_pool *jp, jp_work_func func, void *ctx)
{
  struct jp_thread_pool *th = container_of(jp, struct jp_thread_pool, jp);
  int status = thpool_add_work(th->thpool, func, ctx);
  return status == 0 ? true : false;
}

static void
thread_pool_wait(const struct gkyl_job_pool *jp)
{
  struct jp_thread_pool *th = container_of(jp, struct jp_thread_pool, jp);
  thpool_wait(th->thpool);
}

struct gkyl_job_pool*
gkyl_thread_pool_new(int nthreads)
{
  struct jp_thread_pool *th = gkyl_malloc(sizeof(struct jp_thread_pool));
  // initialize the actual pool object
  th->thpool = thpool_init(nthreads);  

  th->jp.pool_size = nthreads;
  th->jp.add_work = thread_pool_add_work;
  th->jp.wait = thread_pool_wait;

  // set reference counter
  th->jp.ref_count = gkyl_ref_count_init(thread_pool_free);
    
  return &th->jp;
}
// ended inlining thread_pool.c 
// start inlining util.c 
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

// skipping file: gkyl_util.h 
// skipping file: gkyl_alloc.h 

// skipping file: mpack.h 
#include <assert.h>

int
gkyl_search_str_int_pair_by_str(const struct gkyl_str_int_pair pairs[], const char *str, int def)
{
  for (int i=0; pairs[i].str != 0; ++i) {
    if (strcmp(pairs[i].str, str) == 0)
      return pairs[i].val;
  }
  return def;  
}

const char *
gkyl_search_str_int_pair_by_int(const struct gkyl_str_int_pair pairs[], int val, const char *def)
{
  for (int i=0; pairs[i].str != 0; ++i) {
    if (pairs[i].val == val)
      return pairs[i].str;
  }
  return def;  
}

int
gkyl_tm_trigger_check_and_bump(struct gkyl_tm_trigger *tmt, double tcurr)
{
  int status = 0;
  if (tcurr >= tmt->tcurr) {
    status = 1;
    tmt->curr += 1;
    tmt->tcurr += tmt->dt;
  }
  return status;
}

void
gkyl_exit(const char* msg)
{
  fprintf(stderr, "Error: %s\n", msg);
  exit(EXIT_FAILURE);
}

int
gkyl_compare_float(float a, float b, float eps)
{
  //if (isnanf(a) || isnanf(b)) return 0;
  
  float absa = fabs(a), absb = fabs(b), diff = fabs(a-b);

  if (a == b) return 1;
  if (a == 0 || b == 0 || (absa+absb < FLT_MIN)) return diff < eps;
  if (absa < eps) return diff < eps;
  if (absb < eps) return diff < eps;
  return diff/fminf(absa+absb, FLT_MAX) < eps;
}

int
gkyl_compare_double(double a, double b, double eps)
{
  if (isnan(a) || isnan(b)) return 0;
  
  double absa = fabs(a), absb = fabs(b), diff = fabs(a-b);
  if (a == b) return 1;
  if (a == 0 || b == 0 || (absa+absb < DBL_MIN)) return diff < eps;
  if (absa < eps) return diff < eps;
  if (absb < eps) return diff < eps;
  return diff/fmin(absa+absb, DBL_MAX) < eps;
}

struct timespec
gkyl_wall_clock(void)
{
  struct timespec tm = { 0 };
#ifdef GKYL_HAVE_CUDA
  cudaDeviceSynchronize();
#endif
  // we were using CLOCK_REALTIME here
  clock_gettime(CLOCK_MONOTONIC, &tm);
  return tm;
}

struct timespec
gkyl_time_diff(struct timespec start, struct timespec end)
{
  struct timespec tm;
  if ((end.tv_nsec-start.tv_nsec)<0) {
    tm.tv_sec = end.tv_sec-start.tv_sec-1;
    tm.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
  }
  else {
    tm.tv_sec = end.tv_sec-start.tv_sec;
    tm.tv_nsec = end.tv_nsec-start.tv_nsec;
  }
  return tm;  
}

double
gkyl_time_diff_now_sec(struct timespec tm)
{
  return gkyl_time_sec(gkyl_time_diff(tm, gkyl_wall_clock()));
}
   
double
gkyl_time_sec(struct timespec tm)
{
  return tm.tv_sec + 1e-9*tm.tv_nsec;
}

double
gkyl_time_now(void)
{
  return gkyl_time_sec( gkyl_wall_clock() );
}

pcg32_random_t
gkyl_pcg32_init(bool nd_seed)
{
  pcg32_random_t rng;
  int rounds = 5;

  if (nd_seed)
    // seed with external entropy -- the time and some program addresses
    // (which will actually be somewhat random on most modern systems).
    pcg32_srandom_r(&rng, time(NULL) ^ (intptr_t)&printf, 
      (intptr_t)&rounds);
  else
    // seed with a fixed constant
    pcg32_srandom_r(&rng, 42u, 54u);

  return rng;
}

uint32_t
gkyl_pcg32_rand_uint32(pcg32_random_t* rng)
{
  return pcg32_random_r(rng);
}

double
gkyl_pcg32_rand_double(pcg32_random_t* rng)
{
  return ldexp(pcg32_random_r(rng), -32);
}

static void
pcg64_srandom_r(pcg64_random_t* rng, uint64_t seed1, uint64_t seed2,
  uint64_t seq1,  uint64_t seq2)
{
  uint64_t mask = ~0ull >> 1;
  // stream for each generators *must* be distinct
  if ((seq1 & mask) == (seq2 & mask)) 
    seq2 = ~seq2;
  pcg32_srandom_r(rng->gen,   seed1, seq1);
  pcg32_srandom_r(rng->gen+1, seed2, seq2);
}

static int _dummy_global = 0; // just to provide address for use in seed

pcg64_random_t
gkyl_pcg64_init(bool nd_seed)
{
  pcg64_random_t rng;
  int rounds = 5;

  if (nd_seed)
    pcg64_srandom_r(&rng,
      time(NULL) ^ (intptr_t)&printf, ~time(NULL) ^ (intptr_t)&pcg32_random_r,
      (intptr_t)&rounds, (intptr_t)&_dummy_global);
  else
    pcg64_srandom_r(&rng, 42u, 42u, 54u, 54u);

  return rng;
}

uint64_t
gkyl_pcg64_rand_uint64(pcg64_random_t* rng)
{
  return ((uint64_t)(pcg32_random_r(rng->gen)) << 32) | pcg32_random_r(rng->gen+1);
}

double
gkyl_pcg64_rand_double(pcg64_random_t* rng)
{
  return ldexp(gkyl_pcg64_rand_uint64(rng), -64);
}

bool
gkyl_check_file_exists(const char *fname)
{
  return access(fname, F_OK) == 0;
}

int64_t
gkyl_file_size(const char *fname)
{
  struct stat st;
  stat(fname, &st);
  return st.st_size;
}

char*
gkyl_load_file(const char *fname, int64_t *sz)
{
  int64_t msz = gkyl_file_size(fname);
  char *buff = gkyl_malloc(msz);
  FILE *fp = fopen(fname, "r");
  int n = fread(buff, msz, 1, fp);
  *sz = msz;
  fclose(fp);
  return buff;
}

bool
gkyl_msgpack_map_elem_has_key(int nvals, const struct gkyl_msgpack_map_elem *elist,
  const char *key)
{
  bool has_key = false;
  for (int i=0; i<nvals; ++i) {
    if (strcmp(key, elist[i].key) == 0) {
      has_key = true;
      break;
    }
  }
  return has_key;
}

struct gkyl_msgpack_map_elem *
gkyl_msgpack_map_elem_clone(int nvals, const struct gkyl_msgpack_map_elem *elist_in)
{
  struct gkyl_msgpack_map_elem *elist_out = gkyl_malloc(nvals*sizeof(struct gkyl_msgpack_map_elem));

  for (int i=0; i<nvals; ++i) {
    // Copy key name.
    size_t key_len = strlen(elist_in[i].key) + 1;
    elist_out[i].key = gkyl_malloc(key_len);
    strcpy(elist_out[i].key, elist_in[i].key);

    // Copy type.
    elist_out[i].elem_type = elist_in[i].elem_type;

    // Copy value.
    switch (elist_in[i].elem_type) {
      case GKYL_MP_BOOL:
        elist_out[i].bval = elist_in[i].bval;
        break;

      case GKYL_MP_UNSIGNED_INT:
        elist_out[i].uval = elist_in[i].uval;
        break;

      case GKYL_MP_INT:
        elist_out[i].ival = elist_in[i].ival;
        break;

      case GKYL_MP_FLOAT:
        elist_out[i].fval = elist_in[i].fval;
        break;

      case GKYL_MP_DOUBLE:
        elist_out[i].dval = elist_in[i].dval;
        break;

      case GKYL_MP_STRING:
        key_len = strlen(elist_in[i].cval) + 1;
        elist_out[i].cval = gkyl_malloc(key_len);
        strcpy(elist_out[i].cval, elist_in[i].cval);
        break;

      default:
        assert(false); // NYI.
        break;
    }
  }

  return elist_out;
}

struct gkyl_msgpack_map_elem *
gkyl_msgpack_map_elem_union(int numlist_union, int *nvals_union,
  const struct gkyl_msgpack_map_elem **elist_union, int *elist_out_len)
{
  int nvals_tot = 0; // Total number of elements.
  for (int j=0; j<numlist_union; ++j)
    nvals_tot += nvals_union[j];

  assert(nvals_tot > 0);

  elist_out_len[0] = nvals_tot;
  struct gkyl_msgpack_map_elem *elist_out = gkyl_malloc(elist_out_len[0]*sizeof(struct gkyl_msgpack_map_elem));

  int eidx = 0;
  for (int j=0; j<numlist_union; ++j) {

    int nvals_curr = nvals_union[j];
    const struct gkyl_msgpack_map_elem *elist_curr = elist_union[j];

    for (int i=0; i<nvals_curr; ++i) {
      // Copy key name.
      size_t slen = strlen(elist_curr[i].key) + 1;
      elist_out[eidx].key = gkyl_malloc(slen);
      strcpy(elist_out[eidx].key, elist_curr[i].key);

      // Copy type.
      elist_out[eidx].elem_type = elist_curr[i].elem_type;

      // Copy value.
      switch (elist_curr[i].elem_type) {
        case GKYL_MP_BOOL:
          elist_out[eidx].bval = elist_curr[i].bval;
          break;

        case GKYL_MP_UNSIGNED_INT:
          elist_out[eidx].uval = elist_curr[i].uval;
          break;

        case GKYL_MP_INT:
          elist_out[eidx].ival = elist_curr[i].ival;
          break;

        case GKYL_MP_FLOAT:
          elist_out[eidx].fval = elist_curr[i].fval;
          break;

        case GKYL_MP_DOUBLE:
          elist_out[eidx].dval = elist_curr[i].dval;
          break;

        case GKYL_MP_STRING:
          slen = strlen(elist_curr[i].cval) + 1;
          elist_out[eidx].cval = gkyl_malloc(slen);
          strcpy(elist_out[eidx].cval, elist_curr[i].cval);
          break;

	default:
	  assert(false); // NYI.
          break;
      }

      eidx++;
    }
  }

  return elist_out;
}

void
gkyl_msgpack_map_elem_set_double(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key, double value)
{
  for (int i=0; i<nvals; ++i) {
    if (strcmp(key, elist[i].key) == 0) {
      assert(elist[i].elem_type == GKYL_MP_DOUBLE);
      elist[i].dval = value;
      break;
    }
  }
}

void
gkyl_msgpack_map_elem_set_uint(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key, unsigned int value)
{
  for (int i=0; i<nvals; ++i) {
    if (strcmp(key, elist[i].key) == 0) {
      assert(elist[i].elem_type == GKYL_MP_UNSIGNED_INT);
      elist[i].uval = value;
      break;
    }
  }
}

double
gkyl_msgpack_map_elem_get_double(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key)
{
  for (int i=0; i<nvals; ++i) {
    if (strcmp(key, elist[i].key) == 0) {
      assert(elist[i].elem_type == GKYL_MP_DOUBLE);
      return elist[i].dval;
    }
  }
  return 0;
}

unsigned int
gkyl_msgpack_map_elem_get_uint(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key)
{
  for (int i=0; i<nvals; ++i) {
    if (strcmp(key, elist[i].key) == 0) {
      assert(elist[i].elem_type == GKYL_MP_UNSIGNED_INT);
      return elist[i].uval;
    }
  }
  return 0;
}

char *
gkyl_msgpack_map_elem_get_string(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key)
{
  for (int i=0; i<nvals; ++i) {
    if (strcmp(key, elist[i].key) == 0) {
      assert(elist[i].elem_type == GKYL_MP_STRING);
      return elist[i].cval;
    }
  }
  return 0;
}

void
gkyl_msgpack_map_elem_release_string(int nvals, struct gkyl_msgpack_map_elem *elist,
  const char *key)
{
  for (int i=0; i<nvals; ++i) {
    if (strcmp(key, elist[i].key) == 0) {
      assert(elist[i].elem_type == GKYL_MP_STRING);
      MPACK_FREE(elist[i].cval);
      break;
    }
  }
}

void
gkyl_msgpack_map_elem_release(int nvals, struct gkyl_msgpack_map_elem *elist_in)
{
  for (int i=0; i<nvals; ++i) {
    gkyl_free(elist_in[i].key);

    if (elist_in[i].elem_type == GKYL_MP_STRING)
      gkyl_free(elist_in[i].cval);
  }

  gkyl_free(elist_in);
}

struct gkyl_msgpack_data *
gkyl_msgpack_create(int nvals, const struct gkyl_msgpack_map_elem *elist)
{
  struct gkyl_msgpack_data *mdata = gkyl_malloc(sizeof *mdata);
  mdata->meta_sz = 0;
  mdata->meta = 0;

  mpack_writer_t writer;
  mpack_writer_init_growable(&writer, &mdata->meta, &mdata->meta_sz);

  mpack_build_map(&writer);  
  for (int i=0; i<nvals; ++i) {
    mpack_write_cstr(&writer, elist[i].key);

    switch (elist[i].elem_type) {
      case GKYL_MP_BOOL:
        mpack_write_bool(&writer, elist[i].bval);
        break;

      case GKYL_MP_UNSIGNED_INT:
        mpack_write_u64(&writer, elist[i].uval);
        break;

      case GKYL_MP_INT:
        mpack_write_i64(&writer, elist[i].ival);
        break;

      case GKYL_MP_FLOAT:
        mpack_write_float(&writer, elist[i].fval);
        break;

      case GKYL_MP_DOUBLE:
        mpack_write_double(&writer, elist[i].dval);
        break;

      case GKYL_MP_STRING:
        mpack_write_cstr(&writer, elist[i].cval);
        break;

      default:
        assert(false); // NYI.
        break;
    }
  }

  mpack_complete_map(&writer);

  int status = mpack_writer_destroy(&writer);

  if (status != mpack_ok) {
    MPACK_FREE(mdata->meta); // we need to use free here as mpack does its own malloc
    gkyl_free(mdata);
    mdata = 0;
  }

  return mdata;
}

struct gkyl_msgpack_data *
gkyl_msgpack_create_union(int numlist_union, int *nvals_union, const struct gkyl_msgpack_map_elem **elist_union)
{
  struct gkyl_msgpack_data *mdata = gkyl_malloc(sizeof *mdata);
  mdata->meta_sz = 0;
  mdata->meta = 0;

  mpack_writer_t writer;
  mpack_writer_init_growable(&writer, &mdata->meta, &mdata->meta_sz);

  mpack_build_map(&writer);  

  for (int j=0; j<numlist_union; ++j) {

    int nvals = nvals_union[j];
    const struct gkyl_msgpack_map_elem *elist = elist_union[j];

    for (int i=0; i<nvals; ++i) {
      mpack_write_cstr(&writer, elist[i].key);
  
      switch (elist[i].elem_type) {
        case GKYL_MP_BOOL:
          mpack_write_bool(&writer, elist[i].bval);
          break;
  
        case GKYL_MP_UNSIGNED_INT:
          mpack_write_u64(&writer, elist[i].uval);
          break;
  
        case GKYL_MP_INT:
          mpack_write_i64(&writer, elist[i].ival);
          break;
  
        case GKYL_MP_FLOAT:
          mpack_write_float(&writer, elist[i].fval);
          break;
  
        case GKYL_MP_DOUBLE:
          mpack_write_double(&writer, elist[i].dval);
          break;
  
        case GKYL_MP_STRING:
          mpack_write_cstr(&writer, elist[i].cval);
          break;

        default:
          assert(false); // NYI.
          break;
      }
    }
  }

  mpack_complete_map(&writer);

  int status = mpack_writer_destroy(&writer);

  if (status != mpack_ok) {
    MPACK_FREE(mdata->meta); // we need to use free here as mpack does its own malloc
    gkyl_free(mdata);
    mdata = 0;
  }

  return mdata;
}

static void
msgpack_copy_value(mpack_reader_t* r, mpack_writer_t* w)
{
  // Copy the current entry at the top of the reader stack into the writer
  // stack. If it's a map (most cases in Gkeyll so far), loop over the
  // elements.
  mpack_tag_t tag = mpack_read_tag(r);

  switch (mpack_tag_type(&tag)) {
    case mpack_type_nil:
      mpack_write_nil(w);
      break;
    case mpack_type_bool:
      mpack_write_bool(w, mpack_tag_bool_value(&tag));
      break;
    case mpack_type_int:
      mpack_write_int(w, mpack_tag_int_value(&tag));
      break;
    case mpack_type_uint:
      mpack_write_uint(w, mpack_tag_uint_value(&tag));
      break;
    case mpack_type_float:
      mpack_write_float(w, mpack_tag_float_value(&tag));
      break;
    case mpack_type_double:
      mpack_write_double(w, mpack_tag_double_value(&tag));
      break;
    case mpack_type_str: {
      uint32_t len = mpack_tag_str_length(&tag);
      mpack_start_str(w, len);

      char buffer[64];
      uint32_t remaining = len;
      while (remaining > 0) {
        uint32_t chunk = remaining > sizeof(buffer)? sizeof(buffer) : remaining;
        mpack_read_bytes(r, buffer, chunk);
        mpack_write_bytes(w, buffer, chunk);
        remaining -= chunk;
      }

      mpack_done_str(r);
      mpack_finish_str(w);
      break;
    }
    case mpack_type_map: {
      uint32_t count = mpack_tag_map_count(&tag);
      mpack_start_map(w, count);
      for (uint32_t i=0; i<count; i++) {
        msgpack_copy_value(r, w); // Copy key.
        msgpack_copy_value(r, w); // Copy value.
      }
      mpack_done_map(r);
      mpack_finish_map(w);
      break;
    }
    case mpack_type_array: {
      uint32_t count = mpack_tag_array_count(&tag);
      mpack_start_array(w, count);
      for (uint32_t i=0; i<count; i++) {
        msgpack_copy_value(r, w);
      }
      mpack_done_array(r);
      mpack_finish_array(w);
      break;
    }
    default:
      assert(false); // NYI.
      break;
  }
}

struct gkyl_msgpack_data *
gkyl_msgpack_clone(struct gkyl_msgpack_data *mdata_in)
{
  struct gkyl_msgpack_data *mdata_out = gkyl_malloc(sizeof *mdata_out);
  mdata_out->meta_sz = 0;
  mdata_out->meta = 0;

  mpack_reader_t reader;
  mpack_reader_init_data(&reader, mdata_in->meta, mdata_in->meta_sz);
  mpack_writer_t writer;
  mpack_writer_init_growable(&writer, &mdata_out->meta, &mdata_out->meta_sz);

  msgpack_copy_value(&reader, &writer);

  // Check copy was successful.
  if (
      (mpack_reader_destroy(&reader) != mpack_ok || mpack_writer_destroy(&writer) != mpack_ok) ||
      (!(mdata_in->meta_sz == mdata_out->meta_sz && memcmp(mdata_in->meta, mdata_out->meta, mdata_in->meta_sz) == 0))
     ) {
    fprintf(stderr, "gkyl_msgpack_clone: error copying MessagePack.\n");
    MPACK_FREE(mdata_out->meta); // we need to use free here as mpack does its own malloc
    gkyl_free(mdata_out);
    mdata_out = 0;
  }

  return mdata_out;
}

void
gkyl_msgpack_to_map_elem_list(struct gkyl_msgpack_data* mpack_in, int nvals,
  struct gkyl_msgpack_map_elem *elist)
{
  mpack_tree_t tree;
  mpack_tree_init_data(&tree, mpack_in->meta, mpack_in->meta_sz);
  mpack_tree_parse(&tree);
  mpack_node_t root = mpack_tree_root(&tree);

  for (int i=0; i<nvals; ++i) {
    mpack_node_t node = mpack_node_map_cstr(root, elist[i].key);
  
    switch (elist[i].elem_type) {
      case GKYL_MP_BOOL:
        elist[i].bval = mpack_node_bool(node);
        break;
  
      case GKYL_MP_UNSIGNED_INT:
        elist[i].uval = mpack_node_uint(node);
        break;
  
      case GKYL_MP_INT:
        elist[i].ival = mpack_node_int(node);
        break;
  
      case GKYL_MP_FLOAT:
        elist[i].fval = mpack_node_float(node);
        break;
  
      case GKYL_MP_DOUBLE:
        elist[i].dval = mpack_node_double(node);
        break;
  
      case GKYL_MP_STRING:
        elist[i].cval = mpack_node_cstr_alloc(node, mpack_node_strlen(node)+1);
        break;

      default:
        assert(false); // NYI.
        break;
    }
  }

  mpack_tree_destroy(&tree);
}

void
gkyl_msgpack_data_release(struct gkyl_msgpack_data *mdata)
{
  if (!mdata) return;
  if (mdata->meta_sz > 0)
    MPACK_FREE(mdata->meta);
  gkyl_free(mdata);
}

// ended inlining util.c 
