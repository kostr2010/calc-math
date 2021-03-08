#ifndef __DSL_H_INCLUDED__
#define __DSL_H_INCLUDED__

// ====================
// DSL

// destructor
#define DEFAULT_DTOR(class) ~class() = default;

// construcor
#define DEFAULT_CTOR(class) class() = default;

#define NO_DEFAULT_CTOR(class) class() = delete;

// copy semantics
#define COPY_SEMANTICS(class)                                                                      \
    class(const class& r);                                                                         \
    class& operator=(const class& r);

#define DEFAULT_COPY_SEMANTICS(class)                                                              \
    class(const class& r) = default;                                                               \
    class& operator=(const class& r) = default;

#define NO_COPY_SEMANTICS(class)                                                                   \
    class(const class& r) = delete;                                                                \
    class& operator=(const class& r) = delete;

// move semantics
#define MOVE_SEMANTICS(class)                                                                      \
    class(class && r);                                                                             \
    class& operator=(class&& r);

#define DEFAULT_MOVE_SEMANTICS(class)                                                              \
    class(class && r) = default;                                                                   \
    class& operator=(class&& r) = default;

#define NO_MOVE_SEMANTICS(class)                                                                   \
    class(class && r) = delete;                                                                    \
    class& operator=(class&& r) = delete;

#endif