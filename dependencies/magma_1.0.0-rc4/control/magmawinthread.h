/**
 *
 * @file magmawinthread.h
 *
 *  This file handles the mapping from pthreads calls to windows
 *  threads MAGMA is a software package provided by Univ. of
 *  Tennessee, Univ. of California Berkeley and Univ. of Colorado
 *  Denver
 *
 * @version 2.3.1
 * @author Piotr Luszczek
 * @date 2010-11-15
 *
 * This file is originally from PLASMA project, where plasma has been
 * replaced by MAGMA.
 *
 **/
#ifndef MAGMAWINTHREAD_H
#define MAGMAWINTHREAD_H

#include <windows.h>

/*
typedef struct pthread_s {
  HANDLE Hth;
  unsigned IDth;
  void *(*Fth) (void *);
} pthread_t;
*/
typedef struct pthread_s {
  HANDLE hThread;
  unsigned int uThId;
} pthread_t;

typedef HANDLE pthread_mutex_t;
typedef int pthread_mutexattr_t;
typedef int pthread_attr_t;
typedef int pthread_condattr_t;

typedef struct pthread_cond_s {
  HANDLE hSem;
  HANDLE hEvt;
  CRITICAL_SECTION cs;
  int waitCount; /* waiting thread counter */
} pthread_cond_t;

typedef int pthread_attr_t;

#define PTHREAD_MUTEX_INITIALIZER ((pthread_mutex_t) -1)

#define PTHREAD_SCOPE_SYSTEM 1

#define MAGMA_DLLPORT
#define MAGMA_CDECL __cdecl

MAGMA_DLLPORT pthread_t MAGMA_CDECL pthread_self(void);
MAGMA_DLLPORT int MAGMA_CDECL pthread_mutex_init(pthread_mutex_t *mutex, const pthread_mutexattr_t * attr);
MAGMA_DLLPORT int MAGMA_CDECL pthread_mutex_destroy(pthread_mutex_t *mutex);
MAGMA_DLLPORT int MAGMA_CDECL pthread_mutex_lock(pthread_mutex_t *mutex);
MAGMA_DLLPORT int MAGMA_CDECL pthread_mutex_trylock(pthread_mutex_t *mutex);
MAGMA_DLLPORT int MAGMA_CDECL pthread_mutex_unlock(pthread_mutex_t *mutex);
MAGMA_DLLPORT int MAGMA_CDECL pthread_attr_init(pthread_attr_t *attr);
MAGMA_DLLPORT int MAGMA_CDECL pthread_attr_destroy(pthread_attr_t *attr);
MAGMA_DLLPORT int MAGMA_CDECL pthread_attr_setscope(pthread_attr_t *attr, int scope);
MAGMA_DLLPORT int MAGMA_CDECL pthread_create(pthread_t *tid, const pthread_attr_t *attr, void *(*start) (void *), void *arg);
MAGMA_DLLPORT int MAGMA_CDECL pthread_cond_init(pthread_cond_t *cond, const pthread_condattr_t *attr);
MAGMA_DLLPORT int MAGMA_CDECL pthread_cond_destroy(pthread_cond_t *cond);
MAGMA_DLLPORT int MAGMA_CDECL pthread_cond_wait(pthread_cond_t *cond, pthread_mutex_t *mutex);
MAGMA_DLLPORT int MAGMA_CDECL pthread_cond_broadcast(pthread_cond_t *cond);
MAGMA_DLLPORT int MAGMA_CDECL pthread_join(pthread_t thread, void **value_ptr);
MAGMA_DLLPORT int MAGMA_CDECL pthread_equal(pthread_t thread1, pthread_t thread2);

MAGMA_DLLPORT int MAGMA_CDECL pthread_setconcurrency (int);

MAGMA_DLLPORT unsigned int MAGMA_CDECL pthread_self_id(void);

#endif
