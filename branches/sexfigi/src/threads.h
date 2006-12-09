 /*
 				threads.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	A program that uses POSIX threads
*
*	Author:		E.BERTIN
*
*	Contents:	Definitions and shortcuts for POSIX threads.
*
*	Last modify:	03/07/2002
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include <pthread.h>
#include <signal.h>

/*---- Set defines according to machine's specificities and customizing -----*/
/*--------------------------- Technical constants ---------------------------*/
/*---------------------------- Synchro messages -----------------------------*/

#define	STATE_FREE		0
#define	STATE_READY		1
#define	STATE_BUSY		2

/*------------------------------- Other Macros ------------------------------*/

#define QPTHREAD_ATTR_INIT(pthread_attr) \
	{if (pthread_attr_init(pthread_attr)) \
		error(EXIT_FAILURE, \
		"*Error*: pthread_attr_init() failed for ", #pthread_attr );;}

#define QPTHREAD_ATTR_SETDETACHSTATE(pthread_attr, attr) \
	{if (pthread_attr_setdetachstate(pthread_attr, attr)) \
		error(EXIT_FAILURE, \
		"*Error*: pthread_attr_setdetachstate() failed for ", \
		#pthread_attr );;}

#define QPTHREAD_ATTR_DESTROY(pthread_attr) \
	{if (pthread_attr_destroy(pthread_attr)) \
		error(EXIT_FAILURE, \
		"*Error*: pthread_attr_destroy() failed for ",#pthread_attr);;}

#define QPTHREAD_CREATE(pthread, attr, func, arg) \
	{if (pthread_create(pthread, attr, func, arg)) \
		error(EXIT_FAILURE, \
		"*Error*: pthread_create() failed for ", #pthread );;}

#define QPTHREAD_CANCEL(pthread) \
	{if (pthread_cancel(pthread)) \
		warning( \
		"failed to cancel ", #pthread );;}

#define QPTHREAD_JOIN(pthread, ret) \
	{if (pthread_join(pthread, ret)) \
		error(EXIT_FAILURE, \
		"*Error*: pthread_join() failed for ", #pthread );;}

#define QPTHREAD_MUTEX_INIT(mutex, attr) \
	{if (pthread_mutex_init(mutex, attr)) \
		error(EXIT_FAILURE, \
		"*Error*: pthread_mutex_init() failed for ", #mutex );;}

#define QPTHREAD_MUTEX_LOCK(mutex) \
	{if (pthread_mutex_lock(mutex)) \
		error(EXIT_FAILURE, \
		"*Error*: pthread_mutex_lock() failed for ", #mutex );;}

#define QPTHREAD_MUTEX_UNLOCK(mutex) \
	{if (pthread_mutex_unlock(mutex)) \
		error(EXIT_FAILURE, \
		"*Error*: pthread_mutex_unlock() failed for ", #mutex );;}

#define QPTHREAD_MUTEX_DESTROY(mutex) \
	{if (pthread_mutex_destroy(mutex)) \
		error(EXIT_FAILURE, \
		"*Error*: pthread_mutex_destroy() failed for ", #mutex );;}

#define QPTHREAD_COND_INIT(cond, attr) \
	{if (pthread_cond_init(cond, attr)) \
		error(EXIT_FAILURE, \
		"*Error*: pthread_cond_init() failed for ", #cond );;}

#define QPTHREAD_COND_WAIT(cond, mutex) \
	{if (pthread_cond_wait(cond, mutex)) \
		error(EXIT_FAILURE, \
		"*Error*: pthread_cond_wait() failed for ", #cond );;}

#define QPTHREAD_COND_BROADCAST(cond) \
	{if (pthread_cond_broadcast(cond)) \
		error(EXIT_FAILURE, \
		"*Error*: pthread_cond_broadcast() failed for ", #cond );;}

#define QPTHREAD_COND_SIGNAL(cond) \
	{if (pthread_cond_signal(cond)) \
		error(EXIT_FAILURE, \
		"*Error*: pthread_cond_signal() failed for ", #cond );;}

#define QPTHREAD_COND_DESTROY(cond) \
	{if (pthread_cond_destroy(cond)) \
		error(EXIT_FAILURE, \
		"*Error*: pthread_cond_destroy() failed for ", #cond );;}

/*------------------------------- Structures --------------------------------*/
typedef struct _threads_gate_t
  {
  int			ngate;		/* Gate counter */
  int			nthreads;	/* Number of threads to manage */
  void			(*func)(void);	/* Function to execute at wakeup */
  pthread_mutex_t	mutex;		/* Main MutEx */
  pthread_mutex_t	block;		/* Safety Mutex (avoid "rebound") */
  pthread_cond_t	condvar;	/* Main condition variable */
  pthread_cond_t	last;		/* To wake the remaining thread up */
  } threads_gate_t;

/*----------------------------- Global variables ----------------------------*/
 int		nproc;	/* Number of child threads */

/*--------------------------------- Functions -------------------------------*/
threads_gate_t	*threads_gate_init(int nthreads, void (*func)(void));

void		threads_gate_end(threads_gate_t *gate),
		threads_gate_sync(threads_gate_t *gate);

