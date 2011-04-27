#ifndef __ERROR_hh__

#define __ERROR_hh__

#define _ERROR_(a) error(__FILE__,__FUNCTION__,__LINE__,a)
#define _WARNING_(a) warning(__FILE__,__FUNCTION__,__LINE__,a)
#define _ERROR_TRACKING_(a) do {			\
    if (ierr > 0){					\
      error_tracking(__FILE__,__FUNCTION__,__LINE__);	\
      return a;						\
    }							\
  } while(0)


extern int ierr;

void error(const char* file,const char* function ,const int line,const char* message);

void warning(const char* file,const char* function ,const int line,const char* message);

void error_tracking(const char* file,const char* function ,const int line);


#endif
