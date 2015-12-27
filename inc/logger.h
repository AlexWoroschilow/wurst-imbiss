#ifndef SENSEY_LOGGER
#define SENSEY_LOGGER

#include <zlog.h>

void logger_init(const char * filename, const char * category);
zlog_category_t * logger_instance(void);
void logger_destroy();

#define logger_info(...) \
	if(logger_instance()) {\
		zlog_info(logger_instance(), __VA_ARGS__);\
	}

#define logger_debug(...) \
	if(logger_instance()) {\
		zlog_debug(logger_instance(), __VA_ARGS__);\
	}

#define logger_warn(...) \
	if(logger_instance()) {\
		zlog_warn(logger_instance(), __VA_ARGS__);\
	}

#define logger_notice(...) \
	if(logger_instance()) {\
		zlog_notice(logger_instance(), __VA_ARGS__);\
	}

#define logger_error(...) \
	if(logger_instance()) {\
		zlog_error(logger_instance(), __VA_ARGS__);\
	}

#define logger_fatal(...) \
	if(logger_instance()) {\
		zlog_fatal(logger_instance(), __VA_ARGS__);\
	}

#endif
