#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <massert.h>
#include <zlog.h>
#include <logger.h>

zlog_category_t *sensey_logger = NULL;

void logger_init(const char * filename, const char * category) {
	if (strlen(filename) && strlen(category)) {
		massert(!zlog_init(filename), "Can not initialize logger");
		sensey_logger = zlog_get_category(category);
	}
}

zlog_category_t * logger_instance(void) {
	return sensey_logger;
}

void logger_destroy() {
	zlog_fini();
}

