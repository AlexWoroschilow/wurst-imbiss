#ifndef IMBISS_ASSERT
#define IMBISS_ASSERT

#define massert(conditions, message) \
    if (!(conditions)) { \
        fprintf(stderr, "%s - Error in file %s at line %lu\n", \
                        message, __FILE__, (unsigned long) __LINE__); \
        exit(EXIT_FAILURE); \
    }

#endif
