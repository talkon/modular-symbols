#ifndef DEBUG_PRINT_H
#define DEBUG_PRINT_H

// color macros, used to color debugging messages
#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RESET "\x1B[0m"

// macro to probe how far a program was executed before breaking
#define DEBUG_PRINT printf("function " MAG"%s() " RESET "| line " CYN "%d\n" RESET, __func__, __LINE__);

#endif // DEBUG_PRINT_H