#pragma once

#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <map>
#include <string>
#include <stack>
#include <vector>

#ifdef UTIL_CL
#include <CL/cl.h>
#endif

#define ep(x) fprintf(stderr, "%s\n", x)
#define epf(x...) fprintf(stderr, x);  fprintf(stderr, "\n")
#define pf(x...) printf(x); printf("\n")
#define db(x) dbf("%s", x)
#define dbf(x...) printf(x); printf("\n")
#define err(msg...) {fprintf(stderr, msg); fprintf(stderr, "\n"); exit(1);}
#define errif(expr, msg...) if(expr) { err(msg); }

#define panic() {                               \
        fprintf(stderr, "PANIC!\n");            \
        abort();                                \
    }

#ifdef UTIL_CL
#define clx(expr) {cl_int __err = expr; if(__err != CL_SUCCESS) {std::cerr << #expr << " = " << __err << std::endl; exit(1);}}
#endif

bool equals(float a, float b, float epsilon);

std::vector<std::string> split(const std::string &str,
                               const std::string &delims = " \t\n\r");


template<typename Iterator>
std::string join(Iterator begin, Iterator end, const std::string &glue = " ") {
    if(begin == end)
        return "";
    std::string result = *begin;
    ++begin;
    for(; begin != end; ++begin) {
        result += glue;
        result += *begin;
    }
    return result;
}

template<typename Container>
std::string join(Container cont, const std::string &glue = " ") {
    return join(cont.begin(), cont.end(), glue);
}

class Timer {
    class EntryReport {
    public:
        std::string desc;
        size_t n; 
        double min;
        double max;
        double median;
        double mean;
        double total;

        void echo();
    };

    class Entry {
    public:

        Entry(const std::string &desc);
        
        void start();
        void end(const std::string &desc);
#ifdef UTIL_CL
        void add(const cl_event &event);
#endif

        EntryReport report();

    private:
        void add(double time);

        double _start = 0.0;
        std::string _desc;
        std::vector<double> _times;
    };
public:
    ~Timer();

    void push(const std::string &desc);
    void pop(const std::string &desc);
#ifdef UTIL_CL
    void add(const std::string &desc, const cl_event &event);
#endif
    void report();

private:
    Entry *findEntry(const std::string &desc);
    static double getTime_msec();

    std::stack<Entry *> _stack;
    std::map<std::string, Entry *> _entries;
};

std::string pathcat(const std::string &a, const std::string &b);
char *load_file(const char *path);
void load_file(const char *path, void **ret_buf, size_t *ret_n);
