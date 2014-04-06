#include "util.h"

#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <iostream>

using namespace std;

bool equals(float a, float b, float epsilon) {
    return fabs(a - b) < epsilon;
}

vector<string> split(const string &str, const string &delims) {
    vector<string> result;
    char *save;
    char buf[str.length() + 1];

    strcpy(buf, str.c_str());

    for(char *tok = strtok_r(buf, delims.c_str(), &save);
	tok != nullptr;
	tok = strtok_r(nullptr, delims.c_str(), &save))
    {
	result.push_back(tok);
    }

    return result;
}

void Timer::EntryReport::echo() {
    cout << desc;
    if(n == 1) {
	cout << ": t=" << total << endl;
    } else {
	cout << ": n=" << n;
	cout << ", min=" << min;
	cout << ", median=" << median;
	cout << ", mean=" << mean;
	cout << ", max=" << max;
	cout << ", total=" << total;
	cout << endl;
    }
}

Timer::Entry::Entry(const string &desc) {
    _desc = desc;
}

void Timer::Entry::start() {
    errif(_start != 0.0, "entry not started");

    _start = getTime_msec();
}

void Timer::Entry::end(const std::string &desc) {
    errif(_start == 0.0, "entry not started");
    errif(_desc != desc, "timer end out of sequence");

    add( getTime_msec() - _start );
    _start = 0.0;
}

#ifdef UTIL_CL
void Timer::Entry::add(const cl_event &event) {
    cl_ulong time_start, time_end, time;
    clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START,
			    sizeof(time_start), &time_start, NULL);

    clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END,
			    sizeof(time_end), &time_end, NULL);
    // duration in nanoseconds
    time = time_end - time_start;

    // convert to milliseconds
    add( time / 1e6 );
}
#endif

Timer::EntryReport Timer::Entry::report() {
    EntryReport r;
    
    r.desc = _desc;
    
    sort( _times.begin(), _times.end() );

    r.n = _times.size();
    r.min = _times[0];
    r.max = _times[r.n - 1];
    r.median = _times[r.n / 2];

    r.total = 0.0;
    for(auto t: _times) {
	r.total += t;
    }
    r.mean = r.total / r.n;

    return r;
}

void Timer::Entry::add(double time) {
    _times.push_back(time);
}

Timer::~Timer() {
    for(auto entry: _entries) {
	delete entry.second;
    }
}

void Timer::push(const string &desc) {
    Entry *entry = findEntry(desc);
    entry->start();
    _stack.push( entry );
}

void Timer::pop(const string &desc) {
    Entry *entry = _stack.top();
    entry->end(desc);
    _stack.pop();
}

#ifdef UTIL_CL
void Timer::add(const string &desc, const cl_event &event) {
    findEntry(desc)->add(event);
}
#endif

void Timer::report() {
    cout << "--------------------" << endl;
    cout << "--- Timer Report ---" << endl;
    cout << "--------------------" << endl;

    vector<EntryReport> reports;

    for(auto entry: _entries) {
	reports.push_back( entry.second->report() );
    }

    sort( reports.begin(),
	  reports.end(),
	  [](const Timer::EntryReport &a, const Timer::EntryReport &b) {
	      return a.min < b.min;
	  });

    for(auto r: reports) {
	r.echo();
    }
}

Timer::Entry *Timer::findEntry(const std::string &desc) {
    Entry *entry = _entries[desc];
    if(entry == nullptr) {
	entry = new Entry(desc);
	_entries[desc] = entry;
    }
    return entry;
}

double Timer::getTime_msec() {
    timespec curr;
    
    errif(0 != clock_gettime(CLOCK_REALTIME, &curr), "clock_gettime");

    return (curr.tv_sec * 1000) + (curr.tv_nsec / 1e6);
}

// todo: *NIX SPECIFIC
string pathcat(const string &a, const string &b) {
    string result = a;
    if(a[a.length() - 1] != '/')
	result += "/";
    result += b;
    return result;
}

char *load_file(const char *path) {
    void *result;
    size_t n;

    load_file(path, &result, &n);

    return (char*)result;
}

void load_file(const char *path, void **ret_buf, size_t *ret_n) {
    struct stat s;

    int rc = stat(path, &s);
    errif(rc == -1, "failed stating file %s", path);

    size_t n = s.st_size;
    char *buf = (char *)malloc(n + 1);
    errif(!buf, "failed allocating %lu bytes for sequence", *ret_n);

    FILE *f = fopen(path, "r");
    errif(!f, "failed opening %s", path);

    dbf("reading %lu bytes from %s...", n, path);

    size_t nread = fread(buf, n, 1, f);
    errif(nread != 1, "failed reading entire file (%s)", path);
    buf[n] = 0;

    *ret_n = n;
    *ret_buf = buf;

    return;
}
