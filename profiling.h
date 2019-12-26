#pragma once

#include <chrono>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/types.h>
#include <unistd.h>

struct ProfilingPoint {
    rusage usage;
    std::chrono::high_resolution_clock::time_point time;

    ProfilingPoint();
};

float cpu_utilization(const ProfilingPoint &start, const ProfilingPoint &end);
size_t elapsed_time_ms(const ProfilingPoint &start, const ProfilingPoint &end);
