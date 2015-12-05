#!/usr/bin/env python

from __future__ import division
import os,sys,math

original = "../build/original/hcd_s"
parallel = "../build/parallel/hcd_p"

scales = [ 240, 360, 480, 600, 720, 840, 960 ]
pictures = [0]
threads = [1, 2, 4, 8, 10, 12, 14, 16]

N = 100
minr = 32
maxr = 96

def run(*args):
    lines = os.popen(' '.join(map(str, args))).readlines()
    return lines

def extract(*args):
    lines = run(*args)
    time = filter(lambda x: x.isdigit(), lines[0].split())
    return int(time[0])

def average(*args):
    times = [ extract(*args) for i in range(N) ]
    avg = sum(times) / float(N)
    sd = math.sqrt(sum(map(lambda x: (x-avg)**2, times)) / float(N))
    CI = 1.96 * sd / math.sqrt(N)
    return avg, CI

if __name__ == '__main__':

    for pic in pictures:
        for scale in scales:
            with open("data/%d-eye%d.csv" % (scale, pic), "w") as f:
                print >> f, "size,threads,sequential_time,parallel_time,speedup,speedup_ci_low,speedup_ci_high"
                image = "../scales/%d-eye%d.jpg" % (scale, pic)
                for thread in threads:
                    print "Run HCD on <%s> with %d threads" % (image, thread)
                    os.environ["OMP_NUM_THREADS"] = str(thread)
                    s_time, s_CI = average(original, "--source=%s" % image, "--minr=%d" % minr, "--maxr=%d" % maxr)
                    p_time, p_CI = average(parallel, "--source=%s" % image, "--minr=%d" % minr, "--maxr=%d" % maxr)
                    print >> f, "%d,%d,%.2f,%.2f,%.2f,%.2f,%.2f" % (scale, thread, s_time, p_time, s_time / p_time, s_time / (p_time + p_CI), s_time / (p_time - p_CI))
