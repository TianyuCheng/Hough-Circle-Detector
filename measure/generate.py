#!/usr/bin/env python

from __future__ import division
import os,sys,math

parallel = "../build/optimized/hcd_p"
serial   = "../build/original/hcd_s"

scales = [ 240, 360, 480, 600, 720, 840, 960 ]
pictures = [0, 1, 2]

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
        s_times = []
        p_times = []
        for scale in scales:
            print "Run HCD on %d-eye%d.jpg" % (scale, pic)
            image = "../scales/%d-eye%d.jpg" % (scale, pic)
            s_time, s_CI = average(serial, "--source=%s" % image, "--minr=%d" % minr, "--maxr=%d" % maxr)
            p_time, p_CI = average(parallel, "--source=%s" % image, "--minr=%d" % minr, "--maxr=%d" % maxr)
            s_times.append(s_time)
            p_times.append(p_time)

        csv = "eye%d.csv" % pic

        with open("data/%s" % csv, "w") as f:
            print >> f, "size,s_time,p_time,speedup"
            for i, scale in enumerate(scales):
                print >> f, "%d,%d,%d,%f" % (scale, s_times[i], p_times[i], s_times[i] / p_times[i])
