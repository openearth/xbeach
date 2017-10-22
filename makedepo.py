#!/usr/bin/env python
import sys
import re
import getopt
import os


def sortu(a):
    if len(a) == 0:
        return a
    x = a.split()
    x.sort()
    l = [x[0]]
    u = x[0]
    for p in x[1:]:
        if l[-1] != p:
            l = l + [p]
            u = u + " " + p
    return u

def removenil(x):
    u = {}
    for i in x:
        if x[i]:
            u[i] = x[i]
    return u

def deps(mode):
    if mode == "M:F":                    # m1.mod m2.mod: sub.f90
        provn = removenil(prov1)
        for f in provn:
            print(provn[f] + ": " + f)
        return

    for f in src:
        if mode == "O:FHO":            # sub.o: sub.f90 inc.h sub1.o
            rr = req[f].split()
            d = ""
            for r in rr:
                try:
                    if prov[r] != f:
                        d = d + " " + prov[r]
                except KeyError:
                    pass
            print(f + ": " + src[f] + " " + inc[f] + " " + sortu(d))
        if mode == "O:FH":             # sub.o: sub.f90 inc.h
            print(f + ": " + " " + inc[f])
        if mode == "O:M":                # sub.o: m1.mod m2.mod
            print(f + ": " + req[f])

# start of main program
if __name__ == '__main__':
    prov = {}
    prov1= {}
    req = {}
    inc = {}
    src = {}
    findsuffix = re.compile(r'.[^.]*$')
    findquote = re.compile(r'[\'"]')

    optlist, args = getopt.getopt(sys.argv[1:],'m:s:p:')

    mode = "O:FHO"
    prefix=""
    suffix=".o"

    for o,a in optlist:
        if o == "-m":
            mode = a
        if o == "-p":
            prefix = a
        if o == "-s":
            suffix = a



    for f in args :
        try:
            srcfile=open(f)
            basename = re.sub(findsuffix,'',f)
            oname = prefix+basename+suffix
            src[oname] = f
            req[oname] = ""
            prov1[f] = ""
            inc[oname] = ""
            for line in srcfile:
                line = line.strip().split()
                word1 = ""
                if len(line):
                    word1 = line[0]
                word2 = ""
                if len(line) > 1:
                    word2 = line[1]
                if word1.lower() == "module":
                    if word2.lower() == "procedure":
                        continue
                    prov[word2.lower()+".mod"] = oname
                    if prov1[f]:
                        prov1[f] = prov1[f] + " " + word2.lower() + ".mod"
                    else:
                        prov1[f] = word2.lower() + ".mod"
                if word1.lower() == "use":
                    cp = word2.rfind(",")
                    if cp != -1 :
                        word2 = word2[:cp]
                    req[oname] = req[oname] + " " + word2.lower() + ".mod"
                if word1 == "include":
                    # do not use dependencies of files that are not visible
                    word2 = re.sub(findquote,'',word2)
                    #if os.path.exists(word2):
                    if True:
                        inc[oname] = inc[oname] + " " + word2
            srcfile.close()
        except IOError:
            pass


    #inc    = removenil(inc)

    for i in inc:
        inc[i] = sortu(inc[i])

    for i in prov:
        prov[i] = sortu(prov[i])

    for i in req:
        req[i] = sortu(req[i])

    deps(mode)
    #for mode in ["O:FHO" , "O:FH", "O:M", "M:F" ]:
    #    print "testing mode",mode
    #    deps(mode)


