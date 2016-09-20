#-------------------------------------------------
#
# Project created by QtCreator 2016-06-10T10:30:45
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = RandomAggLargeWithSnap
CONFIG   += console
CONFIG   += c++11
CONFIG   -= app_bundle

QMAKE_LFLAGS_WINDOWS += -Wl,--stack,100000000

TEMPLATE = app

DEFINES -= UNICODE

SOURCES += main.cpp \
    vertex.cpp \
    edge.cpp \
    mygraph.cpp \
    include/snap/glib-core/app.cpp \
    include/snap/glib-core/base.cpp \
    include/snap/glib-core/bd.cpp \
    include/snap/glib-core/bits.cpp \
    include/snap/glib-core/blobbs.cpp \
    include/snap/glib-core/console.cpp \
    include/snap/glib-core/dt.cpp \
    include/snap/glib-core/env.cpp \
    include/snap/glib-core/exp.cpp \
    include/snap/glib-core/fl.cpp \
    include/snap/glib-core/gnuplot.cpp \
    include/snap/glib-core/hash.cpp \
    include/snap/glib-core/html.cpp \
    include/snap/glib-core/http.cpp \
    include/snap/glib-core/json.cpp \
    include/snap/glib-core/linalg.cpp \
    include/snap/glib-core/lx.cpp \
    include/snap/glib-core/macro.cpp \
    include/snap/glib-core/md5.cpp \
    include/snap/glib-core/os.cpp \
    include/snap/glib-core/pp.cpp \
    include/snap/glib-core/ss.cpp \
    include/snap/glib-core/tm.cpp \
    include/snap/glib-core/unicode.cpp \
    include/snap/glib-core/unicodestring.cpp \
    include/snap/glib-core/url.cpp \
    include/snap/glib-core/ut.cpp \
    include/snap/glib-core/wch.cpp \
    include/snap/glib-core/xdt.cpp \
    include/snap/glib-core/xfl.cpp \
    include/snap/glib-core/xmath.cpp \
    include/snap/glib-core/xml.cpp \
    include/snap/glib-core/zipfl.cpp \
    include/snap/snap-adv/agm.cpp \
    include/snap/snap-adv/agmfast.cpp \
    include/snap/snap-adv/agmfit.cpp \
    include/snap/snap-adv/cascdynetinf.cpp \
    include/snap/snap-adv/cascnetinf.cpp \
    include/snap/snap-adv/cliques.cpp \
    include/snap/snap-adv/graphcounter.cpp \
    include/snap/snap-adv/kronecker.cpp \
    include/snap/snap-adv/mag.cpp \
    include/snap/snap-adv/ncp.cpp \
    include/snap/snap-adv/subgraphenum.cpp \
    include/snap/snap-core/alg.cpp \
    include/snap/snap-core/anf.cpp \
    include/snap/snap-core/centr.cpp \
    include/snap/snap-core/cmty.cpp \
    include/snap/snap-core/cncom.cpp \
    include/snap/snap-core/ff.cpp \
    include/snap/snap-core/gbase.cpp \
    include/snap/snap-core/ggen.cpp \
    include/snap/snap-core/ghash.cpp \
    include/snap/snap-core/gio.cpp \
    include/snap/snap-core/graph.cpp \
    include/snap/snap-core/gstat.cpp \
    include/snap/snap-core/gsvd.cpp \
    include/snap/snap-core/gviz.cpp \
    include/snap/snap-core/Snap.cpp \
    include/snap/snap-core/statplot.cpp \
    include/snap/snap-core/subgraph.cpp \
    include/snap/snap-core/timenet.cpp \
    include/snap/snap-core/util.cpp

HEADERS += \
    vertex.h \
    edge.h \
    mygraph.h \
    include/snap/glib-core/app.h \
    include/snap/glib-core/base.h \
    include/snap/glib-core/bd.h \
    include/snap/glib-core/bits.h \
    include/snap/glib-core/blobbs.h \
    include/snap/glib-core/console.h \
    include/snap/glib-core/ds.h \
    include/snap/glib-core/dt.h \
    include/snap/glib-core/env.h \
    include/snap/glib-core/exp.h \
    include/snap/glib-core/fds.h \
    include/snap/glib-core/fl.h \
    include/snap/glib-core/gnuplot.h \
    include/snap/glib-core/hash.h \
    include/snap/glib-core/html.h \
    include/snap/glib-core/http.h \
    include/snap/glib-core/json.h \
    include/snap/glib-core/linalg.h \
    include/snap/glib-core/lx.h \
    include/snap/glib-core/macro.h \
    include/snap/glib-core/md5.h \
    include/snap/glib-core/os.h \
    include/snap/glib-core/pp.h \
    include/snap/glib-core/shash.h \
    include/snap/glib-core/ss.h \
    include/snap/glib-core/stdafx.h \
    include/snap/glib-core/tm.h \
    include/snap/glib-core/unicode.h \
    include/snap/glib-core/unicodestring.h \
    include/snap/glib-core/url.h \
    include/snap/glib-core/ut.h \
    include/snap/glib-core/wch.h \
    include/snap/glib-core/xdt.h \
    include/snap/glib-core/xfl.h \
    include/snap/glib-core/xmath.h \
    include/snap/glib-core/xml.h \
    include/snap/glib-core/xmlser.h \
    include/snap/glib-core/zipfl.h \
    include/snap/snap-adv/agm.h \
    include/snap/snap-adv/agmfast.h \
    include/snap/snap-adv/agmfit.h \
    include/snap/snap-adv/cascdynetinf.h \
    include/snap/snap-adv/cascnetinf.h \
    include/snap/snap-adv/circles.h \
    include/snap/snap-adv/cliques.h \
    include/snap/snap-adv/graphcounter.h \
    include/snap/snap-adv/kronecker.h \
    include/snap/snap-adv/mag.h \
    include/snap/snap-adv/ncp.h \
    include/snap/snap-adv/subgraphenum.h \
    include/snap/snap-core/alg.h \
    include/snap/snap-core/anf.h \
    include/snap/snap-core/bfsdfs.h \
    include/snap/snap-core/bignet.h \
    include/snap/snap-core/centr.h \
    include/snap/snap-core/cmty.h \
    include/snap/snap-core/cncom.h \
    include/snap/snap-core/ff.h \
    include/snap/snap-core/gbase.h \
    include/snap/snap-core/ggen.h \
    include/snap/snap-core/ghash.h \
    include/snap/snap-core/gio.h \
    include/snap/snap-core/graph.h \
    include/snap/snap-core/gstat.h \
    include/snap/snap-core/gsvd.h \
    include/snap/snap-core/gviz.h \
    include/snap/snap-core/kcore.h \
    include/snap/snap-core/network.h \
    include/snap/snap-core/Snap.h \
    include/snap/snap-core/statplot.h \
    include/snap/snap-core/stdafx.h \
    include/snap/snap-core/subgraph.h \
    include/snap/snap-core/timenet.h \
    include/snap/snap-core/triad.h \
    include/snap/snap-core/util.h

INCLUDEPATH += "C:\Boost\boost_1_56_0"
INCLUDEPATH += C:\Users\Dumex\Qt_WorkingSpace\RandomAggLargeWithSnap\include\snap\snap-core
INCLUDEPATH += C:\Users\Dumex\Qt_WorkingSpace\RandomAggLargeWithSnap\include\snap\glib-core
INCLUDEPATH += C:\Users\Dumex\Qt_WorkingSpace\RandomAggLargeWithSnap\include\snap\snap-adv
