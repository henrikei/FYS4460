TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    atom.cpp \
    system.cpp \
    lib.cpp

HEADERS += \
    atom.h \
    system.h \
    lib.h

LIBS += -larmadillo -lblas -llapack

COMMON_FLAGS = -std=c++0x
QMAKE_CXXFLAGS += $$COMMON_FLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_FLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_FLAGS
