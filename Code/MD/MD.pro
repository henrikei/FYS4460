TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    atom.cpp \
    system.cpp \
    lib.cpp \
    cell.cpp \
    Modifier/modifier.cpp \
    Modifier/berendsenthermostat.cpp \
    Modifier/andersenthermostat.cpp

HEADERS += \
    atom.h \
    system.h \
    lib.h \
    cell.h \
    Modifier/modifier.h \
    Modifier/berendsenthermostat.h \
    Modifier/andersenthermostat.h

LIBS += -larmadillo -lblas -llapack

COMMON_FLAGS = -std=c++0x
QMAKE_CXXFLAGS += $$COMMON_FLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_FLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_FLAGS
