TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    CylindricalPore/cylindricalpore.cpp

HEADERS += \
    CylindricalPore/cylindricalpore.h

LIBS += -larmadillo -lblas -llapack
