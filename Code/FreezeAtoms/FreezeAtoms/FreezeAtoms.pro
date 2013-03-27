TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    Pores/pores.cpp \
    Pores/cylindricalpore.cpp \
    Pores/circularpores.cpp \
    Pores/reducedensity.cpp

HEADERS += \
    Pores/pores.h \
    Pores/cylindricalpore.h \
    Pores/circularpores.h \
    Pores/reducedensity.h

LIBS += -larmadillo -lblas -llapack
