TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    Pores/pores.cpp \
    Pores/cylindricalpore.cpp \
    Pores/circularpores.cpp

HEADERS += \
    Pores/pores.h \
    Pores/cylindricalpore.h \
    Pores/circularpores.h

LIBS += -larmadillo -lblas -llapack
