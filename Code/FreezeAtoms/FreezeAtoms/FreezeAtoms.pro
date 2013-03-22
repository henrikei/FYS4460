TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    CylindricalPore/pores.cpp \
    Pores/pores.cpp

HEADERS += \
    CylindricalPore/pores.h \
    Pores/pores.h

LIBS += -larmadillo -lblas -llapack
