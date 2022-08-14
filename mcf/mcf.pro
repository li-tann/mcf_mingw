QT -= gui

TARGET = mcf

CONFIG += c++11 console
CONFIG -= app_bundle

CONFIG(debug, debug|release){
    DESTDIR     =   $$PWD/../../bin/debug/plugin
    OBJECTS_DIR =   ./debug/OBJ
}else{
    DESTDIR     =   $$PWD/../../bin/release/plugin
    OBJECTS_DIR =   ./release/OBJ
}

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    mcf.c \
    SWAP_io.c \
    DISP_lib.c \
    mssptr.c \
    m_tri_init2.c \
    m_ssp2_d_heap.c \
    triangle.c


# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

HEADERS += \
    bmp_image.h \
    display.h \
    mcf2.h \
    palsar.h \
    rasterfile.h \
    triangle.h
