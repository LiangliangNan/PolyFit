CONFIG -= qt
TARGET = 3rd_soplex
TEMPLATE = lib
CONFIG += staticlib

win32 { DEFINES += WIN32 WIN64 _WINDOWS }


CONFIG(debug, debug|release) { DEFINES += _DEBUG   NO_SIGACTION NO_STRTOK_R _CRT_SECURE_NO_WARNINGS TPI_NONE NPARASCIP WITH_SCIPDEF ROUNDING_FE}
CONFIG(release, debug|release) { DEFINES += NDEBUG NO_SIGACTION NO_STRTOK_R _CRT_SECURE_NO_WARNINGS TPI_NONE NPARASCIP WITH_SCIPDEF ROUNDING_FE}

INCLUDEPATH += \


SOURCES +=  \
    src/changesoplex.cpp \
    src/clufactor.cpp \
    src/clufactor_rational.cpp \
    src/didxset.cpp \
    src/enter.cpp \
    src/git_hash.cpp \
    src/gzstream.cpp \
    src/idxset.cpp \
    src/leave.cpp \
    src/mpsinput.cpp \
    src/nameset.cpp \
    src/rational.cpp \
    src/ratrecon.cpp \
    src/slufactor.cpp \
    src/slufactor_rational.cpp \
    src/solvedbds.cpp \
    src/solverational.cpp \
    src/solvereal.cpp \
    src/soplex.cpp \
    src/soplexlegacy.cpp \
    src/spxautopr.cpp \
    src/spxbasis.cpp \
    src/spxboundflippingrt.cpp \
    src/spxbounds.cpp \
    src/spxchangebasis.cpp \
    src/spxdantzigpr.cpp \
    src/spxdefaultrt.cpp \
    src/spxdefines.cpp \
    src/spxdesc.cpp \
    src/spxdevexpr.cpp \
    src/spxequilisc.cpp \
    src/spxfastrt.cpp \
    src/spxfileio.cpp \
    src/spxgeometsc.cpp \
    src/spxgithash.cpp \
    src/spxharrisrt.cpp \
    src/spxhybridpr.cpp \
    src/spxid.cpp \
    src/spxleastsqsc.cpp \
    src/spxlpbase_rational.cpp \
    src/spxlpbase_real.cpp \
    src/spxmainsm.cpp \
    src/spxout.cpp \
    src/spxparmultpr.cpp \
    src/spxquality.cpp \
    src/spxscaler.cpp \
    src/spxshift.cpp \
    src/spxsolve.cpp \
    src/spxsolver.cpp \
    src/spxstarter.cpp \
    src/spxsteeppr.cpp \
    src/spxsumst.cpp \
    src/spxvecs.cpp \
    src/spxvectorst.cpp \
    src/spxweightpr.cpp \
    src/spxweightst.cpp \
    src/spxwritestate.cpp \
    src/statistics.cpp \
    src/testsoplex.cpp \
    src/updatevector.cpp \
    src/usertimer.cpp \
    src/validation.cpp \
    src/wallclocktimer.cpp

HEADERS +=  \
    src/array.h \
    src/basevectors.h \
    src/classarray.h \
    src/clufactor.h \
    src/clufactor_rational.h \
    src/cring.h \
    src/dataarray.h \
    src/datahashtable.h \
    src/datakey.h \
    src/dataset.h \
    src/didxset.h \
    src/dsvector.h \
    src/dsvectorbase.h \
    src/dvector.h \
    src/dvectorbase.h \
    src/exceptions.h \
    src/gzstream.h \
    src/idlist.h \
    src/idxset.h \
    src/islist.h \
    src/lpcol.h \
    src/lpcolbase.h \
    src/lpcolset.h \
    src/lpcolsetbase.h \
    src/lprow.h \
    src/lprowbase.h \
    src/lprowset.h \
    src/lprowsetbase.h \
    src/mpsinput.h \
    src/nameset.h \
    src/notimer.h \
    src/random.h \
    src/rational.h \
    src/ratrecon.h \
    src/slinsolver.h \
    src/slinsolver_rational.h \
    src/slufactor.h \
    src/slufactor_rational.h \
    src/sol.h \
    src/solbase.h \
    src/soplex.h \
    src/soplexlegacy.h \
    src/sorter.h \
    src/spxalloc.h \
    src/spxautopr.h \
    src/spxbasis.h \
    src/spxboundflippingrt.h \
    src/spxdantzigpr.h \
    src/spxdefaultrt.h \
    src/spxdefines.h \
    src/spxdevexpr.h \
    src/spxequilisc.h \
    src/spxfastrt.h \
    src/spxfileio.h \
    src/spxgeometsc.h \
    src/spxgithash.h \
    src/spxharrisrt.h \
    src/spxhybridpr.h \
    src/spxid.h \
    src/spxleastsqsc.h \
    src/spxlp.h \
    src/spxlpbase.h \
    src/spxmainsm.h \
    src/spxout.h \
    src/spxparmultpr.h \
    src/spxpricer.h \
    src/spxratiotester.h \
    src/spxscaler.h \
    src/spxsimplifier.h \
    src/spxsolver.h \
    src/spxstarter.h \
    src/spxsteepexpr.h \
    src/spxsteeppr.h \
    src/spxsumst.h \
    src/spxvectorst.h \
    src/spxweightpr.h \
    src/spxweightst.h \
    src/ssvector.h \
    src/ssvectorbase.h \
    src/statistics.h \
    src/svector.h \
    src/svectorbase.h \
    src/svset.h \
    src/svsetbase.h \
    src/timer.h \
    src/timerfactory.h \
    src/unitvector.h \
    src/unitvectorbase.h \
    src/updatevector.h \
    src/usertimer.h \
    src/validation.h \
    src/vector.h \
    src/vectorbase.h \
    src/wallclocktimer.h

