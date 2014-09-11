TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    src/AssertException.cpp \
    src/Checks.cpp \
    src/CurrentTest.cpp \
    src/DeferredTestReporter.cpp \
    src/DeferredTestResult.cpp \
    src/MemoryOutStream.cpp \
    src/ReportAssert.cpp \
    src/Test.cpp \
    src/TestDetails.cpp \
    src/TestList.cpp \
    src/TestReporter.cpp \
    src/TestReporterStdout.cpp \
    src/TestResults.cpp \
    src/TestRunner.cpp \
    src/TimeConstraint.cpp \
    src/Win32/TimeHelpers.cpp \
    src/XmlTestReporter.cpp \
    src/lib/lib.cpp


HEADERS += \
    src/AssertException.h \
    src/CheckMacros.h \
    src/Checks.h \
    src/Config.h \
    src/CurrentTest.h \
    src/DeferredTestReporter.h \
    src/DeferredTestResult.h \
    src/ExecuteTest.h \
    src/MemoryOutStream.h \
    src/ReportAssert.h \
    src/Test.h \
    src/TestDetails.h \
    src/TestList.h \
    src/TestMacros.h \
    src/TestReporter.h \
    src/TestReporterStdout.h \
    src/TestResults.h \
    src/TestRunner.h \
    src/TestSuite.h \
    src/TimeConstraint.h \
    src/Win32/TimeHelpers.h \
    src/UnitTest++.h \
    src/XmlTestReporter.h \
    src/lib/lib.hpp
