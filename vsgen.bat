MKDIR build\windows
CD build\windows
cmake -DCMAKE_CONFIGURATION_TYPES="Debug;Release" ..\.. %*
ECHO VS files have been generated in build\windows
CD ..\..\
