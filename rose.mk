
# Example Makefile for ROSE users
# Modified by Didem Unat 
#
# If ChiLL is available, set path for CHILL
CHILL = /opt/chill_rose
OMEGA = /opt/omega_rose

# Location of Rose Install Dir 
#ROSE_INSTALL_DIR = /opt/rose/compileTree/current
ROSE_INSTALL_DIR = /opt/rose/compileTree/rose-edg4x-June-22-2015
# Location of Rose Src dir 
#ROSE_SRC_DIR = /opt/rose/src/current
ROSE_SRC_DIR = /opt/rose/src/rose-edg4x-June-22-2015
SSA_ext = /opt/rose/SSA_ext

# Location of include directory after "make install"
ROSE_INCLUDE_DIR = -I$(ROSE_INSTALL_DIR)/include -I$(ROSE_INSTALL_DIR)/include/rose -I/usr/include/python2.6 -I$(CHILL) -I$(CHILL)/include -I$(OMEGA) -I$(OMEGA)/include -I$(OMEGA)/include/code_gen/include -I$(OMEGA)/include/omega/include -I$(OMEGA)/include/basic/include

# Location of Boost Install Dir 
BOOST_INSTALL_DIR = /opt/boost/compileTree/boost-1-50-0/
#BOOST_INSTALL_DIR = /home/dunat/software/boost_1_45_0/installTree

# Location of Libtool
LIB_TOOL = $(ROSE_INSTALL_DIR)/libtool

# Location of Boost include directory
BOOST_INCLUDES =  -I$(BOOST_INSTALL_DIR)/include

JAVA_INCLUDE	=-I/opt/java/jdk1.7.0_79//jre/lib/amd64/server
#/opt/java/jdk1.7.0_79/ 
#/opt/java/jdk1.7.0_79/jre/lib/amd64/server
#JAVA_INCLUDE	=-I/home/dunat/software/java-6-sun/include/ -I/usr/lib/jvm/java-6-sun/include/linux

# Location of Dwarf include and lib (if ROSE is configured to use Dwarf)
ROSE_DWARF_INCLUDES = 
ROSE_DWARF_LIBS_WITH_PATH = 
ROSE_INCLUDE_DIR += $(ROSE_DWARF_INCLUDES)
ROSE_LIBS += $(ROSE_DWARF_LIBS_WITH_PATH)



CC                    = gcc
CXX                   = g++
ROSE_INCLUDE               += $(BOOST_INCLUDES) $(JAVA_INCLUDE) $(ROSE_INCLUDE_DIR) $(STDDEF_INCLUDE) -I$(SSA_ext)/src/arraySSA/
CXXFLAGS              += -w -fopenmp -pthread
LDFLAGS               =  -fopenmp -pthread 

# Location of library directory after "make install"
ROSE_LIB_DIR = $(ROSE_INSTALL_DIR)/lib
ROSE_LIBS = $(ROSE_LIB_DIR)/librose.la 
#$(CHILL)/libchill_xform.a
ROSE_LIBS_so = $(ROSE_LIB_DIR)/librose.so

