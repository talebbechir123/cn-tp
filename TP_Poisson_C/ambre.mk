#######################################
# ambre.mk
# Default options for ambre computer
#######################################
CC=gcc
LIBSLOCAL=-L/usr/include/aarch64-linux-gnu -llapack -lblas -lm
INCLUDEBLASLOCAL=-I/usr/include/aarch64-linux-gnu
OPTCLOCAL=-fPIC -march=native