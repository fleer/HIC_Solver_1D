#!/bin/sh
set -ex
wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.4.tar.gz
tar -zxvf gsl-2.4.tar.gz
cd gsl-2.4 && ./configure --prefix=/usr && make && sudo make install
