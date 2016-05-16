#!/bin/bash
#This is a build script that helps with the initial build as described in README file
aclocal
autoheader
automake -ac
autoconf
./configure
make
make install

