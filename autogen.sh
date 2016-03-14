#! /bin/sh

autoreconf -vi

(cd isl; ./autogen.sh)

(cd pet; ./autogen.sh)
