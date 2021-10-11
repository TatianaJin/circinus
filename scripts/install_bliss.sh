#! /usr/bin/env bash

ROOT=`dirname $(dirname $(realpath ${BASH_SOURCE}))`
VERSION=0.73
THIRDPARTY_DIR=$ROOT/third_party

download() {
  if wget http://www.tcs.hut.fi/Software/bliss/bliss-${VERSION}.zip; then
    unzip bliss-${VERSION}.zip
  else
    echo Error downloading bliss from 'http://www.tcs.hut.fi/Software/bliss/bliss-'${VERSION}'.zip'
  fi
}

install() {
  mv bliss-${VERSION} $THIRDPARTY_DIR/bliss && cd $THIRDPARTY_DIR/bliss && make -j10
}

echo Installing bliss-${VERSION} to $THIRDPARTY_DIR/bliss
download && rm bliss-${VERSION}.zip && install && echo Install succeeded || echo Install Failure
