#!/usr/bin/bash

export CAFFE_DIR=$HOME/Projects/CAFFE
export PATH=$PATH:$CAFFE_DIR/bin

# Compilers

export CAFFE_CC=gcc
export CAFFE_CXX=g++

# Aliases

alias caffe='cd $CAFFE_DIR'
