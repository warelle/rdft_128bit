#!/bin/bash

MATRIX_SIZE=2
LIMIT_SIZE=512

function binary_make(){
  make clean
  make
  mv main "./binary/main_$MATRIX_SIZE"
  touch main
}

function exec_program(){
  "./binary/main_$MATRIX_SIZE" > "./result/$MATRIX_SIZE.dat" &
}

while [ $MATRIX_SIZE -le $LIMIT_SIZE ]
do
  cat "test.h" | sed '4d'| sed "4i\#define MATRIX_SIZE $MATRIX_SIZE" > "test_.h"
  rm "test.h"
  mv "test_.h" "test.h"

  binary_make
  exec_program

  ((MATRIX_SIZE += 2))
done
