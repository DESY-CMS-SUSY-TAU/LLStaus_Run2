#!/bin/bash

FILE=$1

root -l $FILE <<< $'.ls\nEvents->Print()\n'
