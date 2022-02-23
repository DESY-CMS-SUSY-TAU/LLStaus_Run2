#!/bin/bash

DIR=$1

echo $DIR

astyle \
--recursive \
--mode=c \
--style=allman \
--indent=spaces=4 \
--indent-classes \
--indent-modifiers \
--indent-switches \
--indent-cases \
--indent-namespaces \
--indent-after-parens \
--add-braces \
--break-blocks \
--pad-comma \
--fill-empty-lines \
$DIR/*.cc \
$DIR/*.h \
