#!/bin/sh
TO=util/buildno.c
echo '#include <stdio.h>' > $TO
echo '#include "util.h"' >> $TO
echo 'int buildno() { return ' >> $TO
NO=`cat buildno | awk '{ print($1+1) }'`
echo $NO > buildno
cat buildno >> $TO
echo '; }' >> $TO
