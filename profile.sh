############################################
# ./profile <experiment name> <S> <M> <cmd>
############################################

expName=$1
S=$2
M=$3

r=$RANDOM

timeLog="tmp.$r.time"
out="tmp.$r.out"


/usr/bin/time -v  -o $timeLog ${@:4}  > $out

InsertionTime=$( grep "InsertionTime" $out | sed -e 's/InsertionTime\t\([^ ]*\) seconds/\1/')
QueryTime=$( grep "QueryTime" $out | sed -e 's/QueryTime\t\([^ ]*\) seconds/\1/')
nMerges=$( grep "nMerges" $out | sed -e 's/nMerges\t\([^ ]*\)/\1/')
sysTime=$(grep "System time"  $timeLog | cut -f2 -d:|tr -d ' ')
usrTime=$(grep "User time"  $timeLog | cut -f2 -d:|tr -d ' ')
wallClock=$(grep "Elapsed (wall clock) time"  $timeLog| sed -e 's/h:mm:ss//'  |sed -e 's/m:ss):/;/'|cut -f2 -d';' |tr -d $' ')
mem=$(grep "Maximum resident set size"  $timeLog | cut -f2 -d:|tr -d ' ')
pgFault=$(grep "Major (requiring I/O) page faults"  $timeLog | cut -f2 -d:|tr -d ' ')
cqfSize=$(ls -lsah ${@: -1}.ser |cut -f6 -d' ')

rm tmp.$r.*

echo -e "$expName\t$S\t$M\t$InsertionTime\t$QueryTime\t$nMerges\t$sysTime\t$usrTime\t$wallClock\t$mem\t$pgFault\t$cqfSize"
