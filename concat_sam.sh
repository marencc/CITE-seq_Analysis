#!/bin/bash
set -e

d=$1
o=$2
n_header_lines="6"
if [ -z "$d" ];
then
    d="./BAM/sjrAlign/"
fi
if [ ! -d "$d" ];
then
    echo "$d directory not found" && exit 1
fi
if [ -z "$o" ];
then
    o="./BAM/sjrAligned.sam"
fi

function error_exit
{
        echo "$1" 1>&2
        exit 1
}
function v_exe
{
    echo "$1"
    eval "$1" || error_exit "Cannot execute command: $1"
}

files=$( ls $d/*.sam )

for f in $files; do # f1= "${files[0]}"
    v_exe "head -n $n_header_lines $f > $o"
    break    
done

n_skip_lines=$(($n_header_lines + 1))
# echo $n_skip_lines
for f in $files; do
    v_exe "tail -n +$n_skip_lines $f >> $o"    
done

echo "Concatenated sam file written to $o"
