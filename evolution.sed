#!/usr/bin/sed -f
/^$/d
s/^Gen \([0-9]\+\)[^{]*{\([^}]*\)} \[g\([0-9]*\)\], \([0-9]*\) .*/%\1 \2 \3 \4/
/^[^%].*/d
s/^%//
s/,/ /g