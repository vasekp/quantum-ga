#!/usr/bin/sed -f
/^$/d
/^[^{]/d
s/{\([^}]*\)}.*/\1/
s/,/ /g
