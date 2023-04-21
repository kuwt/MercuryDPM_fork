#!/bin/sh
# Updates all Mercury licenses (in files) for 2023
for i in $(find Configuration Driver Tools Kernel -depth -type f);
    do
        sed -e 's/Copyright (c) \(20..\)-20../Copyright (c) \1-2023/g'  -i $i;
        sed -e 's/Copyright (c) \(2013\) /Copyright (c) \1-2023/g' -i "" $i;
done
