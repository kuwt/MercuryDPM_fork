#!/bin/sh
# Updates all Mercury licenses (in files) for 2017
for i in `find Configuration Drivers Kernel Tools -type f`; 
    do 
        sed -e 's/Copyright (c) \(20..\)-20../Copyright (c) \1-2020/g' -i "" $i;
        sed -e 's/Copyright (c) \(20..\) /Copyright (c) \1-2020/g' -i "" $i;
        sed -i '' -e '$a\' $i
done
