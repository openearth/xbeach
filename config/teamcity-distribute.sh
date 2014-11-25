#!/bin/bash
. /etc/profile

chmod 775 /opt/teamcity/work/XBeach_unix/install

cd /opt/teamcity/work/XBeach_unix/install

scp -r -i $HOME/xbeach/xbeach-key bin lib geer@h5:/opt/xbeach/trunk_gcc_4.9.1_1.8.1_HEAD
# share trunk 