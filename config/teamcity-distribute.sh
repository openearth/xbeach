#!/bin/bash
. /etc/profile

printenv | grep xbeach_ssh_key_public | tr -d "'" | sed 's/\\n/\n/g'

mkdir -p $HOME/xbeach
printenv | grep xbeach_ssh_key_private | tr -d "'" | sed 's/\\n/\n/g' > $HOME/xbeach/xbeach-key
printenv | grep xbeach_ssh_key_public | tr -d "'" | sed 's/\\n/\n/g' > $HOME/xbeach/xbeach-key.pub

chmod 775 /opt/teamcity/work/XBeach_unix/install

cd /opt/teamcity/work/XBeach_unix/install

scp -r -i $HOME/xbeach/xbeach-key bin lib geer@h5:/opt/xbeach/trunk_gcc_4.9.2_1.8.3_HEAD

# share trunk 