#!/bin/bash
. /etc/profile

chmod 775 /opt/teamcity/work/XBeach_unix/install

cd /opt/teamcity/work/XBeach_unix/install
scp -r -i $HOME/xbeach/xbeach-key bin lib geer@h5:tcxbeach
# share trunk 