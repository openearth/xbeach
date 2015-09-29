#!/bin/bash
. /etc/profile

chmod 775 /opt/teamcity/work/XBeach_unix_trunk/install

cd /opt/teamcity/work/XBeach_unix_trunk/install

copypath="/opt/xbeach/"$XBEACH_PROJECT_ID"_gcc_4.9.2_1.8.3_HEAD"

scp -r -i $HOME/xbeach/xbeach-key bin lib geer@h5:$copypath

cd /opt/teamcity/work/XBeach_unix_trunk/trunk/config

modulepath="/opt/xbeach/"$XBEACH_PROJECT_ID"_gcc_4.9.2_1.8.3_HEAD"

scp -r -i $HOME/xbeach/xbeach-key "xbeach-"$XBEACH_PROJECT_ID"_gcc_4.9.2_1.8.3_HEAD" geer@h5:/opt/xbeach/modules/

# share trunk 