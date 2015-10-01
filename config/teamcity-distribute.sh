#!/bin/bash
. /etc/profile

export SVNPATH=$(pwd)

chmod 775 $SVNPATH"/install"

cd $SVNPATH"/install"

copypath="/opt/xbeach/"$XBEACH_PROJECT_ID"_gcc_4.9.2_1.8.3_HEAD"

scp -r -i $HOME/xbeach/xbeach-key bin lib geer@h5:$copypath

cd $SVNPATH"/trunk/config"

scp -r -i $HOME/xbeach/xbeach-key "xbeach-"$XBEACH_PROJECT_ID"_gcc_4.9.2_1.8.3_HEAD" geer@h5:/opt/xbeach/modules/