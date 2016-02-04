#!/bin/bash
. /etc/profile

export SVNPATH=$(pwd)

chmod 775 $SVNPATH"/install"

cd $SVNPATH"/install"

copypath="/opt/xbeach/"$XBEACH_PROJECT_ID"_gcc_4.9.2_1.8.3_HEAD"

scp -r -i $HOME/.ssh/id_rsa bin lib xbeach@$CLUSTER_HOSTNAME:$copypath

cd $SVNPATH"/trunk/config"

scp -r -i $HOME/.ssh/id_rsa "xbeach-"$XBEACH_PROJECT_ID"_gcc_4.9.2_1.8.3_HEAD" xbeach@$CLUSTER_HOSTNAME:/opt/xbeach/modules/