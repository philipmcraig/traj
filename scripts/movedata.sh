#!/bin/bash

CLUSDIR=/glusterfs/scenario/users/np838619 # met-cluster directory
TJDIR=$CLUSDIR/traj/2014

cd $TJDIR

TJCODE1=$TJDIR/amr/TCb
TJCODE2=$TJDIR/afr/TCb
TJCODE3=$TJDIR/eaa/TCb
TJCODE4=$TJDIR/ara/TCb
TJCODE5=$TJDIR/ari/TCb
TJCODE6=$TJDIR/arp/TCb
TJCODE7=$TJDIR/soa/TCb
TJCODE8=$TJDIR/soi/TCb
TJCODE9=$TJDIR/sop/TCb

TJDATA1=$TJDIR/amr/TD
TJDATA2=$TJDIR/afr/TD
TJDATA3=$TJDIR/eaa/TD
TJDATA4=$TJDIR/ara/TD
TJDATA5=$TJDIR/ari/TD
TJDATA6=$TJDIR/arp/TD
TJDATA7=$TJDIR/soa/TD
TJDATA8=$TJDIR/soi/TD
TJDATA9=$TJDIR/sop/TD

mv $TJCODE1/utraj* $TJDATA1
mv $TJCODE2/utraj* $TJDATA2
mv $TJCODE3/utraj* $TJDATA3
mv $TJCODE4/utraj* $TJDATA4
mv $TJCODE5/utraj* $TJDATA5
mv $TJCODE6/utraj* $TJDATA6
mv $TJCODE7/utraj* $TJDATA7
mv $TJCODE8/utraj* $TJDATA8
mv $TJCODE9/utraj* $TJDATA9
