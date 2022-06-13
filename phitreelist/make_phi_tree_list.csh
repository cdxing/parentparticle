#!/bin/csh
rm Phi_SE_tree.list
rm Phi_ME_tree.list
#cd /star/u/slan/pwg/fastoffline/7p7gev/phitree/out_Phi
#cd /star/data01/pwg/dchen/Ana/19p6GeV/parentparticle/production
cd /star/data05/scratch/dchen/Ana/19p6GeV/parentparticle/production
#ls -d $PWD/*SE.root > /star/u/dchen/ana/19gev_2019/parentparticle/phitreelist/Phi_SE_tree.list
find -type f -name '*SE.root' > /star/u/dchen/ana/19gev_2019/parentparticle/phitreelist/Phi_SE_tree.list
cd -

#cd /star/u/slan/pwg/fastoffline/7p7gev/phitree/out_Phi
#cd /star/data01/pwg/dchen/Ana/19p6GeV/parentparticle/production
cd /star/data05/scratch/dchen/Ana/19p6GeV/parentparticle/production
#ls -d *_ME.root > /star/u/slan/pwg/fastoffline/7p7gev/phitree/List/Phi_ME_tree.list
#ls -d $PWD/*ME.root > /star/u/dchen/ana/19gev_2019/parentparticle/phitreelist/Phi_ME_tree.list
find -type f -name '*ME.root' > /star/u/dchen/ana/19gev_2019/parentparticle/phitreelist/Phi_ME_tree.list
cd -
