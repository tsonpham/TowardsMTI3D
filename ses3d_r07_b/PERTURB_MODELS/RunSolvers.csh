#! /bin/csh

set homedir = `pwd`
foreach model_dname (GAUSSIAN_3D/[5]/MODEL0?1)
    cd $homedir/$model_dname
    set num = `ls MAIN | grep _forward. | wc -l`
    if ($num > 0) then 
        continue 
    endif

    echo $model_dname
    python3 /home/158/tp1732/SES3D/ses3d_r07_b/RunForwardGF.py

    cd $homedir/$model_dname/MAIN
    qsub run_ses3d.pbs
end
